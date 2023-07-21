//  Time-stamp: "Last modified 2023-04-24 16:23:42 mluebke"

/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2021 Marco De Lucia (GFZ Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#ifndef DHT_WRAPPER_H
#define DHT_WRAPPER_H

#include "poet/DHT_Types.hpp"
#include <array>
#include <cstdint>
#include <string>
#include <vector>

extern "C" {
#include "DHT.h"
}

#include <mpi.h>

namespace poet {

struct DHT_SCNotation {
  std::int8_t exp : 8;
  std::int64_t significant : 56;
};

union DHT_Keyelement {
  double fp_elemet;
  DHT_SCNotation sc_notation;
};

using DHT_ResultObject = struct DHTResobj {
  uint32_t length;
  std::vector<std::vector<DHT_Keyelement>> keys;
  std::vector<std::vector<double>> results;
  std::vector<bool> needPhreeqc;
};

/**
 * @brief Return user-defined md5sum
 *
 * This function will calculate a hashsum with the help of md5sum. Therefore the
 * md5sum for a given key is calculated and divided into two 64-bit parts. These
 * will be XORed and returned as the hash.
 *
 * @param key_size Size of key in bytes
 * @param key Pointer to key
 * @return uint64_t Hashsum as an unsigned 64-bit integer
 */
static uint64_t get_md5(int key_size, void *key);

/**
 * @brief C++-Wrapper around DHT implementation
 *
 * Provides an API to interact with the current DHT implentation. This class is
 * POET specific and can't be used outside the POET application.
 *
 */
class DHT_Wrapper {
public:
  /**
   * @brief Construct a new dht wrapper object
   *
   * The constructor will initialize the private dht_object of this class by
   * calling DHT_create with all given parameters. Also the fuzzing buffer will
   * be allocated and all needed parameters extracted from simparams struct.
   *
   * @param dht_comm Communicator which addresses all participating DHT
   * processes
   * @param buckets_per_process Count of buckets to allocate for each
   * process
   * @param key_indices Vector indexing elements of one grid cell used
   * for key creation.
   * @param data_count Count of data entries
   */
  DHT_Wrapper(MPI_Comm dht_comm, uint32_t dht_size,
              const std::vector<std::uint32_t> &key_indices,
              uint32_t data_count);
  /**
   * @brief Destroy the dht wrapper object
   *
   * By destroying this object the DHT will also be freed. Since all statistics
   * are stored inside this object, no statistics will be retrieved during the
   * call of DHT_free. After freeing the DHT the fuzzing buffer will be also
   * freed.
   *
   */
  ~DHT_Wrapper();

  /**
   * @brief Check if values of workpackage are stored in DHT
   *
   * Call DHT_read for all grid cells of the given workpackage and if a
   * previously simulated grid cell was found mark this grid cell as 'not be
   * simulated'. Therefore all values of a grid cell are fuzzed by fuzzForDHT
   * and used as input key. The correspondending retrieved value might be stored
   * directly into the memory area of the work_package and out_result_index is
   * marked with false ('not to be simulated').
   *
   * @param length Count of grid cells inside work package
   * @param[out] out_result_index Indexing work packages which should be
   * simulated
   * @param[in,out] work_package Pointer to current work package
   * @param dt Current timestep of simulation
   */
  auto checkDHT(int length, double dt, const std::vector<double> &work_package,
                std::vector<std::uint32_t> &curr_mapping)
      -> const poet::DHT_ResultObject &;

  /**
   * @brief Write simulated values into DHT
   *
   * Call DHT_write for all grid cells of the given workpackage which was
   * simulated shortly before by the worker. Whether the grid cell was simulated
   * is given by result_index. For every grid cell indicated with true inside
   * result_index write the simulated value into the DHT.
   *
   * @param length Count of grid cells inside work package
   * @param result_index Indexing work packages which was simulated
   * @param work_package Pointer to current work package which was used as input
   * of PHREEQC
   * @param results Pointer to current work package which are the resulting
   * outputs of the PHREEQC simulation
   * @param dt Current timestep of simulation
   */
  void fillDHT(int length, const std::vector<double> &work_package);

  void resultsToWP(std::vector<double> &work_package);

  /**
   * @brief Dump current DHT state into file.
   *
   * This function will simply execute DHT_to_file with given file name (see
   * DHT.h for more info).
   *
   * @param filename Name of the dump file
   * @return int Returns 0 on success, otherwise an error value
   */
  int tableToFile(const char *filename);

  /**
   * @brief Load dump file into DHT.
   *
   * This function will simply execute DHT_from_file with given file name (see
   * DHT.h for more info).
   *
   * @param filename Name of the dump file
   * @return int Returns 0 on success, otherwise an error value
   */
  int fileToTable(const char *filename);

  /**
   * @brief Print a detailed statistic of DHT usage.
   *
   * This function will simply execute DHT_print_statistics with given file name
   * (see DHT.h for more info).
   *
   */
  void printStatistics();

  /**
   * @brief Get the Hits object
   *
   * @return uint64_t Count of hits
   */
  auto getHits() { return this->dht_hits; };

  /**
   * @brief Get the Misses object
   *
   * @return uint64_t Count of read misses
   */
  auto getMisses() { return this->dht_miss; };

  /**
   * @brief Get the Evictions object
   *
   * @return uint64_t Count of evictions
   */
  auto getEvictions() { return this->dht_evictions; };

  void SetSignifVector(std::vector<uint32_t> signif_vec);

  void setBaseTotals(const std::array<double, 2> &bt) {
    this->base_totals = bt;
  }

private:
  uint32_t key_count;
  uint32_t data_count;

  DHT *dht_object;

  std::vector<DHT_Keyelement> fuzzForDHT(int var_count, void *key, double dt);

  uint32_t dht_hits = 0;
  uint32_t dht_miss = 0;
  uint32_t dht_evictions = 0;

  std::vector<uint32_t> dht_signif_vector;
  std::vector<std::uint32_t> dht_prop_type_vector;
  std::vector<std::uint32_t> input_key_elements;

  static constexpr int DHT_KEY_SIGNIF_DEFAULT = 5;
  static constexpr int DHT_KEY_SIGNIF_TOTALS = 10;
  static constexpr int DHT_KEY_SIGNIF_CHARGE = 3;

  DHT_ResultObject dht_results;

  std::array<double, 2> base_totals{0};
};
} // namespace poet

#endif // DHT_WRAPPER_H
