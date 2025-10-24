
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

#include "Base/RInsidePOET.hpp"
#include "DataStructures/NamedVector.hpp"

#include "Chemistry/ChemistryDefs.hpp"

#include "Init/InitialList.hpp"
#include "LookupKey.hpp"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include <LUCX/DHT.h>

#include <mpi.h>

namespace poet {

using DHT_Location = std::pair<std::uint32_t, std::uint32_t>;

/**
 * @brief C++-Wrapper around DHT implementation
 *
 * Provides an API to interact with the current DHT implentation. This class is
 * POET specific and can't be used outside the POET application.
 *
 */
class DHT_Wrapper {
public:
  using DHT_ResultObject = struct DHTResobj {
    std::vector<LookupKey> keys;
    std::vector<DHT_Location> locations;
    std::vector<bool> filledDHT;
  };

  static constexpr std::int32_t DHT_KEY_INPUT_CUSTOM =
      std::numeric_limits<std::int32_t>::min();

  static constexpr int DHT_KEY_SIGNIF_DEFAULT = 5;

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
  DHT_Wrapper(MPI_Comm dht_comm, std::uint64_t dht_size,
              const NamedVector<std::uint32_t> &key_species,
              const std::vector<std::int32_t> &key_indices,
              const std::vector<std::string> &output_names,
              const InitialList::ChemistryHookFunctions &hooks,
              uint32_t data_count, bool with_interp, bool has_het_ids);
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

  DHT_Wrapper &operator=(const DHT_Wrapper &) = delete;
  DHT_Wrapper(const DHT_Wrapper &) = delete;

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
  auto checkDHT(WorkPackage &work_package) -> const DHT_ResultObject &;

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
  void fillDHT(const WorkPackage &work_package);

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
   * @brief Get the Evictions object
   *
   * @return uint64_t Count of evictions
   */
  auto getEvictions() { return this->dht_evictions; };

  void resetCounter() {
    this->dht_hits = 0;
    this->dht_evictions = 0;
  }

  void SetSignifVector(std::vector<uint32_t> signif_vec);

  auto getDataCount() { return this->data_count; }
  auto getCommunicator() { return this->communicator; }
  DHT *getDHT() { return this->dht_object.get(); };

  DHT_ResultObject &getDHTResults() { return this->dht_results; }

  const auto &getKeyElements() const { return this->input_key_elements; }
  const auto &getKeySpecies() const { return this->key_species; }

  void setBaseTotals(double tot_h, double tot_o) {
    this->base_totals = {tot_h, tot_o};
  }

  std::uint32_t getInputCount() const {
    return this->input_key_elements.size();
  }

  std::uint32_t getOutputCount() const { return this->data_count; }

  inline void prepareKeys(const std::vector<std::vector<double>> &input_values,
                          double dt) {
    dht_results.keys.resize(input_values.size());
    for (std::size_t i = 0; i < input_values.size(); i++) {
      if (this->hooks.dht_fuzz.isValid()) {
        dht_results.keys[i] = fuzzForDHT_R(input_values[i], dt);
      } else {
        dht_results.keys[i] = fuzzForDHT(input_values[i], dt);
      }
    }
  }

private:
  uint32_t key_count;
  uint32_t data_count;

  std::unique_ptr<DHT> dht_object;

  MPI_Comm communicator;

  LookupKey fuzzForDHT(const std::vector<double> &cell, double dt);
  LookupKey fuzzForDHT_R(const std::vector<double> &cell, double dt);

  std::vector<double>
  outputToInputAndRates(const std::vector<double> &old_results,
                        const std::vector<double> &new_results);

  std::vector<double>
  inputAndRatesToOutput(const std::vector<double> &dht_data,
                        const std::vector<double> &input_values);

  std::vector<double> outputToRates(const std::vector<double> &old_results,
                                    const std::vector<double> &new_results);

  std::vector<double> ratesToOutput(const std::vector<double> &dht_data,
                                    const std::vector<double> &input_values);

  uint32_t dht_hits = 0;
  uint32_t dht_evictions = 0;

  NamedVector<std::uint32_t> key_species;

  std::vector<std::uint32_t> dht_signif_vector;
  std::vector<std::uint32_t> dht_prop_type_vector;
  std::vector<std::int32_t> input_key_elements;

  const std::vector<std::string> &output_names;

  const InitialList::ChemistryHookFunctions &hooks;
  const bool with_interp;

  DHT_ResultObject dht_results;

  std::array<double, 2> base_totals{0};
  bool has_het_ids{false};
};
} // namespace poet

#endif // DHT_WRAPPER_H
