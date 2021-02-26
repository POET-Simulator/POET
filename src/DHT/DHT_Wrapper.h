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

#include <SimParams.h>

#include <string>
#include <vector>

extern "C" {
#include <DHT.h>
}

#include <mpi.h>

/**
 * @brief Cut double value after signif digit
 *
 * Macro to round a double value by cutting every digit after significant digit
 *
 */
#define ROUND(value, signif) \
  (((int)(pow(10.0, (double)signif) * value)) * pow(10.0, (double)-signif))

namespace poet {
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
   * @param params Simulation parameter object
   * @param dht_comm Communicator which addresses all participating DHT
   * processes
   * @param buckets_per_process Count of buckets to allocate for each process
   * @param data_size Size of data in bytes
   * @param key_size Size of key in bytes
   */
  DHT_Wrapper(SimParams &params, MPI_Comm dht_comm, int buckets_per_process,
              int data_size, int key_size);
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
  void checkDHT(int length, std::vector<bool> &out_result_index,
                double *work_package, double dt);

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
  void fillDHT(int length, std::vector<bool> &result_index,
               double *work_package, double *results, double dt);

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
  uint64_t getHits();

  /**
   * @brief Get the Misses object
   *
   * @return uint64_t Count of read misses
   */
  uint64_t getMisses();

  /**
   * @brief Get the Evictions object
   *
   * @return uint64_t Count of evictions
   */
  uint64_t getEvictions();

 private:
  /**
   * @brief Transform given workpackage into DHT key
   *
   * A given workpackage will be transformed into a DHT key by rounding each
   * value of a workpackage to a given significant digit. Three different types
   * of variables 'act', 'logact' and 'ignore' are used. Those types are given
   * via the dht_signif_vector.
   *
   * If a variable is defined as 'act', dht_log is true and non-negative, the
   * logarithm with base 10 will be applied. After that the value is negated. In
   * case the value is 0 the fuzzing_buffer is also set to 0 at this position.
   * If the value is negative a correspondending warning will be printed to
   * stderr and the fuzzing buffer will be set to 0 at this index.
   *
   * If a variable is defined as 'logact' the value will be cut after the
   * significant digit.
   *
   * If a variable ist defined as 'ignore' the fuzzing_buffer will be set to 0
   * at the index of the variable.
   *
   * If dt_differ is true the current time step of the simulation will be set at
   * the end of the fuzzing_buffer.
   *
   * @param var_count Count of variables for the current work package
   * @param key Pointer to work package handled as the key
   * @param dt Current time step of the simulation
   */
  void fuzzForDHT(int var_count, void *key, double dt);

  /**
   * @brief DHT handle
   *
   * Stores information about the DHT. Will be used as a handle for each DHT
   * library call.
   *
   */
  DHT *dht_object;

  /**
   * @brief Count of hits
   *
   * The counter will be incremented if a previously simulated workpackage can
   * be retrieved with a given key.
   *
   */
  uint64_t dht_hits = 0;

  /**
   * @brief Count of read misses
   *
   * The counter will be incremented if a given key doesn't retrieve a value
   * from the DHT.
   *
   */
  uint64_t dht_miss = 0;

  /**
   * @brief Count of evictions
   *
   * If a value in the DHT must be evicted because of lack of space/reaching the
   * last index etc., this counter will be incremented.
   *
   */
  uint64_t dht_evictions = 0;

  /**
   * @brief Rounded work package values
   *
   * Stores rounded work package values and serves as the DHT key pointer.
   *
   */
  double *fuzzing_buffer;

  /**
   * @brief Indicates change in time step during simulation
   *
   * If set to true, the time step of simulation will differ between iterations,
   * so the current time step must be stored inside the DHT key. Otherwise wrong
   * values would be obtained.
   *
   * If set to false the time step doesn't need to be stored in the DHT key.
   *
   */
  bool dt_differ;

  /**
   * @brief Logarithm before rounding
   *
   * Indicates if the logarithm with base 10 will be applied to a variable
   * before rounding.
   *
   * Defaults to true.
   *
   */
  bool dht_log;

  /**
   * @brief Significant digits for each variable
   *
   * Stores the rounding/significant digits for each variable of the work
   * package.
   *
   */
  std::vector<int> dht_signif_vector;

  /**
   * @brief Type of each variable
   *
   * Defines the type of each variable of the work package.
   *
   */
  std::vector<std::string> dht_prop_type_vector;
};
}  // namespace poet

#endif  // DHT_WRAPPER_H
