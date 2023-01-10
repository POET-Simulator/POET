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

#include "poet/HashFunctions.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <openssl/evp.h>
#include <poet/DHT_Wrapper.hpp>

#include <math.h>

#include <iostream>
#include <vector>

using namespace poet;
using namespace std;

inline DHT_Keyelement round_key_element(double value, std::uint32_t signif) {
  std::int8_t exp = (std::int8_t)std::floor(std::log10(std::fabs(value)));
  std::int64_t significant = value * std::pow(10, signif - exp);

  DHT_Keyelement elem;
  elem.exp = exp;
  elem.significant = significant;

  return elem;
}

DHT_Wrapper::DHT_Wrapper(const poet::SimParams &params, MPI_Comm dht_comm,
                         uint32_t dht_size, uint32_t key_count,
                         uint32_t data_count) {
  poet::initHashCtx(EVP_md5());
  // initialize DHT object
  uint32_t key_size = key_count * sizeof(DHT_Keyelement);
  uint32_t data_size = data_count * sizeof(double);
  uint32_t buckets_per_process = dht_size / (1 + data_size + key_size);
  dht_object = DHT_create(dht_comm, buckets_per_process, data_size, key_size,
                          &poet::hashDHT);

  // extract needed values from sim_param struct
  t_simparams tmp = params.getNumParams();

  this->dt_differ = tmp.dt_differ;
  this->dht_log = tmp.dht_log;

  this->dht_signif_vector = params.getDHTSignifVector();
  this->dht_prop_type_vector = params.getDHTPropTypeVector();
}

DHT_Wrapper::~DHT_Wrapper() {
  // free DHT
  DHT_free(dht_object, NULL, NULL);
  // free fuzzing buffer
  free(fuzzing_buffer);
  poet::freeHashCtx();
}

auto DHT_Wrapper::checkDHT(int length, std::vector<bool> &out_result_index,
                           double *work_package, double dt)
    -> std::vector<std::vector<double>> {
  void *key;
  int res;
  // var count -> count of variables per grid cell
  int var_count = dht_prop_type_vector.size() - 1;
  std::vector<std::vector<double>> data(length);
  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // point to current grid cell
    key = (void *)&(work_package[i * var_count]);
    data[i].resize(dht_prop_type_vector.size());

    // fuzz data (round, logarithm etc.)
    auto vecFuzz = fuzzForDHT(var_count, key, dt);

    // overwrite input with data from DHT, IF value is found in DHT
    res = DHT_read(dht_object, vecFuzz.data(), data[i].data());

    // if DHT_SUCCESS value was found ...
    if (res == DHT_SUCCESS) {
      // ... and grid cell will be marked as 'not to be simulating'
      out_result_index[i] = false;
      dht_hits++;
    }
    // ... otherwise ...
    else if (res == DHT_READ_MISS) {
      // grid cell needs to be simulated by PHREEQC
      dht_miss++;
    } else {
      // MPI ERROR ... WHAT TO DO NOW?
      // RUNNING CIRCLES WHILE SCREAMING
    }
  }

  return data;
}

void DHT_Wrapper::fillDHT(int length, std::vector<bool> &result_index,
                          double *work_package, double *results, double dt) {
  void *key;
  void *data;
  int res;
  // var count -> count of variables per grid cell
  int var_count = dht_prop_type_vector.size();
  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    key = (void *)&(work_package[i * (var_count - 1)]);
    data = (void *)&(results[i * var_count]);

    // If true grid cell was simulated, needs to be inserted into dht
    if (result_index[i]) {
      // fuzz data (round, logarithm etc.)
      auto vecFuzz = fuzzForDHT(var_count - 1, key, dt);

      // insert simulated data with fuzzed key into DHT
      res = DHT_write(dht_object, vecFuzz.data(), data);

      // if data was successfully written ...
      if (res != DHT_SUCCESS) {
        // ... also check if a previously written value was evicted
        if (res == DHT_WRITE_SUCCESS_WITH_EVICTION) {
          // and increment internal eviciton counter
          dht_evictions++;
        } else {
          // MPI ERROR ... WHAT TO DO NOW?
          // RUNNING CIRCLES WHILE SCREAMING
        }
      }
    }
  }
}

int DHT_Wrapper::tableToFile(const char *filename) {
  int res = DHT_to_file(dht_object, filename);
  return res;
}

int DHT_Wrapper::fileToTable(const char *filename) {
  int res = DHT_from_file(dht_object, filename);
  if (res != DHT_SUCCESS)
    return res;

#ifdef DHT_STATISTICS
  DHT_print_statistics(dht_object);
#endif

  return DHT_SUCCESS;
}

void DHT_Wrapper::printStatistics() {
  int res;

  res = DHT_print_statistics(dht_object);

  if (res != DHT_SUCCESS) {
    // MPI ERROR ... WHAT TO DO NOW?
    // RUNNING CIRCLES WHILE SCREAMING
  }
}

uint64_t DHT_Wrapper::getHits() { return this->dht_hits; }

uint64_t DHT_Wrapper::getMisses() { return this->dht_miss; }

uint64_t DHT_Wrapper::getEvictions() { return this->dht_evictions; }

std::vector<DHT_Keyelement> DHT_Wrapper::fuzzForDHT(int var_count, void *key,
                                                    double dt) {

  std::vector<DHT_Keyelement> vecFuzz(var_count);
  std::memset(&vecFuzz[0], 0, sizeof(DHT_Keyelement) * var_count);

  unsigned int i = 0;
  // introduce fuzzing to allow more hits in DHT
  // loop over every variable of grid cell
  for (i = 0; i < (unsigned int)var_count; i++) {
    double &curr_key = ((double *)key)[i];
    if (curr_key != 0) {
      vecFuzz[i] = round_key_element(curr_key, dht_signif_vector[i]);
    }
  }
  // if timestep differs over iterations set current current time step at the
  // end of fuzzing buffer
  if (dt_differ) {
    vecFuzz[var_count] = round_key_element(dt, 55);
  }

  return vecFuzz;
}
