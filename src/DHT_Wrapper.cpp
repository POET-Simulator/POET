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

void poet::DHT_ResultObject::ResultsToMapping(
    std::vector<int32_t> &curr_mapping) {
  uint32_t iMappingIndex = 0;
  for (uint32_t i = 0; i < this->length; i++) {
    if (curr_mapping[i] == -1) {
      continue;
    }
    curr_mapping[i] = (this->needPhreeqc[i] ? iMappingIndex++ : -1);
  }
}

void poet::DHT_ResultObject::ResultsToWP(std::vector<double> &curr_wp) {
  for (uint32_t i = 0; i < this->length; i++) {
    if (!this->needPhreeqc[i]) {
      const uint32_t length = this->results[i].end() - this->results[i].begin();
      std::copy(this->results[i].begin(), this->results[i].end(),
                curr_wp.begin() + (length * i));
    }
  }
}

DHT_Wrapper::DHT_Wrapper(const poet::SimParams &params, MPI_Comm dht_comm,
                         uint32_t dht_size, uint32_t key_count,
                         uint32_t data_count)
    : key_count(key_count), data_count(data_count) {
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
auto DHT_Wrapper::checkDHT(int length, double dt,
                           const std::vector<double> &work_package)
    -> poet::DHT_ResultObject {

  DHT_ResultObject check_data;

  check_data.length = length;
  check_data.keys.resize(length);
  check_data.results.resize(length);
  check_data.needPhreeqc.resize(length);

  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // point to current grid cell
    void *key = (void *)&(work_package[i * this->key_count]);
    auto &data = check_data.results[i];
    auto &key_vector = check_data.keys[i];

    data.resize(this->data_count);
    key_vector = fuzzForDHT(this->key_count, key, dt);

    // overwrite input with data from DHT, IF value is found in DHT
    int res = DHT_read(this->dht_object, key_vector.data(), data.data());

    switch (res) {
    case DHT_SUCCESS:
      check_data.needPhreeqc[i] = false;
      this->dht_hits++;
      break;
    case DHT_READ_MISS:
      check_data.needPhreeqc[i] = true;
      this->dht_miss++;
      break;
    }
  }

  return check_data;
}

void DHT_Wrapper::fillDHT(int length, const DHT_ResultObject &dht_check_data,
                          const std::vector<double> &results) {
  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // If true grid cell was simulated, needs to be inserted into dht
    if (dht_check_data.needPhreeqc[i]) {
      const auto &key = dht_check_data.keys[i];
      void *data = (void *)&(results[i * this->data_count]);
      // fuzz data (round, logarithm etc.)

      // insert simulated data with fuzzed key into DHT
      int res = DHT_write(this->dht_object, (void *)key.data(), data);

      // if data was successfully written ...
      if ((res != DHT_SUCCESS) && (res == DHT_WRITE_SUCCESS_WITH_EVICTION)) {
        dht_evictions++;
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
  constexpr double zero_val = 10E-14;

  std::vector<DHT_Keyelement> vecFuzz(var_count);
  std::memset(&vecFuzz[0], 0, sizeof(DHT_Keyelement) * var_count);

  unsigned int i = 0;
  // introduce fuzzing to allow more hits in DHT
  // loop over every variable of grid cell
  for (i = 0; i < (unsigned int)var_count; i++) {
    double &curr_key = ((double *)key)[i];
    if (curr_key != 0) {
      if (curr_key < zero_val && this->dht_prop_type_vector[i] == "act") {
        continue;
      }
      if (this->dht_prop_type_vector[i] == "ignore") {
        continue;
      }
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
