//  Time-stamp: "Last modified 2023-07-21 17:22:00 mluebke"

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

#include "poet/DHT_Wrapper.hpp"
#include "poet/DHT_Types.hpp"
#include "poet/HashFunctions.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <math.h>
#include <stdexcept>
#include <vector>

using namespace poet;
using namespace std;

inline DHT_SCNotation round_key_element(double value, std::uint32_t signif) {
  std::int8_t exp =
      static_cast<std::int8_t>(std::floor(std::log10(std::fabs(value))));
  std::int64_t significant =
      static_cast<std::int64_t>(value * std::pow(10, signif - exp - 1));

  DHT_SCNotation elem;
  elem.exp = exp;
  elem.significant = significant;

  return elem;
}

DHT_Wrapper::DHT_Wrapper(MPI_Comm dht_comm, uint32_t dht_size,
                         const std::vector<std::uint32_t> &key_indices,
                         uint32_t data_count)
    : key_count(key_indices.size()), data_count(data_count),
      input_key_elements(key_indices) {
  // initialize DHT object
  uint32_t key_size = (key_count + 1) * sizeof(DHT_Keyelement);
  uint32_t data_size = data_count * sizeof(double);
  uint32_t buckets_per_process = dht_size / (1 + data_size + key_size);
  dht_object = DHT_create(dht_comm, buckets_per_process, data_size, key_size,
                          &poet::Murmur2_64A);

  this->dht_signif_vector.resize(key_size, DHT_KEY_SIGNIF_DEFAULT);

  this->dht_prop_type_vector.resize(key_count, DHT_TYPE_DEFAULT);
  this->dht_prop_type_vector[0] = DHT_TYPE_TOTAL;
  this->dht_prop_type_vector[1] = DHT_TYPE_TOTAL;
  this->dht_prop_type_vector[2] = DHT_TYPE_CHARGE;
}

DHT_Wrapper::~DHT_Wrapper() {
  // free DHT
  DHT_free(dht_object, NULL, NULL);
}
auto DHT_Wrapper::checkDHT(int length, double dt,
                           const std::vector<double> &work_package,
                           std::vector<std::uint32_t> &curr_mapping)
    -> const poet::DHT_ResultObject & {

  dht_results.length = length;
  dht_results.keys.resize(length);
  dht_results.results.resize(length);
  dht_results.needPhreeqc.resize(length);

  std::vector<std::uint32_t> new_mapping;

  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // point to current grid cell
    void *key = (void *)&(work_package[i * this->data_count]);
    auto &data = dht_results.results[i];
    auto &key_vector = dht_results.keys[i];

    data.resize(this->data_count);
    key_vector = fuzzForDHT(this->key_count, key, dt);

    // overwrite input with data from DHT, IF value is found in DHT
    int res = DHT_read(this->dht_object, key_vector.data(), data.data());

    switch (res) {
    case DHT_SUCCESS:
      dht_results.needPhreeqc[i] = false;
      this->dht_hits++;
      break;
    case DHT_READ_MISS:
      dht_results.needPhreeqc[i] = true;
      new_mapping.push_back(curr_mapping[i]);
      this->dht_miss++;
      break;
    }
  }

  curr_mapping = std::move(new_mapping);

  return dht_results;
}

void DHT_Wrapper::fillDHT(int length, const std::vector<double> &work_package) {
  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // If true grid cell was simulated, needs to be inserted into dht
    if (dht_results.needPhreeqc[i]) {
      const auto &key = dht_results.keys[i];
      void *data = (void *)&(work_package[i * this->data_count]);
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

void DHT_Wrapper::resultsToWP(std::vector<double> &work_package) {
  for (int i = 0; i < dht_results.length; i++) {
    if (!dht_results.needPhreeqc[i]) {
      std::copy(dht_results.results[i].begin(), dht_results.results[i].end(),
                work_package.begin() + (data_count * i));
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

std::vector<DHT_Keyelement> DHT_Wrapper::fuzzForDHT(int var_count, void *key,
                                                    double dt) {
  constexpr double zero_val = 10E-14;

  std::vector<DHT_Keyelement> vecFuzz(var_count + 1);
  std::memset(&vecFuzz[0], 0, sizeof(DHT_Keyelement) * var_count);

  int totals_i = 0;
  // introduce fuzzing to allow more hits in DHT
  // loop over every variable of grid cell
  for (std::uint32_t i = 0; i < input_key_elements.size(); i++) {
    double curr_key = ((double *)key)[input_key_elements[i]];
    if (curr_key != 0) {
      if (curr_key < zero_val &&
          this->dht_prop_type_vector[i] == DHT_TYPE_DEFAULT) {
        continue;
      }
      if (this->dht_prop_type_vector[i] == DHT_TYPE_TOTAL) {
        curr_key -= base_totals[totals_i++];
      }
      vecFuzz[i].sc_notation =
          round_key_element(curr_key, dht_signif_vector[i]);
    }
  }
  // if timestep differs over iterations set current current time step at the
  // end of fuzzing buffer
  vecFuzz[var_count].fp_elemet = dt;

  return vecFuzz;
}

void poet::DHT_Wrapper::SetSignifVector(std::vector<uint32_t> signif_vec) {
  if (signif_vec.size() != this->key_count) {
    throw std::runtime_error(
        "Significant vector size mismatches count of key elements.");
  }

  this->dht_signif_vector = signif_vec;
}
