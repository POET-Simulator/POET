//  Time-stamp: "Last modified 2023-08-16 16:44:17 mluebke"

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
#include "poet/HashFunctions.hpp"
#include "poet/Interpolation.hpp"
#include "poet/LookupKey.hpp"
#include "poet/Rounding.hpp"
#include "poet/enums.hpp"

#include "poet/RInsidePOET.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

namespace poet {

DHT_Wrapper::DHT_Wrapper(MPI_Comm dht_comm, std::uint64_t dht_size,
                         const NamedVector<std::uint32_t> &key_species,
                         const std::vector<std::int32_t> &key_indices,
                         const std::vector<std::string> &_output_names,
                         const ChemistryParams::Chem_Hook_Functions &_hooks,
                         uint32_t data_count)
    : key_count(key_indices.size()), data_count(data_count),
      input_key_elements(key_indices), communicator(dht_comm),
      key_species(key_species), output_names(_output_names), hooks(_hooks) {
  // initialize DHT object
  // key size = count of key elements + timestep
  uint32_t key_size = (key_count + 1) * sizeof(Lookup_Keyelement);
  uint32_t data_size =
      (data_count + input_key_elements.size()) * sizeof(double);
  uint32_t buckets_per_process =
      static_cast<std::uint32_t>(dht_size / (data_size + key_size));
  dht_object = DHT_create(dht_comm, buckets_per_process, data_size, key_size,
                          &poet::Murmur2_64A);

  dht_signif_vector = key_species.getValues();

  // this->dht_signif_vector.resize(key_size, DHT_KEY_SIGNIF_DEFAULT);

  this->dht_prop_type_vector.resize(key_count, DHT_TYPE_DEFAULT);

  const auto key_names = key_species.getNames();

  auto tot_h = std::find(key_names.begin(), key_names.end(), "H");
  if (tot_h != key_names.end()) {
    this->dht_prop_type_vector[tot_h - key_names.begin()] = DHT_TYPE_TOTAL;
  }
  auto tot_o = std::find(key_names.begin(), key_names.end(), "O");
  if (tot_o != key_names.end()) {
    this->dht_prop_type_vector[tot_o - key_names.begin()] = DHT_TYPE_TOTAL;
  }
  auto charge = std::find(key_names.begin(), key_names.end(), "Charge");
  if (charge != key_names.end()) {
    this->dht_prop_type_vector[charge - key_names.begin()] = DHT_TYPE_CHARGE;
  }
}

DHT_Wrapper::~DHT_Wrapper() {
  // free DHT
  DHT_free(dht_object, NULL, NULL);
}
auto DHT_Wrapper::checkDHT(WorkPackage &work_package)
    -> const DHT_ResultObject & {

  const auto length = work_package.size;

  std::vector<double> bucket_writer(this->data_count +
                                    input_key_elements.size());

  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // point to current grid cell
    auto &key_vector = dht_results.keys[i];

    // overwrite input with data from DHT, IF value is found in DHT
    int res =
        DHT_read(this->dht_object, key_vector.data(), bucket_writer.data());

    switch (res) {
    case DHT_SUCCESS:
      work_package.output[i] = inputAndRatesToOutput(bucket_writer);
      work_package.mapping[i] = CHEM_DHT;
      this->dht_hits++;
      break;
    case DHT_READ_MISS:
      break;
    }
  }

  return dht_results;
}

void DHT_Wrapper::fillDHT(const WorkPackage &work_package) {

  const auto length = work_package.size;

  // loop over every grid cell contained in work package
  dht_results.locations.resize(length);
  dht_results.filledDHT = std::vector<bool>(length, false);
  for (int i = 0; i < length; i++) {
    // If true grid cell was simulated, needs to be inserted into dht
    if (work_package.mapping[i] == CHEM_PQC) {

      // check if calcite or dolomite is absent and present, resp.n and vice
      // versa in input/output. If this is the case -> Do not write to DHT!
      // HACK: hardcoded, should be fixed!
      if (hooks.dht_fill.isValid()) {
        NamedVector<double> old_values(output_names, work_package.input[i]);
        NamedVector<double> new_values(output_names, work_package.output[i]);

        if (hooks.dht_fill(old_values, new_values)) {
          continue;
        }
      }

      uint32_t proc, index;
      auto &key = dht_results.keys[i];
      const auto data =
          outputToInputAndRates(work_package.input[i], work_package.output[i]);
      // void *data = (void *)&(work_package[i * this->data_count]);
      // fuzz data (round, logarithm etc.)

      // insert simulated data with fuzzed key into DHT
      int res = DHT_write(this->dht_object, key.data(),
                          const_cast<double *>(data.data()), &proc, &index);

      dht_results.locations[i] = {proc, index};

      // if data was successfully written ...
      if ((res != DHT_SUCCESS) && (res == DHT_WRITE_SUCCESS_WITH_EVICTION)) {
        dht_evictions++;
      }

      dht_results.filledDHT[i] = true;
    }
  }
}

std::vector<double>
DHT_Wrapper::outputToInputAndRates(const std::vector<double> &old_results,
                                   const std::vector<double> &new_results) {
  const int prefix_size = this->input_key_elements.size();

  std::vector<double> output(prefix_size + this->data_count);
  std::copy(new_results.begin(), new_results.end(),
            output.begin() + prefix_size);

  for (int i = 0; i < prefix_size; i++) {
    const int data_elem_i = input_key_elements[i];
    output[i] = old_results[data_elem_i];
    output[prefix_size + data_elem_i] -= old_results[data_elem_i];
  }

  return output;
}

std::vector<double>
DHT_Wrapper::inputAndRatesToOutput(const std::vector<double> &dht_data) {
  const int prefix_size = this->input_key_elements.size();

  std::vector<double> output{dht_data.begin() + prefix_size, dht_data.end()};

  for (int i = 0; i < prefix_size; i++) {
    const int data_elem_i = input_key_elements[i];
    output[data_elem_i] += dht_data[i];
  }

  return output;
}

// void DHT_Wrapper::resultsToWP(std::vector<double> &work_package) {
//   for (int i = 0; i < dht_results.length; i++) {
//     if (!dht_results.needPhreeqc[i]) {
//       std::copy(dht_results.results[i].begin(), dht_results.results[i].end(),
//                 work_package.begin() + (data_count * i));
//     }
//   }
// }

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

LookupKey DHT_Wrapper::fuzzForDHT_R(const std::vector<double> &cell,
                                    double dt) {
  const auto c_zero_val = std::pow(10, AQUEOUS_EXP);

  NamedVector<double> input_nv(this->output_names, cell);

  const std::vector<double> eval_vec = hooks.dht_fuzz(input_nv);
  assert(eval_vec.size() == this->key_count);
  LookupKey vecFuzz(this->key_count + 1, {.0});

  DHT_Rounder rounder;

  int totals_i = 0;
  // introduce fuzzing to allow more hits in DHT
  // loop over every variable of grid cell
  for (std::uint32_t i = 0; i < eval_vec.size(); i++) {
    double curr_key = eval_vec[i];
    if (curr_key != 0) {
      if (this->dht_prop_type_vector[i] == DHT_TYPE_TOTAL) {
        curr_key -= base_totals[totals_i++];
      }
      vecFuzz[i] =
          rounder.round(curr_key, dht_signif_vector[i],
                        this->dht_prop_type_vector[i] == DHT_TYPE_TOTAL);
    }
  }
  // add timestep to the end of the key as double value
  vecFuzz[this->key_count].fp_element = dt;

  return vecFuzz;
}

LookupKey DHT_Wrapper::fuzzForDHT(const std::vector<double> &cell, double dt) {
  const auto c_zero_val = std::pow(10, AQUEOUS_EXP);

  LookupKey vecFuzz(this->key_count + 1, {.0});
  DHT_Rounder rounder;

  int totals_i = 0;
  // introduce fuzzing to allow more hits in DHT
  // loop over every variable of grid cell
  for (std::uint32_t i = 0; i < input_key_elements.size(); i++) {
    if (input_key_elements[i] == DHT_KEY_INPUT_CUSTOM) {
      continue;
    }
    double curr_key = cell[input_key_elements[i]];
    if (curr_key != 0) {
      if (curr_key < c_zero_val &&
          this->dht_prop_type_vector[i] == DHT_TYPE_DEFAULT) {
        continue;
      }
      if (this->dht_prop_type_vector[i] == DHT_TYPE_TOTAL) {
        curr_key -= base_totals[totals_i++];
      }
      vecFuzz[i] =
          rounder.round(curr_key, dht_signif_vector[i],
                        this->dht_prop_type_vector[i] == DHT_TYPE_TOTAL);
    }
  }
  // add timestep to the end of the key as double value
  vecFuzz[this->key_count].fp_element = dt;

  return vecFuzz;
}

void poet::DHT_Wrapper::SetSignifVector(std::vector<uint32_t> signif_vec) {
  if (signif_vec.size() != this->key_count) {
    throw std::runtime_error(
        "Significant vector size mismatches count of key elements.");
  }

  this->dht_signif_vector = signif_vec;
}
} // namespace poet
