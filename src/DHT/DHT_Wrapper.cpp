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

#include "DHT_Wrapper.h"

#include <math.h>
#include <openssl/md5.h>

#include <iostream>

using namespace poet;
using namespace std;

uint64_t poet::get_md5(int key_size, void *key) {
  MD5_CTX ctx;
  unsigned char sum[MD5_DIGEST_LENGTH];
  uint64_t retval, *v1, *v2;

  // calculate md5 using MD5 functions
  MD5_Init(&ctx);
  MD5_Update(&ctx, key, key_size);
  MD5_Final(sum, &ctx);

  // divide hash in 2 64 bit parts and XOR them
  v1 = (uint64_t *)&sum[0];
  v2 = (uint64_t *)&sum[8];
  retval = *v1 ^ *v2;

  return retval;
}

DHT_Wrapper::DHT_Wrapper(SimParams &params, MPI_Comm dht_comm,
                         int buckets_per_process, int data_size, int key_size) {
  // initialize DHT object
  dht_object =
      DHT_create(dht_comm, buckets_per_process, data_size, key_size, &get_md5);
  // allocate memory for fuzzing buffer
  fuzzing_buffer = (double *)malloc(key_size);

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
}

void DHT_Wrapper::checkDHT(int length, std::vector<bool> &out_result_index,
                           double *work_package, double dt) {
  void *key;
  int res;
  // var count -> count of variables per grid cell
  int var_count = dht_prop_type_vector.size();
  // loop over every grid cell contained in work package
  for (int i = 0; i < length; i++) {
    // point to current grid cell
    key = (void *)&(work_package[i * var_count]);

    // fuzz data (round, logarithm etc.)
    fuzzForDHT(var_count, key, dt);

    // overwrite input with data from DHT, IF value is found in DHT
    res = DHT_read(dht_object, fuzzing_buffer, key);

    // if DHT_SUCCESS value was found ...
    if (res == DHT_SUCCESS) {
      // ... and grid cell will be marked as 'not to be simulating'
      out_result_index[i] = false;
      dht_hits++;

    }
    // ... otherwise ...
    else if (res == DHT_READ_MISS) {
      // grid cell needs to be simulated by PHREEQC
      out_result_index[i] = true;
      dht_miss++;
    } else {
      // MPI ERROR ... WHAT TO DO NOW?
      // RUNNING CIRCLES WHILE SCREAMING
    }
  }
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
    key = (void *)&(work_package[i * var_count]);
    data = (void *)&(results[i * var_count]);

    // If true grid cell was simulated, needs to be inserted into dht
    if (result_index[i]) {
      // fuzz data (round, logarithm etc.)
      fuzzForDHT(var_count, key, dt);

      // insert simulated data with fuzzed key into DHT
      res = DHT_write(dht_object, fuzzing_buffer, data);

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
  if (res != DHT_SUCCESS) return res;

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

void DHT_Wrapper::fuzzForDHT(int var_count, void *key, double dt) {
  unsigned int i = 0;
  // introduce fuzzing to allow more hits in DHT
  // loop over every variable of grid cell
  for (i = 0; i < (unsigned int)var_count; i++) {
    // check if variable is defined as 'act'
    if (dht_prop_type_vector[i] == "act") {
      // if log is enabled (default)
      if (dht_log) {
        // if variable is smaller than 0, which would be a strange result,
        // warn the user and set fuzzing_buffer to 0 at this index
        if (((double *)key)[i] < 0) {
          cerr << "dht_wrapper.cpp::fuzz_for_dht(): Warning! Negative value in "
                  "key!"
               << endl;
          fuzzing_buffer[i] = 0;
        }
        // if variable is 0 set fuzzing buffer to 0
        else if (((double *)key)[i] == 0)
          fuzzing_buffer[i] = 0;
        // otherwise ...
        else
          // round current variable value by applying log with base 10, negate
          // (since the actual values will be between 0 and 1) and cut result
          // after significant digit
          fuzzing_buffer[i] =
              ROUND(-(std::log10(((double *)key)[i])), dht_signif_vector[i]);
      }
      // if log is disabled
      else {
        // just round by cutting after signifanct digit
        fuzzing_buffer[i] = ROUND((((double *)key)[i]), dht_signif_vector[i]);
      }
    }
    // if variable is defined as 'logact' (log was already applied e.g. pH)
    else if (dht_prop_type_vector[i] == "logact") {
      // just round by cutting after signifanct digit
      fuzzing_buffer[i] = ROUND((((double *)key)[i]), dht_signif_vector[i]);
    }
    // if defined ass 'ignore' ...
    else if (dht_prop_type_vector[i] == "ignore") {
      // ... just set fuzzing buffer to 0
      fuzzing_buffer[i] = 0;
    }
    // and finally, if type is not defined, print error message
    else {
      cerr << "dht_wrapper.cpp::fuzz_for_dht(): Warning! Probably wrong "
              "prop_type!"
           << endl;
    }
  }
  // if timestep differs over iterations set current current time step at the
  // end of fuzzing buffer
  if (dt_differ) fuzzing_buffer[var_count] = dt;
}
