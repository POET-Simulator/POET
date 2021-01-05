#include "DHT_Wrapper.h"
#include "DHT.h"
#include <openssl/md5.h>
#include <iostream>

using namespace poet;
using namespace std;
uint64_t get_md5(int key_size, void *key) {
  MD5_CTX ctx;
  unsigned char sum[MD5_DIGEST_LENGTH];
  uint64_t retval, *v1, *v2;

  MD5_Init(&ctx);
  MD5_Update(&ctx, key, key_size);
  MD5_Final(sum, &ctx);

  v1 = (uint64_t *)&sum[0];
  v2 = (uint64_t *)&sum[8];
  retval = *v1 ^ *v2;

  return retval;
}

DHT_Wrapper::DHT_Wrapper(t_simparams *params, MPI_Comm dht_comm,
                         int buckets_per_process, int data_size, int key_size) {
  dht_object =
      DHT_create(dht_comm, buckets_per_process, data_size, key_size, &get_md5);
  fuzzing_buffer = (double *)malloc(key_size);

  this->dt_differ = params->dt_differ;
  this->dht_log = params->dht_log;
  this->dht_signif_vector = params->dht_signif_vector;
  this->dht_prop_type_vector = params->dht_prop_type_vector;
}

DHT_Wrapper::~DHT_Wrapper() {
  DHT_free(dht_object, NULL, NULL);
  free(fuzzing_buffer);
}

void DHT_Wrapper::checkDHT(int length, std::vector<bool> &out_result_index,
                           double *work_package, double dt) {
  void *key;
  int res;
  int var_count = dht_prop_type_vector.size();
  for (int i = 0; i < length; i++) {
    key = (void *)&(work_package[i * var_count]);

    // fuzz data (round, logarithm etc.)
    fuzzForDHT(var_count, key, dt);

    // overwrite input with data from DHT, IF value is found in DHT
    res = DHT_read(dht_object, fuzzing_buffer, key);

    if (res == DHT_SUCCESS) {
      // flag that this line is replaced by DHT-value, do not simulate!!
      out_result_index[i] = false;
      dht_hits++;
    } else if (res == DHT_READ_ERROR) {
      // this line is untouched, simulation is needed
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
  int var_count = dht_prop_type_vector.size();

  for (int i = 0; i < length; i++) {
    key = (void *)&(work_package[i * var_count]);
    data = (void *)&(results[i * var_count]);

    if (result_index[i]) {
      // If true -> was simulated, needs to be inserted into dht

      // fuzz data (round, logarithm etc.)
      fuzzForDHT(var_count, key, dt);

      res = DHT_write(dht_object, fuzzing_buffer, data);

      if (res != DHT_SUCCESS) {
        if (res == DHT_WRITE_SUCCESS_WITH_COLLISION) {
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

  DHT_print_statistics(dht_object);

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
  for (i = 0; i < (unsigned int)var_count; i++) {
    if (dht_prop_type_vector[i] == "act") {
      // with log10
      if (dht_log) {
        if (((double *)key)[i] < 0)
          cerr << "dht_wrapper.cpp::fuzz_for_dht(): Warning! Negative value in "
                  "key!"
               << endl;
        else if (((double *)key)[i] == 0)
          fuzzing_buffer[i] = 0;
        else
          fuzzing_buffer[i] =
              ROUND(-(std::log10(((double *)key)[i])), dht_signif_vector[i]);
      } else {
        // without log10
        fuzzing_buffer[i] = ROUND((((double *)key)[i]), dht_signif_vector[i]);
      }
    } else if (dht_prop_type_vector[i] == "logact") {
      fuzzing_buffer[i] = ROUND((((double *)key)[i]), dht_signif_vector[i]);
    } else if (dht_prop_type_vector[i] == "ignore") {
      fuzzing_buffer[i] = 0;
    } else {
      cerr << "dht_wrapper.cpp::fuzz_for_dht(): Warning! Probably wrong "
              "prop_type!"
           << endl;
    }
  }
  if (dt_differ)
    fuzzing_buffer[var_count] = dt;
}
