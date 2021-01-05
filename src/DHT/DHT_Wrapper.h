#ifndef DHT_WRAPPER_H
#define DHT_WRAPPER_H

#include <mpi.h>

#include <string>
#include <vector>

#include "../util/SimParams.h"
#include "DHT.h"

#define ROUND(value, signif) \
  (((int)(pow(10.0, (double)signif) * value)) * pow(10.0, (double)-signif))

uint64_t get_md5(int key_size, void *key);

namespace poet {
class DHT_Wrapper {
 public:
  DHT_Wrapper(t_simparams *params, MPI_Comm dht_comm, int buckets_per_process,
              int data_size, int key_size);
  ~DHT_Wrapper();

  void checkDHT(int length, std::vector<bool> &out_result_index,
                double *work_package, double dt);
  void fillDHT(int length, std::vector<bool> &result_index,
               double *work_package, double *results, double dt);

  int tableToFile(const char *filename);
  int fileToTable(const char *filename);

  void printStatistics();

  uint64_t getHits();
  uint64_t getMisses();
  uint64_t getEvictions();

 private:
  void fuzzForDHT(int var_count, void *key, double dt);

  DHT *dht_object;

  uint64_t dht_hits = 0;
  uint64_t dht_miss = 0;
  uint64_t dht_evictions = 0;

  double *fuzzing_buffer;

  bool dt_differ;
  bool dht_log;
  std::vector<int> dht_signif_vector;
  std::vector<std::string> dht_prop_type_vector;
};
}  // namespace poet

#endif  // DHT_WRAPPER_H
