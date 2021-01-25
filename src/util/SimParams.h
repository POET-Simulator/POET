#ifndef SIMPARAMS_H
#define SIMPARAMS_H

#include <string>
#include <vector>

typedef struct {
  int world_size;
  int world_rank;

  bool dht_enabled;
  bool dht_log;
  bool dt_differ;
  int dht_snaps;
  int dht_strategy;
  unsigned int dht_size_per_process;
  std::vector<int> dht_signif_vector;
  std::vector<std::string> dht_prop_type_vector;
  std::string dht_file;

  unsigned int wp_size;

  std::string filesim;
  std::string out_dir;

  bool store_result;

  // void* R;
  // void* grid;
} t_simparams;

#endif  // SIMPARAMS_H
