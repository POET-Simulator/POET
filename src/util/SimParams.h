#ifndef PARSER_H
#define PARSER_H

#include "RRuntime.h"

#include <string>

#include "argh.h"

#define PARSER_OK 0
#define PARSER_ERROR 1
#define PARSER_HELP 2

#define DHT_SIZE_PER_PROCESS 1073741824
#define WORK_PACKAGE_SIZE_DEFAULT 5

typedef struct {
  int world_size;
  int world_rank;

  bool dht_enabled;
  bool dht_log;
  bool dt_differ;
  int dht_snaps;
  int dht_strategy;
  unsigned int dht_size_per_process;
  int dht_significant_digits;

  unsigned int wp_size;

  bool store_result;
} t_simparams;

namespace poet {
class SimParams {
 public:
  SimParams(int world_rank, int world_size);
  int parseFromCmdl(char *argv[], RRuntime &R);
  // void parseR(RRuntime &R);

  void initVectorParams(RRuntime &R, int col_count);

  void setDtDiffer(bool dt_differ);

  t_simparams getNumParams();
  
  std::vector<int> getDHTSignifVector();
  std::vector<std::string> getDHTPropTypeVector();
  std::string getDHTFile();

  std::string getFilesim();
  std::string getOutDir();

 private:
  std::list<std::string> checkOptions(argh::parser cmdl);
  std::set<std::string> flaglist{"ignore-result", "dht", "dht-nolog"};
  std::set<std::string> paramlist{"work-package-size", "dht-signif",
                                  "dht-strategy",      "dht-size",
                                  "dht-snaps",         "dht-file"};

  t_simparams simparams;

  int world_rank;
  int world_size;

  std::vector<int> dht_signif_vector;
  std::vector<std::string> dht_prop_type_vector;
  std::string dht_file;

  std::string filesim;
  std::string out_dir;
};
}  // namespace poet
#endif  // PARSER_H