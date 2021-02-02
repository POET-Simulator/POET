#ifndef PARSER_H
#define PARSER_H

#include <RRuntime.h>

#include <string>

#include "SimParams.h"
#include "argh.h"

#define PARSER_OK 0
#define PARSER_ERROR 1
#define PARSER_HELP 2

#define DHT_SIZE_PER_PROCESS 1073741824
#define WORK_PACKAGE_SIZE_DEFAULT 5

namespace poet {
class Parser {
 public:
  Parser(char *argv[], int world_rank, int world_size);
  int parseCmdl();
  void parseR(RRuntime &R);

  t_simparams getParams();

 private:
  std::list<std::string> checkOptions(argh::parser cmdl);
  std::set<std::string> flaglist{"ignore-result", "dht", "dht-nolog"};
  std::set<std::string> paramlist{"work-package-size", "dht-signif",
                                  "dht-strategy",      "dht-size",
                                  "dht-snaps",         "dht-file"};
  argh::parser cmdl;

  t_simparams simparams;

  int world_rank;
  int world_size;
  
};
}  // namespace poet
#endif  // PARSER_H