#ifndef PARSER_H
#define PARSER_H

#include <string>

#include "RRuntime.h"
#include "argh.h"

/** Return value if no error occured */
#define PARSER_OK 0
/** Return value if error occured during parsing of program arguments */
#define PARSER_ERROR -1
/** Return value if user asked for help message with program parameter */
#define PARSER_HELP -2

/** Standard DHT Size (Defaults to 1 GiB) */
#define DHT_SIZE_PER_PROCESS 1073741824
/** Standard work package size */
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
/**
 * @brief Reads information from program arguments and R runtime
 *
 * Providing functions to initialize parameters of the simulation using command
 * line parameters and parameters from the R runtime. This class will also parse
 * arguments from the commandline and decides if argument is known or unknown.
 *
 * Stores and distribute current simulation parameters at any time.
 *
 */
class SimParams {
 public:
  /**
   * @brief Construct a new SimParams object
   *
   * With all given parameters a new instance of this class will be created.
   *
   * @param world_rank Rank of process inside MPI_COMM_WORLD
   * @param world_size Size of communicator MPI_COMM_WORLD
   */
  SimParams(int world_rank, int world_size);

  /**
   * @brief Parse program arguments
   *
   * This is done by the argh.h library.
   *
   * First, the function will check if there is a flag 'help' or 'h'. If this is
   * the case a help message is printed and the function will return with
   * PARSER_HELP.
   *
   * Second, if there are not 2 positional arguments an error will be printed to
   * stderr and the function returns with PARSER_ERROR.
   *
   * Then all given program parameters and flags will be read and checked, if
   * there are known by validateOptions. A list of all unknown options might be
   * returned, printed out and the function will return with PARSER_ERROR.
   * Oterhwise the function continuos.
   *
   * Now all program arguments will be stored inside t_simparams struct, printed
   * out and the function returns with PARSER_OK.
   *
   * Also, all parsed agruments are distributed to the R runtime.
   *
   * @param argv Argument value of the program
   * @param R Instantiated R runtime
   * @return int Returns with 0 if no error occured, otherwise value less than 0
   * is returned.
   */
  int parseFromCmdl(char *argv[], RRuntime &R);

  void initVectorParams(RRuntime &R, int col_count);

  void setDtDiffer(bool dt_differ);

  t_simparams getNumParams();

  std::vector<int> getDHTSignifVector();
  std::vector<std::string> getDHTPropTypeVector();
  std::string getDHTFile();

  std::string getFilesim();
  std::string getOutDir();

 private:
  /**
   * @brief Validate program parameters and flags
   *
   * Therefore this function iterates over the list of flags and parameters and
   * compare them to the class member flagList and paramList. If a program
   * argument is not included it is put to a list. This list will be returned.
   *
   * @return std::list<std::string> List with all unknown parameters. Might be
   * empty.
   */
  std::list<std::string> validateOptions(argh::parser cmdl);

  /**
   * @brief Contains all valid program flags.
   *
   */
  std::set<std::string> flaglist{"ignore-result", "dht", "dht-nolog"};

  /**
   * @brief Contains all valid program parameters.
   *
   */
  std::set<std::string> paramlist{"work-package-size", "dht-signif",
                                  "dht-strategy",      "dht-size",
                                  "dht-snaps",         "dht-file"};

  /**
   * @brief Struct containing all simulation parameters
   *
   * Contains only those values which are standard arithmetic C types.
   *
   */
  t_simparams simparams;

  std::vector<int> dht_signif_vector;
  std::vector<std::string> dht_prop_type_vector;
  std::string dht_file;

  std::string filesim;
  std::string out_dir;
};
}  // namespace poet
#endif  // PARSER_H