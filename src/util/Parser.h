#ifndef PARSER_H
#define PARSER_H

#include <RRuntime.h>

#include <string>

#include "SimParams.h"
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

namespace poet {
/**
 * @brief Reads information from program arguments and R runtime
 *
 * Providing functions to initialize parameters of the simulation using command
 * line parameters and parameters from the R runtime. This class will also parse
 * arguments from the commandline and decides if argument is known or unknown.
 *
 */
class Parser {
 public:
  /**
   * @brief Construct a new Parser object
   *
   * With all given parameters a new instance of this class will be created.
   *
   * @param argv Argument value of the program
   * @param world_rank Rank of process inside MPI_COMM_WORLD
   * @param world_size Size of communicator MPI_COMM_WORLD
   */
  Parser(char *argv[], int world_rank, int world_size);

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
   * @return int Returns with 0 if no error occured, otherwise value less than 0
   * is returned.
   */
  int parseCmdl();

  /**
   * @brief Distribute all known parameters to R runtime.
   *
   * All stored parameters are distributed to the R runtime.
   *
   * @todo This function might be redundant and can be put into parseCmdl.
   *
   * @param R Instance of RRuntime
   */
  void parseR(RRuntime &R);

  /**
   * @brief Get the Params object
   *
   * Return all parsed simulation parameters. Should be called after parseCmdl.
   *
   * @return t_simparams Struct of t_simparams
   */
  t_simparams getParams();

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
  std::list<std::string> validateOptions();

  /**
   * @brief Contains all valid program flags.
   *
   */
  const std::set<std::string> flaglist{"ignore-result", "dht", "dht-nolog"};

  /**
   * @brief Contains all valid program parameters.
   *
   */
  const std::set<std::string> paramlist{"work-package-size", "dht-signif",
                                        "dht-strategy",      "dht-size",
                                        "dht-snaps",         "dht-file"};

  /**
   * @brief Instance of argh class
   *
   * This class will be instantiate inside constructor of this class object.
   *
   */
  argh::parser cmdl;

  /**
   * @brief Struct containing all simulation parameters
   * 
   */
  t_simparams simparams;

  /**
   * @brief Rank of process inside MPI_COMM_WORLD
   * 
   */
  int world_rank;

  /**
   * @brief Size of MPI_COMM_WORLD
   * 
   */
  int world_size;

  int dht_significant_digits;
};
}  // namespace poet
#endif  // PARSER_H