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

namespace poet {

/**
 * @brief Defining all simulation parameters
 *
 */
typedef struct {
  /** Count of processes in MPI_COMM_WORLD */
  int world_size;
  /** rank of proces in MPI_COMM_WORLD */
  int world_rank;
  /** indicates if DHT should be used */
  bool dht_enabled;
  /** apply logarithm to key before rounding */
  bool dht_log;
  /** indicates if timestep dt differs between iterations */
  bool dt_differ;
  /** Indicates, when a DHT snapshot should be written */
  int dht_snaps;
  /** <b>not implemented</b>: How a DHT is distributed over processes */
  int dht_strategy;
  /** Size of DHt per process in byter */
  unsigned int dht_size_per_process;
  /** Default significant digit for rounding */
  int dht_significant_digits;
  /** Default work package size */
  unsigned int wp_size;
  /** indicates if resulting grid should be stored after every iteration */
  bool store_result;
} t_simparams;

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

  /**
   * @brief Init std::vector values
   *
   * This will initialize dht_signif_vector and dht_prop_type_vector internally
   * depending on whether vectors are defined by R-Simulation file or not.
   * 'init_chemistry' must be run beforehand.
   *
   * @param R R runtime
   * @param col_count Count of variables per grid cell (typically the count of
   * columns of each grid cell)
   */
  void initVectorParams(RRuntime &R, int col_count);

  /**
   * @brief Set if dt differs
   *
   * Set a boolean variable if the timestep differs between iterations of
   * simulation.
   *
   * @param dt_differ Boolean value, if dt differs
   */
  void setDtDiffer(bool dt_differ);

  /**
   * @brief Get the numerical params struct
   *
   * Returns a struct which contains all numerical or boolean simulation
   * parameters.
   *
   * @return t_simparams Parameter struct
   */
  t_simparams getNumParams();

  /**
   * @brief Get the DHT_Signif_Vector
   *
   * Returns a vector indicating which significant values are used for each
   * variable of a grid cell.
   *
   * @return std::vector<int> Vector of integers containing information about
   * significant digits
   */
  std::vector<int> getDHTSignifVector();

  /**
   * @brief Get the DHT_Prop_Type_Vector
   *
   * Returns a vector indicating of which type a variable of a grid cell is.
   *
   * @return std::vector<std::string> Vector if strings defining a type of a
   * variable
   */
  std::vector<std::string> getDHTPropTypeVector();

  /**
   * @brief Return name of DHT snapshot.
   *
   * Name of the DHT file which is used to initialize the DHT with a previously
   * written snapshot.
   *
   * @return std::string Absolute paht to the DHT snapshot
   */
  std::string getDHTFile();

  /**
   * @brief Get the filesim name
   *
   * Returns a string containing the absolute path to a R file defining the
   * simulation.
   *
   * @return std::string Absolute path to R file
   */
  std::string getFilesim();

  /**
   * @brief Get the output directory
   *
   * Returns the name of an absolute path where all output files should be
   * stored.
   *
   * @return std::string Absolute path to output path
   */
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

  /**
   * @brief Defines significant digits for each variable of a grid cell
   *
   */
  std::vector<int> dht_signif_vector;

  /**
   * @brief Defines the type of a variable
   *
   */
  std::vector<std::string> dht_prop_type_vector;

  /**
   * @brief Absolute path to a DHT snapshot
   *
   */
  std::string dht_file;

  /**
   * @brief Absolute path to R file containing simulation definitions
   *
   */
  std::string filesim;

  /**
   * @brief Absolute path to output dir
   *
   */
  std::string out_dir;
};
}  // namespace poet
#endif  // PARSER_H