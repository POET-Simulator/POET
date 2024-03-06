/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
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

#ifndef PARSER_H
#define PARSER_H

#include "../DataStructures/DataStructures.hpp"
#include "Macros.hpp"
#include "RInsidePOET.hpp"

#include "argh.hpp" // Argument handler https://github.com/adishavit/argh

#include <cstdint>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include <RInside.h>
#include <Rcpp.h>
// BSD-licenced

/** Standard DHT Size. Defaults to 1 GB (1000 MB) */
constexpr uint32_t DHT_SIZE_PER_PROCESS_MB = 1.5E3;
/** Standard work package size */
#define WORK_PACKAGE_SIZE_DEFAULT 32

namespace poet {

enum { PARSER_OK, PARSER_ERROR, PARSER_HELP };

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
  /** indicating whether the progress bar during chemistry simulation should be
   * printed or not */
  bool print_progressbar;

  bool interp_enabled;
} t_simparams;

using GridParams = struct s_GridParams {
  std::array<uint32_t, 2> n_cells;
  std::array<double, 2> s_cells;
  std::uint8_t dim;
  std::uint32_t total_n;

  std::string type;

  Rcpp::DataFrame init_df;
  std::string input_script;
  std::string database_path;

  std::vector<std::string> props;

  s_GridParams(RInside &R);
};

using DiffusionParams = struct s_DiffusionParams {
  Rcpp::DataFrame initial_t;

  Rcpp::NumericVector alpha;
  Rcpp::List vecinj_inner;

  Rcpp::DataFrame vecinj;
  Rcpp::DataFrame vecinj_index;

  s_DiffusionParams(RInside &R);
};

struct ChemistryParams {
  std::string database_path;
  std::string input_script;

  bool use_dht;
  std::uint64_t dht_size;
  int dht_snaps;
  std::string dht_file;
  std::string dht_outdir;
  NamedVector<std::uint32_t> dht_signifs;

  bool use_interp;
  std::uint64_t pht_size;
  std::uint32_t pht_max_entries;
  std::uint32_t interp_min_entries;
  NamedVector<std::uint32_t> pht_signifs;

  struct Chem_Hook_Functions {
    RHookFunction<bool> dht_fill;
    RHookFunction<std::vector<double>> dht_fuzz;
    RHookFunction<std::vector<std::size_t>> interp_pre;
    RHookFunction<bool> interp_post;
  } hooks;

  void initFromR(RInsidePOET &R);
};

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
  int parseFromCmdl(char *argv[], RInsidePOET &R);

  /**
   * @brief Init std::vector values
   *
   * This will initialize dht_signif_vector and dht_prop_type_vector internally
   * depending on whether vectors are defined by R-Simulation file or not.
   * 'init_chemistry' must be run beforehand.
   *
   * @param R R runtime
   */
  void initVectorParams(RInside &R);

  /**
   * @brief Get the numerical params struct
   *
   * Returns a struct which contains all numerical or boolean simulation
   * parameters.
   *
   * @return t_simparams Parameter struct
   */
  auto getNumParams() const { return this->simparams; };

  /**
   * @brief Get the DHT_Signif_Vector
   *
   * Returns a vector indicating which significant values are used for each
   * variable of a grid cell.
   *
   * @return std::vector<int> Vector of integers containing information about
   * significant digits
   */
  auto getDHTSignifVector() const { return this->dht_signif_vector; };

  auto getPHTSignifVector() const { return this->pht_signif_vector; };

  auto getPHTBucketSize() const { return this->pht_bucket_size; };
  auto getPHTMinEntriesNeeded() const { return this->pht_min_entries_needed; };

  /**
   * @brief Get the filesim name
   *
   * Returns a string containing the absolute path to a R file defining the
   * simulation.
   *
   * @return std::string Absolute path to R file
   */
  auto getFilesim() const { return this->filesim; };

  /**
   * @brief Get the output directory
   *
   * Returns the name of an absolute path where all output files should be
   * stored.
   *
   * @return std::string Absolute path to output path
   */
  auto getOutDir() const { return this->out_dir; };

  const auto &getChemParams() const { return chem_params; }

private:
  std::list<std::string> validateOptions(argh::parser cmdl);

  const std::set<std::string> flaglist{"ignore-result", "dht", "P", "progress",
                                       "interp"};
  const std::set<std::string> paramlist{
      "work-package-size", "dht-strategy",
      "dht-size",          "dht-snaps",
      "dht-file",          "interp-size",
      "interp-min",        "interp-bucket-entries"};

  t_simparams simparams;

  std::vector<uint32_t> dht_signif_vector;
  std::vector<uint32_t> pht_signif_vector;

  uint32_t pht_bucket_size;
  uint32_t pht_min_entries_needed;

  std::string filesim;
  std::string out_dir;

  ChemistryParams chem_params;
};
} // namespace poet
#endif // PARSER_H
