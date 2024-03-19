#pragma once

#include <RInside.h>
#include <Rcpp.h>
#include <Rcpp/vector/instantiation.h>
#include <array>
#include <cstddef>
#include <cstdint>

#include <IPhreeqcPOET.hpp>
#include <map>
#include <utility>
#include <vector>

namespace poet {

using TugType = double;

class InitialList {
public:
  InitialList(RInside &R) : R(R){};

  void initializeFromList(const Rcpp::List &setup);

  void importList(const std::string &file_name);
  Rcpp::List exportList(const std::string &file_name);

private:
  RInside &R;

  // Grid members
  static constexpr const char *grid_key = "Grid";

  enum class GridMembers {
    PQC_SCRIPT_STRING,
    PQC_SCRIPT_FILE,
    PQC_DB_STRING,
    PQC_DB_FILE,
    GRID_DEF,
    GRID_SIZE,
    CONSTANT_CELLS,
    POROSITY,
    ENUM_SIZE
  };

  static constexpr std::size_t size_GridMembers =
      static_cast<std::size_t>(InitialList::GridMembers::ENUM_SIZE);

  static constexpr std::array<const char *, size_GridMembers>
      GridMembersString = {"pqc_in_string",  "pqc_in_file", "pqc_db_string",
                           "pqc_db_file",    "grid_def",    "grid_size",
                           "constant_cells", "porosity"};

  constexpr const char *getGridMemberString(GridMembers member) const {
    return GridMembersString[static_cast<std::size_t>(member)];
  }

  void initGrid(const Rcpp::List &grid_input);
  std::uint8_t dim;

  std::uint32_t n_cols;
  std::uint32_t n_rows;

  double s_cols;
  double s_rows;

  std::vector<std::uint32_t> constant_cells;

  Rcpp::List initial_grid;

  // No export
  Rcpp::NumericMatrix phreeqc_mat;

  // Initialized by grid, modified by chemistry
  std::map<int, std::string> pqc_raw_dumps;

  // Chemistry members
  IPhreeqcPOET::ModulesArray module_sizes;

  // Diffusion members

  static constexpr const char *diffusion_key = "Diffusion";

  void initDiffusion(const Rcpp::List &diffusion_input);
  Rcpp::List boundaries;

  Rcpp::List alpha_x;
  Rcpp::List alpha_y;
};
} // namespace poet