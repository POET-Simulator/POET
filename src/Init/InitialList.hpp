#pragma once

#include <RInside.h>
#include <Rcpp.h>
#include <Rcpp/vector/instantiation.h>
#include <array>
#include <cstddef>
#include <cstdint>

#include <IPhreeqcPOET.hpp>
#include <string>
#include <vector>

#include <DataStructures/Field.hpp>

namespace poet {

using TugType = double;

class InitialList {
public:
  InitialList(RInside &R) : R(R){};

  void initializeFromList(const Rcpp::List &setup);

  void importList(const Rcpp::List &setup);
  Rcpp::List exportList();

private:
  RInside &R;

  enum class ExportList {
    GRID_DIM,
    GRID_SPECS,
    GRID_SPATIAL,
    GRID_CONSTANT,
    GRID_POROSITY,
    DIFFU_TRANSPORT,
    DIFFU_BOUNDARIES,
    DIFFU_ALPHA_X,
    DIFFU_ALPHA_Y,
    CHEM_DATABASE,
    CHEM_PQC_SCRIPTS,
    CHEM_PQC_IDS,
    CHEM_PQC_SOL_ORDER,
    ENUM_SIZE
  };

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

  constexpr const char *GRID_MEMBER_STR(GridMembers member) const {
    return GridMembersString[static_cast<std::size_t>(member)];
  }

  void initGrid(const Rcpp::List &grid_input);
  std::uint8_t dim;

  std::uint32_t n_cols;
  std::uint32_t n_rows;

  double s_cols;
  double s_rows;

  std::vector<std::uint32_t> constant_cells;
  std::vector<double> porosity;

  Rcpp::List initial_grid;

  // No export
  Rcpp::NumericMatrix phreeqc_mat;

public:
  struct DiffusionInit {
    std::uint32_t n_cols;
    std::uint32_t n_rows;

    double s_cols;
    double s_rows;

    std::vector<std::uint32_t> constant_cells;

    std::vector<std::string> transport_names;

    Field initial_grid;
    Field boundaries;
    Field alpha_x;
    Field alpha_y;
  };

  DiffusionInit getDiffusionInit() const;

private:
  // Diffusion members
  static constexpr const char *diffusion_key = "Diffusion";

  enum class DiffusionMembers { BOUNDARIES, ALPHA_X, ALPHA_Y, ENUM_SIZE };

  static constexpr std::size_t size_DiffusionMembers =
      static_cast<std::size_t>(InitialList::DiffusionMembers::ENUM_SIZE);

  static constexpr std::array<const char *, size_DiffusionMembers>
      DiffusionMembersString = {"boundaries", "alpha_x", "alpha_y"};

  constexpr const char *DIFFU_MEMBER_STR(DiffusionMembers member) const {
    return DiffusionMembersString[static_cast<std::size_t>(member)];
  }

  void initDiffusion(const Rcpp::List &diffusion_input);
  Rcpp::List resolveBoundaries(const Rcpp::List &boundaries_list);

  Rcpp::List boundaries;

  Rcpp::List alpha_x;
  Rcpp::List alpha_y;

  std::vector<std::string> transport_names;

  // Chemistry Members
  static constexpr const char *chemistry_key = "Chemistry";

  void initChemistry();

  std::string database;
  std::vector<std::string> pqc_scripts;
  std::vector<int> pqc_ids;

  std::vector<std::string> pqc_sol_order;

public:
  struct ChemistryInit {
    Field initial_grid;

    std::string database;
    std::vector<std::string> pqc_scripts;
    std::vector<int> pqc_ids;
    std::vector<std::string> pqc_sol_order;
  };

  ChemistryInit getChemistryInit() const;
};
} // namespace poet