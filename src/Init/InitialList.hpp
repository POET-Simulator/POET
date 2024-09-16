#pragma once

#include "Base/RInsidePOET.hpp"
#include "DataStructures/NamedVector.hpp"
#include "POETInit.hpp"
#include <Rcpp/vector/instantiation.h>
#include <set>
#include <tug/Boundary.hpp>

#include <RInside.h>
#include <Rcpp.h>
#include <array>
#include <cstddef>
#include <cstdint>

#include <PhreeqcInit.hpp>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <DataStructures/Field.hpp>

namespace poet {

using TugType = double;

class InitialList {
public:
  InitialList(RInside &R) : R(R){};

  void initializeFromList(const Rcpp::List &setup);

  void importList(const Rcpp::List &setup, bool minimal = false);
  Rcpp::List exportList();

  Field getInitialGrid() const { return Field(this->initial_grid); }  

private:
  RInside &R;

  enum class ExportList {
    GRID_DIM,
    GRID_SPECS,
    GRID_SPATIAL,
    GRID_CONSTANT,
    GRID_POROSITY,
    GRID_INITIAL,
    DIFFU_TRANSPORT,
    DIFFU_BOUNDARIES,
    DIFFU_INNER_BOUNDARIES,
    DIFFU_ALPHA_X,
    DIFFU_ALPHA_Y,
    CHEM_DATABASE,
    CHEM_FIELD_HEADER,
    CHEM_PQC_SCRIPTS,
    CHEM_PQC_IDS,
    CHEM_PQC_SOLUTIONS,
    CHEM_PQC_SOLUTION_PRIMARY,
    CHEM_PQC_EXCHANGER,
    CHEM_PQC_KINETICS,
    CHEM_PQC_EQUILIBRIUM,
    CHEM_PQC_SURFACE_COMPS,
    CHEM_PQC_SURFACE_CHARGES,
    CHEM_DHT_SPECIES,
    CHEM_INTERP_SPECIES,
    CHEM_HOOKS,
    AI_SURROGATE_INPUT_SCRIPT,
    ENUM_SIZE // Hack: Last element of the enum to show enum size
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

  std::unique_ptr<PhreeqcInit> phreeqc;

  void prepareGrid(const Rcpp::List &grid_input);

  std::uint8_t dim{0};

  std::uint32_t n_cols{0};
  std::uint32_t n_rows{0};

  double s_cols{0};
  double s_rows{0};

  std::vector<std::uint32_t> constant_cells;
  std::vector<double> porosity;

  Rcpp::List initial_grid;

  // No export
  Rcpp::NumericMatrix phreeqc_mat;

public:
  struct DiffusionInit {
    using BoundaryMap =
        std::map<std::string, std::vector<tug::BoundaryElement<TugType>>>;

    using InnerBoundaryMap =
        std::map<std::string,
                 std::map<std::pair<std::uint32_t, std::uint32_t>, TugType>>;

    uint8_t dim;
    std::uint32_t n_cols;
    std::uint32_t n_rows;

    double s_cols;
    double s_rows;

    std::vector<std::uint32_t> constant_cells;

    std::vector<std::string> transport_names;

    BoundaryMap boundaries;
    InnerBoundaryMap inner_boundaries;

    Field alpha_x;
    Field alpha_y;
  };

  DiffusionInit getDiffusionInit() const;

private:
  // Diffusion members
  static constexpr const char *diffusion_key = "Diffusion";

  enum class DiffusionMembers {
    BOUNDARIES,
    INNER_BOUNDARIES,
    ALPHA_X,
    ALPHA_Y,
    ENUM_SIZE
  };

  static constexpr std::size_t size_DiffusionMembers =
      static_cast<std::size_t>(InitialList::DiffusionMembers::ENUM_SIZE);

  static constexpr std::array<const char *, size_DiffusionMembers>
      DiffusionMembersString = {"boundaries", "inner_boundaries", "alpha_x",
                                "alpha_y"};

  constexpr const char *DIFFU_MEMBER_STR(DiffusionMembers member) const {
    return DiffusionMembersString[static_cast<std::size_t>(member)];
  }

  void initDiffusion(const Rcpp::List &diffusion_input);
  std::pair<Rcpp::List, Rcpp::List>
  resolveBoundaries(const Rcpp::List &boundaries_list,
                    const Rcpp::List &inner_boundaries);

  Rcpp::List boundaries;
  Rcpp::List inner_boundaries;

  Rcpp::List alpha_x;
  Rcpp::List alpha_y;

  std::vector<std::string> transport_names;

  // Chemistry Members
  static constexpr const char *chemistry_key = "Chemistry";

  void initChemistry(const Rcpp::List &chem_input);

  std::vector<std::string> field_header;

  std::string database;
  std::vector<std::string> pqc_scripts;
  std::vector<int> pqc_ids;

  std::vector<std::string> pqc_solutions;
  std::vector<std::string> pqc_solution_primaries;
  Rcpp::List pqc_exchanger;
  Rcpp::List pqc_kinetics;
  Rcpp::List pqc_equilibrium;
  Rcpp::List pqc_surface_comps;
  Rcpp::List pqc_surface_charges;

  NamedVector<std::uint32_t> dht_species;

  NamedVector<std::uint32_t> interp_species;
  
  // Path to R script that the user defines in the input file
  std::string ai_surrogate_input_script;

  Rcpp::List chem_hooks;

  const std::set<std::string> hook_name_list{"dht_fill", "dht_fuzz",
                                             "interp_pre", "interp_post"};

public:
  struct ChemistryHookFunctions {
    poet::DEFunc dht_fill;
    poet::DEFunc dht_fuzz;
    poet::DEFunc interp_pre;
    poet::DEFunc interp_post;
  };

  struct ChemistryInit {
    uint32_t total_grid_cells;

    // std::vector<std::string> field_header;

    std::string database;
    // std::vector<std::string> pqc_scripts;
    // std::vector<int> pqc_ids;

    std::map<int, POETConfig> pqc_config;

    // std::vector<std::string> pqc_sol_order;

    NamedVector<std::uint32_t> dht_species;
    NamedVector<std::uint32_t> interp_species;
    ChemistryHookFunctions hooks;

    std::string ai_surrogate_input_script;
  };

  ChemistryInit getChemistryInit() const;
};
} // namespace poet