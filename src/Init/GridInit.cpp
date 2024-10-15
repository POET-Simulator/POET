#include "InitialList.hpp"

#include "PhreeqcMatrix.hpp"

#include <RInside.h>
#include <Rcpp/Function.h>
#include <Rcpp/vector/Matrix.h>
#include <Rcpp/vector/instantiation.h>
#include <cstdint>
#include <fstream>
#include <map>
#include <memory>
#include <regex>
#include <sstream>
#include <string>
#include <vector>

namespace poet {

// static Rcpp::NumericMatrix pqcMatToR(const PhreeqcMatrix &phreeqc, RInside
// &R) {

//   PhreeqcMatrix::STLExport phreeqc_mat = phreeqc.get();

//   // PhreeqcInit::PhreeqcMat phreeqc_mat = phreeqc->getPhreeqcMat();

//   // add "id" to the front of the column names

//   const std::size_t ncols = phreeqc_mat.names.size();
//   const std::size_t nrows = phreeqc_mat.values.size();

//   phreeqc_mat.names.insert(phreeqc_mat.names.begin(), "ID");

//   Rcpp::NumericMatrix mat(nrows, ncols + 1);

//   for (std::size_t i = 0; i < nrows; i++) {
//     mat(i, 0) = phreeqc_mat.ids[i];
//     for (std::size_t j = 0; j < ncols; ++j) {
//       mat(i, j + 1) = phreeqc_mat.values[i][j];
//     }
//   }

//   Rcpp::colnames(mat) = Rcpp::wrap(phreeqc_mat.names);

//   return mat;
// }

// static inline Rcpp::List matToGrid(RInside &R, const Rcpp::NumericMatrix
// &mat,
//                                    const Rcpp::NumericMatrix &grid) {
//   Rcpp::Function pqc_to_grid_R("pqc_to_grid");

//   return pqc_to_grid_R(mat, grid);
// }

static inline std::map<int, std::string>
replaceRawKeywordIDs(std::map<int, std::string> raws) {
  for (auto &raw : raws) {
    std::string &s = raw.second;
    // find at beginning of line '*_RAW' followed by a number and change this
    // number to 1
    std::regex re(R"((RAW\s+)(\d+))");
    s = std::regex_replace(s, re, "RAW 1");
  }

  return raws;
}

// static inline uint32_t getSolutionCount(std::unique_ptr<PhreeqcInit>
// &phreeqc,
//                                         const Rcpp::List &initial_grid) {
//   PhreeqcInit::ModulesArray mod_array;
//   Rcpp::Function unique_R("unique");

//   std::vector<int> row_ids =
//       Rcpp::as<std::vector<int>>(unique_R(initial_grid["ID"]));

//   // std::vector<std::uint32_t> sizes_vec(sizes.begin(), sizes.end());

//   // Rcpp::Function modify_sizes("modify_module_sizes");
//   // sizes_vec = Rcpp::as<std::vector<std::uint32_t>>(
//   //     modify_sizes(sizes_vec, phreeqc_mat, initial_grid));

//   // std::copy(sizes_vec.begin(), sizes_vec.end(), sizes.begin());

//   return phreeqc->getModuleSizes(row_ids)[POET_SOL];
// }

static std::string readFile(const std::string &path) {
  std::string string_rpath(PATH_MAX, '\0');

  if (realpath(path.c_str(), string_rpath.data()) == nullptr) {
    throw std::runtime_error("Failed to resolve the realpath to file " + path);
  }

  std::ifstream file(string_rpath);

  if (!file.is_open()) {
    throw std::runtime_error("Failed to open file: " + path);
  }

  std::stringstream buffer;
  buffer << file.rdbuf();

  return buffer.str();
}

static Rcpp::List expandGrid(const PhreeqcMatrix &pqc_mat,
                             const std::vector<int> unique_ids,
                             const Rcpp::IntegerMatrix &grid_def) {

  PhreeqcMatrix subset_pqc_mat = pqc_mat.subset(unique_ids);

  PhreeqcMatrix::STLExport phreeqc_mat =
      subset_pqc_mat.get(PhreeqcMatrix::VectorExportType::COLUMN_MAJOR);

  const std::size_t ncols = phreeqc_mat.names.size();
  const std::size_t nrows = phreeqc_mat.values.size() / ncols;

  Rcpp::NumericMatrix pqc_mat_R(nrows, ncols, phreeqc_mat.values.begin());
  Rcpp::colnames(pqc_mat_R) = Rcpp::wrap(phreeqc_mat.names);

  return Rcpp::Function("pqc_to_grid")(pqc_mat_R, grid_def);
}

PhreeqcMatrix InitialList::prepareGrid(const Rcpp::List &grid_input) {
  // parse input values
  std::string script;
  std::string database;

  if (grid_input.containsElementNamed(
          GRID_MEMBER_STR(GridMembers::PQC_SCRIPT_FILE))) {

    script = readFile(Rcpp::as<std::string>(
        grid_input[GRID_MEMBER_STR(GridMembers::PQC_SCRIPT_FILE)]));
  } else {
    script = Rcpp::as<std::string>(
        grid_input[GRID_MEMBER_STR(GridMembers::PQC_SCRIPT_STRING)]);
  }

  if (grid_input.containsElementNamed(
          GRID_MEMBER_STR(GridMembers::PQC_DB_FILE))) {

    database = readFile(Rcpp::as<std::string>(
        grid_input[GRID_MEMBER_STR(GridMembers::PQC_DB_FILE)]));
  } else {
    database = Rcpp::as<std::string>(
        grid_input[GRID_MEMBER_STR(GridMembers::PQC_DB_STRING)]);
  }

  this->database = database;
  this->pqc_script = script;

  Rcpp::IntegerMatrix grid_def =
      grid_input[GRID_MEMBER_STR(GridMembers::GRID_DEF)];
  Rcpp::NumericVector grid_size =
      grid_input[GRID_MEMBER_STR(GridMembers::GRID_SIZE)];
  // Rcpp::NumericVector constant_cells = grid["constant_cells"].;

  this->n_rows = grid_def.nrow();
  this->n_cols = grid_def.ncol();

  this->s_cols = grid_size[0];
  this->s_rows = grid_size[1];

  this->dim = n_cols == 1 ? 1 : 2;

  if (this->n_cols > 1 && this->n_rows == 1) {
    throw std::runtime_error(
        "Dimensions of grid definition does not match the expected format "
        "(n_rows > 1, n_cols >= 1).");
  }

  if (this->dim != grid_size.size()) {
    throw std::runtime_error(
        "Dimensions of grid does not match the dimension of grid definition.");
  }

  if (this->s_rows <= 0 || this->s_cols <= 0) {
    throw std::runtime_error("Grid size must be positive.");
  }

  PhreeqcMatrix pqc_mat = PhreeqcMatrix(database, script);

  this->transport_names = pqc_mat.getSolutionNames(true);

  Rcpp::Function unique_R("unique");
  Rcpp::Function sort_R("sort");

  std::vector<int> unique_ids = Rcpp::as<std::vector<int>>(
      unique_R(Rcpp::IntegerVector(grid_def.begin(), grid_def.end())));

  this->initial_grid = expandGrid(pqc_mat, unique_ids, grid_def);

  const auto pqc_raw_dumps = replaceRawKeywordIDs(pqc_mat.getDumpStringsPQI());

  for (const auto &id : unique_ids) {
    this->pqc_ids.push_back(id);
  }

  return pqc_mat;
  // this->phreeqc_mat = pqcScriptToGrid(phreeqc, R);
  // this->initial_grid = matToGrid(R, this->phreeqc_mat, grid_def);

  // const uint32_t solution_count = getSolutionCount(phreeqc,
  // this->initial_grid);

  // std::vector<std::string> colnames =
  //     Rcpp::as<std::vector<std::string>>(this->initial_grid.names());

  // this->transport_names = std::vector<std::string>(
  //     colnames.begin() + 1,
  //     colnames.begin() + 1 + solution_count); // skip ID

  // std::map<int, std::string> pqc_raw_dumps;

  // pqc_raw_dumps = replaceRawKeywordIDs(phreeqc->raw_dumps());

  // this->pqc_ids =
  //     Rcpp::as<std::vector<int>>(unique_R(this->initial_grid["ID"]));

  // for (const auto &id : this->pqc_ids) {
  //   this->pqc_scripts.push_back(pqc_raw_dumps[id]);
  //   // this->pqc_exchanger.push_back(phreeqc->getExchanger(id));
  //   // this->pqc_kinetics.push_back(phreeqc->getKineticsNames(id));
  //   // this->pqc_equilibrium.push_back(phreeqc->getEquilibriumNames(id));
  //   // this->pqc_surface_comps.push_back(phreeqc->getSurfaceCompNames(id));
  //   //
  //   this->pqc_surface_charges.push_back(phreeqc->getSurfaceChargeNames(id));
  // }
}

} // namespace poet