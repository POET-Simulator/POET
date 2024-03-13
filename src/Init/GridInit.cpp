#include "InitialList.hpp"

#include <IPhreeqcPOET.hpp>
#include <Rcpp/Function.h>
#include <Rcpp/vector/instantiation.h>
#include <map>
#include <regex>
#include <string>

namespace poet {
static Rcpp::NumericMatrix pqcScriptToGrid(IPhreeqcPOET &phreeqc, RInside &R) {

  const auto solution_ids = phreeqc.getSolutionIds();
  auto col_names = phreeqc.getInitNames();
  const auto concentrations = phreeqc.getInitValues();

  // add "id" to the front of the column names

  const std::size_t ncols = col_names.size();
  const std::size_t nrows = solution_ids.size();

  col_names.insert(col_names.begin(), "id");

  Rcpp::NumericMatrix mat(nrows, ncols + 1);

  for (std::size_t i = 0; i < nrows; i++) {
    mat(i, 0) = solution_ids[i];
    for (std::size_t j = 0; j < ncols; ++j) {
      mat(i, j + 1) = concentrations[i][j];
    }
  }

  Rcpp::colnames(mat) = Rcpp::wrap(col_names);

  return mat;
}

static inline Rcpp::List matToGrid(RInside &R, const Rcpp::NumericMatrix &mat,
                                   const Rcpp::NumericMatrix &grid) {
  Rcpp::Function pqc_to_grid("pqc_to_grid");

  return pqc_to_grid(mat, grid);
}

static inline std::map<int, std::string>
replaceRawSolutionsIDs(std::map<int, std::string> raws) {
  for (auto &raw : raws) {
    std::string &s = raw.second;
    // find at beginning of line '*_RAW' followed by a number and change this
    // number to 1
    std::regex re(R"((RAW\s+)(\d+))");
    s = std::regex_replace(s, re, "RAW 1");
  }

  return raws;
}

void InitialList::initGrid(const Rcpp::List &grid_input) {
  // parse input values
  const std::string script = Rcpp::as<std::string>(grid_input["pqc_in"]);
  const std::string database = Rcpp::as<std::string>(grid_input["pqc_db"]);

  std::string script_rp(PATH_MAX, '\0');
  std::string database_rp(PATH_MAX, '\0');

  if (realpath(script.c_str(), script_rp.data()) == nullptr) {
    throw std::runtime_error("Failed to resolve the path of the Phreeqc input" +
                             script);
  }

  if (realpath(database.c_str(), database_rp.data()) == nullptr) {
    throw std::runtime_error("Failed to resolve the path of the database" +
                             database);
  }

  Rcpp::NumericMatrix grid_def = grid_input["grid_def"];
  Rcpp::NumericVector grid_size = grid_input["grid_size"];
  // Rcpp::NumericVector constant_cells = grid["constant_cells"].;

  this->n_rows = grid_def.nrow();
  this->n_cols = grid_def.ncol();

  this->s_rows = grid_size[0];
  this->s_cols = grid_size[1];

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

  IPhreeqcPOET phreeqc(database_rp, script_rp);

  this->phreeqc_mat = pqcScriptToGrid(phreeqc, R);
  this->initial_grid = matToGrid(R, this->phreeqc_mat, grid_def);

  this->pqc_raw_dumps = replaceRawSolutionsIDs(phreeqc.raw_dumps());

  R["pqc_mat"] = this->phreeqc_mat;
  R["grid_def"] = initial_grid;

  R.parseEval("print(pqc_mat)");
  R.parseEval("print(grid_def)");
}
} // namespace poet