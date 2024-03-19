#include <tug/Boundary.hpp>
// leave above Rcpp includes, as eigen seem to have problems with a preceding
// Rcpp include

#include <map>
#include <set>

#include "InitialList.hpp"

#include <Rcpp/proxy/ProtectedProxy.h>
#include <Rcpp/vector/instantiation.h>
#include <Rinternals.h>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace poet {

enum SEXP_TYPE { SEXP_IS_LIST, SEXP_IS_MAT, SEXP_IS_NUM };

const std::map<std::uint8_t, std::string> tug_side_mapping = {
    {tug::BC_SIDE_RIGHT, "E"},
    {tug::BC_SIDE_LEFT, "W"},
    {tug::BC_SIDE_TOP, "N"},
    {tug::BC_SIDE_BOTTOM, "S"}};

static Rcpp::List parseBoundaries2D(const Rcpp::List &boundaries_list,
                                    std::uint32_t n_cols,
                                    std::uint32_t n_rows) {
  Rcpp::List out_list;

  for (const auto &side : tug_side_mapping) {
    Rcpp::NumericVector type =
        (side.first == tug::BC_SIDE_RIGHT || side.first == tug::BC_SIDE_LEFT)
            ? Rcpp::NumericVector(n_rows, tug::BC_TYPE_CLOSED)
            : Rcpp::NumericVector(n_cols, tug::BC_TYPE_CLOSED);

    Rcpp::NumericVector value =
        (side.first == tug::BC_SIDE_RIGHT || side.first == tug::BC_SIDE_LEFT)
            ? Rcpp::NumericVector(n_rows, 0)
            : Rcpp::NumericVector(n_cols, 0);

    if (boundaries_list.containsElementNamed(side.second.c_str())) {
      const Rcpp::List mapping = boundaries_list[side.second];

      const Rcpp::NumericVector cells = mapping["cells"];
      const Rcpp::NumericVector values = mapping["sol_id"];

      if (cells.size() != values.size()) {
        throw std::runtime_error("Boundary [" + side.second +
                                 "] cells and values are not the same length");
      }
    }

    out_list[side.second] = Rcpp::List::create(Rcpp::Named("type") = type,
                                               Rcpp::Named("value") = value);
  }
}

static inline SEXP_TYPE get_datatype(const SEXP &input) {
  Rcpp::Function check_list("is.list");
  Rcpp::Function check_mat("is.matrix");

  if (Rcpp::as<bool>(check_list(input))) {
    return SEXP_IS_LIST;
  } else if (Rcpp::as<bool>(check_mat(input))) {
    return SEXP_IS_MAT;
  } else {
    return SEXP_IS_NUM;
  }
}

static Rcpp::List parseAlphas2D(const SEXP &input,
                                const std::vector<std::string> &transport_names,
                                std::uint32_t n_cols, std::uint32_t n_rows) {
  Rcpp::List out_list;

  SEXP_TYPE input_type = get_datatype(input);

  switch (input_type) {
  case SEXP_IS_LIST: {
    const Rcpp::List input_list(input);

    if (input_list.size() != transport_names.size()) {
      throw std::runtime_error(
          "Length of alphas list does not match number of transport names");
    }

    for (const auto &name : transport_names) {
      if (!input_list.containsElementNamed(name.c_str())) {
        throw std::runtime_error("Alphas list does not contain transport name");
      }

      const Rcpp::NumericMatrix alpha_mat(input_list[name]);

      if (alpha_mat.size() == 1) {
        Rcpp::NumericVector alpha(n_cols * n_rows, alpha_mat(0, 0));
        out_list[name] = Rcpp::wrap(alpha);
      } else {
        if (alpha_mat.nrow() != n_rows || alpha_mat.ncol() != n_cols) {
          throw std::runtime_error(
              "Alpha matrix does not match grid dimensions");
        }

        out_list[name] = alpha_mat;
      }
    }
    break;
  }
  case SEXP_IS_MAT: {
    Rcpp::NumericMatrix input_mat(input);

    Rcpp::NumericVector alpha(n_cols * n_rows, input_mat(0, 0));

    if (input_mat.size() != 1) {
      if (input_mat.nrow() != n_rows || input_mat.ncol() != n_cols) {
        throw std::runtime_error("Alpha matrix does not match grid dimensions");
      }

      for (std::size_t i = 0; i < n_rows; i++) {
        for (std::size_t j = 0; j < n_cols; j++) {
          alpha[i * n_cols + j] = input_mat(i, j);
        }
      }
    }

    for (const auto &name : transport_names) {
      out_list[name] = alpha;
    }

    break;
  }
  case SEXP_IS_NUM: {
    Rcpp::NumericVector alpha(n_cols * n_rows, Rcpp::as<TugType>(input));
    for (const auto &name : transport_names) {
      out_list[name] = alpha;
    }
    break;
  }
  default:
    throw std::runtime_error("Unknown data type");
  }

  return out_list;
}
void InitialList::initDiffusion(const Rcpp::List &diffusion_input) {
  const Rcpp::List &boundaries =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::BOUNDARIES)];
  const Rcpp::NumericVector &alpha_x =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::ALPHA_X)];
  const Rcpp::NumericVector &alpha_y =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::ALPHA_Y)];

  std::vector<std::string> colnames =
      Rcpp::as<std::vector<std::string>>(this->initial_grid.names());

  std::vector<std::string> transport_names(colnames.begin() + 1,
                                           colnames.begin() + 1 +
                                               this->module_sizes[POET_SOL]);

  this->alpha_x =
      parseAlphas2D(alpha_x, transport_names, this->n_cols, this->n_rows);

  this->alpha_y =
      parseAlphas2D(alpha_y, transport_names, this->n_cols, this->n_rows);
}
} // namespace poet