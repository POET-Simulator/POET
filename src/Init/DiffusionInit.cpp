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

enum SEXP_TYPE { SEXP_IS_LIST, SEXP_IS_VEC };

const std::map<std::uint8_t, std::string> tug_side_mapping = {
    {tug::BC_SIDE_RIGHT, "E"},
    {tug::BC_SIDE_LEFT, "W"},
    {tug::BC_SIDE_TOP, "N"},
    {tug::BC_SIDE_BOTTOM, "S"}};

static std::vector<TugType> colMajToRowMaj(const Rcpp::NumericVector &vec,
                                           std::uint32_t n_cols,
                                           std::uint32_t n_rows) {
  if (vec.size() == 1) {
    return std::vector<TugType>(n_cols * n_rows, vec[0]);
  } else {
    if (vec.size() != n_cols * n_rows) {
      throw std::runtime_error("Alpha matrix does not match grid dimensions");
    }

    std::vector<TugType> alpha(n_cols * n_rows);

    for (std::uint32_t i = 0; i < n_cols; i++) {
      for (std::uint32_t j = 0; j < n_rows; j++) {
        alpha[i * n_rows + j] = vec[j * n_cols + i];
      }
    }
    return alpha;
  }
}

Rcpp::List InitialList::resolveBoundaries(const Rcpp::List &boundaries_list) {
  Rcpp::List bound_list;
  Rcpp::Function resolve_R("resolvePqcBound");

  for (const auto &species : this->transport_names) {
    Rcpp::List spec_list;

    for (const auto &side : tug_side_mapping) {
      Rcpp::NumericVector c_type =
          (side.first == tug::BC_SIDE_RIGHT || side.first == tug::BC_SIDE_LEFT)
              ? Rcpp::NumericVector(n_rows, tug::BC_TYPE_CLOSED)
              : Rcpp::NumericVector(n_cols, tug::BC_TYPE_CLOSED);

      Rcpp::NumericVector c_value =
          (side.first == tug::BC_SIDE_RIGHT || side.first == tug::BC_SIDE_LEFT)
              ? Rcpp::NumericVector(n_rows, 0)
              : Rcpp::NumericVector(n_cols, 0);

      if (boundaries_list.containsElementNamed(side.second.c_str())) {
        const Rcpp::List mapping = boundaries_list[side.second];

        const Rcpp::NumericVector cells = mapping["cell"];
        const Rcpp::NumericVector values = mapping["sol_id"];
        const Rcpp::CharacterVector type_str = mapping["type"];

        if (cells.size() != values.size()) {
          throw std::runtime_error(
              "Boundary [" + side.second +
              "] cells and values are not the same length");
        }

        for (std::size_t i = 0; i < cells.size(); i++) {
          if (type_str[i] == "closed") {
            c_type[cells[i] - 1] = tug::BC_TYPE_CLOSED;
          } else if (type_str[i] == "constant") {
            c_type[cells[i] - 1] = tug::BC_TYPE_CONSTANT;
            c_value[cells[i] - 1] = Rcpp::as<TugType>(
                resolve_R(this->phreeqc_mat, Rcpp::wrap(species), values[i]));
          } else {
            throw std::runtime_error("Unknown boundary type");
          }
        }
      }

      spec_list[side.second] = Rcpp::List::create(
          Rcpp::Named("type") = c_type, Rcpp::Named("value") = c_value);
    }

    bound_list[species] = spec_list;
  }

  return bound_list;
}

static inline SEXP_TYPE get_datatype(const SEXP &input) {
  Rcpp::Function check_list("is.list");

  if (Rcpp::as<bool>(check_list(input))) {
    return SEXP_IS_LIST;
  } else {
    return SEXP_IS_VEC;
  }
}

static Rcpp::List parseAlphas(const SEXP &input,
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

      const Rcpp::NumericVector &alpha_col_order_vec = input_list[name];

      out_list[name] =
          Rcpp::wrap(colMajToRowMaj(alpha_col_order_vec, n_cols, n_rows));
    }
    break;
  }
  case SEXP_IS_VEC: {
    const Rcpp::NumericVector alpha(input);
    for (const auto &name : transport_names) {
      out_list[name] = Rcpp::wrap(colMajToRowMaj(alpha, n_cols, n_rows));
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
  const SEXP &alpha_x =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::ALPHA_X)];
  const SEXP &alpha_y =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::ALPHA_Y)];

  this->boundaries = resolveBoundaries(boundaries);

  this->alpha_x =
      parseAlphas(alpha_x, this->transport_names, this->n_cols, this->n_rows);

  this->alpha_y =
      parseAlphas(alpha_y, this->transport_names, this->n_cols, this->n_rows);

  // R["alpha_x"] = this->alpha_x;
  // R["alpha_y"] = this->alpha_y;

  // R.parseEval("print(alpha_x)");
  // R.parseEval("print(alpha_y)");
}
} // namespace poet