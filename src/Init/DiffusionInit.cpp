#include <tug/Boundary.hpp>
// leave above Rcpp includes, as eigen seem to have problems with a preceding
// Rcpp include

#include <map>
#include <set>

#include "DataStructures/Field.hpp"
#include "IPhreeqcPOET.hpp"
#include "InitialList.hpp"

#include <Rcpp/Function.h>
#include <Rcpp/proxy/ProtectedProxy.h>
#include <Rcpp/vector/instantiation.h>
#include <Rinternals.h>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
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

static std::vector<std::string> extend_transport_names(
    std::unique_ptr<IPhreeqcPOET> &phreeqc, const Rcpp::List &boundaries_list,
    const Rcpp::List &inner_boundaries,
    const std::vector<std::string> &old_trans_names, Rcpp::List &initial_grid) {

  std::vector<std::string> transport_names = old_trans_names;
  std::set<int> constant_pqc_ids;

  for (const auto &side : tug_side_mapping) {
    if (!boundaries_list.containsElementNamed(side.second.c_str())) {
      continue;
    }

    Rcpp::List mapping = boundaries_list[side.second];

    const Rcpp::NumericVector cells = mapping["cell"];
    const Rcpp::NumericVector values = mapping["sol_id"];
    const Rcpp::CharacterVector type_str = mapping["type"];

    if (cells.size() != values.size()) {
      throw std::runtime_error("Boundary [" + side.second +
                               "] cells and values are not the same "
                               "length");
    }

    for (std::size_t i = 0; i < cells.size(); i++) {
      if (type_str[i] == "constant") {
        constant_pqc_ids.insert(values[i]);
      }
    }
  }

  if (inner_boundaries.size() > 0) {
    const Rcpp::NumericVector values = inner_boundaries["sol_id"];
    for (std::size_t i = 0; i < values.size(); i++) {
      constant_pqc_ids.insert(values[i]);
    }
  }

  if (!constant_pqc_ids.empty()) {
    for (const auto &pqc_id : constant_pqc_ids) {
      const auto solution_names = phreeqc->getSolutionNames(pqc_id);

      // add those strings which are not already in the transport_names
      for (const auto &name : solution_names) {
        if (std::find(transport_names.begin(), transport_names.end(), name) ==
            transport_names.end()) {
          transport_names.push_back(name);
        }
      }
    }
  }

  return transport_names;
}

static Rcpp::List extend_initial_grid(const Rcpp::List &initial_grid,
                                      std::vector<std::string> transport_names,
                                      std::size_t old_size) {
  std::vector<std::string> names_to_add(transport_names.begin() + old_size,
                                        transport_names.end());

  Rcpp::Function extend_grid_R("add_missing_transport_species");

  return extend_grid_R(initial_grid, Rcpp::wrap(names_to_add), old_size);
}

std::pair<Rcpp::List, Rcpp::List>
InitialList::resolveBoundaries(const Rcpp::List &boundaries_list,
                               const Rcpp::List &inner_boundaries) {
  Rcpp::List bound_list;
  Rcpp::List inner_bound;
  Rcpp::Function resolve_R("resolve_pqc_bound");

  const std::size_t old_transport_size = this->transport_names.size();

  this->transport_names =
      extend_transport_names(this->phreeqc, boundaries_list, inner_boundaries,
                             this->transport_names, this->initial_grid);

  const std::size_t new_transport_size = this->transport_names.size();

  if (old_transport_size != new_transport_size) {
    this->initial_grid = extend_initial_grid(
        this->initial_grid, this->transport_names, old_transport_size);
  }

  const Rcpp::NumericVector &inner_row_vec = inner_boundaries["row"];
  const Rcpp::NumericVector &inner_col_vec = inner_boundaries["col"];
  const Rcpp::NumericVector &inner_pqc_id_vec = inner_boundaries["sol_id"];

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

    if (inner_boundaries.size() > 0) {

      std::vector<std::uint32_t> rows;
      std::vector<std::uint32_t> cols;
      std::vector<TugType> c_value;

      if (inner_row_vec.size() != inner_col_vec.size() ||
          inner_row_vec.size() != inner_pqc_id_vec.size()) {
        throw std::runtime_error(
            "Inner boundary vectors are not the same length");
      }

      for (std::size_t i = 0; i < inner_row_vec.size(); i++) {
        rows.push_back(inner_row_vec[i] - 1);
        cols.push_back(inner_col_vec[i] - 1);
        c_value.push_back(Rcpp::as<TugType>(resolve_R(
            this->phreeqc_mat, Rcpp::wrap(species), inner_pqc_id_vec[i])));
      }

      inner_bound[species] =
          Rcpp::List::create(Rcpp::Named("row") = Rcpp::wrap(rows),
                             Rcpp::Named("col") = Rcpp::wrap(cols),
                             Rcpp::Named("value") = Rcpp::wrap(c_value));
    }
  }

  return std::make_pair(bound_list, inner_bound);
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
  Rcpp::List boundaries;
  Rcpp::List inner_boundaries;

  if (diffusion_input.containsElementNamed(
          DIFFU_MEMBER_STR(DiffusionMembers::BOUNDARIES))) {
    boundaries =
        diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::BOUNDARIES)];
  }

  if (diffusion_input.containsElementNamed(
          DIFFU_MEMBER_STR(DiffusionMembers::INNER_BOUNDARIES))) {
    inner_boundaries =
        diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::INNER_BOUNDARIES)];
  }

  const SEXP &alpha_x =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::ALPHA_X)];
  const SEXP &alpha_y =
      diffusion_input[DIFFU_MEMBER_STR(DiffusionMembers::ALPHA_Y)];

  const auto resolved_boundaries =
      resolveBoundaries(boundaries, inner_boundaries);
  this->boundaries = resolved_boundaries.first;
  this->inner_boundaries = resolved_boundaries.second;

  this->alpha_x =
      parseAlphas(alpha_x, this->transport_names, this->n_cols, this->n_rows);

  this->alpha_y =
      parseAlphas(alpha_y, this->transport_names, this->n_cols, this->n_rows);
}

InitialList::DiffusionInit::BoundaryMap
RcppListToBoundaryMap(const std::vector<std::string> &trans_names,
                      const Rcpp::List &bound_list, uint32_t n_cols,
                      uint32_t n_rows) {
  InitialList::DiffusionInit::BoundaryMap map;

  for (const auto &name : trans_names) {
    const Rcpp::List &conc_list = bound_list[name];

    tug::Boundary<TugType> bc(n_rows, n_cols);

    for (const auto &[tug_index, r_init_name] : tug_side_mapping) {
      const Rcpp::List &side_list = conc_list[tug_index];

      const Rcpp::NumericVector &type = side_list["type"];
      const Rcpp::NumericVector &value = side_list["value"];

      if (type.size() != value.size()) {
        throw std::runtime_error(
            "Boundary type and value are not the same length");
      }

      for (std::size_t i = 0; i < type.size(); i++) {
        if (type[i] == tug::BC_TYPE_CONSTANT) {
          bc.setBoundaryElementConstant(static_cast<tug::BC_SIDE>(tug_index), i,
                                        value[i]);
        }
      }
    }

    map[name] = bc.serialize();
  }

  return map;
}
static InitialList::DiffusionInit::InnerBoundaryMap
RcppListToInnerBoundaryMap(const std::vector<std::string> &trans_names,
                           const Rcpp::List &inner_bound_list,
                           std::uint32_t n_cols, std::uint32_t n_rows) {
  InitialList::DiffusionInit::InnerBoundaryMap map;

  if (inner_bound_list.size() == 0) {
    return map;
  }

  for (const auto &name : trans_names) {
    const Rcpp::List &conc_list = inner_bound_list[name];

    std::map<std::pair<std::uint32_t, std::uint32_t>, TugType> inner_bc;

    const Rcpp::NumericVector &row = conc_list["row"];
    const Rcpp::NumericVector &col = conc_list["col"];
    const Rcpp::NumericVector &value = conc_list["value"];

    if (row.size() != col.size() || row.size() != value.size()) {
      throw std::runtime_error(
          "Inner boundary vectors are not the same length");
    }

    for (std::size_t i = 0; i < row.size(); i++) {
      inner_bc[std::make_pair(row[i], col[i])] = value[i];
    }

    map[name] = inner_bc;
  }

  return map;
}

InitialList::DiffusionInit InitialList::getDiffusionInit() const {
  DiffusionInit diff_init;

  diff_init.dim = this->dim;

  diff_init.n_cols = this->n_cols;
  diff_init.n_rows = this->n_rows;

  diff_init.s_cols = this->s_cols;
  diff_init.s_rows = this->s_rows;

  diff_init.constant_cells = this->constant_cells;
  diff_init.transport_names = this->transport_names;

  diff_init.boundaries = RcppListToBoundaryMap(
      this->transport_names, this->boundaries, this->n_cols, this->n_rows);
  diff_init.inner_boundaries =
      RcppListToInnerBoundaryMap(this->transport_names, this->inner_boundaries,
                                 this->n_cols, this->n_rows);
  diff_init.alpha_x = Field(this->alpha_x);
  diff_init.alpha_y = Field(this->alpha_y);

  return diff_init;
}

} // namespace poet