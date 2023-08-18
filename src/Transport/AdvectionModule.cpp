#include "AdvectionModule.hpp"

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include <Rcpp.h>

namespace poet {

// inline std::vector<std::vector<double>> getFlux(const Rcpp::List &flux_list)
// {
//   std::vector<std::vector<double>> fluxes(flux_lisVt.size());
//   for (std::size_t i = 0; i < flux_list.size(); i++) {
//     const auto flux_2d = Rcpp::as<std::vector<double>>(flux_list[i + 1]);
//     fluxes[i] = flux_2d;
//   }

//   return fluxes;
// }

inline std::array<std::int32_t, 4> getIndices(std::int32_t index,
                                              std::int32_t n, std::int32_t m) {
  std::array<std::int32_t, 4> indices;

  // east index
  indices[0] = (index % n == n - 1 ? -1 : index + 1);
  // south index
  indices[1] = (index + n >= m * n ? -1 : index + n);
  // west index
  indices[2] = (index % n == 0 ? -1 : index - 1);
  // north index
  indices[3] = (index - n < 0 ? -1 : index - n);

  return indices;
}

inline double getFluxApplyConc(bool inbound, std::uint32_t curr_index,
                               std::int32_t neighbor_index,
                               const std::vector<double> &conc) {
  if (inbound && neighbor_index >= 0) {
    return conc[neighbor_index];
  }

  // outbound flow only affects current cell with given index
  return conc[curr_index];
}

void AdvectionModule::simulate(double dt) {
  const auto &n = this->n_cells[0];
  const auto &m = this->n_cells[1];

  // HACK: constant flux for this moment imported from R runtime

  RInsidePOET &R = RInsidePOET::getInstance();

  const auto flux_list =
      Rcpp::as<Rcpp::DataFrame>(R.parseEval("mysetup$advection$const_flux"));
  std::vector<std::vector<double>> flux(flux_list.size());
  for (std::size_t i = 0; i < flux_list.size(); i++) {
    const auto flux_2d = Rcpp::as<std::vector<double>>(flux_list[i]);
    flux[i] = flux_2d;
  }

  auto field_vec = t_field.As2DVector();
  for (auto &species : field_vec) {
    std::vector<double> spec_copy = species;
    for (std::size_t cell_i = 0; cell_i < n * m; cell_i++) {
      if (this->inactive_cells.find(cell_i) != this->inactive_cells.end()) {
        continue;
      }
      const auto indices = getIndices(cell_i, n, m);
      double ds{.0};
      for (std::size_t neighbor_i = 0; neighbor_i < indices.size();
           neighbor_i++) {
        const double &curr_flux = flux[cell_i][neighbor_i];
        const double flux_apply_val = getFluxApplyConc(
            curr_flux > 0, cell_i, indices[neighbor_i], species);
        ds += curr_flux * flux_apply_val;
      }
      spec_copy[cell_i] += ds;
    }
    species = spec_copy;
  }
  t_field = field_vec;
}

void AdvectionModule::initializeParams(RInsidePOET &R) {
  const std::uint32_t n = this->n_cells[0] =
      R.parseEval("mysetup$grid$n_cells[1]");
  const std::uint32_t m = this->n_cells[1] =
      R.parseEval("mysetup$grid$n_cells[2]");

  const std::uint32_t field_size{n * m};

  const auto init_vec =
      Rcpp::as<Rcpp::NumericVector>(R.parseEval("mysetup$advection$init"));

  // Initialize prop names and according values
  const auto prop_names = Rcpp::as<std::vector<std::string>>(init_vec.names());
  std::vector<std::vector<double>> init_field;
  init_field.reserve(prop_names.size());

  for (const auto &val : init_vec) {
    init_field.push_back(std::vector<double>(field_size, val));
  }

  // Set inner constant cells. Is it needed?
  const Rcpp::DataFrame vecinj = R.parseEval("mysetup$advection$vecinj");
  const Rcpp::DataFrame inner_vecinj =
      R.parseEval("mysetup$advection$vecinj_inner");

  // same logic as for diffusion module applied
  for (std::size_t i = 0; i < inner_vecinj.size(); i++) {
    const Rcpp::NumericVector tuple = inner_vecinj[i];
    const std::uint32_t cell_index = (tuple[2] - 1) + ((tuple[1] - 1) * n);
    for (std::size_t prop_i = 0; prop_i < prop_names.size(); prop_i++) {
      const Rcpp::NumericVector curr_prop_vec_inj = vecinj[prop_names[prop_i]];
      init_field[prop_i][cell_index] = curr_prop_vec_inj[tuple[0] - 1];
    }
    this->inactive_cells.insert(cell_index);
  }

  this->t_field = Field(field_size, init_field, prop_names);
}

} // namespace poet
