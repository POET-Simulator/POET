#include "AdvectionModule.hpp"

#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <limits>
#include <string>
#include <vector>

#include <Rcpp.h>

namespace poet {

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
                               const std::vector<double> &conc, double bc) {
  // On inbound flux and non-boundary condition
  if (inbound) {
    if (neighbor_index == -1) {
      return bc;
    }
    return conc[neighbor_index];
  }

  // inbound flow from boundary condition (open) or outbound flow
  return conc[curr_index];
}

inline double
AdvectionModule::calcDeltaConc(std::size_t cell_index, double bc_val,
                               const std::vector<double> &spec_vec,
                               const std::vector<double> &flux) {
  const auto neighbor_indices =
      getIndices(cell_index, this->n_cells[0], this->n_cells[1]);
  double ds{.0};
  for (std::size_t neighbor_i = 0; neighbor_i < neighbor_indices.size();
       neighbor_i++) {
    const double &curr_flux = flux[neighbor_i];
    const bool inbound = curr_flux > 0;
    const double flux_apply_val = getFluxApplyConc(
        inbound, cell_index, neighbor_indices[neighbor_i], spec_vec, bc_val);
    ds += curr_flux * flux_apply_val;
  }

  return ds;
}

std::vector<double>
AdvectionModule::CFLTimeVec(double req_dt,
                            const std::vector<double> &max_fluxes) {
  const auto field_size = this->n_cells[0] * this->n_cells[1];

  double max_dt = std::numeric_limits<double>::max();

  for (std::size_t i = 0; i < field_size; i++) {
    const double dt =
        (cell_volume * density[i] * water_saturation[i] * porosity[i]) /
        max_fluxes[i];
    if (dt < max_dt) {
      max_dt = dt;
    }
  }

  if (max_dt < req_dt) {
    std::uint32_t min_floor = static_cast<std::uint32_t>(req_dt / max_dt);
    std::vector<double> time_vec(min_floor, max_dt);

    double diff = req_dt - (max_dt * min_floor);
    time_vec.push_back(req_dt);
    return time_vec;
  }

  return std::vector<double>(1, req_dt);
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

  MSG("Advection time step requested: " + std::to_string(dt));

  std::vector<double> max_fluxes(flux.size());
  for (std::size_t i = 0; i < max_fluxes.size(); i++) {
    std::array<double, 4> abs_flux;
    for (std::size_t j = 0; j < abs_flux.size(); j++) {
      abs_flux[j] = std::abs(flux[i][j]);
    }
    max_fluxes[i] = *std::max_element(abs_flux.begin(), abs_flux.end());
  }

  const auto time_vec = CFLTimeVec(dt, max_fluxes);

  MSG("CFL yielding " + std::to_string(time_vec.size()) + " inner iterations");

  auto field_vec = t_field.As2DVector();

  // iterate over all inactive cells and set defined values, as chemistry module
  // won't skip those cells until now

  for (const auto &inac_cell : this->inactive_cells) {
    for (std::size_t species_i = 0; species_i < field_vec.size(); species_i++) {
      field_vec[species_i][inac_cell.first] = inac_cell.second[species_i];
    }
  }

  for (std::size_t species_i = 0; species_i < field_vec.size(); species_i++) {
    auto &species = field_vec[species_i];
    std::vector<double> spec_copy = species;
    for (const double &dt : time_vec) {
      for (std::size_t cell_i = 0; cell_i < n * m; cell_i++) {

        // if inactive cell -> skip
        const auto inactive_cell_it = this->inactive_cells.find(cell_i);
        if (this->inactive_cells.find(cell_i) != this->inactive_cells.end()) {
          continue;
        }

        double delta_conc =
            calcDeltaConc(cell_i, this->boundary_condition[species_i], species,
                          flux[species_i]);
        spec_copy[cell_i] +=
            (dt * delta_conc) / (porosity[cell_i] * cell_volume);
      }
      species = spec_copy;
    }
  }

  t_field = field_vec;
}

void AdvectionModule::initializeParams(RInsidePOET &R) {
  const std::uint32_t n = this->n_cells[0] =
      R.parseEval("mysetup$grid$n_cells[1]");
  const std::uint32_t m = this->n_cells[1] =
      R.parseEval("mysetup$grid$n_cells[2]");

  const std::uint32_t x = this->s_cells[0] =
      R.parseEval("mysetup$grid$s_cells[1]");
  const std::uint32_t y = this->s_cells[1] =
      R.parseEval("mysetup$grid$s_cells[2]");

  this->field_size = n * m;
  this->cell_volume = (s_cells[0] * s_cells[1]) / this->field_size;

  // HACK: For now, we neglect porisity, density and water saturation of the
  // cell
  this->porosity = std::vector<double>(field_size, 1.);
  this->density = std::vector<double>(field_size, 1.);
  this->water_saturation = std::vector<double>(field_size, 1.);

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

  for (std::size_t i = 0; i < prop_names.size(); i++) {
    const Rcpp::NumericVector species_vecinj = vecinj[prop_names[i]];
    this->boundary_condition.push_back(species_vecinj[0]);
  }

  // same logic as for diffusion module applied
  for (std::size_t i = 0; i < inner_vecinj.size(); i++) {
    const Rcpp::NumericVector tuple = inner_vecinj[i];
    const std::uint32_t cell_index = (tuple[2] - 1) + ((tuple[1] - 1) * n);
    std::vector<double> curr_cell_state(prop_names.size());
    for (std::size_t prop_i = 0; prop_i < prop_names.size(); prop_i++) {
      const Rcpp::NumericVector curr_prop_vec_inj = vecinj[prop_names[prop_i]];
      init_field[prop_i][cell_index] = curr_prop_vec_inj[tuple[0] - 1];
      curr_cell_state[prop_i] = curr_prop_vec_inj[tuple[0] - 1];
    }
    this->inactive_cells.insert({cell_index, std::move(curr_cell_state)});
  }

  this->t_field = Field(field_size, init_field, prop_names);
}

} // namespace poet
