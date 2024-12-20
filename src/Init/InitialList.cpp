#include "InitialList.hpp"
#include "DataStructures/NamedVector.hpp"
#include "PhreeqcMatrix.hpp"
#include <Rcpp/internal/wrap.h>
#include <Rcpp/iostream/Rstreambuf.h>
#include <Rcpp/proxy/ProtectedProxy.h>
#include <Rcpp/vector/instantiation.h>
#include <cstdint>
#include <string>
#include <vector>

namespace poet {
void InitialList::initializeFromList(const Rcpp::List &setup) {
  PhreeqcMatrix phreeqc = prepareGrid(setup[grid_key]);
  initDiffusion(setup[diffusion_key], phreeqc);
  initChemistry(setup[chemistry_key]);
}

void InitialList::importList(const Rcpp::List &setup, bool minimal) {
  this->dim = Rcpp::as<uint8_t>(setup[static_cast<int>(ExportList::GRID_DIM)]);

  Rcpp::NumericVector grid_specs =
      setup[static_cast<int>(ExportList::GRID_SPECS)];
  this->n_rows = static_cast<std::uint32_t>(grid_specs[0]);
  this->n_cols = static_cast<std::uint32_t>(grid_specs[1]);

  Rcpp::NumericVector spatial =
      setup[static_cast<int>(ExportList::GRID_SPATIAL)];
  this->s_rows = spatial[0];
  this->s_cols = spatial[1];

  this->constant_cells = Rcpp::as<std::vector<uint32_t>>(
      setup[static_cast<int>(ExportList::GRID_CONSTANT)]);

  this->porosity = Rcpp::as<std::vector<double>>(
      setup[static_cast<int>(ExportList::GRID_POROSITY)]);

  if (!minimal) {
    this->initial_grid =
        Rcpp::List(setup[static_cast<int>(ExportList::GRID_INITIAL)]);
  }

  this->transport_names = Rcpp::as<std::vector<std::string>>(
      setup[static_cast<int>(ExportList::DIFFU_TRANSPORT)]);
  if (!minimal) {
    this->boundaries =
        Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_BOUNDARIES)]);
    this->inner_boundaries =
        Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_INNER_BOUNDARIES)]);
    this->alpha_x =
        Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_ALPHA_X)]);
    this->alpha_y =
        Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_ALPHA_Y)]);
  }
  this->database =
      Rcpp::as<std::string>(setup[static_cast<int>(ExportList::CHEM_DATABASE)]);
  this->pqc_script = Rcpp::as<std::string>(
      setup[static_cast<int>(ExportList::CHEM_PQC_SCRIPT)]);
  this->with_h0_o0 =
      Rcpp::as<bool>(setup[static_cast<int>(ExportList::CHEM_PQC_WITH_H0_O0)]);
  this->with_redox =
      Rcpp::as<bool>(setup[static_cast<int>(ExportList::CHEM_PQC_WITH_REDOX)]);
  this->field_header = Rcpp::as<std::vector<std::string>>(
      setup[static_cast<int>(ExportList::CHEM_FIELD_HEADER)]);
  this->pqc_ids = Rcpp::as<std::vector<int>>(
      setup[static_cast<int>(ExportList::CHEM_PQC_IDS)]);
  //   this->pqc_solutions = Rcpp::as<std::vector<std::string>>(
  //       setup[static_cast<int>(ExportList::CHEM_PQC_SOLUTIONS)]);
  //   this->pqc_solution_primaries = Rcpp::as<std::vector<std::string>>(
  //       setup[static_cast<int>(ExportList::CHEM_PQC_SOLUTION_PRIMARY)]);
  //   this->pqc_exchanger =
  //       Rcpp::List(setup[static_cast<int>(ExportList::CHEM_PQC_EXCHANGER)]);
  //   this->pqc_kinetics =
  //       Rcpp::List(setup[static_cast<int>(ExportList::CHEM_PQC_KINETICS)]);
  //   this->pqc_equilibrium =
  //       Rcpp::List(setup[static_cast<int>(ExportList::CHEM_PQC_EQUILIBRIUM)]);
  //   this->pqc_surface_comps =
  //       Rcpp::List(setup[static_cast<int>(ExportList::CHEM_PQC_SURFACE_COMPS)]);
  //   this->pqc_surface_charges =
  //       Rcpp::List(setup[static_cast<int>(ExportList::CHEM_PQC_SURFACE_CHARGES)]);

  this->dht_species = NamedVector<uint32_t>(
      setup[static_cast<int>(ExportList::CHEM_DHT_SPECIES)]);

  this->interp_species = Rcpp::as<NamedVector<uint32_t>>(
      setup[static_cast<int>(ExportList::CHEM_INTERP_SPECIES)]);

  this->chem_hooks =
      Rcpp::as<Rcpp::List>(setup[static_cast<int>(ExportList::CHEM_HOOKS)]);

  this->ai_surrogate_input_script = Rcpp::as<std::string>(setup[static_cast<int>(ExportList::AI_SURROGATE_INPUT_SCRIPT)]);
}

Rcpp::List InitialList::exportList() {
  Rcpp::List out(static_cast<int>(ExportList::ENUM_SIZE));

  out[static_cast<int>(ExportList::GRID_DIM)] = this->dim;
  out[static_cast<int>(ExportList::GRID_SPECS)] =
      Rcpp::NumericVector::create(this->n_rows, this->n_cols);
  out[static_cast<int>(ExportList::GRID_SPATIAL)] =
      Rcpp::NumericVector::create(this->s_rows, this->s_cols);
  out[static_cast<int>(ExportList::GRID_CONSTANT)] =
      Rcpp::wrap(this->constant_cells);
  out[static_cast<int>(ExportList::GRID_POROSITY)] = Rcpp::wrap(this->porosity);
  out[static_cast<int>(ExportList::GRID_INITIAL)] = this->initial_grid;

  out[static_cast<int>(ExportList::DIFFU_TRANSPORT)] =
      Rcpp::wrap(this->transport_names);
  out[static_cast<int>(ExportList::DIFFU_BOUNDARIES)] = this->boundaries;
  out[static_cast<int>(ExportList::DIFFU_INNER_BOUNDARIES)] =
      this->inner_boundaries;
  out[static_cast<int>(ExportList::DIFFU_ALPHA_X)] = this->alpha_x;
  out[static_cast<int>(ExportList::DIFFU_ALPHA_Y)] = this->alpha_y;

  out[static_cast<int>(ExportList::CHEM_DATABASE)] = Rcpp::wrap(this->database);
  out[static_cast<int>(ExportList::CHEM_PQC_SCRIPT)] =
      Rcpp::wrap(this->pqc_script);
  out[static_cast<int>(ExportList::CHEM_PQC_WITH_H0_O0)] =
      Rcpp::wrap(this->with_h0_o0);
  out[static_cast<int>(ExportList::CHEM_PQC_WITH_REDOX)] =
      Rcpp::wrap(this->with_redox);
  out[static_cast<int>(ExportList::CHEM_FIELD_HEADER)] =
      Rcpp::wrap(this->field_header);
  out[static_cast<int>(ExportList::CHEM_PQC_IDS)] = Rcpp::wrap(this->pqc_ids);
  //   out[static_cast<int>(ExportList::CHEM_PQC_SOLUTIONS)] =
  //       Rcpp::wrap(this->pqc_solutions);
  //   out[static_cast<int>(ExportList::CHEM_PQC_SOLUTION_PRIMARY)] =
  //       Rcpp::wrap(this->pqc_solution_primaries);
  //   out[static_cast<int>(ExportList::CHEM_PQC_EXCHANGER)] =
  //       Rcpp::wrap(this->pqc_exchanger);
  //   out[static_cast<int>(ExportList::CHEM_PQC_KINETICS)] =
  //       Rcpp::wrap(this->pqc_kinetics);
  //   out[static_cast<int>(ExportList::CHEM_PQC_EQUILIBRIUM)] =
  //       Rcpp::wrap(this->pqc_equilibrium);
  //   out[static_cast<int>(ExportList::CHEM_PQC_SURFACE_COMPS)] =
  //       Rcpp::wrap(this->pqc_surface_comps);
  //   out[static_cast<int>(ExportList::CHEM_PQC_SURFACE_CHARGES)] =
  //       Rcpp::wrap(this->pqc_surface_charges);
  out[static_cast<int>(ExportList::CHEM_DHT_SPECIES)] = this->dht_species;
  out[static_cast<int>(ExportList::CHEM_INTERP_SPECIES)] =
      Rcpp::wrap(this->interp_species);
  out[static_cast<int>(ExportList::CHEM_HOOKS)] = this->chem_hooks;
  out[static_cast<int>(ExportList::AI_SURROGATE_INPUT_SCRIPT)] = this->ai_surrogate_input_script;

  return out;
}
} // namespace poet