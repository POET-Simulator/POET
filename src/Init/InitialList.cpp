#include "InitialList.hpp"
#include "DataStructures/NamedVector.hpp"
#include <Rcpp/internal/wrap.h>
#include <Rcpp/iostream/Rstreambuf.h>
#include <Rcpp/proxy/ProtectedProxy.h>
#include <Rcpp/vector/instantiation.h>
#include <cstdint>
#include <string>
#include <vector>

namespace poet {
void InitialList::initializeFromList(const Rcpp::List &setup) {
  prepareGrid(setup[grid_key]);
  initDiffusion(setup[diffusion_key]);
  initChemistry(setup[chemistry_key]);
}

void InitialList::importList(const Rcpp::List &setup) {
  this->dim = Rcpp::as<uint8_t>(setup[static_cast<int>(ExportList::GRID_DIM)]);

  Rcpp::NumericVector grid_specs =
      setup[static_cast<int>(ExportList::GRID_SPECS)];
  this->n_rows = grid_specs[0];
  this->n_cols = grid_specs[1];

  Rcpp::NumericVector spatial =
      setup[static_cast<int>(ExportList::GRID_SPATIAL)];
  this->s_rows = spatial[0];
  this->s_cols = spatial[1];

  this->constant_cells = Rcpp::as<std::vector<uint32_t>>(
      setup[static_cast<int>(ExportList::GRID_CONSTANT)]);

  this->porosity = Rcpp::as<std::vector<double>>(
      setup[static_cast<int>(ExportList::GRID_POROSITY)]);

  this->initial_grid =
      Rcpp::List(setup[static_cast<int>(ExportList::GRID_INITIAL)]);

  this->transport_names = Rcpp::as<std::vector<std::string>>(
      setup[static_cast<int>(ExportList::DIFFU_TRANSPORT)]);
  this->boundaries =
      Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_BOUNDARIES)]);
  this->alpha_x =
      Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_ALPHA_X)]);
  this->alpha_y =
      Rcpp::List(setup[static_cast<int>(ExportList::DIFFU_ALPHA_Y)]);

  this->database =
      Rcpp::as<std::string>(setup[static_cast<int>(ExportList::CHEM_DATABASE)]);
  this->pqc_scripts = Rcpp::as<std::vector<std::string>>(
      setup[static_cast<int>(ExportList::CHEM_PQC_SCRIPTS)]);
  this->pqc_ids = Rcpp::as<std::vector<int>>(
      setup[static_cast<int>(ExportList::CHEM_PQC_IDS)]);
  this->pqc_sol_order = Rcpp::as<std::vector<std::string>>(
      setup[static_cast<int>(ExportList::CHEM_PQC_SOL_ORDER)]);

  this->dht_species = NamedVector<uint32_t>(
      setup[static_cast<int>(ExportList::CHEM_DHT_SPECIES)]);

  this->interp_species = Rcpp::as<NamedVector<uint32_t>>(
      setup[static_cast<int>(ExportList::CHEM_INTERP_SPECIES)]);

  this->chem_hooks =
      Rcpp::as<Rcpp::List>(setup[static_cast<int>(ExportList::CHEM_HOOKS)]);
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
  out[static_cast<int>(ExportList::DIFFU_ALPHA_X)] = this->alpha_x;
  out[static_cast<int>(ExportList::DIFFU_ALPHA_Y)] = this->alpha_y;

  out[static_cast<int>(ExportList::CHEM_DATABASE)] = Rcpp::wrap(this->database);
  out[static_cast<int>(ExportList::CHEM_PQC_SCRIPTS)] =
      Rcpp::wrap(this->pqc_scripts);
  out[static_cast<int>(ExportList::CHEM_PQC_IDS)] = Rcpp::wrap(this->pqc_ids);
  out[static_cast<int>(ExportList::CHEM_PQC_SOL_ORDER)] =
      Rcpp::wrap(this->pqc_sol_order);
  out[static_cast<int>(ExportList::CHEM_DHT_SPECIES)] = this->dht_species;
  out[static_cast<int>(ExportList::CHEM_INTERP_SPECIES)] =
      Rcpp::wrap(this->interp_species);
  out[static_cast<int>(ExportList::CHEM_HOOKS)] = this->chem_hooks;

  return out;
}
} // namespace poet