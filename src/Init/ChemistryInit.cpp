#include "InitialList.hpp"

#include <Rcpp.h>
#include <cstddef>
#include <vector>

namespace poet {

void InitialList::initChemistry(const Rcpp::List &chem) {
  this->pqc_solutions = std::vector<std::string>(
      this->transport_names.begin() + 3, this->transport_names.end());

  this->pqc_solution_primaries = this->phreeqc->getSolutionPrimaries();

  if (chem.containsElementNamed("dht_species")) {
    this->dht_species = Rcpp::as<NamedVector<uint32_t>>(chem["dht_species"]);
  }

  if (chem.containsElementNamed("pht_species")) {
    this->interp_species = Rcpp::as<NamedVector<uint32_t>>(chem["pht_species"]);
  }

  if (chem.containsElementNamed("hooks")) {
    this->chem_hooks = Rcpp::List(chem["hooks"]);

    std::vector<std::string> hook_names = this->chem_hooks.names();

    for (const auto &name : hook_names) {
      if (this->hook_name_list.find(name) == this->hook_name_list.end()) {
        Rcpp::Rcerr << "Unknown chemistry hook: " << name << std::endl;
      }
    }
  }

  this->field_header =
      Rcpp::as<std::vector<std::string>>(this->initial_grid.names());
  this->field_header.erase(this->field_header.begin());
}

InitialList::ChemistryInit InitialList::getChemistryInit() const {
  ChemistryInit chem_init;

  chem_init.total_grid_cells = this->n_cols * this->n_rows;

  // chem_init.field_header = this->field_header;

  chem_init.database = database;
  // chem_init.pqc_scripts = pqc_scripts;
  // chem_init.pqc_ids = pqc_ids;

  for (std::size_t i = 0; i < pqc_scripts.size(); i++) {
    POETInitCell cell = {
        pqc_solutions,
        pqc_solution_primaries,
        Rcpp::as<std::vector<std::string>>(pqc_exchanger[i]),
        Rcpp::as<std::vector<std::string>>(pqc_kinetics[i]),
        Rcpp::as<std::vector<std::string>>(pqc_equilibrium[i]),
        Rcpp::as<std::vector<std::string>>(pqc_surface_comps[i]),
        Rcpp::as<std::vector<std::string>>(pqc_surface_charges[i])};

    chem_init.pqc_config[pqc_ids[i]] = {database, pqc_scripts[i], cell};
  }

  // chem_init.pqc_sol_order = pqc_solutions;

  chem_init.dht_species = dht_species;
  chem_init.interp_species = interp_species;

  if (this->chem_hooks.size() > 0) {
    if (this->chem_hooks.containsElementNamed("dht_fill")) {
      chem_init.hooks.dht_fill = this->chem_hooks["dht_fill"];
    }
    if (this->chem_hooks.containsElementNamed("dht_fuzz")) {
      chem_init.hooks.dht_fuzz = this->chem_hooks["dht_fuzz"];
    }
    if (this->chem_hooks.containsElementNamed("interp_pre")) {
      chem_init.hooks.interp_pre = this->chem_hooks["interp_pre"];
    }
    if (this->chem_hooks.containsElementNamed("interp_post")) {
      chem_init.hooks.interp_post = this->chem_hooks["interp_post"];
    }
  }

  return chem_init;
}
} // namespace poet