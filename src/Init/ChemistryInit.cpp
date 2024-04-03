#include "InitialList.hpp"

#include <Rcpp.h>
#include <vector>

namespace poet {

void InitialList::initChemistry(const Rcpp::List &chem) {
  this->pqc_sol_order = this->transport_names;

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
}

InitialList::ChemistryInit InitialList::getChemistryInit() const {
  ChemistryInit chem_init;

  chem_init.total_grid_cells = this->n_cols * this->n_rows;

  chem_init.database = database;
  chem_init.pqc_scripts = pqc_scripts;
  chem_init.pqc_ids = pqc_ids;

  chem_init.pqc_sol_order = pqc_sol_order;

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