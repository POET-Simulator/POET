#include "InitialList.hpp"

namespace poet {

void InitialList::initChemistry(const Rcpp::List &chem) {
  this->dht_defined = chem.containsElementNamed("dht_species");

  if (this->dht_defined) {
    this->dht_species = Rcpp::as<NamedVector<uint32_t>>(chem["dht_species"]);
  }
}

InitialList::ChemistryInit InitialList::getChemistryInit() const {
  ChemistryInit chem_init;

  chem_init.total_grid_cells = this->n_cols * this->n_rows;

  chem_init.database = database;
  chem_init.pqc_scripts = pqc_scripts;
  chem_init.pqc_ids = pqc_ids;

  chem_init.pqc_sol_order = pqc_sol_order;

  chem_init.dht_defined = dht_defined;
  chem_init.dht_species = dht_species;

  return chem_init;
}
} // namespace poet