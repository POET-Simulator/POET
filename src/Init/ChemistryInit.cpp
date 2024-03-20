#include "InitialList.hpp"

namespace poet {
InitialList::ChemistryInit InitialList::getChemistryInit() const {
  ChemistryInit chem_init;

  chem_init.initial_grid = Field(initial_grid);

  chem_init.database = database;
  chem_init.pqc_scripts = pqc_scripts;
  chem_init.pqc_ids = pqc_ids;

  chem_init.pqc_sol_order = pqc_sol_order;

  return chem_init;
}
} // namespace poet