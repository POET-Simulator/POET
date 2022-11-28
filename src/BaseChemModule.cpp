#include <poet/ChemSimSeq.hpp>

poet::BaseChemModule::BaseChemModule(SimParams &params, RInside &R_,
                                     Grid &grid_)
    : R(R_), grid(grid_) {
  t_simparams tmp = params.getNumParams();
  this->world_rank = tmp.world_rank;
  this->world_size = tmp.world_size;
  this->prop_names = this->grid.getPropNames();

  this->n_cells_per_prop = this->grid.getTotalCellCount();
}
