#include "ChemSim.h"

#include <Rcpp.h>

#include <iostream>

using namespace Rcpp;
using namespace poet;

ChemSim::ChemSim(SimParams &params, RRuntime &R_, Grid &grid_)
    : R(R_), grid(grid_) {
  t_simparams tmp = params.getNumParams();
  this->world_rank = tmp.world_rank;
  this->world_size = tmp.world_size;
  this->wp_size = tmp.wp_size;
  this->out_dir = params.getOutDir();
}

void ChemSim::run() {
  double chem_a, chem_b;

  /* start time measuring */
  chem_a = MPI_Wtime();

  R.parseEvalQ(
      "result <- slave_chemistry(setup=mysetup, data=mysetup$state_T)");
  R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");

  /* end time measuring */
  chem_b = MPI_Wtime();

  chem_t += chem_b - chem_a;
}

void ChemSim::end() {
  R["simtime_chemistry"] = chem_t;
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
}

double ChemSim::getChemistryTime() { return this->chem_t; }
