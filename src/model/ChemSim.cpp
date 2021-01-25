#include "ChemSim.h"

#include <Rcpp.h>

#include <iostream>

using namespace Rcpp;
using namespace poet;

ChemSim::ChemSim(t_simparams *params, RRuntime &R_, Grid &grid_)
    : R(R_), grid(grid_) {
  this->world_rank = params->world_rank;
  this->world_size = params->world_size;
  this->wp_size = params->wp_size;
  this->out_dir = params->out_dir;
}

void ChemSim::run() {
  double chem_a, chem_b;

  chem_a = MPI_Wtime();

  R.parseEvalQ(
      "result <- slave_chemistry(setup=mysetup, data=mysetup$state_T)");
  R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");

  chem_b = MPI_Wtime();
  chem_t += chem_b - chem_a;
}

void ChemSim::end() {
  R["simtime_chemistry"] = chem_t;
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
}

double ChemSim::getChemistryTime() { return this->chem_t; }
