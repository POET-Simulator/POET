/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2021 Marco De Lucia (GFZ Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include <poet/ChemSim.hpp>

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
