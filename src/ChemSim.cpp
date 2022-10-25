/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
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

#include "poet/DiffusionModule.hpp"
#include <poet/ChemSim.hpp>
#include <poet/Grid.hpp>

#include <Rcpp.h>

#include <algorithm>
#include <bits/stdint-uintn.h>
#include <iostream>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace poet;

ChemSim::ChemSim(SimParams &params, RRuntime &R_, Grid &grid_)
    : R(R_), grid(grid_) {
  t_simparams tmp = params.getNumParams();
  this->world_rank = tmp.world_rank;
  this->world_size = tmp.world_size;
  this->wp_size = tmp.wp_size;
  this->out_dir = params.getOutDir();

  this->prop_names = this->grid.getPropNames();

  this->n_cells_per_prop = this->grid.getTotalCellCount();
  this->state =
      this->grid.registerState(poet::CHEMISTRY_MODULE_NAME, this->prop_names);
  auto &field = this->state->mem;

  field.resize(this->n_cells_per_prop * this->prop_names.size());
  for (uint32_t i = 0; i < this->prop_names.size(); i++) {
    std::vector<double> prop_vec =
        this->grid.getSpeciesByName(this->prop_names[i]);

    std::copy(prop_vec.begin(), prop_vec.end(),
              field.begin() + (i * this->n_cells_per_prop));
  }
}

void ChemSim::run() {
  double chem_a, chem_b;

  /* start time measuring */
  chem_a = MPI_Wtime();

  std::vector<double> &field = this->state->mem;

  for (uint32_t i = 0; i < this->prop_names.size(); i++) {
    try {
      std::vector<double> t_prop_vec = this->grid.getSpeciesByName(
          this->prop_names[i], poet::DIFFUSION_MODULE_NAME);

      std::copy(t_prop_vec.begin(), t_prop_vec.end(),
                field.begin() + (i * this->n_cells_per_prop));
    } catch (...) {
      continue;
    }
  }

  // HACK: transfer the field into R data structure serving as input for phreeqc
  R["TMP_T"] = field;

  R.parseEvalQ("mysetup$state_T <- setNames(data.frame(matrix(TMP_T, "
               "ncol=length(mysetup$grid$props), nrow=" +
               std::to_string(this->n_cells_per_prop) +
               ")), mysetup$grid$props)");

  R.parseEvalQ(
      "result <- slave_chemistry(setup=mysetup, data=mysetup$state_T)");
  R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");

  // HACK: copy R data structure back to C++ field
  Rcpp::DataFrame result = R.parseEval("mysetup$state_C");

  for (uint32_t i = 0; i < this->prop_names.size(); i++) {
    std::vector<double> c_prop =
        Rcpp::as<std::vector<double>>(result[this->prop_names[i].c_str()]);

    std::copy(c_prop.begin(), c_prop.end(),
              field.begin() + (i * this->n_cells_per_prop));
  }

  /* end time measuring */
  chem_b = MPI_Wtime();

  chem_t += chem_b - chem_a;
}

void ChemSim::end() {
  R["simtime_chemistry"] = chem_t;
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
}

double ChemSim::getChemistryTime() { return this->chem_t; }
