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
#include "poet/SimParams.hpp"
#include <poet/ChemSimSeq.hpp>
#include <poet/Grid.hpp>

#include <Rcpp.h>

#include <algorithm>
#include <bits/stdint-uintn.h>
#include <iostream>
#include <string>
#include <vector>

using namespace Rcpp;
using namespace poet;

ChemSeq::ChemSeq(SimParams &params, RInside &R_, Grid &grid_)
    : BaseChemModule(params, R_, grid_) {

  this->state = this->grid.RegisterState(
      poet::BaseChemModule::CHEMISTRY_MODULE_NAME, this->prop_names);
  std::vector<double> &field = this->state->mem;

  field.resize(this->n_cells_per_prop * this->prop_names.size());
  for (uint32_t i = 0; i < this->prop_names.size(); i++) {
    std::vector<double> prop_vec =
        this->grid.GetSpeciesByName(this->prop_names[i]);

    std::copy(prop_vec.begin(), prop_vec.end(),
              field.begin() + (i * this->n_cells_per_prop));
  }
}

poet::ChemSeq::~ChemSeq() {
  if (this->phreeqc_rm) {
    delete this->phreeqc_rm;
  }
}

void poet::ChemSeq::InitModule(poet::ChemistryParams &chem_params) {
  this->phreeqc_rm = new PhreeqcWrapper(this->n_cells_per_prop);

  this->phreeqc_rm->SetupAndLoadDB(chem_params);
  this->phreeqc_rm->InitFromFile(chem_params.input_script);
}

void ChemSeq::Simulate(double dt) {
  double chem_a, chem_b;

  /* start time measuring */
  chem_a = MPI_Wtime();

  std::vector<double> &field = this->state->mem;

  // HACK: transfer the field into R data structure serving as input for phreeqc
  R["TMP_T"] = field;

  R.parseEvalQ("mysetup$state_T <- setNames(data.frame(matrix(TMP_T, "
               "ncol=length(mysetup$grid$props), nrow=" +
               std::to_string(this->n_cells_per_prop) +
               ")), mysetup$grid$props)");

  this->phreeqc_rm->SetInternalsFromWP(field, this->n_cells_per_prop);
  this->phreeqc_rm->SetTime(0);
  this->phreeqc_rm->SetTimeStep(dt);
  this->phreeqc_rm->RunCells();

  // HACK: we will copy resulting field into global grid field. Maybe we will
  // find a more performant way here ...
  std::vector<double> vecSimResult;
  this->phreeqc_rm->GetWPFromInternals(vecSimResult, this->n_cells_per_prop);
  std::copy(vecSimResult.begin(), vecSimResult.end(), field.begin());

  R["TMP_C"] = field;

  R.parseEvalQ("mysetup$state_C <- setNames(data.frame(matrix(TMP_C, "
               "ncol=length(mysetup$grid$props), nrow=" +
               std::to_string(this->n_cells_per_prop) +
               ")), mysetup$grid$props)");

  /* end time measuring */
  chem_b = MPI_Wtime();

  chem_t += chem_b - chem_a;
}

void ChemSeq::End() {
  R["simtime_chemistry"] = this->chem_t;
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
}

double ChemSeq::getChemistryTime() { return this->chem_t; }
