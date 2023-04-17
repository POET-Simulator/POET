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

#include "poet/SimParams.hpp"
#include "tug/BoundaryCondition.hpp"
#include "tug/Diffusion.hpp"
#include <Rcpp.h>
#include <algorithm>
#include <cstdint>
#include <poet/DiffusionModule.hpp>
#include <poet/Grid.hpp>

#include <array>
#include <cassert>
#include <mpi.h>
#include <string>
#include <vector>

using namespace poet;

static constexpr double ZERO_MULTIPLICATOR = 10E-14;

constexpr std::array<uint8_t, 4> borders = {
    tug::bc::BC_SIDE_LEFT, tug::bc::BC_SIDE_RIGHT, tug::bc::BC_SIDE_TOP,
    tug::bc::BC_SIDE_BOTTOM};

inline const char *convert_bc_to_config_file(uint8_t in) {
  switch (in) {
  case tug::bc::BC_SIDE_TOP:
    return "N";
  case tug::bc::BC_SIDE_RIGHT:
    return "E";
  case tug::bc::BC_SIDE_BOTTOM:
    return "S";
  case tug::bc::BC_SIDE_LEFT:
    return "W";
  }
  return "";
}

DiffusionModule::DiffusionModule(poet::DiffusionParams diffu_args, Grid &grid_)
    : grid(grid_) {
  this->diff_input.setGridCellN(grid_.GetGridCellsCount(GRID_X_DIR),
                                grid_.GetGridCellsCount(GRID_Y_DIR));
  this->diff_input.setDomainSize(grid_.GetGridSize(GRID_X_DIR),
                                 grid_.GetGridSize(GRID_Y_DIR));

  this->dim = grid_.GetGridDimension();

  this->n_cells_per_prop = grid_.GetTotalCellCount();

  this->initialize(diffu_args);
}

void DiffusionModule::initialize(poet::DiffusionParams args) {
  // const poet::DiffusionParams args(this->R);

  // name of props
  //
  this->prop_names = Rcpp::as<std::vector<std::string>>(args.initial_t.names());
  this->prop_count = this->prop_names.size();

  this->state =
      this->grid.RegisterState(DIFFUSION_MODULE_NAME, this->prop_names);

  std::vector<double> &field = this->state->mem;

  // initialize field and alphas
  field.resize(this->n_cells_per_prop * this->prop_count);
  this->alpha.reserve(this->prop_count);

  for (uint32_t i = 0; i < this->prop_count; i++) {
    // get alphas - we cannot assume alphas are provided in same order as
    // initial input
    this->alpha.push_back(args.alpha[this->prop_names[i]]);

    double val = args.initial_t[prop_names[i]];

    std::vector<double> prop_vec(n_cells_per_prop, val);
    std::copy(prop_vec.begin(), prop_vec.end(),
              field.begin() + (i * this->n_cells_per_prop));
    if (this->dim == this->DIM_2D) {
      tug::bc::BoundaryCondition bc(this->grid.GetGridCellsCount(GRID_X_DIR),
                                    this->grid.GetGridCellsCount(GRID_Y_DIR));
      this->bc_vec.push_back(bc);
    } else {
      tug::bc::BoundaryCondition bc(this->grid.GetGridCellsCount(GRID_X_DIR));
      this->bc_vec.push_back(bc);
    }
  }

  // apply boundary conditions to each ghost node
  uint8_t bc_size = (this->dim == this->DIM_1D ? 2 : 4);
  for (uint8_t i = 0; i < bc_size; i++) {
    const auto &side = borders[i];
    std::vector<uint32_t> vecinj_i = Rcpp::as<std::vector<uint32_t>>(
        args.vecinj_index[convert_bc_to_config_file(side)]);
    for (int i = 0; i < this->prop_count; i++) {
      std::vector<double> bc_vec = args.vecinj[this->prop_names[i]];
      tug::bc::BoundaryCondition &curr_bc = *(this->bc_vec.begin() + i);
      for (int j = 0; j < vecinj_i.size(); j++) {
        if (vecinj_i[j] == 0) {
          continue;
        }
        if (this->dim == this->DIM_2D) {
          curr_bc(side, j) = {tug::bc::BC_TYPE_CONSTANT,
                              bc_vec[vecinj_i[j] - 1]};
        }
        if (this->dim == this->DIM_1D) {
          curr_bc(side) = {tug::bc::BC_TYPE_CONSTANT, bc_vec[vecinj_i[j] - 1]};
        }
      }
    }
  }

  // apply inner grid constant cells
  // NOTE: opening a scope here for distinguish variable names
  if (args.vecinj_inner.size() != 0) {
    // get indices of constant grid cells
    // Rcpp::NumericVector indices_const_cells = args.vecinj_inner(Rcpp::_, 0);
    // this->index_constant_cells =
    //     Rcpp::as<std::vector<uint32_t>>(indices_const_cells);

    // // get indices to vecinj for constant cells
    // Rcpp::NumericVector vecinj_indices = args.vecinj_inner(Rcpp::_, 1);

    // apply inner constant cells for every concentration
    for (int i = 0; i < this->prop_count; i++) {
      std::vector<double> bc_vec = args.vecinj[this->prop_names[i]];
      tug::bc::BoundaryCondition &curr_bc = *(this->bc_vec.begin() + i);
      for (int j = 0; j < args.vecinj_inner.size(); j++) {
        std::vector<double> inner_tuple =
            Rcpp::as<std::vector<double>>(args.vecinj_inner[j]);
        tug::bc::boundary_condition bc = {tug::bc::BC_TYPE_CONSTANT,
                                          bc_vec[inner_tuple[0] - 1]};

        this->index_constant_cells.push_back(inner_tuple[1]);
        uint32_t x = inner_tuple[1];
        uint32_t y = (this->dim == this->DIM_1D ? 0 : inner_tuple[2]);

        curr_bc.setInnerBC(bc, x, y);
      }
    }
  }
}

void DiffusionModule::simulate(double dt) {
  double sim_a_transport, sim_b_transport;

  sim_b_transport = MPI_Wtime();

  std::cout << "DiffusionModule::simulate(): Starting diffusion ..."
            << std::flush;

  std::vector<double> &curr_field = this->state->mem;

  this->diff_input.setTimestep(dt);

  double *field = curr_field.data();
  for (int i = 0; i < this->prop_count; i++) {
    double *in_field = &field[i * this->n_cells_per_prop];
    std::vector<double> in_alpha(this->n_cells_per_prop, this->alpha[i]);
    this->diff_input.setBoundaryCondition(this->bc_vec[i]);

    if (this->dim == this->DIM_1D) {
      tug::diffusion::BTCS_1D(this->diff_input, in_field, in_alpha.data());
    } else {
      tug::diffusion::ADI_2D(this->diff_input, in_field, in_alpha.data());
    }
  }

  std::cout << " done!\n";

  sim_a_transport = MPI_Wtime();

  transport_t += sim_a_transport - sim_b_transport;
}

void DiffusionModule::end() {
  // R["simtime_transport"] = transport_t;
  // R.parseEvalQ("profiling$simtime_transport <- simtime_transport");
}

double DiffusionModule::getTransportTime() { return this->transport_t; }
