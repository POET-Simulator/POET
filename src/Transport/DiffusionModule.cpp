// SPDX-License-Identifier: GPL-2.0-or-later
/*
 * DiffusionModule.cpp - Coupling module betweeen POET and tug
 * Copyright (C) 2018-2025 Max Luebke (University of Potsdam)
 */

#include "DiffusionModule.hpp"

#include "Base/Macros.hpp"
#include "Init/InitialList.hpp"

#include <chrono>
#include <string>
#include <tug/Boundary.hpp>
#include <tug/Diffusion.hpp>

using namespace poet;

void DiffusionModule::simulate(double requested_dt) {
  // MSG("Starting diffusion ...");
  const auto start_diffusion_t = std::chrono::high_resolution_clock::now();

  const auto &n_rows = this->param_list.n_rows;
  const auto &n_cols = this->param_list.n_cols;

  for (const auto &sol_name : this->param_list.transport_names) {
#if defined(POET_TUG_BTCS)
    tug::Diffusion<TugType> diffusion_solver(
        this->transport_field[sol_name].data(), n_rows, n_cols);
#elif defined(POET_TUG_FTCS)
    tug::Diffusion<TugType, tug::FTCS_APPROACH> diffusion_solver(
        this->transport_field[sol_name].data(), n_rows, n_cols);
#else
#error "No valid approach defined"
#endif

    tug::RowMajMatMap<TugType> alpha_x_map(
        this->param_list.alpha_x[sol_name].data(), n_rows, n_cols);
    tug::RowMajMatMap<TugType> alpha_y_map(
        this->param_list.alpha_y[sol_name].data(), n_rows, n_cols);

    diffusion_solver.setAlphaX(alpha_x_map);
    diffusion_solver.setAlphaY(alpha_y_map);

    auto &boundary = diffusion_solver.getBoundaryConditions();

    boundary.deserialize(this->param_list.boundaries[sol_name]);

    if (!this->param_list.inner_boundaries[sol_name].empty()) {
      boundary.setInnerBoundaries(this->param_list.inner_boundaries[sol_name]);
    }

    diffusion_solver.setTimestep(requested_dt);
    diffusion_solver.run();
  }

  const auto end_diffusion_t = std::chrono::high_resolution_clock::now();

  const auto consumed_time_seconds =
      std::chrono::duration_cast<std::chrono::duration<double>>(
          end_diffusion_t - start_diffusion_t);

  MSG("Diffusion done in " + std::to_string(consumed_time_seconds.count()) +
      " [s]");

  transport_t += consumed_time_seconds.count();
}

double DiffusionModule::getTransportTime() { return this->transport_t; }
