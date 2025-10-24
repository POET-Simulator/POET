/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
**
** Copyright (C) 2023-2024 Marco De Lucia (GFZ Potsdam), Max Luebke (University
** of Potsdam)
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
