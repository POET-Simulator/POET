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

#include <Eigen/src/Core/Map.h>
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <chrono>

#include "tug/Boundary.hpp"
#include "tug/Grid.hpp"
#include "tug/Simulation.hpp"

#include <string>
#include <vector>

using namespace poet;

static inline std::vector<TugType>
MatrixToVec(const Eigen::MatrixX<TugType> &mat) {
  std::vector<TugType> vec(mat.rows() * mat.cols());

  for (std::uint32_t i = 0; i < mat.cols(); i++) {
    for (std::uint32_t j = 0; j < mat.rows(); j++) {
      vec[j * mat.cols() + i] = mat(j, i);
    }
  }

  return vec;
}

static inline Eigen::MatrixX<TugType>
VecToMatrix(const std::vector<TugType> &vec, std::uint32_t n_rows,
            std::uint32_t n_cols) {
  Eigen::MatrixX<TugType> mat(n_rows, n_cols);

  for (std::uint32_t i = 0; i < n_cols; i++) {
    for (std::uint32_t j = 0; j < n_rows; j++) {
      mat(j, i) = vec[j * n_cols + i];
    }
  }

  return mat;
}

// static constexpr double ZERO_MULTIPLICATOR = 10E-14;

void DiffusionModule::simulate(double requested_dt) {
  // MSG("Starting diffusion ...");
  const auto start_diffusion_t = std::chrono::high_resolution_clock::now();

  const auto &n_rows = this->param_list.n_rows;
  const auto &n_cols = this->param_list.n_cols;

  tug::Grid<TugType> grid(param_list.n_rows, param_list.n_cols);
  tug::Boundary boundary(grid);

  grid.setDomain(param_list.s_rows, param_list.s_cols);

  tug::Simulation<TugType> sim(grid, boundary);
  sim.setIterations(1);

  for (const auto &sol_name : this->param_list.transport_names) {
    auto &species_conc = this->transport_field[sol_name];

    Eigen::MatrixX<TugType> conc = VecToMatrix(species_conc, n_rows, n_cols);
    Eigen::MatrixX<TugType> alpha_x =
        VecToMatrix(this->param_list.alpha_x[sol_name], n_rows, n_cols);
    Eigen::MatrixX<TugType> alpha_y =
        VecToMatrix(this->param_list.alpha_y[sol_name], n_rows, n_cols);

    boundary.deserialize(this->param_list.boundaries[sol_name]);

    if (!this->param_list.inner_boundaries[sol_name].empty()) {
      boundary.setInnerBoundaries(this->param_list.inner_boundaries[sol_name]);
    }

    grid.setAlpha(alpha_x, alpha_y);
    grid.setConcentrations(conc);

    sim.setTimestep(requested_dt);

    sim.run();

    species_conc = MatrixToVec(grid.getConcentrations());
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
