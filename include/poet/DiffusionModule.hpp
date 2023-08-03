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

#ifndef DIFFUSION_MODULE_H
#define DIFFUSION_MODULE_H

#include "Field.hpp"
#include "SimParams.hpp"
#include "poet/SimParams.hpp"
#include <array>
#include <cmath>
#include <cstdint>
#include <poet/Grid.hpp>
#include <string>
#include <tug/BoundaryCondition.hpp>
#include <tug/Diffusion.hpp>
#include <vector>

namespace poet {
/**
 * @brief Class describing transport simulation
 *
 * Offers simple methods to run an iteration and end the simulation.
 *
 */

constexpr const char *DIFFUSION_MODULE_NAME = "state_T";

class DiffusionModule {
public:
  /**
   * @brief Construct a new TransportSim object
   *
   * The instance will only be initialized with given R object.
   *
   * @param R RRuntime object
   */
  DiffusionModule(const poet::DiffusionParams &diffu_args,
                  const poet::GridParams &grid_params);

  /**
   * @brief Run simulation for one iteration
   *
   * This will simply call the R function 'master_advection'
   *
   */
  void simulate(double dt);

  /**
   * @brief End simulation
   *
   * All measured timings are distributed to the R runtime
   *
   */
  void end();

  /**
   * @brief Get the transport time
   *
   * @return double time spent in transport
   */
  double getTransportTime();

  /**
   * \brief Returns the current diffusion field.
   *
   * \return Reference to the diffusion field.
   */
  Field &getField() { return this->t_field; }

private:
  /**
   * @brief Instance of RRuntime
   *
   */
  // RRuntime &R;

  enum { DIM_1D = 1, DIM_2D };

  void initialize(const poet::DiffusionParams &args);

  uint8_t dim;

  uint32_t prop_count;

  tug::diffusion::TugInput diff_input;
  std::vector<double> alpha;
  std::vector<uint32_t> index_constant_cells;
  std::vector<std::string> prop_names;

  std::vector<tug::bc::BoundaryCondition> bc_vec;
  Field t_field;

  uint32_t n_cells_per_prop;

  /**
   * @brief time spent for transport
   *
   */
  double transport_t = 0.f;
};
} // namespace poet

#endif // DIFFUSION_MODULE_H
