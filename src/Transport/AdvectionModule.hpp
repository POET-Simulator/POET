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

#ifndef ADVECTION_MODULE_H
#define ADVECTION_MODULE_H

#include "../Base/RInsidePOET.hpp"
#include "../DataStructures/DataStructures.hpp"

#include <Rcpp.h>
#include <array>
#include <cstdint>
#include <map>
#include <utility>
#include <vector>

namespace poet {
/**
 * @brief Class describing transport simulation using advection
 *
 * Quick and dirty implementation, assuming homogenous and constant porisity,
 * pressure and temperature. Open boundary conditions are assumed.
 *
 */

class AdvectionModule {
public:
  /**
   * @brief Construct a new TransportSim object
   *
   * The instance will only be initialized with given R object.
   *
   * @param R RRuntime object
   */
  AdvectionModule(RInsidePOET &R) { initializeParams(R); }

  /**
   * @brief Run simulation for one iteration
   *
   * This will simply call the R function 'master_advection'
   *
   */
  void simulate(double dt);

  /**
   * @brief Get the transport time
   *
   * @return double time spent in transport
   */
  double getTransportTime() const { return this->transport_t; }

  /**
   * \brief Returns the current diffusion field.
   *
   * \return Reference to the diffusion field.
   */
  const Field &getField() const { return this->t_field; }

private:
  void initializeParams(RInsidePOET &R);

  std::array<std::uint32_t, 2> n_cells;
  std::set<std::uint32_t> inactive_cells;

  Field t_field;

  /**
   * @brief time spent for transport
   *
   */
  double transport_t = 0.f;
};
} // namespace poet

#endif // ADVECTION_MODULE_H
