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

#include "DataStructures/Field.hpp"
#include "Init/InitialList.hpp"

#include <sys/types.h>

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
  DiffusionModule(const InitialList::DiffusionInit &init_list, Field field)
      : param_list(init_list), transport_field(field){};

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
  double getTransportTime();

  /**
   * \brief Returns the current diffusion field.
   *
   * \return Reference to the diffusion field.
   */
  Field &getField() { return this->transport_field; }

private:
  /**
   * @brief Instance of RRuntime
   *
   */

  InitialList::DiffusionInit param_list;

  Field transport_field;

  /**
   * @brief time spent for transport
   *
   */
  double transport_t = 0.;
};
} // namespace poet

#endif // DIFFUSION_MODULE_H
