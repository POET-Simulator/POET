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

#ifndef TRANSPORT_SIM_H
#define TRANSPORT_SIM_H

#include <RRuntime.h>

namespace poet {
/**
 * @brief Class describing transport simulation
 *
 * Offers simple methods to run an iteration and end the simulation.
 *
 */
class TransportSim {
 public:
  /**
   * @brief Construct a new TransportSim object
   *
   * The instance will only be initialized with given R object.
   *
   * @param R RRuntime object
   */
  TransportSim(RRuntime &R);

  /**
   * @brief Run simulation for one iteration
   *
   * This will simply call the R function 'master_advection'
   *
   */
  void run();

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

 private:
  /**
   * @brief Instance of RRuntime
   *
   */
  RRuntime &R;

  /**
   * @brief time spent for transport
   *
   */
  double transport_t = 0.f;
};
}  // namespace poet

#endif  // TRANSPORT_SIM_H
