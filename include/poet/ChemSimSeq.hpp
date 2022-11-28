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

#ifndef CHEMSIMSEQ_H
#define CHEMSIMSEQ_H

#include "DHT_Wrapper.hpp"
#include "Grid.hpp"
#include "PhreeqcWrapper.hpp"
#include "RInside.h"
#include "SimParams.hpp"
#include <bits/stdint-uintn.h>
#include <cstddef>
#include <cstdint>
#include <mpi.h>

#include <string>
#include <vector>

/** Number of data elements that are kept free at each work package */
#define BUFFER_OFFSET 5

/** Message tag indicating work */
#define TAG_WORK 42
/** Message tag indicating to finish loop */
#define TAG_FINISH 43
/** Message tag indicating timing profiling */
#define TAG_TIMING 44
/** Message tag indicating collecting DHT performance */
#define TAG_DHT_PERF 45
/** Message tag indicating simulation reached the end of an itertation */
#define TAG_DHT_ITER 47

namespace poet {
/**
 * @brief Base class of the chemical simulation
 *
 * Providing member functions to run an iteration and to end a simulation. Also
 * containing basic parameters for simulation.
 *
 */
class BaseChemModule {
public:
  BaseChemModule(SimParams &params, RInside &R_, Grid &grid_);

  virtual void InitModule(poet::ChemistryParams &chem_params){};

  static constexpr const char *CHEMISTRY_MODULE_NAME = "state_c";

protected:
  int world_rank;
  int world_size;
  RInside &R;

  Grid &grid;
  double chem_t = 0.f;
  double current_sim_time = 0;

  uint32_t n_cells_per_prop;
  std::vector<std::string> prop_names;
};

class ChemSeq : public BaseChemModule {

public:
  /**
   * @brief Construct a new ChemSim object
   *
   * Creating a new instance of class ChemSim will just extract simulation
   * parameters from SimParams object.
   *
   * @param params SimParams object
   * @param R_ R runtime
   * @param grid_ Initialized grid
   */
  ChemSeq(SimParams &params, RInside &R_, Grid &grid_);

  ~ChemSeq();

  void InitModule(poet::ChemistryParams &chem_params);
  /**
   * @brief Run iteration of simulation in sequential mode
   *
   * This will call the correspondending R function slave_chemistry, followed by
   * the execution of master_chemistry.
   *
   * @todo change function name. Maybe 'slave' to 'seq'.
   *
   */
  void Simulate(double dt);

  /**
   * @brief End simulation
   *
   * End the simulation by distribute the measured runtime of simulation to the
   * R runtime.
   *
   */
  void End();

  /**
   * @brief Get the Chemistry Time
   *
   * @return double Runtime of sequential chemistry simulation in seconds
   */
  double getChemistryTime();

private:
  poet::StateMemory *state;
  PhreeqcWrapper *phreeqc_rm = std::nullptr_t();
};

} // namespace poet
#endif // CHEMSIMSEQ_H
