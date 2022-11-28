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

#ifndef GRID_H
#define GRID_H

#include "PhreeqcRM.h"
#include "poet/SimParams.hpp"
#include <Rcpp.h>
#include <array>
#include <cstdint>
#include <map>
#include <string>
#include <vector>

#define MAX_DIM 2

namespace poet {

enum { GRID_X_DIR, GRID_Y_DIR, GRID_Z_DIR };

using StateMemory = struct s_StateMemory {
  std::vector<double> mem;
  std::vector<std::string> props;
};

constexpr const char *GRID_MODULE_NAME = "grid_init";

/**
 * @brief Class describing the grid
 *
 * Providing methods to shuffle and unshuffle grid (for the master) as also to
 * import and export a work package (for worker).
 *
 * @todo find better abstraction
 *
 */
class Grid {
private:
  using prop_vec = std::vector<std::string>;

public:
  /**
   * @brief Construct a new Grid object
   *
   * This will call the default constructor and initializes private RRuntime
   * with given R runtime.
   *
   * @param R
   */
  Grid(poet::GridParams grid_args);

  ~Grid();
  /**
   * @brief Init the grid
   *
   * At this moment init will only declare and define a variable inside the R
   * runtime called grid_tmp since the whole Grid initialization and management
   * is done by the R runtime. This may change in the future.
   *
   */
  void init_from_R();

  auto getGridDimension() -> uint8_t;
  auto getTotalCellCount() -> uint32_t;
  auto getGridCellsCount(uint8_t direction) -> uint32_t;
  auto getDomainSize(uint8_t dircetion) -> uint32_t;

  StateMemory *registerState(std::string module_name,
                             std::vector<std::string> props);
  StateMemory *getStatePointer(std::string module_name);
  void deregisterState(std::string module_name);

  auto getSpeciesCount() -> uint32_t;
  auto getPropNames() -> prop_vec;

  auto getSpeciesByName(std::string name,
                        std::string module_name = poet::GRID_MODULE_NAME)
      -> std::vector<double>;

private:
  PhreeqcRM *initPhreeqc(const poet::GridParams &params);

  std::uint8_t dim;
  std::array<std::uint32_t, MAX_DIM> d_spatial;
  std::array<std::uint32_t, MAX_DIM> n_cells;

  std::map<std::string, StateMemory *> state_register;
  std::map<std::string, std::vector<double>> grid_init;

  prop_vec prop_names;

  PhreeqcRM *phreeqc_rm = std::nullptr_t();

  void initFromPhreeqc(const poet::GridParams &grid_args);
  void initFromRDS(const poet::GridParams &grid_args);
  void initFromScratch(const poet::GridParams &grid_args);
};
} // namespace poet
#endif // GRID_H
