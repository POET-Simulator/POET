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

#include "RRuntime.hpp"
#include "poet/SimParams.hpp"
#include <Rcpp.h>
#include <array>
#include <bits/stdint-intn.h>
#include <bits/stdint-uintn.h>
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
  Grid(RRuntime &R, poet::GridParams grid_args);

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

  /**
   * @brief Shuffle the grid and export it to C memory area
   *
   * This will call shuffle_field inside R runtime, set the resulting grid as
   * buffered data frame in RRuntime object and write R grid to continous C
   * memory area.
   *
   * @param[in,out] buffer Pointer to C memory area
   */
  void shuffleAndExport(double *buffer);

  /**
   * @brief Unshuffle the grid and import it from C memory area into R runtime
   *
   * Write C memory area into temporary R grid variable and unshuffle it.
   *
   * @param buffer Pointer to C memory area
   */
  void importAndUnshuffle(double *buffer);

  /**
   * @brief Import a C memory area as a work package into R runtime
   *
   * Get a skeleton from getSkeletonDataFrame inside R runtime and set this as
   * buffer data frame of RRuntime object. Now convert C memory area to R data
   * structure.
   *
   * @param buffer Pointer to C memory area
   * @param wp_size Count of grid cells per work package
   */
  void importWP(double *buffer, unsigned int wp_size);

  /**
   * @brief Export a work package from R runtime into C memory area
   *
   * Set buffer data frame in RRuntime object to data frame holding the results
   * and convert this to C memory area.
   *
   * @param buffer Pointer to C memory area
   */
  void exportWP(double *buffer);

private:
  /**
   * @brief Instance of RRuntime
   *
   */
  RRuntime R;

  // uint32_t species_count;

  std::uint8_t dim;
  std::array<std::uint32_t, MAX_DIM> d_spatial;
  std::array<std::uint32_t, MAX_DIM> n_cells;

  std::map<std::string, StateMemory *> state_register;
  // std::map<std::string, std::vector<std::string>> prop_register;

  prop_vec prop_names;

  std::map<std::string, std::vector<double>> grid_init;

  /**
   * @brief Get a SkeletonDataFrame
   *
   * Return a skeleton with \e n rows of current grid
   *
   * @param rows number of rows to return skeleton
   * @return Rcpp::DataFrame Can be seen as a skeleton. The content of the data
   * frame might be irrelevant.
   */
  Rcpp::DataFrame getSkeletonDataFrame(unsigned int rows);
};
} // namespace poet
#endif // GRID_H
