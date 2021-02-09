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

#ifndef GRID_H
#define GRID_H

#include <RRuntime.h>
#include <Rcpp.h>

namespace poet {
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
 public:
  /**
   * @brief Construct a new Grid object
   *
   * This will call the default constructor and initializes private RRuntime
   * with given R runtime.
   *
   * @param R
   */
  Grid(RRuntime &R) : R(R){};

  /**
   * @brief Init the grid
   *
   * At this moment init will only declare and define a variable inside the R
   * runtime called grid_tmp since the whole Grid initialization and management
   * is done by the R runtime. This may change in the future.
   *
   */
  void init();

  /**
   * @brief Returns the number of elements for each gridcell
   *
   * @return unsigned int Number of elements
   */
  unsigned int getCols();

  /**
   * @brief Returns the number of grid cells
   *
   * @return unsigned int Number of grid cells
   */
  unsigned int getRows();

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

  /**
   * @brief Number of columns of grid
   *
   */
  unsigned int ncol;

  /**
   * @brief Number of rows of grid
   *
   */
  unsigned int nrow;

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
}  // namespace poet
#endif  // GRID_H
