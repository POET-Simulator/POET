/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2023 Marco De Lucia, Max Luebke (GFZ Potsdam)
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

#include "SimParams.hpp"

#include <RInside.h>
#include <Rcpp.h>
#include <array>
#include <cstddef>
#include <cstdint>
#include <list>
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

using FlowInputOutputInfo = struct inOut_info {
  std::string input_field;
  std::string output_field;
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

public:
  Grid();

  ~Grid();

  void InitModuleFromParams(const poet::GridParams &grid_args);

  void SetGridDimension(uint8_t dim);
  void SetGridCellCount(uint32_t n_x, uint32_t n_y = 0, uint32_t n_z = 0);
  void SetGridSize(double s_x, double s_y = 0., double s_z = 0.);
  void SetPropNames(const std::vector<std::string> &prop_names);

  void PushbackModuleFlow(const std::string &input, const std::string &output);
  void PreModuleFieldCopy(uint32_t tick);

  void InitGridFromScratch(const Rcpp::DataFrame &init_cell);
  void InitGridFromRDS();
  void InitGridFromPhreeqc();

  auto GetGridDimension() const -> uint8_t;
  auto GetTotalCellCount() const -> uint32_t;
  auto GetGridCellsCount(uint8_t direction) const -> uint32_t;
  auto GetGridSize(uint8_t direction) const -> uint32_t;

  auto RegisterState(std::string module_name, std::vector<std::string> props)
      -> StateMemory *;
  auto GetStatePointer(std::string module_name) -> StateMemory *;

  auto GetInitialGrid() const -> StateMemory *;

  auto GetSpeciesCount() const -> uint32_t;
  auto GetPropNames() const -> std::vector<std::string>;

  auto GetSpeciesByName(std::string name,
                        std::string module_name = poet::GRID_MODULE_NAME) const
      -> std::vector<double>;

  void WriteFieldsToR(RInside &R);

private:
  std::vector<FlowInputOutputInfo> flow_vec;

  std::uint8_t dim = 0;
  std::array<double, MAX_DIM> grid_size;
  std::array<std::uint32_t, MAX_DIM> n_cells;

  std::map<std::string, StateMemory *> state_register;
  StateMemory *grid_init = std::nullptr_t();

  std::vector<std::string> prop_names;
};
} // namespace poet
#endif // GRID_H
