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

#include "poet/SimParams.hpp"
#include <RInside.h>
#include <Rcpp.h>
#include <Rcpp/internal/wrap.h>
#include <algorithm>
#include <cassert>
#include <cstdint>
#include <new>
#include <numeric>
#include <poet/Grid.hpp>
#include <stdexcept>
#include <string>
#include <vector>

using namespace poet;
using namespace Rcpp;

enum { INIT_SCRATCH, INIT_PHREEQC, INIT_RDS };

static inline int8_t get_type_id(std::string type) {
  if (type == "scratch") {
    return INIT_SCRATCH;
  }
  if (type == "phreeqc") {
    return INIT_PHREEQC;
  }
  if (type == "rds") {
    return INIT_RDS;
  }
  return -1;
}

void poet::Grid::InitModuleFromParams(const poet::GridParams &grid_args) {
  assert(grid_args.n_cells.size() == grid_args.s_cells.size());

  this->SetGridDimension(grid_args.n_cells.size());

  std::copy(grid_args.n_cells.begin(), grid_args.n_cells.end(),
            this->n_cells.begin());
  std::copy(grid_args.s_cells.begin(), grid_args.s_cells.end(),
            this->grid_size.begin());

  auto init_type = get_type_id(grid_args.type);
  assert(init_type >= 0);

  switch (init_type) {
  case INIT_SCRATCH:
    this->SetPropNames(grid_args.props);
    this->InitGridFromScratch(grid_args.init_df);
    break;
  case INIT_RDS:
    this->SetPropNames(grid_args.props);
    this->InitGridFromRDS();
    break;
  case INIT_PHREEQC:
    this->InitGridFromPhreeqc();
  }
}

void poet::Grid::SetGridDimension(uint8_t dim) {
  assert(dim < MAX_DIM + 1);

  this->dim = dim;
}

void poet::Grid::SetGridCellCount(uint32_t n_x, uint32_t n_y, uint32_t n_z) {
  assert(this->dim < MAX_DIM + 1);

  this->n_cells[0] = n_x;
  this->n_cells[1] = n_y;
  // this->n_cells[2] = n_z;
}

void poet::Grid::SetGridSize(double s_x, double s_y, double s_z) {
  assert(this->dim < MAX_DIM + 1);

  this->grid_size[0] = s_x;
  this->grid_size[1] = s_y;
  // this->grid_size[2] = s_z;
}

void poet::Grid::SetPropNames(const std::vector<std::string> &prop_names) {
  this->prop_names = prop_names;
}

void poet::Grid::PushbackModuleFlow(const std::string &input,
                                    const std::string &output) {
  FlowInputOutputInfo element = {input, output};
  this->flow_vec.push_back(element);
}

void poet::Grid::PreModuleFieldCopy(uint32_t tick) {
  FlowInputOutputInfo curr_element = this->flow_vec.at(tick);

  const std::string input_module_name = curr_element.input_field;
  StateMemory *out_state = this->GetStatePointer(curr_element.output_field);

  std::vector<std::string> &mod_props = out_state->props;
  std::vector<double> &mod_field = out_state->mem;

  // copy output of another module  as input for another module
  for (uint32_t i = 0; i < mod_props.size(); i++) {
    try {
      std::vector<double> t_prop_vec =
          this->GetSpeciesByName(mod_props[i], input_module_name);

      std::copy(t_prop_vec.begin(), t_prop_vec.end(),
                mod_field.begin() + (i * this->GetTotalCellCount()));
    } catch (...) {
      // TODO: there might be something wrong ...
      continue;
    }
  }
}

Grid::Grid() {
  this->n_cells.fill(0);
  this->grid_size.fill(0);
};

Grid::~Grid() {
  for (auto &[key, val] : this->state_register) {
    delete val;
  }
}

void poet::Grid::InitGridFromScratch(const Rcpp::DataFrame &init_cell) {
  const uint32_t vec_size = this->GetTotalCellCount();

  StateMemory *curr_state =
      this->RegisterState(poet::GRID_MODULE_NAME, this->prop_names);
  curr_state->props = this->prop_names;

  std::vector<double> &curr_field = curr_state->mem;

  for (auto const &prop : this->prop_names) {
    // size of the vector shall be 1
    double prop_val = ((Rcpp::NumericVector)init_cell[prop.c_str()])[0];
    std::vector<double> prop_vec(vec_size, prop_val);

    curr_field.insert(curr_field.end(), prop_vec.begin(), prop_vec.end());
  }

  this->grid_init = curr_state;
}
void poet::Grid::InitGridFromRDS() {
  // TODO
}
void poet::Grid::InitGridFromPhreeqc() {
  // TODO
}

poet::StateMemory *Grid::RegisterState(std::string module_name,
                                       std::vector<std::string> props) {
  poet::StateMemory *mem = new poet::StateMemory;

  mem->props = props;

  const auto [it, success] = this->state_register.insert({module_name, mem});

  if (!success) {
    delete mem;
    throw std::bad_alloc();
  }

  return mem;
}

poet::StateMemory *Grid::GetStatePointer(std::string module_name) {
  const auto it = this->state_register.find(module_name);

  if (it == this->state_register.end()) {
    throw std::range_error(
        std::string("Requested module " + module_name + " was not found"));
  }

  return it->second;
}

auto Grid::GetGridSize(uint8_t direction) const -> uint32_t {
  return this->grid_size.at(direction);
}

auto Grid::GetGridDimension() const -> uint8_t { return this->dim; }

auto Grid::GetTotalCellCount() const -> uint32_t {
  uint32_t sum = 1;
  for (auto const &dim : this->n_cells) {
    sum *= (dim != 0 ? dim : 1);
  }

  return sum;
}

auto Grid::GetGridCellsCount(uint8_t direction) const -> uint32_t {
  return this->n_cells.at(direction);
}

auto poet::Grid::GetInitialGrid() const -> StateMemory * {
  return this->grid_init;
}

auto Grid::GetSpeciesCount() const -> uint32_t {
  return this->prop_names.size();
}

auto Grid::GetPropNames() const -> std::vector<std::string> {
  return this->prop_names;
}

auto poet::Grid::GetSpeciesByName(std::string name,
                                  std::string module_name) const
    -> std::vector<double> {

  auto const it = this->state_register.find(module_name);

  if (it == this->state_register.end()) {
    throw std::range_error(
        std::string("Requested module " + module_name + " was not found"));
  }

  poet::StateMemory const *module_memory = it->second;

  auto const &prop_vector = module_memory->props;
  auto const find_res = std::find(prop_vector.begin(), prop_vector.end(), name);
  if (find_res == prop_vector.end()) {
    throw std::range_error(std::string(
        "Species " + name + " was not found for module " + module_name));
  }
  uint32_t prop_index = find_res - prop_vector.begin();

  uint32_t begin_vec = prop_index * this->GetTotalCellCount(),
           end_vec = ((prop_index + 1) * this->GetTotalCellCount());

  return std::vector<double>(module_memory->mem.begin() + begin_vec,
                             module_memory->mem.begin() + end_vec);
}

// HACK: Helper function
void poet::Grid::WriteFieldsToR(RInside &R) {

  for (auto const &elem : this->state_register) {
    std::string field_name = elem.first;
    StateMemory *field = elem.second;

    const uint32_t n_cells_per_prop = field->mem.size() / field->props.size();

    R["TMP"] = Rcpp::wrap(field->mem);
    R["TMP_PROPS"] = Rcpp::wrap(field->props);
    R.parseEval(std::string(
        "mysetup$" + field_name + " <- setNames(data.frame(matrix(TMP, nrow=" +
        std::to_string(n_cells_per_prop) + ")), TMP_PROPS)"));
  }
}
