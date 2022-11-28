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

#include "PhreeqcRM.h"
#include "poet/SimParams.hpp"
#include <Rcpp.h>
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

Grid::Grid(poet::GridParams grid_args) {
  // get dimension of grid
  this->dim = grid_args.n_cells.size();

  // assert that dim is 1 or 2
  assert(this->dim <= MAX_DIM && this->dim > 0);
  // also assert that vetor spatial discretization is equal in size of
  // dimensions
  assert(this->dim == grid_args.s_cells.size());

  this->n_cells.fill(0);
  this->d_spatial.fill(0);

  // store grid dimensions
  std::copy(grid_args.n_cells.begin(), grid_args.n_cells.end(),
            this->n_cells.begin());
  std::copy(grid_args.s_cells.begin(), grid_args.s_cells.end(),
            this->d_spatial.begin());

  auto init_type = get_type_id(grid_args.type);
  assert(init_type >= 0);

  uint32_t vec_size = this->getTotalCellCount();

  if (init_type != INIT_PHREEQC) {
    this->prop_names = grid_args.props;
  }

  for (auto const &prop : this->prop_names) {
    std::vector<double> prop_vec(vec_size);
    // size of the vector shall be 1
    if (init_type == INIT_SCRATCH) {
      double prop_val =
          ((Rcpp::NumericVector)grid_args.init_df[prop.c_str()])[0];
      std::fill(prop_vec.begin(), prop_vec.end(), prop_val);
    } else if (init_type == INIT_RDS) {
      prop_vec = Rcpp::as<std::vector<double>>(grid_args.init_df[prop.c_str()]);
    }
    // TODO: define behaviour for phreeqc script input

    // NOTE: if there are props defined more than once the first written value
    // will be used
    this->grid_init.insert(std::pair(prop, prop_vec));
  }
};

Grid::~Grid() {
  for (auto &[key, val] : this->state_register) {
    delete val;
  }
  if (this->phreeqc_rm != nullptr) {
    delete this->phreeqc_rm;
  }
}

poet::StateMemory *Grid::registerState(std::string module_name,
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

poet::StateMemory *Grid::getStatePointer(std::string module_name) {
  const auto it = this->state_register.find(module_name);

  if (it == this->state_register.end()) {
    throw std::range_error(
        std::string("Requested module " + module_name + " was not found"));
  }

  return it->second;
}

void Grid::deregisterState(std::string module_name) {
  const auto it = this->state_register.find(module_name);

  if (it == this->state_register.end()) {
    throw std::range_error(
        std::string("Requested module " + module_name + " was not found"));
  }

  auto *p = it->second;
  delete p;

  this->state_register.erase(it);
}

auto Grid::getDomainSize(uint8_t dircetion) -> uint32_t {
  return this->d_spatial.at(dircetion);
}

auto Grid::getGridDimension() -> uint8_t { return this->dim; }

auto Grid::getTotalCellCount() -> uint32_t {
  uint32_t sum = 1;
  for (auto const &dim : this->n_cells) {
    sum *= (dim != 0 ? dim : 1);
  }

  return sum;
}

auto Grid::getGridCellsCount(uint8_t direction) -> uint32_t {
  return this->n_cells.at(direction);
}

auto Grid::getSpeciesCount() -> uint32_t { return this->prop_names.size(); }

auto Grid::getPropNames() -> prop_vec { return this->prop_names; }

auto poet::Grid::getSpeciesByName(std::string name, std::string module_name)
    -> std::vector<double> {

  // if module name references to initial grid
  if (module_name == poet::GRID_MODULE_NAME) {
    auto const it = this->grid_init.find(name);

    if (it == this->grid_init.end()) {
      throw std::range_error(
          std::string("Species " + name + " was not found for initial grid"));
    }

    return it->second;
  }

  // else find module with name
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

  uint32_t begin_vec = prop_index * this->getTotalCellCount(),
           end_vec = ((prop_index + 1) * this->getTotalCellCount());

  return std::vector<double>(module_memory->mem.begin() + begin_vec,
                             module_memory->mem.begin() + end_vec);
}

PhreeqcRM *poet::Grid::initPhreeqc(const poet::GridParams &params) {
  PhreeqcRM *chem_model = new PhreeqcRM(this->getTotalCellCount(), 1);

  // FIXME: hardcoded options ...
  chem_model->SetErrorHandlerMode(1);
  chem_model->SetComponentH2O(false);
  chem_model->SetRebalanceFraction(0.5);
  chem_model->SetRebalanceByCell(true);
  chem_model->UseSolutionDensityVolume(false);
  chem_model->SetPartitionUZSolids(false);

  // Set concentration units
  // 1, mg/L; 2, mol/L; 3, kg/kgs
  chem_model->SetUnitsSolution(2);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem_model->SetUnitsPPassemblage(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem_model->SetUnitsExchange(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem_model->SetUnitsSurface(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem_model->SetUnitsGasPhase(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem_model->SetUnitsSSassemblage(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem_model->SetUnitsKinetics(1);

  // Set representative volume
  std::vector<double> rv;
  rv.resize(this->getTotalCellCount(), 1.0);
  chem_model->SetRepresentativeVolume(rv);

  // Set initial porosity
  std::vector<double> por;
  por.resize(this->getTotalCellCount(), 0.05);
  chem_model->SetPorosity(por);

  // Set initial saturation
  std::vector<double> sat;
  sat.resize(this->getTotalCellCount(), 1.0);
  chem_model->SetSaturation(sat);

  // Load database
  chem_model->LoadDatabase(params.database_path);
  chem_model->RunFile(true, true, true, params.input_script);
  chem_model->RunString(true, false, true, "DELETE; -all");

  return chem_model;
}
