#include "poet/ChemistryModule.hpp"
#include "PhreeqcRM.h"
#include "poet/DHT_Wrapper.hpp"

#include <algorithm>
#include <bits/stdint-uintn.h>
#include <cassert>
#include <cstdint>
#include <map>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#ifndef POET_USE_PRM
poet::ChemistryModule::ChemistryModule(uint32_t nxyz, uint32_t wp_size,
                                       MPI_Comm communicator)
    : PhreeqcRM(nxyz, 1), group_comm(communicator), wp_size(wp_size) {

  MPI_Comm_size(communicator, &this->comm_size);
  MPI_Comm_rank(communicator, &this->comm_rank);

  this->is_sequential = (this->comm_size == 1);
  this->is_master = (this->comm_rank == 0);

  if (!is_sequential && is_master) {
    MPI_Bcast(&wp_size, 1, MPI_UINT32_T, 0, this->group_comm);
  }
}

poet::ChemistryModule::~ChemistryModule() {
  if (dht) {
    delete dht;
  }
}

poet::ChemistryModule
poet::ChemistryModule::createWorker(MPI_Comm communicator) {
  uint32_t wp_size;
  MPI_Bcast(&wp_size, 1, MPI_UINT32_T, 0, communicator);

  return ChemistryModule(wp_size, wp_size, communicator);
}

#endif

void poet::ChemistryModule::RunInitFile(const std::string &input_script_path) {
#ifndef POET_USE_PRM
  if (this->is_master) {
    int f_type = CHEM_INIT;
    PropagateFunctionType(f_type);

    int count = input_script_path.size();
    ChemBCast(&count, 1, MPI_INT);
    ChemBCast(const_cast<char *>(input_script_path.data()), count, MPI_CHAR);
  }
#endif

  this->RunFile(true, true, false, input_script_path);
  this->RunString(true, false, false, "DELETE; -all; PRINT; -warnings 0;");

  this->FindComponents();

  PhreeqcRM::initializePOET(this->speciesPerModule, this->prop_names);
  this->prop_count = prop_names.size();

  char exchange = (speciesPerModule[1] == 0 ? -1 : 1);
  char kinetics = (speciesPerModule[2] == 0 ? -1 : 1);
  char equilibrium = (speciesPerModule[3] == 0 ? -1 : 1);
  char surface = (speciesPerModule[4] == 0 ? -1 : 1);

#ifdef POET_USE_PRM
  std::vector<int> ic1;
  ic1.resize(this->nxyz * 7, -1);
  // TODO: hardcoded reaction modules
  for (int i = 0; i < nxyz; i++) {
    ic1[i] = 1;                   // Solution 1
    ic1[nxyz + i] = equilibrium;  // Equilibrium 1
    ic1[2 * nxyz + i] = exchange; // Exchange none
    ic1[3 * nxyz + i] = surface;  // Surface none
    ic1[4 * nxyz + i] = -1;       // Gas phase none
    ic1[5 * nxyz + i] = -1;       // Solid solutions none
    ic1[6 * nxyz + i] = kinetics; // Kinetics 1
  }

  this->InitialPhreeqc2Module(ic1);
#else
  std::vector<int> ic1;
  ic1.resize(this->nxyz * 7, -1);
  // TODO: hardcoded reaction modules
  for (int i = 0; i < nxyz; i++) {
    ic1[i] = 1;                   // Solution 1
    ic1[nxyz + i] = equilibrium;  // Equilibrium 1
    ic1[2 * nxyz + i] = exchange; // Exchange none
    ic1[3 * nxyz + i] = surface;  // Surface none
    ic1[4 * nxyz + i] = -1;       // Gas phase none
    ic1[5 * nxyz + i] = -1;       // Solid solutions none
    ic1[6 * nxyz + i] = kinetics; // Kinetics 1
  }

  this->InitialPhreeqc2Module(ic1);
#endif
}

#ifndef POET_USE_PRM
void poet::ChemistryModule::mergeFieldWithModule(const SingleCMap &input_map,
                                                 std::uint32_t n_cells) {

  if (is_master) {
    int f_type = CHEM_INIT_SPECIES;
    PropagateFunctionType(f_type);
  }

  std::vector<std::string> essentials_backup{
      prop_names.begin() + speciesPerModule[0], prop_names.end()};

  std::vector<std::string> new_solution_names{
      this->prop_names.begin(), this->prop_names.begin() + speciesPerModule[0]};

  if (is_master) {
    for (const auto &init_val : input_map) {
      std::string name = init_val.first;
      if (std::find(new_solution_names.begin(), new_solution_names.end(),
                    name) == new_solution_names.end()) {
        int size = name.size();
        ChemBCast(&size, 1, MPI_INT);
        ChemBCast(name.data(), name.size(), MPI_CHAR);
        new_solution_names.push_back(name);
      }
    }
    int end = 0;
    ChemBCast(&end, 1, MPI_INT);
  } else {
    constexpr int MAXSIZE = 128;
    MPI_Status status;
    int recv_size;
    char recv_buffer[MAXSIZE];
    while (1) {
      ChemBCast(&recv_size, 1, MPI_INT);
      if (recv_size == 0) {
        break;
      }
      ChemBCast(recv_buffer, recv_size, MPI_CHAR);
      recv_buffer[recv_size] = '\0';
      new_solution_names.push_back(std::string(recv_buffer));
    }
  }

  // now sort the new values
  std::sort(new_solution_names.begin() + 3, new_solution_names.end());

  // and append other processes than solutions
  std::vector<std::string> new_prop_names = new_solution_names;
  new_prop_names.insert(new_prop_names.end(), essentials_backup.begin(),
                        essentials_backup.end());

  std::vector<std::string> old_prop_names{this->prop_names};
  this->prop_names = std::move(new_prop_names);
  this->prop_count = prop_names.size();

  this->SetPOETSolutionNames(new_solution_names);

  if (is_master) {
    this->n_cells = n_cells;

    this->field.clear();
    this->field.reserve(this->prop_count * n_cells);

    std::vector<double> init_values;
    this->getDumpedField(init_values);

    const std::vector<std::string> ess_names = old_prop_names;

    for (int i = 0; i < prop_names.size(); i++) {
      double value = init_values[i];
      std::vector<double> curr_vec(n_cells, value);
      this->field.insert(field.end(), curr_vec.begin(), curr_vec.end());
    }
  }
}

void poet::ChemistryModule::SetDHTEnabled(bool enable, uint32_t size_mb) {
  constexpr uint32_t MB_FACTOR = 1E6;

  if (this->is_master) {
    int ftype = CHEM_DHT_ENABLE;
    PropagateFunctionType(ftype);
    ChemBCast(&enable, 1, MPI_CXX_BOOL);
    ChemBCast(&size_mb, 1, MPI_UINT32_T);
  }

  this->dht_enabled = enable;

  if (enable) {
    MPI_Comm dht_comm;

    if (this->is_master) {
      MPI_Comm_split(this->group_comm, MPI_UNDEFINED, this->comm_rank,
                     &dht_comm);
      return;
    }

    MPI_Comm_split(this->group_comm, 1, this->comm_rank, &dht_comm);

    if (this->dht) {
      delete this->dht;
    }

    const uint32_t dht_size = size_mb * MB_FACTOR;

    this->dht =
        new DHT_Wrapper(dht_comm, dht_size, this->prop_count, this->prop_count);
  }
}

void poet::ChemistryModule::SetDHTSnaps(int type, const std::string &out_dir) {
  if (this->is_master) {
    int ftype = CHEM_DHT_SNAPS;
    PropagateFunctionType(ftype);

    int str_size = out_dir.size();

    ChemBCast(&type, 1, MPI_INT);
    ChemBCast(&str_size, 1, MPI_INT);
    ChemBCast(const_cast<char *>(out_dir.data()), str_size, MPI_CHAR);
  }

  this->dht_file_out_dir = out_dir;
}

void poet::ChemistryModule::SetDHTSignifVector(
    std::vector<uint32_t> signif_vec) {
  if (this->is_master) {
    if (signif_vec.size() != this->prop_count) {
      throw std::runtime_error(
          "Significant vector sizes mismatches prop count.");
    }

    int ftype = CHEM_DHT_SIGNIF_VEC;
    PropagateFunctionType(ftype);
    ChemBCast(signif_vec.data(), signif_vec.size(), MPI_UINT32_T);

    return;
  }

  this->dht->SetSignifVector(signif_vec);
}

void poet::ChemistryModule::SetDHTPropTypeVector(
    std::vector<uint32_t> proptype_vec) {
  if (this->is_master) {
    if (proptype_vec.size() != this->prop_count) {
      throw std::runtime_error("Prop type vector sizes mismatches prop count.");
    }

    int ftype = CHEM_DHT_PROP_TYPE_VEC;
    PropagateFunctionType(ftype);
    ChemBCast(proptype_vec.data(), proptype_vec.size(), MPI_UINT32_T);

    return;
  }

  this->dht->SetPropTypeVector(proptype_vec);
}

void poet::ChemistryModule::ReadDHTFile(const std::string &input_file) {
  if (this->is_master) {
    int ftype = CHEM_DHT_READ_FILE;
    PropagateFunctionType(ftype);
    int str_size = input_file.size();

    ChemBCast(&str_size, 1, MPI_INT);
    ChemBCast(const_cast<char *>(input_file.data()), str_size, MPI_CHAR);
  }

  if (!this->dht_enabled) {
    throw std::runtime_error("DHT file cannot be read. DHT is not enabled.");
  }

  if (!this->is_master) {
    WorkerReadDHTDump(input_file);
  }
}

#else // POET_USE_PRM

inline void poet::ChemistryModule::PrmToPoetField(std::vector<double> &field) {
  this->getDumpedField(field);
}

void poet::ChemistryModule::RunCells() {
  PhreeqcRM::RunCells();

  std::vector<double> tmp_field;

  PrmToPoetField(tmp_field);
  this->field = tmp_field;

  for (uint32_t i = 0; i < field.size(); i++) {
    uint32_t back_i = static_cast<uint32_t>(i / this->nxyz);
    uint32_t mod_i = i % this->nxyz;

    field[i] = tmp_field[back_i + (mod_i * this->prop_count)];
  }
}

#endif
