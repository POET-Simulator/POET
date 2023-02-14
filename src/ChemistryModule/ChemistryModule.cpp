#include "poet/ChemistryModule.hpp"
#include "PhreeqcRM.h"
#include "poet/DHT_Wrapper.hpp"

#include <cassert>
#include <cstdint>
#include <mpi.h>
#include <stdexcept>
#include <string>
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

  std::vector<std::string> props;
  this->speciesPerModule.reserve(this->MODULE_COUNT);

  std::vector<std::string> curr_prop_names = this->GetComponents();
  props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  this->speciesPerModule.push_back(curr_prop_names.size());

  curr_prop_names = this->GetEquilibriumPhases();
  props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  char equilibrium = (curr_prop_names.empty() ? -1 : 1);
  this->speciesPerModule.push_back(curr_prop_names.size());

  curr_prop_names = this->GetExchangeNames();
  props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  char exchange = (curr_prop_names.empty() ? -1 : 1);
  this->speciesPerModule.push_back(curr_prop_names.size());

  curr_prop_names = this->GetSurfaceNames();
  props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  char surface = (curr_prop_names.empty() ? -1 : 1);
  this->speciesPerModule.push_back(curr_prop_names.size());

  // curr_prop_names = this->GetGasComponents();
  // props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  // char gas = (curr_prop_names.empty() ? -1 : 1);
  // this->speciesPerModule.push_back(curr_prop_names.size());

  // curr_prop_names = this->GetSolidSolutionNames();
  // props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  // char ssolutions = (curr_prop_names.empty() ? -1 : 1);
  // this->speciesPerModule.push_back(curr_prop_names.size());

  curr_prop_names = this->GetKineticReactions();
  props.insert(props.end(), curr_prop_names.begin(), curr_prop_names.end());
  char kinetics = (curr_prop_names.empty() ? -1 : 1);
  this->speciesPerModule.push_back(curr_prop_names.size());

  this->prop_count = props.size();

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
  if (!this->is_master || this->is_sequential) {
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
  }
#endif
  this->prop_names = props;
}

#ifndef POET_USE_PRM
void poet::ChemistryModule::InitializeField(
    uint32_t n_cells, const poet::ChemistryModule::SingleCMap &mapped_values) {
  if (this->is_master) {
    this->field.reserve(this->prop_count * n_cells);
    for (const std::string &prop_name : this->prop_names) {
      const auto m_it = mapped_values.find(prop_name);
      if (m_it == mapped_values.end()) {
        throw std::domain_error(
            "Prop names vector does not match any key in given map.");
      }

      const std::vector<double> field_row(n_cells, m_it->second);
      const auto chem_field_end = this->field.end();
      this->field.insert(chem_field_end, field_row.begin(), field_row.end());
    }
    this->n_cells = n_cells;
  }
}

void poet::ChemistryModule::InitializeField(
    const poet::ChemistryModule::VectorCMap &mapped_values) {
  if (this->is_master) {
    for (const std::string &prop_name : this->prop_names) {
      const auto m_it = mapped_values.find(prop_name);
      if (m_it == mapped_values.end()) {
        throw std::domain_error(
            "Prop names vector does not match any key in given map.");
      }

      const auto field_row = m_it->second;
      const auto chem_field_end = this->field.end();
      this->field.insert(chem_field_end, field_row.begin(), field_row.end());
    }
    this->n_cells = mapped_values.begin()->second.size();
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

void poet::ChemistryModule::GetWPFromInternals(std::vector<double> &vecWP,
                                               uint32_t wp_size) {
  std::vector<double> vecCurrOutput;

  vecWP.clear();
  vecWP.reserve(this->prop_count * wp_size);

  this->GetConcentrations(vecCurrOutput);
  vecWP.insert(vecWP.end(), vecCurrOutput.begin(), vecCurrOutput.end());

  if (this->speciesPerModule[1] != 0) {
    this->GetPPhaseMoles(vecCurrOutput);
    vecWP.insert(vecWP.end(), vecCurrOutput.begin(), vecCurrOutput.end());
  }

  // NOTE: Block for 'Surface' and 'Exchange' is missing because of missing
  // Getters @ PhreeqcRM
  // ...
  // BLOCK_END

  if (this->speciesPerModule[4] != 0) {
    this->GetKineticsMoles(vecCurrOutput);
    vecWP.insert(vecWP.end(), vecCurrOutput.begin(), vecCurrOutput.end());
  }
}
void poet::ChemistryModule::SetInternalsFromWP(const std::vector<double> &vecWP,
                                               uint32_t wp_size) {
  uint32_t iCurrElements;

  auto itStart = vecWP.begin();
  auto itEnd = itStart;

  // this->SetMappingForWP(iCurrWPSize);

  int nchem = this->GetChemistryCellCount();

  iCurrElements = this->speciesPerModule[0];

  itEnd += iCurrElements * wp_size;
  this->SetConcentrations(std::vector<double>(itStart, itEnd));
  itStart = itEnd;

  // Equlibirum Phases
  if ((iCurrElements = this->speciesPerModule[1]) != 0) {
    itEnd += iCurrElements * wp_size;
    this->SetPPhaseMoles(std::vector<double>(itStart, itEnd));
    itStart = itEnd;
  }

  // // NOTE: Block for 'Surface' and 'Exchange' is missing because of missing
  // // setters @ PhreeqcRM
  // // ...
  // // BLOCK_END

  if ((iCurrElements = this->speciesPerModule[4]) != 0) {
    itEnd += iCurrElements * wp_size;
    this->SetKineticsMoles(std::vector<double>(itStart, itEnd));
    itStart = itEnd;
  }
}

#else //POET_USE_PRM

inline void poet::ChemistryModule::PrmToPoetField(std::vector<double> &field) {
  GetConcentrations(field);
  int col = GetSelectedOutputColumnCount();
  int rows = GetSelectedOutputRowCount();

  field.reserve(field.size() + 3 * rows);

  std::vector<double> so;
  GetSelectedOutput(so);

  for (int j = 0; j < col; j += 2) {
    const auto start = so.begin() + (j * rows);
    const auto end = start + rows;
    field.insert(field.end(), start, end);
  }
}

void poet::ChemistryModule::RunCells() {
  PhreeqcRM::RunCells();
  PrmToPoetField(this->field);
}

#endif
