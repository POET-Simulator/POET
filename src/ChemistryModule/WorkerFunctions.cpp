#include "IrmResult.h"
#include "poet/ChemistryModule.hpp"

#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

inline std::string get_string(int root, MPI_Comm communicator) {
  int count;
  MPI_Bcast(&count, 1, MPI_INT, root, communicator);

  char *buffer = new char[count + 1];
  MPI_Bcast(buffer, count, MPI_CHAR, root, communicator);

  buffer[count] = '\0';

  std::string ret_str(buffer);
  delete[] buffer;

  return ret_str;
}

void poet::ChemistryModule::WorkerLoop() {
  struct worker_s timings;

  // HACK: defining the worker iteration count here, which will increment after
  // each CHEM_ITER_END message
  uint32_t iteration = 1;
  bool loop = true;

  while (loop) {
    int func_type;
    PropagateFunctionType(func_type);

    switch (func_type) {
    case CHEM_INIT: {
      RunInitFile(get_string(0, this->group_comm));
      break;
    }
    case CHEM_DHT_ENABLE: {
      bool enable;
      ChemBCast(&enable, 1, MPI_CXX_BOOL);

      uint32_t size_mb;
      ChemBCast(&size_mb, 1, MPI_UINT32_T);

      SetDHTEnabled(enable, size_mb);
      break;
    }
    case CHEM_DHT_SIGNIF_VEC: {
      std::vector<uint32_t> input_vec(this->prop_count);
      ChemBCast(input_vec.data(), this->prop_count, MPI_UINT32_T);

      SetDHTSignifVector(input_vec);
      break;
    }
    case CHEM_DHT_PROP_TYPE_VEC: {
      std::vector<uint32_t> input_vec(this->prop_count);
      ChemBCast(input_vec.data(), this->prop_count, MPI_UINT32_T);

      SetDHTPropTypeVector(input_vec);
      break;
    }
    case CHEM_DHT_SNAPS: {
      int type;
      ChemBCast(&type, 1, MPI_INT);

      SetDHTSnaps(type, get_string(0, this->group_comm));

      break;
    }
    case CHEM_DHT_READ_FILE: {
      ReadDHTFile(get_string(0, this->group_comm));
      break;
    }
    case CHEM_WORK_LOOP: {
      WorkerProcessPkgs(timings, iteration);
      break;
    }
    case CHEM_PERF: {
      int type;
      ChemBCast(&type, 1, MPI_INT);
      if (type < WORKER_DHT_HITS) {
        WorkerPerfToMaster(type, timings);
        break;
      }
      WorkerMetricsToMaster(type);
      break;
    }
    case CHEM_BREAK_MAIN_LOOP: {
      WorkerPostSim(iteration);
      loop = false;
      break;
    }
    default: {
      throw std::runtime_error("Worker received unknown tag from master.");
    }
    }
  }
}

void poet::ChemistryModule::WorkerProcessPkgs(struct worker_s &timings,
                                              uint32_t &iteration) {
  MPI_Status probe_status;
  bool loop = true;

  MPI_Barrier(this->group_comm);

  while (loop) {
    double idle_a = MPI_Wtime();
    MPI_Probe(0, MPI_ANY_TAG, this->group_comm, &probe_status);
    double idle_b = MPI_Wtime();

    switch (probe_status.MPI_TAG) {
    case LOOP_WORK: {
      timings.idle_t += idle_b - idle_a;
      int count;
      MPI_Get_count(&probe_status, MPI_DOUBLE, &count);

      WorkerDoWork(probe_status, count, timings);
      break;
    }
    case LOOP_END: {
      WorkerPostIter(probe_status, iteration);
      iteration++;
      loop = false;
      break;
    }
    }
  }
}

void poet::ChemistryModule::WorkerDoWork(MPI_Status &probe_status,
                                         int double_count,
                                         struct worker_s &timings) {
  int local_work_package_size = 0;

  static int counter = 1;

  double dht_get_start, dht_get_end;
  double phreeqc_time_start, phreeqc_time_end;
  double dht_fill_start, dht_fill_end;

  uint32_t iteration;
  double dt;
  double current_sim_time;

  const uint32_t n_cells_times_props = this->prop_count * this->wp_size;
  std::vector<double> vecCurrWP(n_cells_times_props + BUFFER_OFFSET);
  int count = double_count;

  /* receive */
  MPI_Recv(vecCurrWP.data(), count, MPI_DOUBLE, 0, LOOP_WORK, this->group_comm,
           MPI_STATUS_IGNORE);

  /* decrement count of work_package by BUFFER_OFFSET */
  count -= BUFFER_OFFSET;

  /* check for changes on all additional variables given by the 'header' of
   * mpi_buffer */

  // work_package_size
  local_work_package_size = vecCurrWP[count];

  // current iteration of simulation
  iteration = vecCurrWP[count + 1];

  // current timestep size
  dt = vecCurrWP[count + 2];

  // current simulation time ('age' of simulation)
  current_sim_time = vecCurrWP[count + 3];

  /* 4th double value is currently a placeholder */
  // placeholder = mpi_buffer[count+4];

  // std::vector<double> vecCurrWP(
  //     mpi_buffer,
  //     mpi_buffer + (local_work_package_size * this->prop_names.size()));
  vecCurrWP.resize(n_cells_times_props);
  std::vector<int32_t> vecMappingWP(this->wp_size);

  DHT_ResultObject DHT_Results;

  for (uint32_t i = 0; i < local_work_package_size; i++) {
    vecMappingWP[i] = i;
  }

  if (local_work_package_size != this->wp_size) {
    // std::vector<double> vecFiller(
    //     (this->wp_size - local_work_package_size) * prop_count, 0);
    // vecCurrWP.insert(vecCurrWP.end(), vecFiller.begin(), vecFiller.end());

    // set all remaining cells to inactive
    for (int i = local_work_package_size; i < this->wp_size; i++) {
      vecMappingWP[i] = -1;
    }
  }

  if (dht_enabled) {
    /* check for values in DHT */
    dht_get_start = MPI_Wtime();
    DHT_Results = dht->checkDHT(local_work_package_size, dt, vecCurrWP);
    dht_get_end = MPI_Wtime();

    DHT_Results.ResultsToMapping(vecMappingWP);
  }

  phreeqc_time_start = MPI_Wtime();

  WorkerRunWorkPackage(vecCurrWP, vecMappingWP, current_sim_time, dt);

  phreeqc_time_end = MPI_Wtime();

  if (dht_enabled) {
    DHT_Results.ResultsToWP(vecCurrWP);
  }

  /* send results to master */
  MPI_Request send_req;
  MPI_Isend(vecCurrWP.data(), count, MPI_DOUBLE, 0, LOOP_WORK, MPI_COMM_WORLD,
            &send_req);

  if (dht_enabled) {
    /* write results to DHT */
    dht_fill_start = MPI_Wtime();
    dht->fillDHT(local_work_package_size, DHT_Results, vecCurrWP);
    dht_fill_end = MPI_Wtime();

    timings.dht_get += dht_get_end - dht_get_start;
    timings.dht_fill += dht_fill_end - dht_fill_start;
  }

  timings.phreeqc_t += phreeqc_time_end - phreeqc_time_start;

  MPI_Wait(&send_req, MPI_STATUS_IGNORE);
}

void poet::ChemistryModule::WorkerPostIter(MPI_Status &prope_status,
                                           uint32_t iteration) {
  MPI_Recv(NULL, 0, MPI_DOUBLE, 0, LOOP_END, this->group_comm,
           MPI_STATUS_IGNORE);
  if (this->dht_enabled) {
    this->dht->printStatistics();

    if (this->dht_snaps_type == DHT_FILES_ITEREND) {
      WorkerWriteDHTDump(iteration);
    }
  }
}
void poet::ChemistryModule::WorkerPostSim(uint32_t iteration) {
  /* before death, submit profiling/timings to master*/

  // double timings_serialized[4];
  // timings_serialized[0] = timings.phreeqc_t;
  // timings_serialized[1] = timings.dht_get;
  // timings_serialized[2] = timings.dht_fill;
  // timings_serialized[3] = timings.idle_t;

  // // timings
  // MPI_Send(timings_serialized, 4, MPI_DOUBLE, 0, 0, this->group_comm);

  // // MPI_Send(&phreeqc_count, 1, MPI_INT, 0, TAG_TIMING, MPI_COMM_WORLD);
  // // MPI_Send(&idle_t, 1, MPI_DOUBLE, 0, TAG_TIMING, MPI_COMM_WORLD);

  // if (this->dht_enabled) {
  //   // dht_perf
  //   int dht_perf[3];
  //   dht_perf[0] = dht->getHits();
  //   dht_perf[1] = dht->getMisses();
  //   dht_perf[2] = dht->getEvictions();
  //   MPI_Send(dht_perf, 3, MPI_INT, 0, 0, this->group_comm);
  // }

  if (this->dht_enabled && this->dht_snaps_type > DHT_FILES_SIMEND) {
    WorkerWriteDHTDump(iteration);
  }
}

void poet::ChemistryModule::WorkerWriteDHTDump(uint32_t iteration) {
  std::stringstream out;
  out << this->dht_file_out_dir << "/iter_" << std::setfill('0') << std::setw(4)
      << iteration << ".dht";
  int res = dht->tableToFile(out.str().c_str());
  if (res != DHT_SUCCESS && this->comm_rank == 2)
    std::cerr
        << "CPP: Worker: Error in writing current state of DHT to file.\n";
  else if (this->comm_rank == 2)
    std::cout << "CPP: Worker: Successfully written DHT to file " << out.str()
              << "\n";
}
void poet::ChemistryModule::WorkerReadDHTDump(
    const std::string &dht_input_file) {
  int res = dht->fileToTable((char *)dht_input_file.c_str());
  if (res != DHT_SUCCESS) {
    if (res == DHT_WRONG_FILE) {
      if (this->comm_rank == 1)
        std::cerr
            << "CPP: Worker: Wrong file layout! Continue with empty DHT ...\n";
    } else {
      if (this->comm_rank == 1)
        std::cerr << "CPP: Worker: Error in loading current state of DHT from "
                     "file. Continue with empty DHT ...\n";
    }
  } else {
    if (this->comm_rank == 2)
      std::cout << "CPP: Worker: Successfully loaded state of DHT from file "
                << dht_input_file << "\n";
  }
}

IRM_RESULT
poet::ChemistryModule::WorkerRunWorkPackage(std::vector<double> &vecWP,
                                            std::vector<int32_t> &vecMapping,
                                            double dSimTime, double dTimestep) {
  if (this->wp_size != vecMapping.size()) {
    return IRM_INVALIDARG;
  }

  if ((this->wp_size * this->prop_count) != vecWP.size()) {
    return IRM_INVALIDARG;
  }

  // check if we actually need to start phreeqc
  bool bRunPhreeqc = false;
  for (const auto &aMappingNum : vecMapping) {
    if (aMappingNum != -1) {
      bRunPhreeqc = true;
      break;
    }
  }

  if (!bRunPhreeqc) {
    return IRM_OK;
  }

  std::vector<double> vecCopy;

  vecCopy = vecWP;
  for (uint32_t i = 0; i < this->prop_count; i++) {
    for (uint32_t j = 0; j < this->wp_size; j++) {
      vecWP[(i * this->wp_size) + j] = vecCopy[(j * this->prop_count) + i];
    }
  }

  IRM_RESULT result;
  this->PhreeqcRM::CreateMapping(vecMapping);
  this->SetInternalsFromWP(vecWP, this->wp_size);

  this->PhreeqcRM::SetTime(dSimTime);
  this->PhreeqcRM::SetTimeStep(dTimestep);

  result = this->PhreeqcRM::RunCells();

  this->GetWPFromInternals(vecWP, this->wp_size);

  vecCopy = vecWP;
  for (uint32_t i = 0; i < this->prop_count; i++) {
    for (uint32_t j = 0; j < this->wp_size; j++) {
      vecWP[(j * this->prop_count) + i] = vecCopy[(i * this->wp_size) + j];
    }
  }

  return result;
}

void poet::ChemistryModule::WorkerPerfToMaster(int type,
                                               const struct worker_s &timings) {
  switch (type) {
  case WORKER_PHREEQC: {
    MPI_Gather(&timings.phreeqc_t, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
               this->group_comm);
    break;
  }
  case WORKER_DHT_GET: {
    MPI_Gather(&timings.dht_get, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
               this->group_comm);
    break;
  }
  case WORKER_DHT_FILL: {
    MPI_Gather(&timings.dht_fill, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
               this->group_comm);
    break;
  }
  case WORKER_IDLE: {
    MPI_Gather(&timings.idle_t, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
               this->group_comm);
    break;
  }
  default: {
    throw std::runtime_error("Unknown perf type in master's message.");
  }
  }
}

void poet::ChemistryModule::WorkerMetricsToMaster(int type) {
  uint32_t value;
  switch (type) {
  case WORKER_DHT_HITS: {
    value = dht->getHits();
    break;
  }
  case WORKER_DHT_MISS: {
    value = dht->getMisses();
    break;
  }
  case WORKER_DHT_EVICTIONS: {
    value = dht->getEvictions();
    break;
  }
  default: {
    throw std::runtime_error("Unknown perf type in master's message.");
  }
  }
  MPI_Gather(&value, 1, MPI_UINT32_T, NULL, 1, MPI_UINT32_T, 0,
             this->group_comm);
}
