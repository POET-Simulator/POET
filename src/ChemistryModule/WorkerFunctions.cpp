//  Time-stamp: "Last modified 2023-07-27 15:06:02 mluebke"

#include "IrmResult.h"
#include "poet/ChemistryModule.hpp"
#include "poet/DHT_Wrapper.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace poet {

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
    case CHEM_INIT_SPECIES: {
      Field dummy{0};
      initializeField(dummy);
      break;
    }
    case CHEM_DHT_ENABLE: {
      bool enable;
      ChemBCast(&enable, 1, MPI_CXX_BOOL);

      std::uint64_t size_mb;
      ChemBCast(&size_mb, 1, MPI_UINT64_T);

      std::vector<std::string> name_dummy;

      SetDHTEnabled(enable, size_mb, name_dummy);
      break;
    }
    case CHEM_DHT_SIGNIF_VEC: {
      std::vector<uint32_t> input_vec;

      SetDHTSignifVector(input_vec);
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

  vecCurrWP.resize(n_cells_times_props);
  std::vector<std::uint32_t> vecMappingWP(local_work_package_size);

  {
    std::uint32_t i = 0;
    std::generate(vecMappingWP.begin(), vecMappingWP.end(),
                  [&] { return i++; });
  }

  if (dht_enabled) {
    /* check for values in DHT */
    dht_get_start = MPI_Wtime();
    dht->checkDHT(local_work_package_size, dt, vecCurrWP, vecMappingWP);
    dht_get_end = MPI_Wtime();

    // DHT_Results.ResultsToMapping(vecMappingWP);
  }

  phreeqc_time_start = MPI_Wtime();

  if (WorkerRunWorkPackage(vecCurrWP, vecMappingWP, current_sim_time, dt) !=
      IRM_OK) {
    throw std::runtime_error("Phreeqc threw an error!");
  };

  phreeqc_time_end = MPI_Wtime();

  if (dht_enabled) {
    dht->resultsToWP(vecCurrWP);
  }

  /* send results to master */
  MPI_Request send_req;
  MPI_Isend(vecCurrWP.data(), count, MPI_DOUBLE, 0, LOOP_WORK, MPI_COMM_WORLD,
            &send_req);

  if (dht_enabled) {
    /* write results to DHT */
    dht_fill_start = MPI_Wtime();
    dht->fillDHT(local_work_package_size, vecCurrWP);
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
    dht_hits.push_back(dht->getHits());
    dht_evictions.push_back(dht->getEvictions());
    dht->resetCounter();

    if (this->dht_snaps_type == DHT_SNAPS_ITEREND) {
      WorkerWriteDHTDump(iteration);
    }
  }
}
void poet::ChemistryModule::WorkerPostSim(uint32_t iteration) {
  if (this->dht_enabled && this->dht_snaps_type == DHT_SNAPS_SIMEND) {
    WorkerWriteDHTDump(iteration);
  }
}

void poet::ChemistryModule::WorkerWriteDHTDump(uint32_t iteration) {
  std::stringstream out;
  out << this->dht_file_out_dir << "/iter_" << std::setfill('0')
      << std::setw(this->file_pad) << iteration << ".dht";
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
poet::ChemistryModule::WorkerRunWorkPackage(
    std::vector<double> &vecWP, std::vector<std::uint32_t> &vecMapping,
    double dSimTime, double dTimestep) {
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

  IRM_RESULT result;
  this->PhreeqcRM::setPOETMapping(vecMapping);
  this->setDumpedField(vecWP);

  this->PhreeqcRM::SetTime(dSimTime);
  this->PhreeqcRM::SetTimeStep(dTimestep);

  result = this->PhreeqcRM::RunCells();

  this->getDumpedField(vecWP);

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
  MPI_Comm worker_comm = dht->getCommunicator();
  int worker_rank;
  MPI_Comm_rank(worker_comm, &worker_rank);

  MPI_Comm &group_comm = this->group_comm;

  auto reduce_and_send = [&worker_rank, &worker_comm, &group_comm](
                             std::vector<std::uint32_t> &send_buffer, int tag) {
    std::vector<uint32_t> to_master(send_buffer.size());
    MPI_Reduce(send_buffer.data(), to_master.data(), send_buffer.size(),
               MPI_UINT32_T, MPI_SUM, 0, worker_comm);

    if (worker_rank == 0) {
      MPI_Send(to_master.data(), to_master.size(), MPI_UINT32_T, 0, tag,
               group_comm);
    }
  };

  switch (type) {
  case WORKER_DHT_HITS: {
    reduce_and_send(dht_hits, WORKER_DHT_HITS);
    break;
  }
  case WORKER_DHT_EVICTIONS: {
    reduce_and_send(dht_evictions, WORKER_DHT_EVICTIONS);
    break;
  }
  default: {
    throw std::runtime_error("Unknown perf type in master's message.");
  }
  }
}

} // namespace poet
