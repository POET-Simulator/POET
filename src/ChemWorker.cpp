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

#include "poet/ChemSimPar.hpp"
#include "poet/DHT_Wrapper.hpp"
#include "poet/SimParams.hpp"
#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <mpi.h>
#include <ostream>
#include <string>
#include <vector>

using namespace poet;
using namespace std;
using namespace Rcpp;

ChemWorker::ChemWorker(SimParams &params, RInside &R_, Grid &grid_,
                       MPI_Comm dht_comm)
    : BaseChemModule(params, R_, grid_) {
  t_simparams tmp = params.getNumParams();
  this->wp_size = tmp.wp_size;
  this->out_dir = params.getOutDir();

  this->dt_differ = tmp.dt_differ;
  this->dht_enabled = tmp.dht_enabled;
  this->dht_size_per_process = tmp.dht_size_per_process;
  this->dht_snaps = tmp.dht_snaps;

  this->dht_file = params.getDHTFile();

  this->mpi_buffer =
      new double[(this->wp_size * this->prop_names.size()) + BUFFER_OFFSET];
  this->mpi_buffer_results = new double[wp_size * this->prop_names.size()];

  if (world_rank == 1)
    cout << "CPP: Worker: DHT usage is " << (dht_enabled ? "ON" : "OFF")
         << endl;

  if (this->dht_enabled) {

    uint32_t iKeyCount = this->prop_names.size() + (dt_differ);
    uint32_t iDataCount = this->prop_names.size();

    if (world_rank == 1)
      cout << "CPP: Worker: data count: " << iDataCount << " entries" << endl
           << "CPP: Worker: key count: " << iKeyCount << " entries" << endl
           << "CPP: Worker: memory per process "
           << params.getNumParams().dht_size_per_process / std::pow(10, 6)
           << " MByte" << endl;

    dht = new DHT_Wrapper(params, dht_comm,
                          params.getNumParams().dht_size_per_process, iKeyCount,
                          iDataCount);

    if (world_rank == 1)
      cout << "CPP: Worker: DHT created!" << endl;

    if (!dht_file.empty())
      readFile();
    // set size
    dht_flags.resize(wp_size, true);
  }

  this->timing.fill(0.0);
}

ChemWorker::~ChemWorker() {
  delete this->mpi_buffer;
  delete this->mpi_buffer_results;
  if (dht_enabled)
    delete this->dht;

  if (this->phreeqc_rm) {
    delete this->phreeqc_rm;
  }
}

void ChemWorker::loop() {
  MPI_Status probe_status;
  while (1) {
    double idle_a = MPI_Wtime();
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status);
    double idle_b = MPI_Wtime();

    /* there is a work package to receive */
    if (probe_status.MPI_TAG == TAG_WORK) {
      idle_t += idle_b - idle_a;
      doWork(probe_status);
    }
    /* end of iteration */
    else if (probe_status.MPI_TAG == TAG_DHT_ITER) {
      postIter();
    }
    /* end of simulation */
    else if (probe_status.MPI_TAG == TAG_FINISH) {
      finishWork();
      break;
    }
  }
}

void poet::ChemWorker::InitModule(poet::ChemistryParams &chem_params) {
  this->phreeqc_rm = new PhreeqcWrapper(this->wp_size);

  this->phreeqc_rm->SetupAndLoadDB(chem_params);
  this->phreeqc_rm->InitFromFile(chem_params.input_script);
}

void ChemWorker::doWork(MPI_Status &probe_status) {
  int count;
  int local_work_package_size = 0;

  static int counter = 1;

  double dht_get_start, dht_get_end;
  double phreeqc_time_start, phreeqc_time_end;
  double dht_fill_start, dht_fill_end;

  double dt;

  bool bNoPhreeqc = false;

  /* get number of doubles to be received */
  MPI_Get_count(&probe_status, MPI_DOUBLE, &count);

  /* receive */
  MPI_Recv(mpi_buffer, count, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);

  /* decrement count of work_package by BUFFER_OFFSET */
  count -= BUFFER_OFFSET;

  /* check for changes on all additional variables given by the 'header' of
   * mpi_buffer */

  // work_package_size
  local_work_package_size = mpi_buffer[count];

  // current iteration of simulation
  this->iteration = mpi_buffer[count + 1];

  // current timestep size
  dt = mpi_buffer[count + 2];

  // current simulation time ('age' of simulation)
  current_sim_time = mpi_buffer[count + 3];

  /* 4th double value is currently a placeholder */
  // placeholder = mpi_buffer[count+4];

  std::vector<double> vecCurrWP(
      mpi_buffer,
      mpi_buffer + (local_work_package_size * this->prop_names.size()));
  std::vector<int32_t> vecMappingWP(this->wp_size);

  DHT_ResultObject DHT_Results;

  for (uint32_t i = 0; i < local_work_package_size; i++) {
    vecMappingWP[i] = i;
  }

  if (local_work_package_size != this->wp_size) {
    std::vector<double> vecFiller(
        (this->wp_size - local_work_package_size) * this->prop_names.size(), 0);
    vecCurrWP.insert(vecCurrWP.end(), vecFiller.begin(), vecFiller.end());

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
  this->phreeqc_rm->RunWorkPackage(vecCurrWP, vecMappingWP, current_sim_time,
                                   dt);
  phreeqc_time_end = MPI_Wtime();

  if (dht_enabled) {
    DHT_Results.ResultsToWP(vecCurrWP);
  }

  /* send results to master */
  MPI_Request send_req;
  MPI_Isend(vecCurrWP.data(), count, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD,
            &send_req);

  if (dht_enabled) {
    /* write results to DHT */
    dht_fill_start = MPI_Wtime();
    dht->fillDHT(local_work_package_size, DHT_Results, vecCurrWP);
    dht_fill_end = MPI_Wtime();

    timing[1] += dht_get_end - dht_get_start;
    timing[2] += dht_fill_end - dht_fill_start;
  }

  timing[0] += phreeqc_time_end - phreeqc_time_start;

  MPI_Wait(&send_req, MPI_STATUS_IGNORE);
}

void ChemWorker::postIter() {
  MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_DHT_ITER, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);
  if (dht_enabled) {
    dht->printStatistics();

    if (dht_snaps == 2) {
      writeFile();
    }
  }
  // synchronize all processes
  MPI_Barrier(MPI_COMM_WORLD);
}

void ChemWorker::writeFile() {
  cout.flush();
  std::stringstream out;
  out << out_dir << "/iter_" << setfill('0') << setw(3) << this->iteration
      << ".dht";
  int res = dht->tableToFile(out.str().c_str());
  if (res != DHT_SUCCESS && world_rank == 2)
    cerr << "CPP: Worker: Error in writing current state of DHT to file."
         << endl;
  else if (world_rank == 2)
    cout << "CPP: Worker: Successfully written DHT to file " << out.str()
         << endl;
}

void ChemWorker::readFile() {
  int res = dht->fileToTable((char *)dht_file.c_str());
  if (res != DHT_SUCCESS) {
    if (res == DHT_WRONG_FILE) {
      if (world_rank == 1)
        cerr << "CPP: Worker: Wrong file layout! Continue with empty DHT ..."
             << endl;
    } else {
      if (world_rank == 1)
        cerr << "CPP: Worker: Error in loading current state of DHT from "
                "file. Continue with empty DHT ..."
             << endl;
    }
  } else {
    if (world_rank == 2)
      cout << "CPP: Worker: Successfully loaded state of DHT from file "
           << dht_file << endl;
    std::cout.flush();
  }
}

void ChemWorker::finishWork() {
  /* before death, submit profiling/timings to master*/
  MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_FINISH, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);

  // timings
  MPI_Send(timing.data(), timing.size(), MPI_DOUBLE, 0, TAG_TIMING,
           MPI_COMM_WORLD);

  MPI_Send(&phreeqc_count, 1, MPI_INT, 0, TAG_TIMING, MPI_COMM_WORLD);
  MPI_Send(&idle_t, 1, MPI_DOUBLE, 0, TAG_TIMING, MPI_COMM_WORLD);

  if (dht_enabled) {
    // dht_perf
    int dht_perf[3];
    dht_perf[0] = dht->getHits();
    dht_perf[1] = dht->getMisses();
    dht_perf[2] = dht->getEvictions();
    MPI_Send(dht_perf, 3, MPI_INT, 0, TAG_DHT_PERF, MPI_COMM_WORLD);
  }

  if (dht_enabled && dht_snaps > 0)
    writeFile();
}
