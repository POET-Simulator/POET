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

#include "poet/SimParams.hpp"
#include <Rcpp.h>

#include <iostream>
#include <ostream>
#include <string>

#include <poet/ChemSimPar.hpp>
#include <vector>

using namespace poet;
using namespace std;
using namespace Rcpp;

ChemWorker::ChemWorker(SimParams &params, RInside &R_, Grid &grid_,
                       MPI_Comm dht_comm, poet::ChemistryParams chem_args)
    : ChemSim(params, R_, grid_, chem_args) {
  t_simparams tmp = params.getNumParams();

  this->dt_differ = tmp.dt_differ;
  this->dht_enabled = tmp.dht_enabled;
  this->dht_size_per_process = tmp.dht_size_per_process;
  this->dht_snaps = tmp.dht_snaps;

  this->dht_file = params.getDHTFile();

  mpi_buffer = (double *)calloc(
      (wp_size * (this->prop_names.size())) + BUFFER_OFFSET, sizeof(double));
  mpi_buffer_results =
      (double *)calloc(wp_size * (this->prop_names.size()), sizeof(double));

  if (world_rank == 1)
    cout << "CPP: Worker: DHT usage is " << (dht_enabled ? "ON" : "OFF")
         << endl;

  if (dht_enabled) {
    int data_size = this->prop_names.size() * sizeof(double);
    int key_size =
        this->prop_names.size() * sizeof(double) + (dt_differ * sizeof(double));
    int dht_buckets_per_process =
        dht_size_per_process / (1 + data_size + key_size);

    if (world_rank == 1)
      cout << "CPP: Worker: data size: " << data_size << " bytes" << endl
           << "CPP: Worker: key size: " << key_size << " bytes" << endl
           << "CPP: Worker: buckets per process " << dht_buckets_per_process
           << endl;

    dht = new DHT_Wrapper(params, dht_comm, dht_buckets_per_process, data_size,
                          key_size);

    if (world_rank == 1)
      cout << "CPP: Worker: DHT created!" << endl;

    if (!dht_file.empty())
      readFile();
    // set size
    dht_flags.resize(wp_size, true);
    // assign all elements to true (default)
    dht_flags.assign(wp_size, true);
  }

  timing[0] = 0.0;
  timing[1] = 0.0;
  timing[2] = 0.0;
}

ChemWorker::~ChemWorker() {
  free(mpi_buffer);
  free(mpi_buffer_results);
  if (dht_enabled)
    delete dht;

  delete this->dht;
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

void ChemWorker::doWork(MPI_Status &probe_status) {
  int count;
  int local_work_package_size = 0;

  static int counter = 1;

  double dht_get_start, dht_get_end;
  double phreeqc_time_start, phreeqc_time_end;
  double dht_fill_start, dht_fill_end;

  double dt;

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
  if (mpi_buffer[count] != local_work_package_size) { // work_package_size
    local_work_package_size = mpi_buffer[count];
    R["work_package_size"] = local_work_package_size;
    R.parseEvalQ("mysetup$work_package_size <- work_package_size");
  }

  // current iteration of simulation
  if (mpi_buffer[count + 1] != iteration) {
    iteration = mpi_buffer[count + 1];
    R["iter"] = iteration;
    R.parseEvalQ("mysetup$iter <- iter");
  }

  // current timestep size
  if (mpi_buffer[count + 2] != dt) {
    dt = mpi_buffer[count + 2];
    R["dt"] = dt;
    R.parseEvalQ("mysetup$dt <- dt");
  }

  // current simulation time ('age' of simulation)
  if (mpi_buffer[count + 3] != current_sim_time) {
    current_sim_time = mpi_buffer[count + 3];
    R["simulation_time"] = current_sim_time;
    R.parseEvalQ("mysetup$simulation_time <- simulation_time");
  }
  /* 4th double value is currently a placeholder */
  // if (mpi_buffer[count+4] != placeholder) {
  //     placeholder = mpi_buffer[count+4];
  //     R["mysetup$placeholder"] = placeholder;
  // }

  if (dht_enabled) {
    /* resize helper vector dht_flags of work_package_size changes */
    if ((int)dht_flags.size() != local_work_package_size) {
      dht_flags.resize(local_work_package_size, true); // set size
      dht_flags.assign(local_work_package_size,
                       true); // assign all elements to true (default)
    }

    /* check for values in DHT */
    dht_get_start = MPI_Wtime();
    dht->checkDHT(local_work_package_size, dht_flags, mpi_buffer, dt);
    dht_get_end = MPI_Wtime();

    /* distribute dht_flags to R Runtime */
    R["dht_flags"] = as<LogicalVector>(wrap(dht_flags));
  }

  /* Convert grid to R runtime */
  size_t rowCount = local_work_package_size;
  size_t colCount = this->prop_names.size();

  std::vector<std::vector<double>> input(colCount);

  for (size_t i = 0; i < rowCount; i++) {
    for (size_t j = 0; j < colCount; j++) {
      input[j].push_back(mpi_buffer[i * colCount + j]);
    }
  }

  R["work_package_full"] = Rcpp::as<Rcpp::DataFrame>(Rcpp::wrap(input));

  if (dht_enabled) {
    R.parseEvalQ("work_package <- work_package_full[dht_flags,]");
  } else {
    R.parseEvalQ("work_package <- work_package_full");
  }

  R.parseEvalQ("work_package <- as.matrix(work_package)");

  unsigned int nrows = R.parseEval("nrow(work_package)");

  if (nrows > 0) {
    /*Single Line error Workaround*/
    if (nrows <= 1) {
      // duplicate line to enable correct simmulation
      R.parseEvalQ("work_package <- work_package[rep(1:nrow(work_package), "
                   "times = 2), ]");
    }
    /* Run PHREEQC */
    phreeqc_time_start = MPI_Wtime();
    R.parseEvalQ("result <- as.data.frame(slave_chemistry(setup=mysetup, "
                 "data = work_package))");
    phreeqc_time_end = MPI_Wtime();
  } else {
    // undefined behaviour, isn't it?
  }

  phreeqc_count++;

  if (dht_enabled) {
    R.parseEvalQ("result_full <- work_package_full");
    if (nrows > 0)
      R.parseEvalQ("result_full[dht_flags,] <- result");
  } else {
    R.parseEvalQ("result_full <- result");
  }

  /* convert grid to C domain */

  std::vector<std::vector<double>> output =
      Rcpp::as<std::vector<std::vector<double>>>(R.parseEval("result_full"));

  for (size_t i = 0; i < rowCount; i++) {
    for (size_t j = 0; j < colCount; j++) {
      mpi_buffer_results[i * colCount + j] = output[j][i];
    }
  }

  // grid.exportWP(mpi_buffer_results);
  /* send results to master */
  MPI_Request send_req;
  MPI_Isend(mpi_buffer_results, count, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD,
            &send_req);

  if (dht_enabled) {
    /* write results to DHT */
    dht_fill_start = MPI_Wtime();
    dht->fillDHT(local_work_package_size, dht_flags, mpi_buffer,
                 mpi_buffer_results, dt);
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
  out << out_dir << "/iter_" << setfill('0') << setw(3) << iteration << ".dht";
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
