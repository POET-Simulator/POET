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

#include "poet/DiffusionModule.hpp"
#include <Rcpp.h>

#include <bits/stdint-uintn.h>
#include <iostream>

#include <poet/ChemSim.hpp>
#include <poet/Grid.hpp>
#include <string>
#include <vector>

using namespace poet;
using namespace std;
using namespace Rcpp;

ChemMaster::ChemMaster(SimParams &params, RInside &R_, Grid &grid_)
    : ChemSim(params, R_, grid_) {
  t_simparams tmp = params.getNumParams();

  this->wp_size = tmp.wp_size;
  this->dht_enabled = tmp.dht_enabled;

  this->out_dir = params.getOutDir();

  uint32_t grid_size = grid.getTotalCellCount() * this->prop_names.size();

  /* allocate memory */
  workerlist = (worker_struct *)calloc(world_size - 1, sizeof(worker_struct));
  send_buffer = (double *)calloc(
      (wp_size * this->prop_names.size()) + BUFFER_OFFSET, sizeof(double));
  mpi_buffer = (double *)calloc(grid_size, sizeof(double));

  /* calculate distribution of work packages */
  uint32_t mod_pkgs = grid_size % this->wp_size;
  uint32_t n_packages = (uint32_t)(grid.getTotalCellCount() / this->wp_size) +
                        (mod_pkgs != 0 ? 1 : 0);

  this->wp_sizes_vector =
      std::vector<uint32_t>(n_packages - mod_pkgs, this->wp_size);
  for (uint32_t i = 0; i < mod_pkgs; i++) {
    this->wp_sizes_vector.push_back(this->wp_size - 1);
  }
}

ChemMaster::~ChemMaster() {
  free(mpi_buffer);
  free(workerlist);
}

void ChemMaster::simulate(double dt) {
  /* declare most of the needed variables here */
  double chem_a, chem_b;
  double seq_a, seq_b, seq_c, seq_d;
  double worker_chemistry_a, worker_chemistry_b;
  double sim_e_chemistry, sim_f_chemistry;
  int pkg_to_send, pkg_to_recv;
  int free_workers;
  int i_pkgs;

  /* start time measurement of whole chemistry simulation */
  chem_a = MPI_Wtime();

  /* start time measurement of sequential part */
  seq_a = MPI_Wtime();

  std::vector<double> &field = this->state->mem;

  for (uint32_t i = 0; i < this->prop_names.size(); i++) {
    try {
      std::vector<double> t_prop_vec = this->grid.getSpeciesByName(
          this->prop_names[i], poet::DIFFUSION_MODULE_NAME);

      std::copy(t_prop_vec.begin(), t_prop_vec.end(),
                field.begin() + (i * this->n_cells_per_prop));
    } catch (...) {
      continue;
    }
  }

  // HACK: transfer the field into R data structure serving as input for phreeqc
  R["TMP_T"] = field;

  R.parseEvalQ("mysetup$state_T <- setNames(data.frame(matrix(TMP_T, "
               "ncol=length(mysetup$grid$props), nrow=" +
               std::to_string(this->n_cells_per_prop) +
               ")), mysetup$grid$props)");

  /* shuffle grid */
  // grid.shuffleAndExport(mpi_buffer);
  this->shuffleField(field, this->n_cells_per_prop, this->prop_names.size(),
                     mpi_buffer);

  /* retrieve needed data from R runtime */
  iteration = (int)R.parseEval("mysetup$iter");
  // dt = (double)R.parseEval("mysetup$requested_dt");
  current_sim_time = (double)R.parseEval("mysetup$simulation_time") - dt;

  /* setup local variables */
  pkg_to_send = wp_sizes_vector.size();
  pkg_to_recv = wp_sizes_vector.size();
  work_pointer = mpi_buffer;
  free_workers = world_size - 1;
  i_pkgs = 0;

  /* end time measurement of sequential part */
  seq_b = MPI_Wtime();
  seq_t += seq_b - seq_a;

  /* start time measurement of chemistry time needed for send/recv loop */
  worker_chemistry_a = MPI_Wtime();

  /* start send/recv loop */
  // while there are still packages to recv
  while (pkg_to_recv > 0) {
    // print a progressbar to stdout
    printProgressbar((int)i_pkgs, (int)wp_sizes_vector.size());
    // while there are still packages to send
    if (pkg_to_send > 0) {
      // send packages to all free workers ...
      sendPkgs(pkg_to_send, i_pkgs, free_workers);
    }
    // ... and try to receive them from workers who has finished their work
    recvPkgs(pkg_to_recv, pkg_to_send > 0, free_workers);
  }

  // Just to complete the progressbar
  cout << endl;

  /* stop time measurement of chemistry time needed for send/recv loop */
  worker_chemistry_b = MPI_Wtime();
  worker_t = worker_chemistry_b - worker_chemistry_a;

  /* start time measurement of sequential part */
  seq_c = MPI_Wtime();

  /* unshuffle grid */
  // grid.importAndUnshuffle(mpi_buffer);
  this->unshuffleField(mpi_buffer, this->n_cells_per_prop,
                       this->prop_names.size(), field);

  /* do master stuff */

  /* start time measurement of master chemistry */
  sim_e_chemistry = MPI_Wtime();

  // HACK: We don't need to call master_chemistry here since our result is
  // already written to the memory as a data frame

  // R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");
  R["TMP_T"] = Rcpp::wrap(field);
  R.parseEval(std::string("mysetup$state_C <- setNames(data.frame(matrix(TMP_T, nrow=" +
                          to_string(this->n_cells_per_prop) +
                          ")), mysetup$grid$props)"));

  /* end time measurement of master chemistry */
  sim_f_chemistry = MPI_Wtime();
  chem_master += sim_f_chemistry - sim_e_chemistry;

  /* end time measurement of sequential part */
  seq_d = MPI_Wtime();
  seq_t += seq_d - seq_c;

  /* end time measurement of whole chemistry simulation */
  chem_b = MPI_Wtime();
  chem_t += chem_b - chem_a;

  /* advise workers to end chemistry iteration */
  for (int i = 1; i < world_size; i++) {
    MPI_Send(NULL, 0, MPI_DOUBLE, i, TAG_DHT_ITER, MPI_COMM_WORLD);
  }
}

void ChemMaster::sendPkgs(int &pkg_to_send, int &count_pkgs,
                          int &free_workers) {
  /* declare variables */
  double master_send_a, master_send_b;
  int local_work_package_size;
  int end_of_wp;

  /* start time measurement */
  master_send_a = MPI_Wtime();

  /* search for free workers and send work */
  for (int p = 0; p < world_size - 1; p++) {
    if (workerlist[p].has_work == 0 && pkg_to_send > 0) /* worker is free */ {
      /* to enable different work_package_size, set local copy of
       * work_package_size to pre-calculated work package size vector */

      local_work_package_size = (int)wp_sizes_vector[count_pkgs];
      count_pkgs++;

      /* note current processed work package in workerlist */
      workerlist[p].send_addr = work_pointer;

      /* push work pointer to next work package */
      end_of_wp = local_work_package_size * grid.getSpeciesCount();
      work_pointer = &(work_pointer[end_of_wp]);

      // fill send buffer starting with work_package ...
      std::memcpy(send_buffer, workerlist[p].send_addr,
                  (end_of_wp) * sizeof(double));
      // followed by: work_package_size
      send_buffer[end_of_wp] = (double)local_work_package_size;
      // current iteration of simulation
      send_buffer[end_of_wp + 1] = (double)iteration;
      // size of timestep in seconds
      send_buffer[end_of_wp + 2] = dt;
      // current time of simulation (age) in seconds
      send_buffer[end_of_wp + 3] = current_sim_time;
      // placeholder for work_package_count
      send_buffer[end_of_wp + 4] = 0.;

      /* ATTENTION Worker p has rank p+1 */
      MPI_Send(send_buffer, end_of_wp + BUFFER_OFFSET, MPI_DOUBLE, p + 1,
               TAG_WORK, MPI_COMM_WORLD);

      /* Mark that worker has work to do */
      workerlist[p].has_work = 1;
      free_workers--;
      pkg_to_send -= 1;
    }
  }
  master_send_b = MPI_Wtime();
  send_t += master_send_b - master_send_a;
}

void ChemMaster::recvPkgs(int &pkg_to_recv, bool to_send, int &free_workers) {
  /* declare most of the variables here */
  int need_to_receive = 1;
  double master_recv_a, master_recv_b;
  double idle_a, idle_b;
  int p, size;

  MPI_Status probe_status;
  master_recv_a = MPI_Wtime();
  /* start to loop as long there are packages to recv and the need to receive
   */
  while (need_to_receive && pkg_to_recv > 0) {
    // only of there are still packages to send and free workers are available
    if (to_send && free_workers > 0)
      // non blocking probing
      MPI_Iprobe(MPI_ANY_SOURCE, TAG_WORK, MPI_COMM_WORLD, &need_to_receive,
                 &probe_status);
    else {
      idle_a = MPI_Wtime();
      // blocking probing
      MPI_Probe(MPI_ANY_SOURCE, TAG_WORK, MPI_COMM_WORLD, &probe_status);
      idle_b = MPI_Wtime();
      master_idle += idle_b - idle_a;
    }

    /* if need_to_receive was set to true above, so there is a message to
     * receive */
    if (need_to_receive) {
      p = probe_status.MPI_SOURCE;
      MPI_Get_count(&probe_status, MPI_DOUBLE, &size);
      MPI_Recv(workerlist[p - 1].send_addr, size, MPI_DOUBLE, p, TAG_WORK,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      workerlist[p - 1].has_work = 0;
      pkg_to_recv -= 1;
      free_workers++;
    }
  }
  master_recv_b = MPI_Wtime();
  recv_t += master_recv_b - master_recv_a;
}

void ChemMaster::printProgressbar(int count_pkgs, int n_wp, int barWidth) {
  /* visual progress */
  double progress = (float)(count_pkgs + 1) / n_wp;

  cout << "[";
  int pos = barWidth * progress;
  for (int iprog = 0; iprog < barWidth; ++iprog) {
    if (iprog < pos)
      cout << "=";
    else if (iprog == pos)
      cout << ">";
    else
      cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
  /* end visual progress */
}

void ChemMaster::end() {
  /* call end() from base class */
  ChemSim::end();

  /* now we get to the part of the master */
  double *timings;
  int *dht_perfs;

  Rcpp::NumericVector phreeqc_time;
  Rcpp::NumericVector dht_get_time;
  Rcpp::NumericVector dht_fill_time;
  Rcpp::IntegerVector phreeqc_counts;
  Rcpp::NumericVector idle_worker;

  int phreeqc_tmp;

  timings = (double *)calloc(3, sizeof(double));

  int dht_hits = 0;
  int dht_miss = 0;
  int dht_evictions = 0;

  if (dht_enabled) {
    dht_perfs = (int *)calloc(3, sizeof(int));
  }

  double idle_worker_tmp;

  /* loop over all workers *
   * ATTENTION Worker p has rank p+1 */
  for (int p = 0; p < world_size - 1; p++) {
    /* Send termination message to worker */
    MPI_Send(NULL, 0, MPI_DOUBLE, p + 1, TAG_FINISH, MPI_COMM_WORLD);

    /* ... and receive all timings and metrics from each worker */
    MPI_Recv(timings, 3, MPI_DOUBLE, p + 1, TAG_TIMING, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    phreeqc_time.push_back(timings[0], "w" + to_string(p + 1));

    MPI_Recv(&phreeqc_tmp, 1, MPI_INT, p + 1, TAG_TIMING, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    phreeqc_counts.push_back(phreeqc_tmp, "w" + to_string(p + 1));

    MPI_Recv(&idle_worker_tmp, 1, MPI_DOUBLE, p + 1, TAG_TIMING, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    idle_worker.push_back(idle_worker_tmp, "w" + to_string(p + 1));

    if (dht_enabled) {
      dht_get_time.push_back(timings[1], "w" + to_string(p + 1));
      dht_fill_time.push_back(timings[2], "w" + to_string(p + 1));

      MPI_Recv(dht_perfs, 3, MPI_INT, p + 1, TAG_DHT_PERF, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      dht_hits += dht_perfs[0];
      dht_miss += dht_perfs[1];
      dht_evictions += dht_perfs[2];
    }
  }

  /* distribute all data to the R runtime */
  R["simtime_chemistry"] = chem_t;
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
  R["simtime_workers"] = worker_t;
  R.parseEvalQ("profiling$simtime_workers <- simtime_workers");
  R["simtime_chemistry_master"] = chem_master;
  R.parseEvalQ(
      "profiling$simtime_chemistry_master <- simtime_chemistry_master");

  R["seq_master"] = seq_t;
  R.parseEvalQ("profiling$seq_master <- seq_master");

  R["idle_master"] = master_idle;
  R.parseEvalQ("profiling$idle_master <- idle_master");
  R["idle_worker"] = idle_worker;
  R.parseEvalQ("profiling$idle_worker <- idle_worker");

  R["phreeqc_time"] = phreeqc_time;
  R.parseEvalQ("profiling$phreeqc <- phreeqc_time");

  R["phreeqc_count"] = phreeqc_counts;
  R.parseEvalQ("profiling$phreeqc_count <- phreeqc_count");

  if (dht_enabled) {
    R["dht_hits"] = dht_hits;
    R.parseEvalQ("profiling$dht_hits <- dht_hits");
    R["dht_miss"] = dht_miss;
    R.parseEvalQ("profiling$dht_miss <- dht_miss");
    R["dht_evictions"] = dht_evictions;
    R.parseEvalQ("profiling$dht_evictions <- dht_evictions");
    R["dht_get_time"] = dht_get_time;
    R.parseEvalQ("profiling$dht_get_time <- dht_get_time");
    R["dht_fill_time"] = dht_fill_time;
    R.parseEvalQ("profiling$dht_fill_time <- dht_fill_time");
  }

  /* do some cleanup */
  free(timings);

  if (dht_enabled)
    free(dht_perfs);
}

void ChemMaster::shuffleField(const std::vector<double> &in_field,
                              uint32_t size_per_prop, uint32_t prop_count,
                              double *out_buffer) {
  uint32_t wp_count = this->wp_sizes_vector.size();
  uint32_t write_i = 0;
  for (uint32_t i = 0; i < wp_count; i++) {
    for (uint32_t j = i; j < size_per_prop; j += wp_count) {
      for (uint32_t k = 0; k < prop_count; k++) {
        out_buffer[(write_i * prop_count) + k] =
            in_field[(k * size_per_prop) + j];
      }
      write_i++;
    }
  }
}

void ChemMaster::unshuffleField(const double *in_buffer, uint32_t size_per_prop,
                                uint32_t prop_count,
                                std::vector<double> &out_field) {
  uint32_t wp_count = this->wp_sizes_vector.size();
  uint32_t read_i = 0;

  for (uint32_t i = 0; i < wp_count; i++) {
    for (uint32_t j = i; j < size_per_prop; j += wp_count) {
      for (uint32_t k = 0; k < prop_count; k++) {
        out_field[(k * size_per_prop) + j] =
            in_buffer[(read_i * prop_count) + k];
      }
      read_i++;
    }
  }
}

double ChemMaster::getSendTime() { return this->send_t; }

double ChemMaster::getRecvTime() { return this->recv_t; }

double ChemMaster::getIdleTime() { return this->master_idle; }

double ChemMaster::getWorkerTime() { return this->worker_t; }

double ChemMaster::getChemMasterTime() { return this->chem_master; }

double ChemMaster::getSeqTime() { return this->seq_t; }
