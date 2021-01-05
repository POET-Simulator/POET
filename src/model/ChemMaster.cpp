#include <Rcpp.h>
#include <mpi.h>

#include <iostream>

#include "ChemSim.h"

using namespace poet;
using namespace std;
using namespace Rcpp;

#define TAG_WORK 42

ChemMaster::ChemMaster(t_simparams *params, RRuntime &R_, Grid &grid_)
    : ChemSim(params, R_, grid_) {
  this->wp_size = params->wp_size;
  this->out_dir = params->out_dir;
}

void ChemMaster::runIteration() {
  double seq_a, seq_b, seq_c, seq_d;
  double worker_chemistry_a, worker_chemistry_b;
  double sim_e_chemistry, sim_f_chemistry;
  int pkg_to_send, pkg_to_recv;
  int free_workers;
  int i_pkgs;

  seq_a = MPI_Wtime();
  grid.shuffleAndExport(mpi_buffer);
  // retrieve data from R runtime
  iteration = (int)R.parseEval("mysetup$iter");
  dt = (double)R.parseEval("mysetup$requested_dt");
  current_sim_time =
      (double)R.parseEval("mysetup$simulation_time-mysetup$requested_dt");

  // setup local variables
  pkg_to_send = wp_sizes_vector.size();
  pkg_to_recv = wp_sizes_vector.size();
  work_pointer = mpi_buffer;
  free_workers = world_size - 1;
  i_pkgs = 0;

  seq_b = MPI_Wtime();
  seq_t += seq_b - seq_a;

  worker_chemistry_a = MPI_Wtime();
  while (pkg_to_recv > 0) {
    // TODO: Progressbar into IO instance.
    printProgressbar((int)i_pkgs, (int)wp_sizes_vector.size());
    if (pkg_to_send > 0) {
      sendPkgs(pkg_to_send, i_pkgs, free_workers);
    }
    recvPkgs(pkg_to_recv, pkg_to_send > 0, free_workers);
  }

  // Just to complete the progressbar
  cout << endl;

  worker_chemistry_b = MPI_Wtime();
  worker_t = worker_chemistry_b - worker_chemistry_a;

  seq_c = MPI_Wtime();
  grid.importAndUnshuffle(mpi_buffer);
  /* do master stuff */
  sim_e_chemistry = MPI_Wtime();
  R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");
  sim_f_chemistry = MPI_Wtime();
  chem_master += sim_f_chemistry - sim_e_chemistry;
  seq_d = MPI_Wtime();
  seq_t += seq_d - seq_c;
}

void ChemMaster::sendPkgs(int &pkg_to_send, int &count_pkgs,
                          int &free_workers) {
  double master_send_a, master_send_b;
  int local_work_package_size;
  int end_of_wp;

  // start time measurement
  master_send_a = MPI_Wtime();
  /*search for free workers and send work*/
  for (int p = 0; p < world_size - 1; p++) {
    if (workerlist[p].has_work == 0 && pkg_to_send > 0) /* worker is free */ {
      // to enable different work_package_size, set local copy of
      // work_package_size to either global work_package size or
      // remaining 'to_send' packages to_send >= work_package_size ?
      // local_work_package_size = work_package_size :
      // local_work_package_size = to_send;

      local_work_package_size = (int)wp_sizes_vector[count_pkgs];
      count_pkgs++;

      // cout << "CPP: sending pkg n. " << count_pkgs << " with size "
      // << local_work_package_size << endl;

      /*push pointer forward to next work package, after taking the
       * current one*/
      workerlist[p].send_addr = work_pointer;

      end_of_wp = local_work_package_size * grid.getCols();
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

      workerlist[p].has_work = 1;
      free_workers--;
      pkg_to_send -= 1;
    }
  }
  master_send_b = MPI_Wtime();
  send_t += master_send_b - master_send_a;
}

void ChemMaster::recvPkgs(int &pkg_to_recv, bool to_send, int &free_workers) {
  int need_to_receive = 1;
  double master_recv_a, master_recv_b;
  double idle_a, idle_b;
  int p, size;

  MPI_Status probe_status;
  master_recv_a = MPI_Wtime();
  while (need_to_receive && pkg_to_recv > 0) {
    if (to_send && free_workers > 0)
      MPI_Iprobe(MPI_ANY_SOURCE, TAG_WORK, MPI_COMM_WORLD, &need_to_receive,
                 &probe_status);
    else {
      idle_a = MPI_Wtime();
      MPI_Probe(MPI_ANY_SOURCE, TAG_WORK, MPI_COMM_WORLD, &probe_status);
      idle_b = MPI_Wtime();
      master_idle += idle_b - idle_a;
    }

    if (need_to_receive) {
      p = probe_status.MPI_SOURCE;
      size;
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

void ChemMaster::prepareSimulation() {
  workerlist = (worker_struct *)calloc(world_size - 1, sizeof(worker_struct));
  send_buffer = (double *)calloc((wp_size * (grid.getCols())) + BUFFER_OFFSET,
                                 sizeof(double));

  R.parseEvalQ(
      "wp_ids <- distribute_work_packages(len=nrow(mysetup$state_T), "
      "package_size=work_package_size)");

  // we only sort once the vector
  R.parseEvalQ("ordered_ids <- order(wp_ids)");
  R.parseEvalQ("wp_sizes_vector <- compute_wp_sizes(wp_ids)");
  R.parseEval("stat_wp_sizes(wp_sizes_vector)");
  wp_sizes_vector = as<std::vector<int>>(R["wp_sizes_vector"]);

  mpi_buffer =
      (double *)calloc(grid.getRows() * grid.getCols(), sizeof(double));
}

void ChemMaster::finishSimulation() {
  free(mpi_buffer);
  free(workerlist);
}

double ChemMaster::getSendTime() { return this->send_t; }

double ChemMaster::getRecvTime() { return this->recv_t; }

double ChemMaster::getIdleTime() { return this->master_idle; }

double ChemMaster::getWorkerTime() { return this->worker_t; }

double ChemMaster::getChemMasterTime() { return this->chem_master; }

double ChemMaster::getSeqTime() { return this->seq_t; }
