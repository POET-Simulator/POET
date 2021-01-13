#include "Profiler.h"

#include <Rcpp.h>
#include <mpi.h>

#include <string>
#include <iostream>

using namespace Rcpp;
using namespace std;

int poet::Profiler::startProfiling(t_simparams &params, ChemMaster &chem,
                                   TransportSim &trans, RRuntime &R,
                                   double simtime) {
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
  int dht_collision = 0;

  if (params.dht_enabled) {
    dht_hits = 0;
    dht_miss = 0;
    dht_collision = 0;
    dht_perfs = (int *)calloc(3, sizeof(int));
  }

  double idle_worker_tmp;

  for (int p = 0; p < params.world_size - 1; p++) {
    /* ATTENTION Worker p has rank p+1 */
    /* Send termination message to worker */
    MPI_Send(NULL, 0, MPI_DOUBLE, p + 1, TAG_FINISH, MPI_COMM_WORLD);

    MPI_Recv(timings, 3, MPI_DOUBLE, p + 1, TAG_TIMING, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    phreeqc_time.push_back(timings[0], "w" + to_string(p + 1));

    MPI_Recv(&phreeqc_tmp, 1, MPI_INT, p + 1, TAG_TIMING, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    phreeqc_counts.push_back(phreeqc_tmp, "w" + to_string(p + 1));

    MPI_Recv(&idle_worker_tmp, 1, MPI_DOUBLE, p + 1, TAG_TIMING, MPI_COMM_WORLD,
             MPI_STATUS_IGNORE);
    idle_worker.push_back(idle_worker_tmp, "w" + to_string(p + 1));

    if (params.dht_enabled) {
      dht_get_time.push_back(timings[1], "w" + to_string(p + 1));
      dht_fill_time.push_back(timings[2], "w" + to_string(p + 1));

      MPI_Recv(dht_perfs, 3, MPI_INT, p + 1, TAG_DHT_PERF,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      dht_hits += dht_perfs[0];
      dht_miss += dht_perfs[1];
      cout << "profiler miss = " << dht_miss << endl;
      dht_collision += dht_perfs[2];
    }
  }

  R.parseEvalQ("profiling <- list()");

  R["simtime"] = simtime;
  R.parseEvalQ("profiling$simtime <- simtime");
  R["simtime_transport"] = trans.getTransportTime();
  R.parseEvalQ("profiling$simtime_transport <- simtime_transport");
  R["simtime_chemistry"] = chem.getChemistryTime();
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
  R["simtime_workers"] = chem.getWorkerTime();
  R.parseEvalQ("profiling$simtime_workers <- simtime_workers");
  R["simtime_chemistry_master"] = chem.getChemMasterTime();
  R.parseEvalQ(
      "profiling$simtime_chemistry_master <- simtime_chemistry_master");

  R["seq_master"] = chem.getSeqTime();
  R.parseEvalQ("profiling$seq_master <- seq_master");

  // R["master_send"] = master_send;
  // R.parseEvalQ("profiling$master_send <- master_send");
  // R["master_recv"] = master_recv;
  // R.parseEvalQ("profiling$master_recv <- master_recv");

  R["idle_master"] = chem.getIdleTime();
  R.parseEvalQ("profiling$idle_master <- idle_master");
  R["idle_worker"] = idle_worker;
  R.parseEvalQ("profiling$idle_worker <- idle_worker");

  R["phreeqc_time"] = phreeqc_time;
  R.parseEvalQ("profiling$phreeqc <- phreeqc_time");

  R["phreeqc_count"] = phreeqc_counts;
  R.parseEvalQ("profiling$phreeqc_count <- phreeqc_count");

  if (params.dht_enabled) {
    R["dht_hits"] = dht_hits;
    R.parseEvalQ("profiling$dht_hits <- dht_hits");
    R["dht_miss"] = dht_miss;
    R.parseEvalQ("profiling$dht_miss <- dht_miss");
    R["dht_collision"] = dht_collision;
    R.parseEvalQ("profiling$dht_collisions <- dht_collision");
    R["dht_get_time"] = dht_get_time;
    R.parseEvalQ("profiling$dht_get_time <- dht_get_time");
    R["dht_fill_time"] = dht_fill_time;
    R.parseEvalQ("profiling$dht_fill_time <- dht_fill_time");
  }

  free(timings);

  if (params.dht_enabled) free(dht_perfs);

  string r_vis_code;
  r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
  R.parseEval(r_vis_code);

  return 0;
}
