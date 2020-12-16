#include "worker.h"
#include "dht_wrapper.h"
#include "global_buffer.h"
#include "model/Grid.h"
#include "util/RRuntime.h"
#include <Rcpp.h>
#include <iostream>
#include <mpi.h>

using namespace poet;
using namespace Rcpp;

void worker_function(RRuntime &R, Grid &grid) {
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  MPI_Status probe_status;
  int count;

  int local_work_package_size;
  int iteration;
  double dt, current_sim_time;

  double idle_a, idle_b;
  double cummul_idle = 0.f;

  double dht_get_start = 0, dht_get_end = 0;
  double dht_fill_start = 0, dht_fill_end = 0;
  double phreeqc_time_start = 0, phreeqc_time_end = 0;
  int phreeqc_count = 0;

  // timing[0] -> phreeqc
  // timing[1] -> dht_get
  // timing[2] -> dht_fill
  double timing[3];
  timing[0] = 0.0;
  timing[1] = 0.0;
  timing[2] = 0.0;

  // dht_perf[0] -> hits
  // dht_perf[1] -> miss
  // dht_perf[2] -> collisions
  uint64_t dht_perf[3];

  if (dht_enabled) {
    dht_flags.resize(work_package_size, true); // set size
    dht_flags.assign(work_package_size,
                     true); // assign all elements to true (default)
    dht_hits = 0;
    dht_miss = 0;
    dht_collision = 0;

    // MDL: This code has now been moved to kin.cpp
    // /*Load significance vector from R setup file (or set default)*/
    // bool signif_vector_exists = R.parseEval("exists('signif_vector')");
    // if (signif_vector_exists)
    // {
    //     dht_significant_digits_vector =
    //     as<std::vector<int>>(R["signif_vector"]);
    // } else
    // {
    //     dht_significant_digits_vector.assign(dht_object->key_size /
    //     sizeof(double), dht_significant_digits);
    // }

    // /*Load property type vector from R setup file (or set default)*/
    // bool prop_type_vector_exists = R.parseEval("exists('prop_type')");
    // if (prop_type_vector_exists)
    // {
    //     prop_type_vector = as<std::vector<string>>(R["prop_type"]);
    // } else
    // {
    //     prop_type_vector.assign(dht_object->key_size / sizeof(double),
    //     "normal");
    // }
  }

  // initialization of helper variables
  iteration = 0;
  dt = 0;
  current_sim_time = 0;
  local_work_package_size = 0;

  /*worker loop*/
  while (1) {
    /*Wait for Message*/
    idle_a = MPI_Wtime();
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status);
    idle_b = MPI_Wtime();

    if (probe_status.MPI_TAG == TAG_WORK) { /* do work */

      cummul_idle += idle_b - idle_a;

      /* get number of doubles sent */
      MPI_Get_count(&probe_status, MPI_DOUBLE, &count);

      /* receive */
      MPI_Recv(mpi_buffer, count, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      // decrement count of work_package by BUFFER_OFFSET
      count -= BUFFER_OFFSET;
      // check for changes on all additional variables given by the 'header' of
      // mpi_buffer
      if (mpi_buffer[count] != local_work_package_size) { // work_package_size
        local_work_package_size = mpi_buffer[count];
        R["work_package_size"] = local_work_package_size;
        R.parseEvalQ("mysetup$work_package_size <- work_package_size");
      }
      if (mpi_buffer[count + 1] !=
          iteration) { // current iteration of simulation
        iteration = mpi_buffer[count + 1];
        R["iter"] = iteration;
        R.parseEvalQ("mysetup$iter <- iter");
      }
      if (mpi_buffer[count + 2] != dt) { // current timestep size
        dt = mpi_buffer[count + 2];
        R["dt"] = dt;
        R.parseEvalQ("mysetup$dt <- dt");
      }
      if (mpi_buffer[count + 3] !=
          current_sim_time) { // current simulation time ('age' of simulation)
        current_sim_time = mpi_buffer[count + 3];
        R["simulation_time"] = current_sim_time;
        R.parseEvalQ("mysetup$simulation_time <- simulation_time");
      }
      /* 4th double value is currently a placeholder */
      // if (mpi_buffer[count+4] != placeholder) {
      //     placeholder = mpi_buffer[count+4];
      //     R["mysetup$placeholder"] = placeholder;
      // }

      /* get df with right structure to fill in work package */
      // R.parseEvalQ("skeleton <- head(mysetup$state_C, work_package_size)");
      // R["skeleton"] = grid.buildDataFrame(work_package_size);
      //// R.parseEval("print(rownames(tmp2)[1:5])");
      //// R.parseEval("print(head(tmp2, 2))");
      //// R.parseEvalQ("tmp2$id <- as.double(rownames(tmp2))");

      ////Rcpp::DataFrame buffer = R.parseEval("tmp2");
      // R.setBufferDataFrame("skeleton");

      if (dht_enabled) {
        // DEBUG
        // cout << "RANK " << world_rank << " start checking DHT\n";

        // resize helper vector dht_flags of work_package_size changes
        if ((int)dht_flags.size() != local_work_package_size) {
          dht_flags.resize(local_work_package_size, true); // set size
          dht_flags.assign(local_work_package_size,
                           true); // assign all elements to true (default)
        }

        dht_get_start = MPI_Wtime();
        check_dht(local_work_package_size, dht_flags, mpi_buffer, dt);
        dht_get_end = MPI_Wtime();

        // DEBUG
        // cout << "RANK " << world_rank << " checking DHT complete \n";

        R["dht_flags"] = as<LogicalVector>(wrap(dht_flags));
        // R.parseEvalQ("print(head(dht_flags))");
      }


      /* work */
      // R.from_C_domain(mpi_buffer);
      ////convert_C_buffer_2_R_Dataframe(mpi_buffer, buffer);
      // R["work_package_full"] = R.getBufferDataFrame();
      // R["work_package"] = buffer;

      // DEBUG
      // R.parseEvalQ("print(head(work_package_full))");
      // R.parseEvalQ("print( c(length(dht_flags), nrow(work_package_full)) )");

      grid.importWP(mpi_buffer, work_package_size);

      if (dht_enabled) {
        R.parseEvalQ("work_package <- work_package_full[dht_flags,]");
      } else {
        R.parseEvalQ("work_package <- work_package_full");
      }

      // DEBUG
      // R.parseEvalQ("print(head(work_package),2)");

      // R.parseEvalQ("rownames(work_package) <- work_package$id");
      // R.parseEval("print(paste('id %in% colnames(work_package)', 'id' %in%
      // colnames(work_package)"); R.parseEvalQ("id_store <-
      // rownames(work_package)"); //"[, ncol(work_package)]");
      // R.parseEvalQ("work_package$id <- NULL");
      R.parseEvalQ("work_package <- as.matrix(work_package)");

      unsigned int nrows = R.parseEval("nrow(work_package)");

      if (nrows > 0) {
        /*Single Line error Workaround*/
        if (nrows <= 1) {
          // duplicate line to enable correct simmulation
          R.parseEvalQ("work_package <- work_package[rep(1:nrow(work_package), "
                       "times = 2), ]");
        }

        phreeqc_count++;

        phreeqc_time_start = MPI_Wtime();
        // MDL
        // R.parseEvalQ("print('Work_package:\n'); print(head(work_package ,
        // 2)); cat('RCpp: worker_function:', local_rank, ' \n')");
        R.parseEvalQ("result <- as.data.frame(slave_chemistry(setup=mysetup, "
                     "data = work_package))");
        phreeqc_time_end = MPI_Wtime();
        // R.parseEvalQ("result$id <- id_store");
      } else {
        // cout << "Work-Package is empty, skipping phreeqc!" << endl;
      }

      if (dht_enabled) {
        R.parseEvalQ("result_full <- work_package_full");
        if (nrows > 0)
          R.parseEvalQ("result_full[dht_flags,] <- result");
      } else {
        R.parseEvalQ("result_full <- result");
      }

      // R.setBufferDataFrame("result_full");
      ////Rcpp::DataFrame result = R.parseEval("result_full");
      ////convert_R_Dataframe_2_C_buffer(mpi_buffer_results, result);
      // R.to_C_domain(mpi_buffer_results);

      grid.exportWP(mpi_buffer_results);
      /* send results to master */
      MPI_Request send_req;
      MPI_Isend(mpi_buffer_results, count, MPI_DOUBLE, 0, TAG_WORK,
                MPI_COMM_WORLD, &send_req);

      if (dht_enabled) {
        dht_fill_start = MPI_Wtime();
        fill_dht(local_work_package_size, dht_flags, mpi_buffer,
                 mpi_buffer_results, dt);
        dht_fill_end = MPI_Wtime();

        timing[1] += dht_get_end - dht_get_start;
        timing[2] += dht_fill_end - dht_fill_start;
      }

      timing[0] += phreeqc_time_end - phreeqc_time_start;

      MPI_Wait(&send_req, MPI_STATUS_IGNORE);


    } else if (probe_status.MPI_TAG == TAG_FINISH) { /* recv and die */
      /* before death, submit profiling/timings to master*/

      MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_FINISH, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);

      // timings
      MPI_Send(timing, 3, MPI_DOUBLE, 0, TAG_TIMING, MPI_COMM_WORLD);

      MPI_Send(&phreeqc_count, 1, MPI_INT, 0, TAG_TIMING, MPI_COMM_WORLD);
      MPI_Send(&cummul_idle, 1, MPI_DOUBLE, 0, TAG_TIMING, MPI_COMM_WORLD);

      if (dht_enabled) {
        // dht_perf
        dht_perf[0] = dht_hits;
        dht_perf[1] = dht_miss;
        dht_perf[2] = dht_collision;
        MPI_Send(dht_perf, 3, MPI_UNSIGNED_LONG_LONG, 0, TAG_DHT_PERF,
                 MPI_COMM_WORLD);
      }
      break;

    } else if ((probe_status.MPI_TAG == TAG_DHT_STATS)) {
      MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_DHT_STATS, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      print_statistics();
      MPI_Barrier(MPI_COMM_WORLD);

    } else if ((probe_status.MPI_TAG == TAG_DHT_STORE)) {
      char *outdir;
      MPI_Get_count(&probe_status, MPI_CHAR, &count);
      outdir = (char *)calloc(count + 1, sizeof(char));
      MPI_Recv(outdir, count, MPI_CHAR, 0, TAG_DHT_STORE, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      int res = table_to_file((char *)outdir);
      if (res != DHT_SUCCESS) {
        if (world_rank == 2)
          cerr << "CPP: Worker: Error in writing current state of DHT to file "
                  "(TAG_DHT_STORE)"
               << endl;
      } else {
        if (world_rank == 2)
          cout << "CPP: Worker: Successfully written DHT to file " << outdir
               << endl;
      }
      free(outdir);
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }
}
