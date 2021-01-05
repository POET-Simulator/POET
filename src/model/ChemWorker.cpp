#include <Rcpp.h>
#include <mpi.h>

#include <iostream>
#include <string>

#include "ChemSim.h"

using namespace poet;
using namespace std;
using namespace Rcpp;

ChemWorker::ChemWorker(t_simparams *params_, RRuntime &R_, Grid &grid_)
    : params(params_), ChemSim(params_, R_, grid_) {
  this->dt_differ = params->dt_differ;
  this->dht_enabled = params->dht_enabled;
  this->dht_size_per_process = params->dht_size_per_process;
  this->dht_file = params->dht_file;
}

void ChemWorker::prepareSimulation(MPI_Comm dht_comm) {
  mpi_buffer = (double *)calloc((wp_size * (grid.getCols())) + BUFFER_OFFSET,
                                sizeof(double));
  mpi_buffer_results =
      (double *)calloc(wp_size * (grid.getCols()), sizeof(double));

  if (world_rank == 1)
    cout << "Worker: DHT usage is " << (dht_enabled ? "ON" : "OFF") << endl;

  if (dht_enabled) {
    int data_size = grid.getCols() * sizeof(double);
    int key_size =
        grid.getCols() * sizeof(double) + (dt_differ * sizeof(double));
    int dht_buckets_per_process =
        dht_size_per_process / (1 + data_size + key_size);

    if (world_rank == 1)
      cout << "CPP: Worker: data size: " << data_size << " bytes" << endl
           << "CPP: Worker: key size: " << key_size << "bytes" << endl
           << "CPP: Worker: buckets per process " << dht_buckets_per_process
           << endl;

    dht = new DHT_Wrapper(params, dht_comm, dht_buckets_per_process, data_size,
                          key_size);

    if (world_rank == 1) cout << "CPP: Worker: DHT created!" << endl;
  }

  if (!dht_file.empty()) readFile();

  // set size
  dht_flags.resize(params->wp_size, true);
  // assign all elements to true (default)
  dht_flags.assign(params->wp_size, true);

  timing[0] = 0.0;
  timing[1] = 0.0;
  timing[2] = 0.0;
}

void ChemWorker::loop() {
  MPI_Status probe_status;
  while (1) {
    double idle_a = MPI_Wtime();
    MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status);
    double idle_b = MPI_Wtime();

    if (probe_status.MPI_TAG == TAG_WORK) {
      idle_t += idle_b - idle_a;
      doWork(probe_status);
    } else if (probe_status.MPI_TAG == TAG_FINISH) {
      finishWork();
      break;
    } else if (probe_status.MPI_TAG == TAG_DHT_ITER) {
      postIter();
    }
  }
}

void ChemWorker::doWork(MPI_Status &probe_status) {
  int count;
  int local_work_package_size = 0;

  double dht_get_start, dht_get_end;
  double phreeqc_time_start, phreeqc_time_end;
  double dht_fill_start, dht_fill_end;

  /* get number of doubles sent */
  MPI_Get_count(&probe_status, MPI_DOUBLE, &count);

  /* receive */
  MPI_Recv(mpi_buffer, count, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);

  // decrement count of work_package by BUFFER_OFFSET
  count -= BUFFER_OFFSET;
  // check for changes on all additional variables given by the 'header' of
  // mpi_buffer

  // work_package_size
  if (mpi_buffer[count] != local_work_package_size) {  // work_package_size
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
      dht_flags.resize(local_work_package_size, true);  // set size
      dht_flags.assign(local_work_package_size,
                       true);  // assign all elements to true (default)
    }

    dht_get_start = MPI_Wtime();
    dht->checkDHT(local_work_package_size, dht_flags, mpi_buffer, dt);
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

  grid.importWP(mpi_buffer, params->wp_size);

  if (params->dht_enabled) {
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
      R.parseEvalQ(
          "work_package <- work_package[rep(1:nrow(work_package), "
          "times = 2), ]");
    }

    phreeqc_count++;

    phreeqc_time_start = MPI_Wtime();
    // MDL
    // R.parseEvalQ("print('Work_package:\n'); print(head(work_package ,
    // 2)); cat('RCpp: worker_function:', local_rank, ' \n')");
    R.parseEvalQ(
        "result <- as.data.frame(slave_chemistry(setup=mysetup, "
        "data = work_package))");
    phreeqc_time_end = MPI_Wtime();
    // R.parseEvalQ("result$id <- id_store");
  } else {
    // cout << "Work-Package is empty, skipping phreeqc!" << endl;
  }

  if (dht_enabled) {
    R.parseEvalQ("result_full <- work_package_full");
    if (nrows > 0) R.parseEvalQ("result_full[dht_flags,] <- result");
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
  MPI_Isend(mpi_buffer_results, count, MPI_DOUBLE, 0, TAG_WORK, MPI_COMM_WORLD,
            &send_req);

  if (dht_enabled) {
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
  std::stringstream out;
  out << out_dir << "/iter_" << setfill('0') << setw(3) << iteration << ".dht";
  int res = dht->tableToFile(out.str().c_str());
  if (res != DHT_SUCCESS && world_rank == 2)
    cerr << "CPP: Worker: Errir in writing current state of DHT to file."
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
  if (dht_enabled && dht_snaps > 0) writeFile();

  double dht_perf[3];
  /* before death, submit profiling/timings to master*/
  MPI_Recv(NULL, 0, MPI_DOUBLE, 0, TAG_FINISH, MPI_COMM_WORLD,
           MPI_STATUS_IGNORE);

  // timings
  MPI_Send(timing, 3, MPI_DOUBLE, 0, TAG_TIMING, MPI_COMM_WORLD);

  MPI_Send(&phreeqc_count, 1, MPI_INT, 0, TAG_TIMING, MPI_COMM_WORLD);
  MPI_Send(&idle_t, 1, MPI_DOUBLE, 0, TAG_TIMING, MPI_COMM_WORLD);

  if (dht_enabled) {
    // dht_perf
    dht_perf[0] = dht->getHits();
    dht_perf[1] = dht->getMisses();
    dht_perf[2] = dht->getEvictions();
    MPI_Send(dht_perf, 3, MPI_UNSIGNED_LONG_LONG, 0, TAG_DHT_PERF,
             MPI_COMM_WORLD);
  }

  free(mpi_buffer);
  free(mpi_buffer_results);
  delete dht;
}
