#include <Rcpp.h>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

// #include "DHT.h"  // MPI-DHT Implementation
// #include "argh.h"  // Argument handler https://github.com/adishavit/argh
// BSD-licenced #include "dht_wrapper.h" #include "global_buffer.h"
#include <ChemSim.h>
#include <Grid.h>
#include <Parser.h>
#include <RRuntime.h>
#include <SimParams.h>
#include <TransportSim.h>
// #include "worker.h"

//#define DHT_SIZE_PER_PROCESS 1073741824

using namespace std;
using namespace poet;
using namespace Rcpp;

// double *mpi_buffer;
// double *mpi_buffer_results;

// uint32_t work_package_size;
// #define WORK_PACKAGE_SIZE_DEFAULT 5

// bool store_result;

// std::set<std::string> paramList() {
//   std::set<std::string> options;
//   // global
//   options.insert("work-package-size");
//   // only DHT
//   options.insert("dht-signif");
//   options.insert("dht-strategy");
//   options.insert("dht-size");
//   options.insert("dht-snaps");
//   options.insert("dht-file");

//   return options;
// }

// std::set<std::string> flagList() {
//   std::set<std::string> options;
//   // global
//   options.insert("ignore-result");
//   // only DHT
//   options.insert("dht");
//   options.insert("dht-log");

//   return options;
// }

// std::list<std::string> checkOptions(argh::parser cmdl) {
//   std::list<std::string> retList;
//   std::set<std::string> flist = flagList();
//   std::set<std::string> plist = paramList();

//   for (auto &flag : cmdl.flags()) {
//     if (!(flist.find(flag) != flist.end())) retList.push_back(flag);
//   }

//   for (auto &param : cmdl.params()) {
//     if (!(plist.find(param.first) != plist.end()))
//       retList.push_back(param.first);
//   }

//   return retList;
// }

// typedef struct {
//   char has_work;
//   double *send_addr;
// } worker_struct;

int main(int argc, char *argv[]) {
  double sim_start, sim_b_transport, sim_a_transport, sim_b_chemistry,
      sim_a_chemistry, sim_end;

  double cummul_transport = 0.f;
  double cummul_chemistry = 0.f;
  double cummul_master_seq = 0.f;

  // argh::parser cmdl(argv);
  // int dht_significant_digits;
  // cout << "CPP: Start Init (MPI)" << endl;

  t_simparams params;
  int world_size, world_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /*Create custom Communicator with all processes except 0 (the master) for DHT
   * storage*/
  MPI_Comm dht_comm;

  if (world_rank == 0) {
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &dht_comm);
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &dht_comm);
  }
  // cout << "Done";

  // if (cmdl[{"help", "h"}]) {
  //   if (params.world_rank == 0) {
  //     cout << "Todo" << endl
  //          << "See README.md for further information." << endl;
  //   }
  //   MPI_Finalize();
  //   return EXIT_SUCCESS;
  // }

  // /*INIT is now done separately in an R file provided here as argument!*/
  // if (!cmdl(2)) {
  //   if (params.world_rank == 0) {
  //     cerr << "ERROR. Kin needs 2 positional arguments: " << endl
  //          << "1) the R script defining your simulation and" << endl
  //          << "2) the directory prefix where to save results and profiling"
  //          << endl;
  //   }
  //   MPI_Finalize();
  //   return EXIT_FAILURE;
  // }

  // std::list<std::string> optionsError = checkOptions(cmdl);
  // if (!optionsError.empty()) {
  //   if (params.world_rank == 0) {
  //     cerr << "Unrecognized option(s):\n" << endl;
  //     for (auto option : optionsError) {
  //       cerr << option << endl;
  //     }
  //     cerr << "\nMake sure to use available options. Exiting!" << endl;
  //   }
  //   MPI_Finalize();
  //   return EXIT_FAILURE;
  // }

  // /*Parse DHT arguments*/
  // params.dht_enabled = cmdl["dht"];
  // // cout << "CPP: DHT is " << ( dht_enabled ? "ON" : "OFF" ) << '\n';

  // if (params.dht_enabled) {
  //   cmdl("dht-strategy", 0) >> params.dht_strategy;
  //   // cout << "CPP: DHT strategy is " << dht_strategy << endl;

  //   cmdl("dht-signif", 5) >> dht_significant_digits;
  //   // cout << "CPP: DHT significant digits = " << dht_significant_digits <<
  //   // endl;

  //   params.dht_log = !(cmdl["dht-nolog"]);
  //   // cout << "CPP: DHT logarithm before rounding: " << ( dht_logarithm ?
  //   "ON"
  //   // : "OFF" ) << endl;

  //   cmdl("dht-size", DHT_SIZE_PER_PROCESS) >> params.dht_size_per_process;
  //   // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process
  //   <<
  //   // endl;

  //   cmdl("dht-snaps", 0) >> params.dht_snaps;

  //   cmdl("dht-file") >> params.dht_file;
  // }

  // /*Parse work package size*/
  // cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> params.wp_size;

  // /*Parse output options*/
  // store_result = !cmdl["ignore-result"];

  RRuntime R(argc, argv);
  cout << "CPP: R Init (RInside) on process " << world_rank << endl;

  // if local_rank == 0 then master else worker
  // R["local_rank"] = params.world_rank;

  /*Loading Dependencies*/
  std::string r_load_dependencies =
      "suppressMessages(library(Rmufits));"
      "suppressMessages(library(RedModRphree));"
      "source('kin_r_library.R');"
      "source('parallel_r_library.R');";
  R.parseEvalQ(r_load_dependencies);

  Parser parser(argv, world_rank, world_size);
  int pret = parser.parseCmdl();

  if (pret == PARSER_ERROR) {
    MPI_Finalize();
    return EXIT_FAILURE;
  } else if (pret == PARSER_HELP) {
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  parser.parseR(R);
  params = parser.getParams();

  // if (params.world_rank == 0) {
  //   cout << "CPP: Complete results storage is " << (params.store_result ?
  //   "ON" : "OFF")
  //        << endl;
  //   cout << "CPP: Work Package Size: " << params.wp_size << endl;
  //   cout << "CPP: DHT is " << (params.dht_enabled ? "ON" : "OFF") << '\n';

  //   if (params.dht_enabled) {
  //     cout << "CPP: DHT strategy is " << params.dht_strategy << endl;
  //     // cout << "CPP: DHT key default digits (ignored if 'signif_vector' is
  //     "
  //     //         "defined) = "
  //     //      << dht_significant_digits << endl;
  //     cout << "CPP: DHT logarithm before rounding: "
  //          << (params.dht_log ? "ON" : "OFF") << endl;
  //     cout << "CPP: DHT size per process (Byte) = "
  //          << params.dht_size_per_process << endl;
  //     cout << "CPP: DHT save snapshots is " << params.dht_snaps << endl;
  //     cout << "CPP: DHT load file is " << params.dht_file << endl;
  //   }
  // }

  // std::string filesim;
  // cmdl(1) >> filesim;            // <- first positional argument
  // R["filesim"] = wrap(filesim);  // assign a char* (string) to 'filesim'
  // R.parseEvalQ(
  //     "source(filesim)");  // eval the init string, ignoring any returns

  // only rank 0 initializes goes through the whole initialization
  if (world_rank == 0) {
    // cmdl(2) >> params.out_dir;  // <- second positional argument
    // R["fileout"] =
    //     wrap(params.out_dir);  // assign a char* (string) to 'fileout'

    // Note: R::sim_init() checks if the directory already exists,
    // if not it makes it

    // pass the boolean "store_result" to the R process
    // R["store_result"] = store_result;

    // get timestep vector from grid_init function ...
    std::string master_init_code = "mysetup <- master_init(setup=setup)";
    R.parseEval(master_init_code);

    params.dt_differ = R.parseEval("mysetup$dt_differ");
    MPI_Bcast(&(params.dt_differ), 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  } else {  // workers will only read the setup DataFrame defined by input file
    R.parseEval("mysetup <- setup");
    MPI_Bcast(&(params.dt_differ), 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  }

  if (world_rank == 0) {
    cout << "CPP: R init done on process with rank " << params.world_rank
         << endl;
  }

  // initialize chemistry on all processes
  std::string init_chemistry_code = "mysetup <- init_chemistry(setup=mysetup)";
  R.parseEval(init_chemistry_code);

  Grid grid(R);
  // params.grid = &grid;
  grid.init();
  /* Retrieve state_C from R context for MPI buffer generation */
  // Rcpp::DataFrame state_C = R.parseEval("mysetup$state_C");

  /* Init Parallel helper functions */
  // R["n_procs"] = params.world_size - 1; /* worker count */
  // R["work_package_size"] = params.wp_size;

  // Removed additional field for ID in previous versions
  // if (params.world_rank == 0) {
  //   mpi_buffer =
  //       (double *)calloc(grid.getRows() * grid.getCols(), sizeof(double));
  // } else {
  //   mpi_buffer = (double *)calloc(
  //       (params.wp_size * (grid.getCols())) + BUFFER_OFFSET, sizeof(double));
  //   mpi_buffer_results =
  //       (double *)calloc(params.wp_size * (grid.getCols()), sizeof(double));
  // }

  // if (params.world_rank == 0) {
  //   cout << "CPP: parallel init completed (buffers allocated)!" << endl;
  // }

  // MDL: pass to R the DHT stuff (basically, only for storing of
  // simulation parameters). These 2 variables are always defined:
  // R["dht_enabled"] = params.dht_enabled;
  // R["dht_log"] = params.dht_log;

  // params.R = &R;

  if (params.dht_enabled) {
    // cout << "\nCreating DHT\n";
    // determine size of dht entries
    // int dht_data_size = grid.getCols() * sizeof(double);
    // int dht_key_size =
    //     grid.getCols() * sizeof(double) + (params.dt_differ *
    //     sizeof(double));

    // // // determine bucket count for preset memory usage
    // // // bucket size is key + value + 1 byte for status
    // int dht_buckets_per_process =
    //     params.dht_size_per_process / (1 + dht_data_size + dht_key_size);

    // MDL : following code moved here from worker.cpp
    /*Load significance vector from R setup file (or set default)*/
    bool signif_vector_exists = R.parseEval("exists('signif_vector')");
    if (signif_vector_exists) {
      params.dht_signif_vector = as<std::vector<int>>(R["signif_vector"]);
    } else {
      params.dht_signif_vector.assign(grid.getCols(), params.dht_significant_digits);
    }

    /*Load property type vector from R setup file (or set default)*/
    bool prop_type_vector_exists = R.parseEval("exists('prop_type')");
    if (prop_type_vector_exists) {
      params.dht_prop_type_vector = as<std::vector<string>>(R["prop_type"]);
    } else {
      params.dht_prop_type_vector.assign(grid.getCols(), "act");
    }

    if (params.world_rank == 0) {
      // // print only on master, values are equal on all workes
      // cout << "CPP: dht_data_size: " << dht_data_size << "\n";
      // cout << "CPP: dht_key_size: " << dht_key_size << "\n";
      // cout << "CPP: dht_buckets_per_process: " << dht_buckets_per_process
      //      << endl;

      // MDL: new output on signif_vector and prop_type
      if (signif_vector_exists) {
        cout << "CPP: using problem-specific rounding digits: " << endl;
        R.parseEval(
            "print(data.frame(prop=prop, type=prop_type, "
            "digits=signif_vector))");
      } else {
        cout << "CPP: using DHT default rounding digits = "
             << params.dht_significant_digits << endl;
      }

      // MDL: pass to R the DHT stuff. These variables exist
      // only if dht_enabled is true
      R["dht_final_signif"] = params.dht_signif_vector;
      R["dht_final_proptype"] = params.dht_prop_type_vector;
    }

    // if (params.dht_strategy == 0) {
    //   if (params.world_rank != 0) {
    //     dht_object = DHT_create(dht_comm, dht_buckets_per_process,
    //                             dht_data_size, dht_key_size, get_md5);

    //     // storing for access from worker and callback functions
    //     fuzzing_buffer = (double *)malloc(dht_key_size);
    //   }
    // } else {
    //   dht_object = DHT_create(MPI_COMM_WORLD, dht_buckets_per_process,
    //                           dht_data_size, dht_key_size, get_md5);
    // }

    // if (params.world_rank == 0) {
    //   cout << "CPP: DHT successfully created!" << endl;
    // }
  }

  // MDL: store all parameters
  if (params.world_rank == 0) {
    cout << "CPP: Calling R Function to store calling parameters" << endl;
    R.parseEvalQ("StoreSetup(setup=mysetup)");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (params.world_rank == 0) { /* This is executed by the master */
    ChemMaster master(&params, R, grid);
    TransportSim trans(R);

    Rcpp::NumericVector master_send;
    Rcpp::NumericVector master_recv;

    double *timings;
    uint64_t *dht_perfs = NULL;

    // a temporary send buffer
    sim_start = MPI_Wtime();

    // Iteration Count is dynamic, retrieving value from R (is only needed by
    // master for the following loop)
    uint32_t maxiter = R.parseEval("mysetup$maxiter");

    /*SIMULATION LOOP*/
    for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
      cout << "CPP: Evaluating next time step" << endl;
      R.parseEvalQ("mysetup <- master_iteration_setup(mysetup)");

      /*displaying iteration number, with C++ and R iterator*/
      cout << "CPP: Going through iteration " << iter << endl;
      cout << "CPP: R's $iter: " << ((uint32_t)(R.parseEval("mysetup$iter")))
           << ". Iteration" << endl;

      cout << "CPP: Calling Advection" << endl;

      trans.run();
      // sim_b_transport = MPI_Wtime();
      // R.parseEvalQ("mysetup <- master_advection(setup=mysetup)");
      // sim_a_transport = MPI_Wtime();

      // if (iter == 1) master.prepareSimulation();

      cout << "CPP: Chemistry" << endl;
      /*Fallback for sequential execution*/

      if (params.world_size == 1) {
        master.ChemSim::run();
      } else {
        master.run();
      }

      // MDL master_iteration_end just writes on disk state_T and
      // state_C after every iteration if the cmdline option
      // --ignore-results is not given (and thus the R variable
      // store_result is TRUE)
      R.parseEvalQ("mysetup <- master_iteration_end(setup=mysetup)");

      // cummul_transport += trans.getTransportTime();
      // cummul_chemistry += master.getChemistryTime();

      cout << endl
           << "CPP: End of *coupling* iteration " << iter << "/" << maxiter
           << endl
           << endl;

      // master_send.push_back(master.getSendTime(), "it_" + to_string(iter));
      // master_recv.push_back(master.getRecvTime(), "it_" + to_string(iter));

      for (int i = 1; i < params.world_size; i++) {
        MPI_Send(NULL, 0, MPI_DOUBLE, i, TAG_DHT_ITER, MPI_COMM_WORLD);
      }

      MPI_Barrier(MPI_COMM_WORLD);

    }  // END SIMULATION LOOP

    cout << "CPP: finished simulation loop" << endl;

    sim_end = MPI_Wtime();

    cout << "CPP: start timing profiling" << endl;

    R.parseEvalQ("profiling <- list()");

    R["simtime"] = sim_end - sim_start;
    R.parseEvalQ("profiling$simtime <- simtime");

    trans.end();

    if (params.world_size == 1) {
      master.ChemSim::end();
    } else {
      master.end();
    }

    // Rcpp::NumericVector phreeqc_time;
    // Rcpp::NumericVector dht_get_time;
    // Rcpp::NumericVector dht_fill_time;
    // Rcpp::IntegerVector phreeqc_counts;
    // Rcpp::NumericVector idle_worker;

    // int phreeqc_tmp;

    // timings = (double *)calloc(3, sizeof(double));

    // int dht_hits = 0;
    // int dht_miss = 0;
    // int dht_collision = 0;

    // if (params.dht_enabled) {
    //   dht_hits = 0;
    //   dht_miss = 0;
    //   dht_collision = 0;
    //   dht_perfs = (uint64_t *)calloc(3, sizeof(uint64_t));
    // }

    // double idle_worker_tmp;

    // for (int p = 0; p < params.world_size - 1; p++) {
    //   /* ATTENTION Worker p has rank p+1 */
    //   /* Send termination message to worker */
    //   MPI_Send(NULL, 0, MPI_DOUBLE, p + 1, TAG_FINISH, MPI_COMM_WORLD);

    //   MPI_Recv(timings, 3, MPI_DOUBLE, p + 1, TAG_TIMING, MPI_COMM_WORLD,
    //            MPI_STATUS_IGNORE);
    //   phreeqc_time.push_back(timings[0], "w" + to_string(p + 1));

    //   MPI_Recv(&phreeqc_tmp, 1, MPI_INT, p + 1, TAG_TIMING, MPI_COMM_WORLD,
    //            MPI_STATUS_IGNORE);
    //   phreeqc_counts.push_back(phreeqc_tmp, "w" + to_string(p + 1));

    //   MPI_Recv(&idle_worker_tmp, 1, MPI_DOUBLE, p + 1, TAG_TIMING,
    //            MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //   idle_worker.push_back(idle_worker_tmp, "w" + to_string(p + 1));

    //   if (params.dht_enabled) {
    //     dht_get_time.push_back(timings[1], "w" + to_string(p + 1));
    //     dht_fill_time.push_back(timings[2], "w" + to_string(p + 1));

    //     MPI_Recv(dht_perfs, 3, MPI_UNSIGNED_LONG_LONG, p + 1, TAG_DHT_PERF,
    //              MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //     dht_hits += dht_perfs[0];
    //     dht_miss += dht_perfs[1];
    //     dht_collision += dht_perfs[2];
    //   }
    // }

    // R.parseEvalQ("profiling <- list()");

    // R["simtime"] = sim_end - sim_start;
    // R.parseEvalQ("profiling$simtime <- simtime");
    // R["simtime_transport"] = cummul_transport;
    // R.parseEvalQ("profiling$simtime_transport <- simtime_transport");
    // R["simtime_chemistry"] = cummul_chemistry;
    // R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
    // R["simtime_workers"] = master.getWorkerTime();
    // R.parseEvalQ("profiling$simtime_workers <- simtime_workers");
    // R["simtime_chemistry_master"] = master.getChemMasterTime();
    // R.parseEvalQ(
    //     "profiling$simtime_chemistry_master <- simtime_chemistry_master");

    // R["seq_master"] = cummul_master_seq;
    // R.parseEvalQ("profiling$seq_master <- seq_master");

    // // R["master_send"] = master_send;
    // // R.parseEvalQ("profiling$master_send <- master_send");
    // // R["master_recv"] = master_recv;
    // // R.parseEvalQ("profiling$master_recv <- master_recv");

    // R["idle_master"] = master.getIdleTime();
    // R.parseEvalQ("profiling$idle_master <- idle_master");
    // R["idle_worker"] = idle_worker;
    // R.parseEvalQ("profiling$idle_worker <- idle_worker");

    // R["phreeqc_time"] = phreeqc_time;
    // R.parseEvalQ("profiling$phreeqc <- phreeqc_time");

    // R["phreeqc_count"] = phreeqc_counts;
    // R.parseEvalQ("profiling$phreeqc_count <- phreeqc_count");

    // if (params.dht_enabled) {
    //   R["dht_hits"] = dht_hits;
    //   R.parseEvalQ("profiling$dht_hits <- dht_hits");
    //   R["dht_miss"] = dht_miss;
    //   R.parseEvalQ("profiling$dht_miss <- dht_miss");
    //   R["dht_collision"] = dht_collision;
    //   R.parseEvalQ("profiling$dht_collisions <- dht_collision");
    //   R["dht_get_time"] = dht_get_time;
    //   R.parseEvalQ("profiling$dht_get_time <- dht_get_time");
    //   R["dht_fill_time"] = dht_fill_time;
    //   R.parseEvalQ("profiling$dht_fill_time <- dht_fill_time");
    // }

    // free(timings);

    // if (params.dht_enabled) free(dht_perfs);

    string r_vis_code;
    r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
    R.parseEval(r_vis_code);

    cout << "CPP: Done! Results are stored as R objects into <"
         << params.out_dir << "/timings.rds>" << endl;
    /*exporting results and profiling data*/

    // std::string r_vis_code;
    // r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
    // R.parseEval(r_vis_code);
  } else { /*This is executed by the workers*/
    ChemWorker worker(&params, R, grid, dht_comm);
    // worker.prepareSimulation(dht_comm);
    worker.loop();
  }

  cout << "CPP: finished, cleanup of process " << params.world_rank << endl;

  // if (params.dht_enabled) {
  //   if (params.dht_strategy == 0) {
  //     if (params.world_rank != 0) {
  //       DHT_free(dht_object, NULL, NULL);
  //     }
  //   } else {
  //     DHT_free(dht_object, NULL, NULL);
  //   }
  // }

  // free(mpi_buffer);
  MPI_Finalize();

  if (params.world_rank == 0) {
    cout << "CPP: done, bye!" << endl;
  }

  exit(0);
}
