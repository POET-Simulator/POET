#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <Rcpp.h>

#include <mpi.h> // mpi header file

#include "DHT.h" // MPI-DHT Implementation
#include "argh.h" // Argument handler https://github.com/adishavit/argh BSD-licenced
#include "dht_wrapper.h"
#include "global_buffer.h"
#include "worker.h"
#include "util/RRuntime.h"

using namespace std;
using namespace poet;
using namespace Rcpp;

double *mpi_buffer;
double *mpi_buffer_results;

uint32_t work_package_size;
#define WORK_PACKAGE_SIZE_DEFAULT 5

bool store_result;

std::set<std::string> paramList() {
  std::set<std::string> options;
  // global
  options.insert("work-package-size");
  // only DHT
  options.insert("dht-signif");
  options.insert("dht-strategy");
  options.insert("dht-size");
  options.insert("dht-snaps");
  options.insert("dht-file");

  return options;
}

std::set<std::string> flagList() {
  std::set<std::string> options;
  // global
  options.insert("ignore-result");
  // only DHT
  options.insert("dht");
  options.insert("dht-log");

  return options;
}

std::list<std::string> checkOptions(argh::parser cmdl) {
  std::list<std::string> retList;
  std::set<std::string> flist = flagList();
  std::set<std::string> plist = paramList();

  for (auto &flag : cmdl.flags()) {
    if (!(flist.find(flag) != flist.end()))
      retList.push_back(flag);
  }

  for (auto &param : cmdl.params()) {
    if (!(plist.find(param.first) != plist.end()))
      retList.push_back(param.first);
  }

  return retList;
}

typedef struct {
  char has_work;
  double *send_addr;
} worker_struct;

int main(int argc, char *argv[]) {
  double sim_start, sim_b_transport, sim_a_transport, sim_b_chemistry,
      sim_a_chemistry, sim_end;

  double cummul_transport = 0.f;
  double cummul_chemistry = 0.f;
  double cummul_workers = 0.f;
  double cummul_chemistry_master = 0.f;

  double cummul_master_seq_pre_loop = 0.f;
  double cummul_master_seq_loop = 0.f;
  double master_idle = 0.f;

  double master_send_a, master_send_b;
  double cummul_master_send = 0.f;
  double master_recv_a, master_recv_b;
  double cummul_master_recv = 0.f;

  double sim_a_seq, sim_b_seq, sim_c_seq, sim_d_seq;
  double idle_a, idle_b;

  double sim_c_chemistry, sim_d_chemistry;
  double sim_e_chemistry, sim_f_chemistry;

  argh::parser cmdl(argv);

  // cout << "CPP: Start Init (MPI)" << endl;

  MPI_Init(&argc, &argv);

  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /*Create custom Communicator with all processes except 0 (the master) for DHT
   * storage*/
  // only needed if strategy == 0, but done anyway
  MPI_Group group_world;
  MPI_Group dht_group;
  MPI_Comm dht_comm;
  int *process_ranks;

  // make a list of processes in the new communicator
  process_ranks = (int *)malloc(world_size * sizeof(int));
  for (int I = 1; I < world_size; I++)
    process_ranks[I - 1] = I;

  // get the group under MPI_COMM_WORLD
  MPI_Comm_group(MPI_COMM_WORLD, &group_world);
  // create the new group
  MPI_Group_incl(group_world, world_size - 1, process_ranks, &dht_group);
  // create the new communicator
  MPI_Comm_create(MPI_COMM_WORLD, dht_group, &dht_comm);
  free(process_ranks); // cleanup
  // cout << "Done";

  if (cmdl[{"help", "h"}]) {
    if (world_rank == 0) {
      cout << "Todo" << endl
           << "See README.md for further information." << endl;
    }
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  /*INIT is now done separately in an R file provided here as argument!*/
  if (!cmdl(2)) {
    if (world_rank == 0) {
      cerr << "ERROR. Kin needs 2 positional arguments: " << endl
           << "1) the R script defining your simulation and" << endl
           << "2) the directory prefix where to save results and profiling"
           << endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  std::list<std::string> optionsError = checkOptions(cmdl);
  if (!optionsError.empty()) {
    if (world_rank == 0) {
      cerr << "Unrecognized option(s):\n" << endl;
      for (auto option : optionsError) {
        cerr << option << endl;
      }
      cerr << "\nMake sure to use available options. Exiting!" << endl;
    }
    MPI_Finalize();
    return EXIT_FAILURE;
  }

  /*Parse DHT arguments*/
  dht_enabled = cmdl["dht"];
  // cout << "CPP: DHT is " << ( dht_enabled ? "ON" : "OFF" ) << '\n';

  if (dht_enabled) {
    cmdl("dht-strategy", 0) >> dht_strategy;
    // cout << "CPP: DHT strategy is " << dht_strategy << endl;

    cmdl("dht-signif", 5) >> dht_significant_digits;
    // cout << "CPP: DHT significant digits = " << dht_significant_digits <<
    // endl;

    dht_logarithm = cmdl["dht-log"];
    // cout << "CPP: DHT logarithm before rounding: " << ( dht_logarithm ? "ON"
    // : "OFF" ) << endl;

    cmdl("dht-size", DHT_SIZE_PER_PROCESS) >> dht_size_per_process;
    // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process <<
    // endl;

    cmdl("dht-snaps", 0) >> dht_snaps;

    cmdl("dht-file") >> dht_file;
  }

  /*Parse work package size*/
  cmdl("work-package-size", WORK_PACKAGE_SIZE_DEFAULT) >> work_package_size;

  /*Parse output options*/
  store_result = !cmdl["ignore-result"];

  if (world_rank == 0) {
    cout << "CPP: Complete results storage is " << (store_result ? "ON" : "OFF")
         << endl;
    cout << "CPP: Work Package Size: " << work_package_size << endl;
    cout << "CPP: DHT is " << (dht_enabled ? "ON" : "OFF") << '\n';

    if (dht_enabled) {
      cout << "CPP: DHT strategy is " << dht_strategy << endl;
      cout << "CPP: DHT key default digits (ignored if 'signif_vector' is "
              "defined) = "
           << dht_significant_digits << endl;
      cout << "CPP: DHT logarithm before rounding: "
           << (dht_logarithm ? "ON" : "OFF") << endl;
      cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process
           << endl;
      cout << "CPP: DHT save snapshots is " << dht_snaps << endl;
      cout << "CPP: DHT load file is " << dht_file << endl;
    }
  }

  cout << "CPP: R Init (RInside) on process " << world_rank << endl;
  RRuntime R(argc, argv);

  // if local_rank == 0 then master else worker
  R["local_rank"] = world_rank;

  /*Loading Dependencies*/
  std::string r_load_dependencies = "suppressMessages(library(Rmufits));"
                                    "suppressMessages(library(RedModRphree));"
                                    "source('kin_r_library.R');"
                                    "source('parallel_r_library.R');";
  R.parseEvalQ(r_load_dependencies);

  std::string filesim;
  cmdl(1) >> filesim;              // <- first positional argument
  R["filesim"] = wrap(filesim);    // assign a char* (string) to 'filesim'
  R.parseEvalQ("source(filesim)"); // eval the init string, ignoring any returns

  std::string out_dir;
  if (world_rank ==
      0) { // only rank 0 initializes goes through the whole initialization
    cmdl(2) >> out_dir;           // <- second positional argument
    R["fileout"] = wrap(out_dir); // assign a char* (string) to 'fileout'

    // Note: R::sim_init() checks if the directory already exists,
    // if not it makes it

    // pass the boolean "store_result" to the R process
    R["store_result"] = store_result;

    // get timestep vector from grid_init function ...
    std::string master_init_code = "mysetup <- master_init(setup=setup)";
    R.parseEval(master_init_code);

    dt_differ = R.parseEval("mysetup$dt_differ");
    MPI_Bcast(&dt_differ, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  } else { // workers will only read the setup DataFrame defined by input file
    R.parseEval("mysetup <- setup");
    MPI_Bcast(&dt_differ, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  }

  if (world_rank == 0) {
    cout << "CPP: R init done on process with rank " << world_rank << endl;
  }

  // initialize chemistry on all processes
  std::string init_chemistry_code = "mysetup <- init_chemistry(setup=mysetup)";
  R.parseEval(init_chemistry_code);

  /* Retrieve state_C from R context for MPI buffer generation */
  Rcpp::DataFrame state_C = R.parseEval("mysetup$state_C");

  /* Init Parallel helper functions */
  R["n_procs"] = world_size - 1; /* worker count */
  R["work_package_size"] = work_package_size;

  // Removed additional field for ID in previous versions
  if (world_rank == 0) {
    mpi_buffer =
        (double *)calloc(state_C.nrow() * (state_C.ncol()), sizeof(double));
  } else {
    mpi_buffer = (double *)calloc(
        (work_package_size * (state_C.ncol())) + BUFFER_OFFSET, sizeof(double));
    mpi_buffer_results =
        (double *)calloc(work_package_size * (state_C.ncol()), sizeof(double));
  }

  if (world_rank == 0) {
    cout << "CPP: parallel init completed (buffers allocated)!" << endl;
  }

  // MDL: pass to R the DHT stuff (basically, only for storing of
  // simulation parameters). These 2 variables are always defined:
  R["dht_enabled"] = dht_enabled;
  R["dht_log"] = dht_logarithm;

  if (dht_enabled) {
    // cout << "\nCreating DHT\n";
    // determine size of dht entries
    int dht_data_size = state_C.ncol() * sizeof(double);
    int dht_key_size =
        state_C.ncol() * sizeof(double) + (dt_differ * sizeof(double));

    // determine bucket count for preset memory usage
    // bucket size is key + value + 1 byte for status
    int dht_buckets_per_process =
        dht_size_per_process / (1 + dht_data_size + dht_key_size);

    // MDL : following code moved here from worker.cpp
    /*Load significance vector from R setup file (or set default)*/
    bool signif_vector_exists = R.parseEval("exists('signif_vector')");
    if (signif_vector_exists) {
      dht_significant_digits_vector = as<std::vector<int>>(R["signif_vector"]);
    } else {
      dht_significant_digits_vector.assign(
          dht_object->key_size / sizeof(double), dht_significant_digits);
    }

    /*Load property type vector from R setup file (or set default)*/
    bool prop_type_vector_exists = R.parseEval("exists('prop_type')");
    if (prop_type_vector_exists) {
      prop_type_vector = as<std::vector<string>>(R["prop_type"]);
    } else {
      prop_type_vector.assign(dht_object->key_size / sizeof(double), "act");
    }

    if (world_rank == 0) {
      // print only on master, values are equal on all workes
      cout << "CPP: dht_data_size: " << dht_data_size << "\n";
      cout << "CPP: dht_key_size: " << dht_key_size << "\n";
      cout << "CPP: dht_buckets_per_process: " << dht_buckets_per_process
           << endl;

      // MDL: new output on signif_vector and prop_type
      if (signif_vector_exists) {
        cout << "CPP: using problem-specific rounding digits: " << endl;
        R.parseEval("print(data.frame(prop=prop, type=prop_type, "
                    "digits=signif_vector))");
      } else {
        cout << "CPP: using DHT default rounding digits = "
             << dht_significant_digits << endl;
      }

      // MDL: pass to R the DHT stuff. These variables exist
      // only if dht_enabled is true
      R["dht_final_signif"] = dht_significant_digits_vector;
      R["dht_final_proptype"] = prop_type_vector;
    }

    if (dht_strategy == 0) {
      if (world_rank != 0) {
        dht_object = DHT_create(dht_comm, dht_buckets_per_process,
                                dht_data_size, dht_key_size, get_md5);

        // storing for access from worker and callback functions
        fuzzing_buffer = (double *)malloc(dht_key_size);
      }
    } else {
      dht_object = DHT_create(MPI_COMM_WORLD, dht_buckets_per_process,
                              dht_data_size, dht_key_size, get_md5);
    }

    if (world_rank == 0) {
      cout << "CPP: DHT successfully created!" << endl;
    }
  }

  // MDL: store all parameters
  if (world_rank == 0) {
    cout << "CPP: Calling R Function to store calling parameters" << endl;
    R.parseEvalQ("StoreSetup(setup=mysetup)");
  }

  MPI_Barrier(MPI_COMM_WORLD);

  if (world_rank == 0) { /* This is executed by the master */

    Rcpp::NumericVector master_send;
    Rcpp::NumericVector master_recv;

    sim_a_seq = MPI_Wtime();

    worker_struct *workerlist =
        (worker_struct *)calloc(world_size - 1, sizeof(worker_struct));
    int need_to_receive;
    MPI_Status probe_status;
    double *timings;
    uint64_t *dht_perfs = NULL;

    int local_work_package_size;

    // a temporary send buffer
    double *send_buffer;
    send_buffer = (double *)calloc(
        (work_package_size * (state_C.ncol())) + BUFFER_OFFSET, sizeof(double));

    // helper variables
    int iteration;
    double dt, current_sim_time;

    int n_wp = 1; // holds the actual number of wp which is
                  // computed later in R::distribute_work_packages()
    std::vector<int> wp_sizes_vector; // vector with the sizes of
                                      // each package

    sim_start = MPI_Wtime();

    // Iteration Count is dynamic, retrieving value from R (is only needed by
    // master for the following loop)
    uint32_t maxiter = R.parseEval("mysetup$maxiter");

    sim_b_seq = MPI_Wtime();

    cummul_master_seq_pre_loop += sim_b_seq - sim_a_seq;

    /*SIMULATION LOOP*/
    for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
      sim_a_seq = MPI_Wtime();

      cummul_master_send = 0.f;
      cummul_master_recv = 0.f;

      cout << "CPP: Evaluating next time step" << endl;
      R.parseEvalQ("mysetup <- master_iteration_setup(mysetup)");

      /*displaying iteration number, with C++ and R iterator*/
      cout << "CPP: Going through iteration " << iter << endl;
      cout << "CPP: R's $iter: " << ((uint32_t)(R.parseEval("mysetup$iter")))
           << ". Iteration" << endl;

      cout << "CPP: Calling Advection" << endl;

      sim_b_transport = MPI_Wtime();
      R.parseEvalQ("mysetup <- master_advection(setup=mysetup)");
      sim_a_transport = MPI_Wtime();

      cout << "CPP: Chemistry" << endl;
      /*Fallback for sequential execution*/
      sim_b_chemistry = MPI_Wtime();

      if (world_size == 1) {
        // MDL : the transformation of values into pH and pe
        // takes now place in master_advection() so the
        // following line is unneeded
        // R.parseEvalQ("mysetup$state_T <-
        // RedModRphree::Act2pH(mysetup$state_T)");
        R.parseEvalQ(
            "result <- slave_chemistry(setup=mysetup, data=mysetup$state_T)");
        R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");
      } else { /*send work to workers*/

        // NEW: only in the first iteration we call
        // R::distribute_work_packages()!!
        if (iter == 1) {
          R.parseEvalQ(
              "wp_ids <- distribute_work_packages(len=nrow(mysetup$state_T), "
              "package_size=work_package_size)");

          // we only sort once the vector
          R.parseEvalQ("ordered_ids <- order(wp_ids)");
          R.parseEvalQ("wp_sizes_vector <- compute_wp_sizes(wp_ids)");
          n_wp = (int)R.parseEval("length(wp_sizes_vector)");
          wp_sizes_vector = as<std::vector<int>>(R["wp_sizes_vector"]);
          cout << "CPP: Total number of work packages: " << n_wp << endl;
          R.parseEval("stat_wp_sizes(wp_sizes_vector)");
        }

        /* shuffle and extract data
           MDL: we now apply :Act2pH directly in master_advection
         */
        // R.parseEval("tmp <-
        // shuffle_field(RedModRphree::Act2pH(mysetup$state_T), ordered_ids)");
        R.parseEval("tmp <- shuffle_field(mysetup$state_T, ordered_ids)");
        R.setBufferDataFrame("tmp");
        R.to_C_domain(mpi_buffer);
        //Rcpp::DataFrame chemistry_data = R.parseEval("tmp");

        //convert_R_Dataframe_2_C_buffer(mpi_buffer, chemistry_data);
        // cout << "CPP: shuffle_field() done" << endl;

        /* send and receive work; this is done by counting
         * the wp */
        int pkg_to_send = n_wp;
        int pkg_to_recv = n_wp;
        size_t colCount = R.getBufferNCol();
        int free_workers = world_size - 1;
        double *work_pointer = mpi_buffer;
        sim_c_chemistry = MPI_Wtime();

        /* visual progress */
        float progress = 0.0;
        int barWidth = 70;

        // retrieve data from R runtime
        iteration = (int)R.parseEval("mysetup$iter");
        dt = (double)R.parseEval("mysetup$requested_dt");
        current_sim_time =
            (double)R.parseEval("mysetup$simulation_time-mysetup$requested_dt");

        int count_pkgs = 0;

        sim_b_seq = MPI_Wtime();

        sim_c_chemistry = MPI_Wtime();

        while (pkg_to_recv > 0) // start dispatching work packages
        {
          /* visual progress */
          progress = (float)(count_pkgs + 1) / n_wp;

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

          if (pkg_to_send > 0) {

            master_send_a = MPI_Wtime();
            /*search for free workers and send work*/
            for (int p = 0; p < world_size - 1; p++) {
              if (workerlist[p].has_work == 0 &&
                  pkg_to_send > 0) /* worker is free */ {

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

                int end_of_wp = local_work_package_size * colCount;
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
                MPI_Send(send_buffer, end_of_wp + BUFFER_OFFSET, MPI_DOUBLE,
                         p + 1, TAG_WORK, MPI_COMM_WORLD);

                workerlist[p].has_work = 1;
                free_workers--;
                pkg_to_send -= 1;
              }
            }
            master_send_b = MPI_Wtime();
            cummul_master_send += master_send_b - master_send_a;
          }

          /*check if there are results to receive and receive them*/
          need_to_receive = 1;
          master_recv_a = MPI_Wtime();
          while (need_to_receive && pkg_to_recv > 0) {

            if (pkg_to_send > 0 && free_workers > 0)
              MPI_Iprobe(MPI_ANY_SOURCE, TAG_WORK, MPI_COMM_WORLD,
                         &need_to_receive, &probe_status);
            else {
              idle_a = MPI_Wtime();
              MPI_Probe(MPI_ANY_SOURCE, TAG_WORK, MPI_COMM_WORLD,
                        &probe_status);
              idle_b = MPI_Wtime();
              master_idle += idle_b - idle_a;
            }

            if (need_to_receive) {
              int p = probe_status.MPI_SOURCE;
              int size;
              MPI_Get_count(&probe_status, MPI_DOUBLE, &size);
              MPI_Recv(workerlist[p - 1].send_addr, size, MPI_DOUBLE, p,
                       TAG_WORK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
              workerlist[p - 1].has_work = 0;
              pkg_to_recv -= 1;
              free_workers++;
            }
          }
          master_recv_b = MPI_Wtime();
          cummul_master_recv += master_recv_b - master_recv_a;
        }

        sim_c_seq = MPI_Wtime();

        // don't overwrite last progress
        cout << endl;

        sim_d_chemistry = MPI_Wtime();
        cummul_workers += sim_d_chemistry - sim_c_chemistry;

        //convert_C_buffer_2_R_Dataframe(mpi_buffer, chemistry_data);
        R.from_C_domain(mpi_buffer);

        R["chemistry_data"] = R.getBufferDataFrame();

        /* unshuffle results */
        R.parseEval("result <- unshuffle_field(chemistry_data, ordered_ids)");

        /* do master stuff */
        sim_e_chemistry = MPI_Wtime();
        R.parseEvalQ("mysetup <- master_chemistry(setup=mysetup, data=result)");
        sim_f_chemistry = MPI_Wtime();
        cummul_chemistry_master += sim_f_chemistry - sim_e_chemistry;
      }
      sim_a_chemistry = MPI_Wtime();

      // MDL master_iteration_end just writes on disk state_T and
      // state_C after every iteration if the cmdline option
      // --ignore-results is not given (and thus the R variable
      // store_result is TRUE)
      R.parseEvalQ("mysetup <- master_iteration_end(setup=mysetup)");

      cummul_transport += sim_a_transport - sim_b_transport;
      cummul_chemistry += sim_a_chemistry - sim_b_chemistry;

      cout << endl
           << "CPP: End of *coupling* iteration " << iter << "/" << maxiter
           << endl
           << endl;

      if (dht_enabled) {
        for (int i = 1; i < world_size; i++) {
          MPI_Send(NULL, 0, MPI_DOUBLE, i, TAG_DHT_STATS, MPI_COMM_WORLD);
        }

        // MPI_Barrier(MPI_COMM_WORLD);

        if (dht_snaps == 2) {
          std::stringstream outfile;
          outfile << out_dir << "/iter_" << std::setfill('0') << std::setw(3)
                  << iter << ".dht";
          for (int i = 1; i < world_size; i++) {
            MPI_Send(outfile.str().c_str(), outfile.str().size(), MPI_CHAR, i,
                     TAG_DHT_STORE, MPI_COMM_WORLD);
          }
          MPI_Barrier(MPI_COMM_WORLD);
        }
      }

      sim_d_seq = MPI_Wtime();

      cummul_master_seq_loop +=
          ((sim_b_seq - sim_a_seq) - (sim_a_transport - sim_b_transport)) +
          (sim_d_seq - sim_c_seq);
      master_send.push_back(cummul_master_send, "it_" + to_string(iter));
      master_recv.push_back(cummul_master_recv, "it_" + to_string(iter));

    } // END SIMULATION LOOP

    sim_end = MPI_Wtime();

    if (dht_enabled && dht_snaps > 0) {
      cout << "CPP: Master: Instruct workers to write DHT to file ..." << endl;
      std::string outfile;
      outfile = out_dir + ".dht";
      for (int i = 1; i < world_size; i++) {
        MPI_Send(outfile.c_str(), outfile.size(), MPI_CHAR, i, TAG_DHT_STORE,
                 MPI_COMM_WORLD);
      }
      MPI_Barrier(MPI_COMM_WORLD);
      cout << "CPP: Master: ... done" << endl;
    }

    Rcpp::NumericVector phreeqc_time;
    Rcpp::NumericVector dht_get_time;
    Rcpp::NumericVector dht_fill_time;
    Rcpp::IntegerVector phreeqc_counts;
    Rcpp::NumericVector idle_worker;

    int phreeqc_tmp;

    timings = (double *)calloc(3, sizeof(double));

    if (dht_enabled) {
      dht_hits = 0;
      dht_miss = 0;
      dht_collision = 0;
      dht_perfs = (uint64_t *)calloc(3, sizeof(uint64_t));
    }

    double idle_worker_tmp;

    for (int p = 0; p < world_size - 1; p++) {
      /* ATTENTION Worker p has rank p+1 */
      /* Send termination message to worker */
      MPI_Send(NULL, 0, MPI_DOUBLE, p + 1, TAG_FINISH, MPI_COMM_WORLD);

      MPI_Recv(timings, 3, MPI_DOUBLE, p + 1, TAG_TIMING, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      phreeqc_time.push_back(timings[0], "w" + to_string(p + 1));

      MPI_Recv(&phreeqc_tmp, 1, MPI_INT, p + 1, TAG_TIMING, MPI_COMM_WORLD,
               MPI_STATUS_IGNORE);
      phreeqc_counts.push_back(phreeqc_tmp, "w" + to_string(p + 1));

      MPI_Recv(&idle_worker_tmp, 1, MPI_DOUBLE, p + 1, TAG_TIMING,
               MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      idle_worker.push_back(idle_worker_tmp, "w" + to_string(p + 1));

      if (dht_enabled) {
        dht_get_time.push_back(timings[1], "w" + to_string(p + 1));
        dht_fill_time.push_back(timings[2], "w" + to_string(p + 1));

        MPI_Recv(dht_perfs, 3, MPI_UNSIGNED_LONG_LONG, p + 1, TAG_DHT_PERF,
                 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        dht_hits += dht_perfs[0];
        dht_miss += dht_perfs[1];
        dht_collision += dht_perfs[2];
      }
    }

    R.parseEvalQ("profiling <- list()");

    R["simtime"] = sim_end - sim_start;
    R.parseEvalQ("profiling$simtime <- simtime");
    R["simtime_transport"] = cummul_transport;
    R.parseEvalQ("profiling$simtime_transport <- simtime_transport");
    R["simtime_chemistry"] = cummul_chemistry;
    R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");
    R["simtime_workers"] = cummul_workers;
    R.parseEvalQ("profiling$simtime_workers <- simtime_workers");
    R["simtime_chemistry_master"] = cummul_chemistry_master;
    R.parseEvalQ(
        "profiling$simtime_chemistry_master <- simtime_chemistry_master");

    R["seq_master_prep"] = cummul_master_seq_pre_loop;
    R.parseEvalQ("profiling$seq_master_prep <- seq_master_prep");
    R["seq_master_loop"] = cummul_master_seq_loop;
    R.parseEvalQ("profiling$seq_master_loop <- seq_master_loop");

    // R["master_send"] = master_send;
    // R.parseEvalQ("profiling$master_send <- master_send");
    // R["master_recv"] = master_recv;
    // R.parseEvalQ("profiling$master_recv <- master_recv");

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
      R["dht_collision"] = dht_collision;
      R.parseEvalQ("profiling$dht_collisions <- dht_collision");
      R["dht_get_time"] = dht_get_time;
      R.parseEvalQ("profiling$dht_get_time <- dht_get_time");
      R["dht_fill_time"] = dht_fill_time;
      R.parseEvalQ("profiling$dht_fill_time <- dht_fill_time");
    }

    free(workerlist);
    free(timings);

    if (dht_enabled)
      free(dht_perfs);

    cout << "CPP: Done! Results are stored as R objects into <" << out_dir
         << "/timings.rds>" << endl;
    /*exporting results and profiling data*/

    std::string r_vis_code;
    r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
    R.parseEval(r_vis_code);
  } else { /*This is executed by the workers*/
    if (!dht_file.empty()) {
      int res = file_to_table((char *)dht_file.c_str());
      if (res != DHT_SUCCESS) {
        if (res == DHT_WRONG_FILE) {
          if (world_rank == 2)
            cerr << "CPP: Worker: Wrong File" << endl;
        } else {
          if (world_rank == 2)
            cerr << "CPP: Worker: Error in loading current state of DHT from "
                    "file"
                 << endl;
        }
        return EXIT_FAILURE;
      } else {
        if (world_rank == 2)
          cout << "CPP: Worker: Successfully loaded state of DHT from file "
               << dht_file << endl;
        std::cout.flush();
      }
    }
    worker_function(R);
    free(mpi_buffer_results);
  }

  cout << "CPP: finished, cleanup of process " << world_rank << endl;

  if (dht_enabled) {

    if (dht_strategy == 0) {
      if (world_rank != 0) {
        DHT_free(dht_object, NULL, NULL);
      }
    } else {
      DHT_free(dht_object, NULL, NULL);
    }
  }

  free(mpi_buffer);
  MPI_Finalize();

  if (world_rank == 0) {
    cout << "CPP: done, bye!" << endl;
  }

  exit(0);
}
