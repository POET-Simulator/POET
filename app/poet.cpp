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

#include <Rcpp.h>
#include <poet/ChemSim.hpp>
#include <poet/Grid.hpp>
#include <poet/RRuntime.hpp>
#include <poet/SimParams.hpp>
#include <poet/DiffusionModule.hpp>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <poet.h>

using namespace std;
using namespace poet;
using namespace Rcpp;

int main(int argc, char *argv[]) {
  double sim_start, sim_end;

  int world_size, world_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  /*Create custom Communicator with all processes except 0 (the master) for DHT
   * storage */
  MPI_Comm dht_comm;

  if (world_rank == 0) {
    MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, world_rank, &dht_comm);
  } else {
    MPI_Comm_split(MPI_COMM_WORLD, 1, world_rank, &dht_comm);
  }

  if (world_rank == 0) {
    cout << "Running POET in version " << poet_version << endl << endl;
  }

  /* initialize R runtime */
  RRuntime R(argc, argv);

  /*Loading Dependencies*/
  // TODO: kann raus
  std::string r_load_dependencies = "suppressMessages(library(Rmufits));"
                                    "suppressMessages(library(RedModRphree));"
                                    "source('../R_lib/kin_r_library.R');"
                                    "source('../R_lib/parallel_r_library.R');";
  R.parseEvalQ(r_load_dependencies);

  SimParams params(world_rank, world_size);
  int pret = params.parseFromCmdl(argv, R);

  if (pret == PARSER_ERROR) {
    MPI_Finalize();
    return EXIT_FAILURE;
  } else if (pret == PARSER_HELP) {
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  cout << "CPP: R Init (RInside) on process " << world_rank << endl;

  // HACK: we disable master_init and dt_differ propagation here for testing
  // purposes
  //
  // bool dt_differ;
  R.parseEvalQ("mysetup <- setup");
  if (world_rank == 0) { // get timestep vector from
                         // grid_init function ... //
    std::string master_init_code = "mysetup <- master_init(setup=setup)";
    R.parseEval(master_init_code);

    //   dt_differ = R.parseEval("mysetup$dt_differ");
    //   // ... and broadcast it to every other rank unequal to 0
    //   MPI_Bcast(&dt_differ, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
    // }
    // /* workers will only read the setup DataFrame defined by input file */
    // else {
    //   R.parseEval("mysetup <- setup");
    //   MPI_Bcast(&dt_differ, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  }

  // params.setDtDiffer(dt_differ);

  // initialize chemistry on all processes
  // TODO: einlesen einer initial matrix (DataFrame)
  // std::string init_chemistry_code = "mysetup <-
  // init_chemistry(setup=mysetup)"; R.parseEval(init_chemistry_code);

  // TODO: Grid anpassen


  Grid grid(R, poet::GridParams(R));
  // grid.init_from_R();

  params.initVectorParams(R, grid.getSpeciesCount());

  // MDL: store all parameters
  if (world_rank == 0) {
    cout << "CPP: Calling R Function to store calling parameters" << endl;
    R.parseEvalQ("StoreSetup(setup=mysetup)");
  }

  if (world_rank == 0) {
    cout << "CPP: Init done on process with rank " << world_rank << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  /* THIS IS EXECUTED BY THE MASTER */
  if (world_rank == 0) {
    ChemMaster master(params, R, grid);
    DiffusionModule diffusion(poet::DiffusionParams(R), grid);

    // diffusion.initialize();

    sim_start = MPI_Wtime();

    /* Iteration Count is dynamic, retrieving value from R (is only needed by
     * master for the following loop) */
    uint32_t maxiter = R.parseEval("mysetup$iterations");

    /* SIMULATION LOOP */
    for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
      // cout << "CPP: Evaluating next time step" << endl;
      // R.parseEvalQ("mysetup <- master_iteration_setup(mysetup)");

      double dt = Rcpp::as<double>(
          R.parseEval("mysetup$timesteps[" + std::to_string(iter) + "]"));
      cout << "CPP: Next time step is " << dt << "[s]" << endl;

      /* displaying iteration number, with C++ and R iterator */
      cout << "CPP: Going through iteration " << iter << endl;
      cout << "CPP: R's $iter: " << ((uint32_t)(R.parseEval("mysetup$iter")))
           << ". Iteration" << endl;

      cout << "CPP: Calling Advection" << endl;

      /* run transport */
      // TODO: transport to diffusion
      diffusion.simulate(dt);

      cout << "CPP: Chemistry" << endl;

      /* Fallback for sequential execution */
      // TODO: use new grid
      if (world_size == 1) {
        master.ChemSim::run();
      }
      /* otherwise run parallel */
      else {
        master.run();
      }

      // MDL master_iteration_end just writes on disk state_T and
      // state_C after every iteration if the cmdline option
      // --ignore-results is not given (and thus the R variable
      // store_result is TRUE)
      R.parseEvalQ("mysetup <- master_iteration_end(setup=mysetup)");

      cout << endl
           << "CPP: End of *coupling* iteration " << iter << "/" << maxiter
           << endl
           << endl;

      MPI_Barrier(MPI_COMM_WORLD);

    } // END SIMULATION LOOP

    cout << "CPP: finished simulation loop" << endl;

    sim_end = MPI_Wtime();

    cout << "CPP: start timing profiling" << endl;

    R.parseEvalQ("profiling <- list()");

    R["simtime"] = sim_end - sim_start;
    R.parseEvalQ("profiling$simtime <- simtime");

    diffusion.end();

    if (world_size == 1) {
      master.ChemSim::end();
    } else {
      master.end();
    }

    string r_vis_code;
    r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
    R.parseEval(r_vis_code);

    cout << "CPP: Done! Results are stored as R objects into <"
         << params.getOutDir() << "/timings.rds>" << endl;
  }
  /* THIS IS EXECUTED BY THE WORKERS */
  else {
    ChemWorker worker(params, R, grid, dht_comm);
    worker.loop();
  }

  cout << "CPP: finished, cleanup of process " << world_rank << endl;
  MPI_Finalize();

  if (world_rank == 0) {
    cout << "CPP: done, bye!" << endl;
  }

  exit(0);
}
