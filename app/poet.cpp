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

#include <RInside.h>
#include <Rcpp.h>
#include <cstdint>
#include <poet/ChemSimPar.hpp>
#include <poet/ChemSimSeq.hpp>
#include <poet/DiffusionModule.hpp>
#include <poet/Grid.hpp>
#include <poet/SimParams.hpp>

#include <cstring>
#include <iostream>
#include <string>
#include <vector>

#include <mpi.h>
#include <poet.h>

using namespace std;
using namespace poet;
using namespace Rcpp;

template <class ChemistryInstance>
inline double RunMasterLoop(SimParams &params, RInside &R, Grid &grid,
                            ChemistryParams &chem_params) {

  DiffusionModule diffusion(poet::DiffusionParams(R), grid);
  /* Iteration Count is dynamic, retrieving value from R (is only needed by
   * master for the following loop) */
  uint32_t maxiter = R.parseEval("mysetup$iterations");

  double sim_time = .0;

  ChemistryInstance C(params, R, grid);
  C.InitModule(chem_params);
  /* SIMULATION LOOP */

  double dStartTime = MPI_Wtime();
  for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
    uint32_t tick = 0;
    // cout << "CPP: Evaluating next time step" << endl;
    // R.parseEvalQ("mysetup <- master_iteration_setup(mysetup)");

    double dt = Rcpp::as<double>(
        R.parseEval("mysetup$timesteps[" + std::to_string(iter) + "]"));
    cout << "CPP: Next time step is " << dt << "[s]" << endl;

    /* displaying iteration number, with C++ and R iterator */
    cout << "CPP: Going through iteration " << iter << endl;
    cout << "CPP: R's $iter: " << ((uint32_t)(R.parseEval("mysetup$iter")))
         << ". Iteration" << endl;

    /* run transport */
    // TODO: transport to diffusion
    diffusion.simulate(dt);

    grid.PreModuleFieldCopy(tick++);

    cout << "CPP: Chemistry" << endl;

    C.Simulate(dt);
    grid.PreModuleFieldCopy(tick++);

    R["req_dt"] = dt;
    R["simtime"] = (sim_time += dt);

    R.parseEval("mysetup$req_dt <- req_dt");
    R.parseEval("mysetup$simtime <- simtime");

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

  R.parseEvalQ("profiling <- list()");
  C.End();
  diffusion.end();

  return MPI_Wtime() - dStartTime;
}

int main(int argc, char *argv[]) {
  double dSimTime, sim_end;

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
  RInside R(argc, argv);

  /*Loading Dependencies*/
  // TODO: kann raus
  std::string r_load_dependencies = "source('../R_lib/kin_r_library.R');";
  R.parseEvalQ(r_load_dependencies);

  SimParams params(world_rank, world_size);
  int pret = params.parseFromCmdl(argv, R);

  if (pret == poet::PARSER_ERROR) {
    MPI_Finalize();
    return EXIT_FAILURE;
  } else if (pret == poet::PARSER_HELP) {
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  cout << "CPP: R Init (RInside) on process " << world_rank << endl;

  // HACK: we disable master_init and dt_differ propagation here for testing
  // purposes
  //
  bool dt_differ = false;
  R.parseEvalQ("mysetup <- setup");
  if (world_rank == 0) { // get timestep vector from
                         // grid_init function ... //
    std::string master_init_code = "mysetup <- master_init(setup=setup)";
    R.parseEval(master_init_code);

    //   dt_differ = R.parseEval("mysetup$dt_differ");
    //   // ... and broadcast it to every other rank unequal to 0
    MPI_Bcast(&dt_differ, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  }
  /* workers will only read the setup DataFrame defined by input file */
  else {
    // R.parseEval("mysetup <- setup");
    MPI_Bcast(&dt_differ, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
  }

  params.setDtDiffer(dt_differ);

  // initialize chemistry on all processes
  // TODO: einlesen einer initial matrix (DataFrame)
  // std::string init_chemistry_code = "mysetup <-
  // init_chemistry(setup=mysetup)"; R.parseEval(init_chemistry_code);

  // TODO: Grid anpassen

  Grid grid;

  grid.InitModuleFromParams(GridParams(R));
  grid.PushbackModuleFlow(poet::DIFFUSION_MODULE_NAME,
                          poet::BaseChemModule::CHEMISTRY_MODULE_NAME);
  grid.PushbackModuleFlow(poet::BaseChemModule::CHEMISTRY_MODULE_NAME,
                          poet::DIFFUSION_MODULE_NAME);

  params.initVectorParams(R, grid.GetSpeciesCount());

  // MDL: store all parameters
  if (world_rank == 0) {
    cout << "CPP: Calling R Function to store calling parameters" << endl;
    R.parseEvalQ("StoreSetup(setup=mysetup)");
  }

  if (world_rank == 0) {
    cout << "CPP: Init done on process with rank " << world_rank << endl;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  poet::ChemistryParams chem_params(R);

  /* THIS IS EXECUTED BY THE MASTER */
  if (world_rank == 0) {
    if (world_size == 1) {
      dSimTime = RunMasterLoop<ChemSeq>(params, R, grid, chem_params);
    } else {
      dSimTime = RunMasterLoop<ChemMaster>(params, R, grid, chem_params);
    }

    cout << "CPP: finished simulation loop" << endl;

    cout << "CPP: start timing profiling" << endl;

    R["simtime"] = dSimTime;
    R.parseEvalQ("profiling$simtime <- simtime");

    string r_vis_code;
    r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
    R.parseEval(r_vis_code);

    cout << "CPP: Done! Results are stored as R objects into <"
         << params.getOutDir() << "/timings.rds>" << endl;
  }
  /* THIS IS EXECUTED BY THE WORKERS */
  else {
    ChemWorker worker(params, R, grid, dht_comm);
    worker.InitModule(chem_params);
    worker.loop();
  }

  cout << "CPP: finished, cleanup of process " << world_rank << endl;
  MPI_Finalize();

  if (world_rank == 0) {
    cout << "CPP: done, bye!" << endl;
  }

  exit(0);
}
