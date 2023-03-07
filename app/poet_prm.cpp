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

#include "poet/ChemistryModule.hpp"
#include <RInside.h>
#include <Rcpp.h>
#include <Rcpp/internal/wrap.h>
#include <cstdint>
#include <cstdlib>
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

void set_chem_parameters(poet::ChemistryModule &chem, uint32_t wp_size,
                         const std::string &database_path) {
  chem.SetErrorHandlerMode(1);
  chem.SetComponentH2O(false);
  chem.SetRebalanceFraction(0.5);
  chem.SetRebalanceByCell(true);
  chem.UseSolutionDensityVolume(false);
  chem.SetPartitionUZSolids(false);

  // Set concentration units
  // 1, mg/L; 2, mol/L; 3, kg/kgs
  chem.SetUnitsSolution(2);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem.SetUnitsPPassemblage(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem.SetUnitsExchange(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem.SetUnitsSurface(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem.SetUnitsGasPhase(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem.SetUnitsSSassemblage(1);
  // 0, mol/L cell; 1, mol/L water; 2 mol/L rock
  chem.SetUnitsKinetics(1);

  // Set representative volume
  std::vector<double> rv;
  rv.resize(wp_size, 1.0);
  chem.SetRepresentativeVolume(rv);

  // Set initial porosity
  std::vector<double> por;
  por.resize(wp_size, 1);
  chem.SetPorosity(por);

  // Set initial saturation
  std::vector<double> sat;
  sat.resize(wp_size, 1.0);
  chem.SetSaturation(sat);

  // Load database
  chem.LoadDatabase(database_path);
}

inline double RunMasterLoop(SimParams &params, RInside &R, Grid &grid,
                            ChemistryParams &chem_params,
                            const GridParams &g_params, uint32_t nxyz_master) {

  DiffusionModule diffusion(poet::DiffusionParams(R), grid);
  /* Iteration Count is dynamic, retrieving value from R (is only needed by
   * master for the following loop) */
  uint32_t maxiter = R.parseEval("mysetup$iterations");

  double sim_time = .0;

  ChemistryModule chem(grid.GetTotalCellCount(), MPI_COMM_WORLD);
  set_chem_parameters(chem, nxyz_master, chem_params.database_path);
  chem.RunInitFile(chem_params.input_script);

  chem.SetSelectedOutputOn(true);
  chem.SetTimeStep(0);
  chem.RunCells();

  StateMemory *chem_state = grid.RegisterState("state_C", chem.GetPropNames());
  auto &chem_field = chem_state->mem;
  chem_field = chem.GetField();
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

    chem.SetTimeStep(dt);

    chem.SetConcentrations(chem_field);
    chem.SetTimeStep(dt);
    chem.RunCells();
    chem_field = chem.GetField();

    grid.WriteFieldsToR(R);
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

    // MPI_Barrier(MPI_COMM_WORLD);

  } // END SIMULATION LOOP

  R.parseEvalQ("profiling <- list()");

  R["simtime_chemistry"] = chem.GetChemistryTime();
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");

  chem.MpiWorkerBreak();
  diffusion.end();

  return MPI_Wtime() - dStartTime;
}

int main(int argc, char *argv[]) {
  double dSimTime, sim_end;

  int world_size, world_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  if (world_rank == 0) {
    cout << "Running POET in version " << poet_version << endl << endl;
  }

  if (world_rank > 0) {
    {
      uint32_t nxyz;
      MPI_Bcast(&nxyz, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
      ChemistryModule worker(nxyz, MPI_COMM_WORLD);
      worker.MpiWorker();
    }

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "CPP: finished, cleanup of process " << world_rank << endl;
    MPI_Finalize();

    return EXIT_SUCCESS;
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

  R.parseEvalQ("mysetup <- setup");
  // if (world_rank == 0) { // get timestep vector from
  // grid_init function ... //
  std::string master_init_code = "mysetup <- master_init(setup=setup)";
  R.parseEval(master_init_code);

  Grid grid;
  GridParams g_params(R);

  grid.InitModuleFromParams(g_params);
  grid.PushbackModuleFlow(poet::DIFFUSION_MODULE_NAME, CHEMISTRY_MODULE_NAME);
  grid.PushbackModuleFlow(CHEMISTRY_MODULE_NAME, poet::DIFFUSION_MODULE_NAME);

  params.initVectorParams(R, grid.GetSpeciesCount());

  // MDL: store all parameters
  if (world_rank == 0) {
    cout << "CPP: Calling R Function to store calling parameters" << endl;
    R.parseEvalQ("StoreSetup(setup=mysetup)");
  }

  if (world_rank == 0) {
    cout << "CPP: Init done on process with rank " << world_rank << endl;
  }

  // MPI_Barrier(MPI_COMM_WORLD);

  poet::ChemistryParams chem_params(R);

  /* THIS IS EXECUTED BY THE MASTER */
  uint32_t nxyz = grid.GetTotalCellCount();
  MPI_Bcast(&nxyz, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);
  uint32_t nxyz_master = grid.GetTotalCellCount();

  dSimTime = RunMasterLoop(params, R, grid, chem_params, g_params, nxyz_master);

  cout << "CPP: finished simulation loop" << endl;

  cout << "CPP: start timing profiling" << endl;

  R["simtime"] = dSimTime;
  R.parseEvalQ("profiling$simtime <- simtime");

  string r_vis_code;
  r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
  R.parseEval(r_vis_code);

  cout << "CPP: Done! Results are stored as R objects into <"
       << params.getOutDir() << "/timings.rds>" << endl;

  MPI_Barrier(MPI_COMM_WORLD);

  cout << "CPP: finished, cleanup of process " << world_rank << endl;
  MPI_Finalize();

  if (world_rank == 0) {
    cout << "CPP: done, bye!" << endl;
  }

  exit(0);
}
