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

#include "Base/Grid.hpp"
#include "Base/Macros.hpp"
#include "Base/RInsidePOET.hpp"
#include "Base/SimParams.hpp"
#include "Chemistry/ChemistryModule.hpp"
#include "Transport/AdvectionModule.hpp"
#include "Transport/DiffusionModule.hpp"

#include <poet.hpp>

#include <Rcpp.h>
#include <cstdint>
#include <cstdlib>

#include <cstring>
#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include <mpi.h>

using namespace std;
using namespace poet;
using namespace Rcpp;

poet::ChemistryModule::SingleCMap DFToHashMap(const Rcpp::DataFrame &df) {
  std::unordered_map<std::string, double> out_map;
  vector<string> col_names = Rcpp::as<vector<string>>(df.names());

  for (const auto &name : col_names) {
    double val = df[name.c_str()];
    out_map.insert({name, val});
  }

  return out_map;
}

// HACK: this is a step back as the order and also the count of fields is
// predefined, but it will change in the future
void writeFieldsToR(RInside &R, const Field &trans, const Field &chem) {
  R["TMP"] = Rcpp::wrap(trans.AsVector());
  R["TMP_PROPS"] = Rcpp::wrap(trans.GetProps());
  R.parseEval(std::string(
      "mysetup$state_T <- setNames(data.frame(matrix(TMP, nrow=" +
      std::to_string(trans.GetRequestedVecSize()) + ")), TMP_PROPS)"));

  R["TMP"] = Rcpp::wrap(chem.AsVector());
  R["TMP_PROPS"] = Rcpp::wrap(chem.GetProps());
  R.parseEval(std::string(
      "mysetup$state_C <- setNames(data.frame(matrix(TMP, nrow=" +
      std::to_string(chem.GetRequestedVecSize()) + ")), TMP_PROPS)"));
}

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

inline double RunMasterLoop(SimParams &params, RInsidePOET &R,
                            const GridParams &g_params, uint32_t nxyz_master) {

  // DiffusionParams d_params{R};
  // DiffusionModule diffusion(d_params, g_params);

  AdvectionModule advection(R);

  /* Iteration Count is dynamic, retrieving value from R (is only needed by
   * master for the following loop) */
  uint32_t maxiter = R.parseEval("mysetup$iterations");

  double sim_time = .0;

  ChemistryModule chem(nxyz_master, params.getNumParams().wp_size, maxiter,
                       params.getChemParams(), MPI_COMM_WORLD);

  set_chem_parameters(chem, nxyz_master, params.getChemParams().database_path);
  chem.RunInitFile(params.getChemParams().input_script);

  chem.initializeField(advection.getField());

  if (params.getNumParams().print_progressbar) {
    chem.setProgressBarPrintout(true);
  }

  /* SIMULATION LOOP */

  double dSimTime{0};
  for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
    double start_t = MPI_Wtime();
    uint32_t tick = 0;
    // cout << "CPP: Evaluating next time step" << endl;
    // R.parseEvalQ("mysetup <- master_iteration_setup(mysetup)");

    double dt = Rcpp::as<double>(
        R.parseEval("mysetup$timesteps[" + std::to_string(iter) + "]"));

    //  cout << "CPP: Next time step is " << dt << "[s]" << endl;
    MSG("Next time step is " + std::to_string(dt) + " [s]");

    /* displaying iteration number, with C++ and R iterator */
    MSG("Going through iteration " + std::to_string(iter));
    MSG("R's $iter: " +
        std::to_string((uint32_t)(R.parseEval("mysetup$iter"))) +
        ". Iteration");

    /* run transport */
    // TODO: transport to diffusion
    advection.simulate(dt);

    chem.getField().update(advection.getField());

    MSG("Chemistry step");

    chem.SetTimeStep(dt);
    chem.RunCells();

    writeFieldsToR(R, advection.getField(), chem.GetField());
    advection.getField().update(chem.GetField());

    R["req_dt"] = dt;
    R["simtime"] = (sim_time += dt);

    R.parseEval("mysetup$req_dt <- req_dt");
    R.parseEval("mysetup$simtime <- simtime");

    // MDL master_iteration_end just writes on disk state_T and
    // state_C after every iteration if the cmdline option
    // --ignore-results is not given (and thus the R variable
    // store_result is TRUE)
    R.parseEvalQ("mysetup <- master_iteration_end(setup=mysetup)");

    MSG("End of *coupling* iteration " + std::to_string(iter) + "/" +
        std::to_string(maxiter));
    MSG();

    // MPI_Barrier(MPI_COMM_WORLD);
    double end_t = MPI_Wtime();
    dSimTime += end_t - start_t;
  } // END SIMULATION LOOP

  R.parseEvalQ("profiling <- list()");

  R["simtime_chemistry"] = chem.GetChemistryTime();
  R.parseEvalQ("profiling$simtime_chemistry <- simtime_chemistry");

  R["chemistry_loop"] = chem.GetMasterLoopTime();
  R.parseEvalQ("profiling$chemistry_loop <- chemistry_loop");

  R["chemistry_sequential"] = chem.GetMasterSequentialTime();
  R.parseEvalQ("profiling$simtime_sequential <- chemistry_sequential");

  R["idle_master"] = chem.GetMasterIdleTime();
  R.parseEvalQ("profiling$idle_master <- idle_master");

  R["idle_worker"] = Rcpp::wrap(chem.GetWorkerIdleTimings());
  R.parseEvalQ("profiling$idle_worker <- idle_worker");

  R["phreeqc_time"] = Rcpp::wrap(chem.GetWorkerPhreeqcTimings());
  R.parseEvalQ("profiling$phreeqc <- phreeqc_time");

  R["simtime_transport"] = advection.getTransportTime();
  R.parseEvalQ("profiling$simtime_transport <- simtime_transport");

  // R["phreeqc_count"] = phreeqc_counts;
  // R.parseEvalQ("profiling$phreeqc_count <- phreeqc_count");

  if (params.getChemParams().use_dht) {
    R["dht_hits"] = Rcpp::wrap(chem.GetWorkerDHTHits());
    R.parseEvalQ("profiling$dht_hits <- dht_hits");
    R["dht_evictions"] = Rcpp::wrap(chem.GetWorkerDHTEvictions());
    R.parseEvalQ("profiling$dht_evictions <- dht_evictions");
    R["dht_corrupt_buckets"] = Rcpp::wrap(chem.GetWorkerDHTCorruptBuckets());
    R.parseEvalQ("profiling$dht_corrupt_buckets <- dht_corrupt_buckets");
    R["dht_get_time"] = Rcpp::wrap(chem.GetWorkerDHTGetTimings());
    R.parseEvalQ("profiling$dht_get_time <- dht_get_time");
    R["dht_fill_time"] = Rcpp::wrap(chem.GetWorkerDHTFillTimings());
    R.parseEvalQ("profiling$dht_fill_time <- dht_fill_time");
  }
  // if (params.getChemParams().use_interp) {
  //   R["interp_w"] = Rcpp::wrap(chem.GetWorkerInterpolationWriteTimings());
  //   R.parseEvalQ("profiling$interp_write <- interp_w");
  //   R["interp_r"] = Rcpp::wrap(chem.GetWorkerInterpolationReadTimings());
  //   R.parseEvalQ("profiling$interp_read <- interp_r");
  //   R["interp_g"] = Rcpp::wrap(chem.GetWorkerInterpolationGatherTimings());
  //   R.parseEvalQ("profiling$interp_gather <- interp_g");
  //   R["interp_fc"] =
  //       Rcpp::wrap(chem.GetWorkerInterpolationFunctionCallTimings());
  //   R.parseEvalQ("profiling$interp_function_calls <- interp_fc");
  //   R["interp_calls"] = Rcpp::wrap(chem.GetWorkerInterpolationCalls());
  //   R.parseEvalQ("profiling$interp_calls <- interp_calls");
  //   R["interp_cached"] = Rcpp::wrap(chem.GetWorkerPHTCacheHits());
  //   R.parseEvalQ("profiling$interp_cached <- interp_cached");
  // }

  chem.MasterLoopBreak();

  return dSimTime;
}

int main(int argc, char *argv[]) {
  double dSimTime, sim_end;

  int world_size, world_rank;

  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  RInsidePOET &R = RInsidePOET::getInstance();

  if (world_rank == 0) {
    MSG("Running POET version " + std::string(poet_version));
  }

  if (world_rank > 0) {
    {
      SimParams params(world_rank, world_size);
      int pret = params.parseFromCmdl(argv, R);

      if (pret == poet::PARSER_ERROR) {
        MPI_Finalize();
        return EXIT_FAILURE;
      } else if (pret == poet::PARSER_HELP) {
        MPI_Finalize();
        return EXIT_SUCCESS;
      }

      // ChemistryModule worker(nxyz, nxyz, MPI_COMM_WORLD);
      ChemistryModule worker = poet::ChemistryModule::createWorker(
          MPI_COMM_WORLD, params.getChemParams());
      set_chem_parameters(worker, worker.GetWPSize(),
                          params.getChemParams().database_path);

      worker.WorkerLoop();
    }

    MPI_Barrier(MPI_COMM_WORLD);

    MSG("finished, cleanup of process " + std::to_string(world_rank));

    MPI_Finalize();

    return EXIT_SUCCESS;
  }

  /*Loading Dependencies*/
  // TODO: kann raus
  R.parseEvalQ(kin_r_library);

  SimParams params(world_rank, world_size);
  int pret = params.parseFromCmdl(argv, R);

  if (pret == poet::PARSER_ERROR) {
    MPI_Finalize();
    return EXIT_FAILURE;
  } else if (pret == poet::PARSER_HELP) {
    MPI_Finalize();
    return EXIT_SUCCESS;
  }

  MSG("RInside initialized on process " + std::to_string(world_rank));

  R.parseEvalQ("mysetup <- setup");
  // if (world_rank == 0) { // get timestep vector from
  // grid_init function ... //
  std::string master_init_code = "mysetup <- master_init(setup=setup)";
  R.parseEval(master_init_code);

  GridParams g_params(R);

  params.initVectorParams(R);

  // MDL: store all parameters
  if (world_rank == 0) {
    MSG("Calling R Function to store calling parameters");
    R.parseEvalQ("StoreSetup(setup=mysetup)");
  }

  if (world_rank == 0) {
    MSG("Init done on process with rank " + std::to_string(world_rank));
  }

  // MPI_Barrier(MPI_COMM_WORLD);

  uint32_t nxyz_master = (world_size == 1 ? g_params.total_n : 1);

  dSimTime = RunMasterLoop(params, R, g_params, nxyz_master);

  MSG("finished simulation loop");

  MSG("start timing profiling");

  R["simtime"] = dSimTime;
  R.parseEvalQ("profiling$simtime <- simtime");

  string r_vis_code;
  r_vis_code = "saveRDS(profiling, file=paste0(fileout,'/timings.rds'));";
  R.parseEval(r_vis_code);

  MSG("Done! Results are stored as R objects into <" + params.getOutDir() +
      "/timings.rds>");

  MPI_Barrier(MPI_COMM_WORLD);

  MSG("finished, cleanup of process " + std::to_string(world_rank));
  MPI_Finalize();

  if (world_rank == 0) {
    MSG("done, bye!");
  }

  exit(0);
}
