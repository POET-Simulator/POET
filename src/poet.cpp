/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
**
** Copyright (C) 2023-2024 Marco De Lucia (GFZ Potsdam), Max Luebke (University
** of Potsdam)
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

#include "Base/Macros.hpp"
#include "Base/RInsidePOET.hpp"
#include "CLI/CLI.hpp"
#include "Chemistry/ChemistryModule.hpp"
#include "DataStructures/Field.hpp"
#include "Init/InitialList.hpp"
#include "Transport/DiffusionModule.hpp"
#include <RInside.h>
#include <Rcpp.h>
#include <Rcpp/DataFrame.h>
#include <Rcpp/Function.h>
#include <Rcpp/vector/instantiation.h>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mpi.h>
#include <string>
#include <mutex>
#include <condition_variable>
#include "Chemistry/SurrogateModels/AI_functions.hpp"
#include <CLI/CLI.hpp>
#include <poet.hpp>
#include <vector>

using namespace std;
using namespace poet;
using namespace Rcpp;

static int MY_RANK = 0;

static std::unique_ptr<Rcpp::List> global_rt_setup;

// we need some lazy evaluation, as we can't define the functions
// before the R runtime is initialized
static poet::DEFunc master_init_R;
static poet::DEFunc master_iteration_end_R;
static poet::DEFunc store_setup_R;
static poet::DEFunc ReadRObj_R;
static poet::DEFunc SaveRObj_R;
static poet::DEFunc source_R;

static void init_global_functions(RInside &R) {
  R.parseEval(kin_r_library);
  master_init_R = DEFunc("master_init");
  master_iteration_end_R = DEFunc("master_iteration_end");
  store_setup_R = DEFunc("StoreSetup");
  source_R = DEFunc("source");
  ReadRObj_R = DEFunc("ReadRObj");
  SaveRObj_R = DEFunc("SaveRObj");
}


// HACK: this is a step back as the order and also the count of fields is
// predefined, but it will change in the future
// static inline void writeFieldsToR(RInside &R, const Field &trans,
//                                   const Field &chem) {

//   Rcpp::DataFrame t_field(trans.asSEXP());
//   R["TMP"] = t_field;
//   R.parseEval("mysetup$state_T <- TMP");

//   R["TMP"] = chem.asSEXP();
//   R.parseEval("mysetup$state_C <- TMP");
// }

enum ParseRet { PARSER_OK, PARSER_ERROR, PARSER_HELP };

int parseInitValues(int argc, char **argv, RuntimeParameters &params) {

  CLI::App app{"POET - Potsdam rEactive Transport simulator"};

  app.add_flag("-P,--progress", params.print_progress,
               "Print progress bar during chemical simulation");

  /*Parse work package size*/
  app.add_option(
         "-w,--work-package-size", params.work_package_size,
         "Work package size to distribute to each worker for chemistry module")
      ->check(CLI::PositiveNumber)
      ->default_val(RuntimeParameters::WORK_PACKAGE_SIZE_DEFAULT);

  /* Parse DHT arguments */
  auto *dht_group = app.add_option_group("DHT", "DHT related options");

  dht_group->add_flag("--dht", params.use_dht, "Enable DHT");

  // cout << "CPP: DHT is " << ( dht_enabled ? "ON" : "OFF" ) << '\n';

  dht_group
      ->add_option("--dht-size", params.dht_size,
                   "DHT size per process in Megabyte")
      ->check(CLI::PositiveNumber)
      ->default_val(RuntimeParameters::DHT_SIZE_DEFAULT);
  // cout << "CPP: DHT size per process (Byte) = " << dht_size_per_process <<
  // endl;

  dht_group->add_option(
      "--dht-snaps", params.dht_snaps,
      "Save snapshots of DHT to disk: \n0 = disabled (default)\n1 = After "
      "simulation\n2 = After each iteration");

  auto *interp_group =
      app.add_option_group("Interpolation", "Interpolation related options");

  interp_group->add_flag("--interp", params.use_interp, "Enable interpolation");
  interp_group
      ->add_option("--interp-size", params.interp_size,
                   "Size of the interpolation table in Megabyte")
      ->check(CLI::PositiveNumber)
      ->default_val(RuntimeParameters::INTERP_SIZE_DEFAULT);
  interp_group
      ->add_option("--interp-min", params.interp_min_entries,
                   "Minimum number of entries in the interpolation table")
      ->check(CLI::PositiveNumber)
      ->default_val(RuntimeParameters::INTERP_MIN_ENTRIES_DEFAULT);
  interp_group
      ->add_option(
          "--interp-bucket-entries", params.interp_bucket_entries,
          "Maximum number of entries in each bucket of the interpolation table")
      ->check(CLI::PositiveNumber)
      ->default_val(RuntimeParameters::INTERP_BUCKET_ENTRIES_DEFAULT);

  app.add_flag("--ai-surrogate", params.use_ai_surrogate,
               "Enable AI surrogate for chemistry module");

  app.add_flag("--rds", params.as_rds,
               "Save output as .rds file instead of .qs");

  std::string init_file;
  std::string runtime_file;

  app.add_option("runtime_file", runtime_file,
                 "Runtime R script defining the simulation")
      ->required()
      ->check(CLI::ExistingFile);

  app.add_option(
         "init_file", init_file,
         "Initial R script defining the simulation, produced by poet_init")
      ->required()
      ->check(CLI::ExistingFile);

  app.add_option("out_dir", params.out_dir,
                 "Output directory of the simulation")
      ->required();

  try {
    app.parse(argc, argv);
  } catch (const CLI::ParseError &e) {
    app.exit(e);
    return -1;
  }

  // set the output extension
  params.out_ext = params.as_rds ? "rds" : "qs";

  if (MY_RANK == 0) {
    // MSG("Complete results storage is " + BOOL_PRINT(simparams.store_result));
    MSG("Output format/extension is " + params.out_ext);
    MSG("Work Package Size: " + std::to_string(params.work_package_size));
    MSG("DHT is " + BOOL_PRINT(params.use_dht));
    MSG("AI Surrogate is " + BOOL_PRINT(params.use_ai_surrogate));
    #ifndef USE_AI_SURROGATE
      if (params.use_ai_surrogate) {
        throw std::runtime_error("AI Surrogate functions can only be used if they are included during compile time.\n \
          Please use the CMake flag -DUSE_AI_SURROGATE=ON.");
      }
    #endif

    if (params.use_dht) {
      // MSG("DHT strategy is " + std::to_string(simparams.dht_strategy));
      // MDL: these should be outdated (?)
      // MSG("DHT key default digits (ignored if 'signif_vector' is "
      // 	"defined) = "
      // 	 << simparams.dht_significant_digits);
      // MSG("DHT logarithm before rounding: "
      // 	 << (simparams.dht_log ? "ON" : "OFF"));
      MSG("DHT size per process (Megabyte) = " +
          std::to_string(params.dht_size));
      MSG("DHT save snapshots is " + BOOL_PRINT(params.dht_snaps));
      // MSG("DHT load file is " + chem_params.dht_file);
    }

    if (params.use_interp) {
      MSG("PHT interpolation enabled: " + BOOL_PRINT(params.use_interp));
      MSG("PHT interp-size = " + std::to_string(params.interp_size));
      MSG("PHT interp-min  = " + std::to_string(params.interp_min_entries));
      MSG("PHT interp-bucket-entries = " +
          std::to_string(params.interp_bucket_entries));
    }
  }
  // chem_params.dht_outdir = out_dir;

  /* distribute information to R runtime */
  // if local_rank == 0 then master else worker
  // R["local_rank"] = MY_RANK;
  // assign a char* (string) to 'filesim'
  // R["filesim"] = wrap(runtime_file);
  // assign a char* (string) to 'fileout'
  // R["fileout"] = wrap(out_dir);
  // pass the boolean "store_result" to the R process
  // R["store_result"] = simparams.store_result;
  // // worker count
  // R["n_procs"] = simparams.world_size - 1;
  // // work package size
  // R["work_package_size"] = simparams.wp_size;
  // // dht enabled?
  // R["dht_enabled"] = chem_params.use_dht;
  // // log before rounding?
  // R["dht_log"] = simparams.dht_log;

  try {

    Rcpp::List init_params_(ReadRObj_R(init_file));
    params.init_params = init_params_;

    global_rt_setup = std::make_unique<Rcpp::List>();
    *global_rt_setup = source_R(runtime_file, Rcpp::Named("local", true));
    *global_rt_setup = global_rt_setup->operator[]("value");

    // MDL add "out_ext" for output format to R setup
    (*global_rt_setup)["out_ext"] = params.out_ext;

    params.timesteps =
        Rcpp::as<std::vector<double>>(global_rt_setup->operator[]("timesteps"));

  } catch (const std::exception &e) {
    ERRMSG("Error while parsing R scripts: " + std::string(e.what()));
    return ParseRet::PARSER_ERROR;
  }

  return ParseRet::PARSER_OK;
}

// HACK: this is a step back as the order and also the count of fields is
// predefined, but it will change in the future
void call_master_iter_end(RInside &R, const Field &trans, const Field &chem) {
  R["TMP"] = Rcpp::wrap(trans.AsVector());
  R["TMP_PROPS"] = Rcpp::wrap(trans.GetProps());
  R.parseEval(std::string("state_T <- setNames(data.frame(matrix(TMP, nrow=" +
                          std::to_string(trans.GetRequestedVecSize()) +
                          ")), TMP_PROPS)"));

  R["TMP"] = Rcpp::wrap(chem.AsVector());
  R["TMP_PROPS"] = Rcpp::wrap(chem.GetProps());
  R.parseEval(std::string("state_C <- setNames(data.frame(matrix(TMP, nrow=" +
                          std::to_string(chem.GetRequestedVecSize()) +
                          ")), TMP_PROPS)"));
  R["setup"] = *global_rt_setup;
  R.parseEval("setup <- master_iteration_end(setup, state_T, state_C)");
  *global_rt_setup = R["setup"];
}

static Rcpp::List RunMasterLoop(RInsidePOET &R, const RuntimeParameters &params,
                                DiffusionModule &diffusion,
                                ChemistryModule &chem) {

  /* Iteration Count is dynamic, retrieving value from R (is only needed by
   * master for the following loop) */
  uint32_t maxiter = params.timesteps.size();
  
  if (params.print_progress) {
    chem.setProgressBarPrintout(true);
  } 

  /* For the weights and biases of the AI surrogate
   * model to use in an inference function with Eigen */
  std::mutex Eigen_model_mutex; 
  static EigenModel Eigen_model;
  /* For the training data */
  std::mutex training_data_buffer_mutex;
  std::condition_variable training_data_buffer_full;
  bool start_training, end_training;
  TrainingData training_data_buffer;
  if (params.use_ai_surrogate) {  
    MSG("AI: Initialize model");
    Python_Keras_load_model(R["model_file_path"], params.cuda_src_dir);
    if (!params.disable_training) {
      MSG("AI: Initialize training thread");
      Python_Keras_training_thread(&Eigen_model, &Eigen_model_mutex,
                                  &training_data_buffer, &training_data_buffer_mutex,
                                  &training_data_buffer_full, &start_training, &end_training,
                                  params);
    }
    if (!params.use_Keras_predictions) {
      // Initialize Eigen model for custom inference function
      MSG("AI: Use custom C++ prediction function");
      // Get Keras weights from Python
      std::vector<std::vector<std::vector<double>>> cpp_weights = Python_Keras_get_weights();
      // Set model size
      size_t num_layers = cpp_weights.size() / 2;
      Eigen_model.weight_matrices.resize(num_layers); 
      Eigen_model.biases.resize(num_layers);
      for (size_t i = 0; i < cpp_weights.size(); i += 2) {
          size_t rows = cpp_weights[i][0].size();
          size_t cols = cpp_weights[i].size();
          Eigen_model.weight_matrices[i / 2].resize(rows, cols);
          size_t bias_size = cpp_weights[i + 1][0].size();
          Eigen_model.biases[i / 2].resize(bias_size);
      }
      // Set initial model weights
      update_weights(&Eigen_model, cpp_weights);
    }
    MSG("AI: Surrogate model initialized");
  }

  R["TMP_PROPS"] = Rcpp::wrap(chem.getField().GetProps());
  R["field_nrow"] = chem.getField().GetRequestedVecSize();

    /* SIMULATION LOOP */
  double dSimTime{0};
  for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
    double start_t = MPI_Wtime();
    
    const double &dt = params.timesteps[iter - 1];

    //  cout << "CPP: Next time step is " << dt << "[s]" << endl;
    MSG("Next time step is " + std::to_string(dt) + " [s]");

    /* displaying iteration number, with C++ and R iterator */
    MSG("Going through iteration " + std::to_string(iter));

    /* run transport */
    diffusion.simulate(dt);

    chem.getField().update(diffusion.getField());

    MSG("Chemistry step");
    
    if (params.use_ai_surrogate) {
      double ai_start_t = MPI_Wtime();
      double ai_start_steps = MPI_Wtime();
      // Get current values from the tug field for the ai predictions
      R["TMP"] = Rcpp::wrap(chem.getField().AsVector());
      R.parseEval(std::string("predictors <- ") + 
        "set_field(TMP, TMP_PROPS, field_nrow, ai_surrogate_species)");


      double ai_end_t = MPI_Wtime();
      R["diff_to_R"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();

      // Apply preprocessing
      MSG("AI Preprocessing");
      R.parseEval("predictors_scaled <- preprocess(predictors)");

      ai_end_t = MPI_Wtime();
      R["R_preprocessing"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();
      

      std::vector<std::vector<double>> x = R["predictors_scaled"];
      ai_end_t = MPI_Wtime();
      R["R_preprocessed_to_cxx"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();

      MSG("AI: Predict");
      if (params.use_Keras_predictions) {  // Predict with Keras default function
        R["TMP"] = Python_Keras_predict(R["predictors_scaled"], params.batch_size);

      } else {  // Predict with custom Eigen function
        R["TMP"] = Eigen_predict(Eigen_model, R["predictors_scaled"], params.batch_size, &Eigen_model_mutex);
      }

      ai_end_t = MPI_Wtime();
      R["cxx_inference"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();
      
      R["xyz_THROWAWAY"] = x;
       ai_end_t = MPI_Wtime();
      std::cout << "C++ predictions back to R: " << ai_end_t - ai_start_steps << std::endl;
      R["cxx_predictions_to_R"] = ai_end_t - ai_start_steps; 
      ai_start_steps = MPI_Wtime();

      // Apply postprocessing
      MSG("AI: Postprocesing");
      R.parseEval(std::string("predictions_scaled <- ") + 
        "set_field(TMP, ai_surrogate_species, field_nrow, ai_surrogate_species, byrow = TRUE)");
      R.parseEval("predictions <- postprocess(predictions_scaled)");

      ai_end_t = MPI_Wtime();
      R["R_postprocessing"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();

      // Validate prediction and write valid predictions to chem field
      MSG("AI: Validate");
      R.parseEval("validity_vector <- validate_predictions(predictors, predictions)");




      ai_end_t = MPI_Wtime();
      R["R_validate"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();




      MSG("AI: Marking valid");
      chem.set_ai_surrogate_validity_vector(R.parseEval("validity_vector"));



      ai_end_t = MPI_Wtime();
      R["validity_to_cxx"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();


      std::vector<std::vector<double>> RTempField =
        R.parseEval("set_valid_predictions(predictors, predictions, validity_vector)");

      Field predictions_field =
          Field(R.parseEval("nrow(predictors)"), RTempField,
                R.parseEval("colnames(predictors)"));

      MSG("AI: Update field with AI predictions");
      chem.getField().update(predictions_field);


      ai_end_t = MPI_Wtime();
      R["update_field"] = ai_end_t - ai_start_steps;
      ai_start_steps = MPI_Wtime();
      
      // store time for output file
      // double ai_end_t = MPI_Wtime();
       ai_end_t = MPI_Wtime();
      R["ai_prediction_time"] = ai_end_t - ai_start_t;

      if (!params.disable_training) {
        // Add to training data buffer:
        // Input values for which the predictions were invalid
        MSG("AI: Add invalid input data to training data buffer");
        std::vector<std::vector<double>> invalid_x = 
          R.parseEval("get_invalid_values(predictors_scaled, validity_vector)");
        training_data_buffer_mutex.lock();
        training_data_buffer_append(training_data_buffer.x, invalid_x);
        training_data_buffer_mutex.unlock();
      }

      ai_end_t = MPI_Wtime();
      R["append_to_training_buffer"] = ai_end_t - ai_start_steps;
    }

    // Run simulation step
    MSG("Simulate chemistry");
    chem.simulate(dt);

    // MPI_Barrier(MPI_COMM_WORLD);
    double end_t = MPI_Wtime();
    dSimTime += end_t - start_t;
    R["totaltime"] = dSimTime;

    // MDL master_iteration_end just writes on disk state_T and
    // state_C after every iteration if the cmdline option
    // --ignore-results is not given (and thus the R variable
    // store_result is TRUE)
    call_master_iter_end(R, diffusion.getField(), chem.getField());

    /* AI surrogate iterative training*/
    if (params.use_ai_surrogate && !params.disable_training) {
      // Add to training data buffer targets:
      // True values for invalid predictions      
      MSG("AI: Add invalid target data to training data buffer");
      R.parseEval("target_scaled <- preprocess(state_C[ai_surrogate_species])");
      std::vector<std::vector<double>> invalid_y = 
        R.parseEval("get_invalid_values(target_scaled, validity_vector)");
      training_data_buffer_mutex.lock();
      training_data_buffer_append(training_data_buffer.y, invalid_y);

      // Signal to training thread if training data buffer is full
      if (training_data_buffer.y[0].size() >= params.training_data_size) {
        start_training = true;
        training_data_buffer_full.notify_one();
      }
      training_data_buffer_mutex.unlock();
      R["n_training_runs"] = training_data_buffer.n_training_runs;
    }


    diffusion.getField().update(chem.getField());

    MSG("End of *coupling* iteration " + std::to_string(iter) + "/" +
        std::to_string(maxiter));
    MSG();
  } // END SIMULATION LOOP

  Rcpp::List chem_profiling;
  chem_profiling["simtime"] = chem.GetChemistryTime();
  chem_profiling["loop"] = chem.GetMasterLoopTime();
  chem_profiling["sequential"] = chem.GetMasterSequentialTime();
  chem_profiling["idle_master"] = chem.GetMasterIdleTime();
  chem_profiling["idle_worker"] = Rcpp::wrap(chem.GetWorkerIdleTimings());
  chem_profiling["phreeqc_time"] = Rcpp::wrap(chem.GetWorkerPhreeqcTimings());

  Rcpp::List diffusion_profiling;
  diffusion_profiling["simtime"] = diffusion.getTransportTime();

  if (params.use_dht) {
    chem_profiling["dht_hits"] = Rcpp::wrap(chem.GetWorkerDHTHits());
    chem_profiling["dht_evictions"] = Rcpp::wrap(chem.GetWorkerDHTEvictions());
    chem_profiling["dht_get_time"] = Rcpp::wrap(chem.GetWorkerDHTGetTimings());
    chem_profiling["dht_fill_time"] =
        Rcpp::wrap(chem.GetWorkerDHTFillTimings());
  }

  if (params.use_interp) {
    chem_profiling["interp_w"] =
        Rcpp::wrap(chem.GetWorkerInterpolationWriteTimings());
    chem_profiling["interp_r"] =
        Rcpp::wrap(chem.GetWorkerInterpolationReadTimings());
    chem_profiling["interp_g"] =
        Rcpp::wrap(chem.GetWorkerInterpolationGatherTimings());
    chem_profiling["interp_fc"] =
        Rcpp::wrap(chem.GetWorkerInterpolationFunctionCallTimings());
    chem_profiling["interp_calls"] =
        Rcpp::wrap(chem.GetWorkerInterpolationCalls());
    chem_profiling["interp_cached"] = Rcpp::wrap(chem.GetWorkerPHTCacheHits());
  }

  Rcpp::List profiling;
  profiling["simtime"] = dSimTime;
  profiling["chemistry"] = chem_profiling;
  profiling["diffusion"] = diffusion_profiling;

  if (params.use_ai_surrogate) {
    MSG("Finalize Python and wind down training thread");
    Python_finalize(&Eigen_model_mutex, &training_data_buffer_mutex,
                    &training_data_buffer_full, &start_training, &end_training);
  }

  chem.MasterLoopBreak();

  return profiling;
}

std::vector<std::string> getSpeciesNames(const Field &&field, int root,
                                         MPI_Comm comm) {
  std::uint32_t n_elements;
  std::uint32_t n_string_size;

  int rank;
  MPI_Comm_rank(comm, &rank);

  const bool is_master = root == rank;

  // first, the master sends all the species names iterative
  if (is_master) {
    n_elements = field.GetProps().size();
    MPI_Bcast(&n_elements, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);

    for (std::uint32_t i = 0; i < n_elements; i++) {
      n_string_size = field.GetProps()[i].size();
      MPI_Bcast(&n_string_size, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);
      MPI_Bcast(const_cast<char *>(field.GetProps()[i].c_str()), n_string_size,
                MPI_CHAR, root, MPI_COMM_WORLD);
    }

    return field.GetProps();
  }

  // now all the worker stuff
  MPI_Bcast(&n_elements, 1, MPI_UINT32_T, root, comm);

  std::vector<std::string> species_names_out(n_elements);

  for (std::uint32_t i = 0; i < n_elements; i++) {
    MPI_Bcast(&n_string_size, 1, MPI_UINT32_T, root, MPI_COMM_WORLD);

    char recv_buf[n_string_size];

    MPI_Bcast(recv_buf, n_string_size, MPI_CHAR, root, MPI_COMM_WORLD);

    species_names_out[i] = std::string(recv_buf, n_string_size);
  }

  return species_names_out;
}

int main(int argc, char *argv[]) {
  int world_size;

  MPI_Init(&argc, &argv);

  {
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &MY_RANK);

    RInsidePOET &R = RInsidePOET::getInstance();

    if (MY_RANK == 0) {
      MSG("Running POET version " + std::string(poet_version));
    }

    init_global_functions(R);

    RuntimeParameters run_params;

    if (parseInitValues(argc, argv, run_params) != 0) {
      MPI_Finalize();
      return 0;
    }

    // switch (parseInitValues(argc, argv, run_params)) {
    // case ParseRet::PARSER_ERROR:
    // case ParseRet::PARSER_HELP:
    //   MPI_Finalize();
    //   return 0;
    // case ParseRet::PARSER_OK:
    //   break;
    // }

    InitialList init_list(R);
    init_list.importList(run_params.init_params, MY_RANK != 0);

    MSG("RInside initialized on process " + std::to_string(MY_RANK));

    std::cout << std::flush;

    MPI_Barrier(MPI_COMM_WORLD);

    ChemistryModule chemistry(run_params.work_package_size,
                              init_list.getChemistryInit(), MPI_COMM_WORLD);

    const ChemistryModule::SurrogateSetup surr_setup = {

        getSpeciesNames(init_list.getInitialGrid(), 0, MPI_COMM_WORLD),
        run_params.use_dht,
        run_params.dht_size,
        run_params.use_interp,
        run_params.interp_bucket_entries,
        run_params.interp_size,
        run_params.interp_min_entries,
        run_params.use_ai_surrogate};

    chemistry.masterEnableSurrogates(surr_setup);

    if (MY_RANK > 0) {
      chemistry.WorkerLoop();
    } else {
      // R.parseEvalQ("mysetup <- setup");
      // // if (MY_RANK == 0) { // get timestep vector from
      // // grid_init function ... //

      *global_rt_setup = master_init_R(*global_rt_setup, run_params.out_dir,
                                       init_list.getInitialGrid().asSEXP());

      // MDL: store all parameters
      // MSG("Calling R Function to store calling parameters");
      // R.parseEvalQ("StoreSetup(setup=mysetup)");
      R["out_ext"] = run_params.out_ext;
      R["out_dir"] = run_params.out_dir;

      // MPI_Barrier(MPI_COMM_WORLD);

      DiffusionModule diffusion(init_list.getDiffusionInit(),
                                init_list.getInitialGrid());

      chemistry.masterSetField(init_list.getInitialGrid());

      if (run_params.use_ai_surrogate) {
        // Load default function implementations 
        R.parseEvalQ(ai_surrogate_r_library);
        /* Use dht species for model input and output */
        R["ai_surrogate_species"] =
            init_list.getChemistryInit().dht_species.getNames();

        const std::string ai_surrogate_input_script =
            init_list.getChemistryInit().ai_surrogate_input_script;

        MSG("AI: Sourcing user-provided script");
        R.parseEvalQ(ai_surrogate_input_script);
        if (!Rcpp::as<bool>(R.parseEval("exists(\"model_file_path\")"))) {
          throw std::runtime_error("AI surrogate input script must contain a value for model_file_path");
        }
        
        /* AI surrogate training and inference parameters. (Can be set by declaring a 
        variable of the same name in one of the the R input scripts)*/
        run_params.use_Keras_predictions = false; // Default inference function is custom C++ / Eigen implementation
        run_params.disable_training = false; // Model will be trained per default
        run_params.batch_size = 2560; // default value determined in test on the UP Turing cluster
        run_params.training_epochs = 20; // 
        run_params.training_data_size = init_list.getDiffusionInit().n_rows *
                                        init_list.getDiffusionInit().n_cols; // Default value is number of cells in field
        run_params.save_model_path = ""; // Model is only saved if a path is set in the input field
        if (Rcpp::as<bool>(R.parseEval("exists(\"batch_size\")"))) {
          run_params.batch_size = R["batch_size"];
        }
        if (Rcpp::as<bool>(R.parseEval("exists(\"training_epochs\")"))) {
          run_params.training_epochs = R["training_epochs"];
        }
        if (Rcpp::as<bool>(R.parseEval("exists(\"training_data_size\")"))) {
          run_params.training_data_size = R["training_data_size"];
        }
        if (Rcpp::as<bool>(R.parseEval("exists(\"use_Keras_predictions\")"))) {
          run_params.use_Keras_predictions = R["use_Keras_predictions"];
        }
        if (Rcpp::as<bool>(R.parseEval("exists(\"disable_training\")"))) {
          run_params.disable_training = R["disable_training"];
        }
        if (Rcpp::as<bool>(R.parseEval("exists(\"save_model_path\")"))) {
          run_params.save_model_path = Rcpp::as<std::string>(R["save_model_path"]);
          std::cout << "AI: Model will be saved as \"" << run_params.save_model_path << "\"" << std::endl;
        }
        
        MSG("AI: Initialize Python for AI surrogate functions");
        std::string python_keras_file = std::string(SRC_DIR) +
          "/src/Chemistry/SurrogateModels/AI_Python_functions/keras_AI_surrogate.py";
        Python_Keras_setup(python_keras_file);
      }

      MSG("Init done on process with rank " + std::to_string(MY_RANK));

      Rcpp::List profiling = RunMasterLoop(R, run_params, diffusion, chemistry);

      MSG("finished simulation loop");

      R["profiling"] = profiling;
      R["setup"] = *global_rt_setup;
      R["setup$out_ext"] = run_params.out_ext;

      string r_vis_code;
      r_vis_code = "SaveRObj(x = profiling, path = paste0(out_dir, "
                   "'/timings.', setup$out_ext));";
      R.parseEval(r_vis_code);

      MSG("Done! Results are stored as R objects into <" + run_params.out_dir +
          "/timings." + run_params.out_ext);
    }
  }

  MSG("finished, cleanup of process " + std::to_string(MY_RANK));

  MPI_Finalize();

  if (MY_RANK == 0) {
    MSG("done, bye!");
  }

  exit(0);
}
