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
#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <memory>
#include <mpi.h>
#include <set>
#include <string>

#include <Model.hpp>
#include <NAABackend.hpp>
#include <TrainingBackend.hpp>
#include <TrainingData.hpp>

#include <CLI/CLI.hpp>

#include <poet.hpp>
#include <vector>

// using namespace std;
using namespace poet;
using namespace Rcpp;

static int MY_RANK = 0;

static std::unique_ptr<Rcpp::List> global_rt_setup;

// we need some lazy evaluation, as we can't define the functions
// before the R runtime is initialized
static poet::DEFunc master_init_R;
static poet::DEFunc master_iteration_end_R;
// MDL: unused -> static poet::DEFunc store_setup_R;
static poet::DEFunc ReadRObj_R;
static poet::DEFunc SaveRObj_R;
static poet::DEFunc source_R;

static void init_global_functions(RInside &R) {
  R.parseEval(kin_r_library);
  master_init_R = DEFunc("master_init");
  master_iteration_end_R = DEFunc("master_iteration_end");
  // MDL: unused -> store_setup_R = DEFunc("StoreSetup");
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

  auto *ai_option_group =
      app.add_option_group("ai_surrogate", "AI Surrogate related options");

  ai_option_group->add_flag("--ai", params.ai,
                            "Enable AI surrogate for chemistry module");
  ai_option_group
      ->add_option("--ai-backend", params.ai_backend,
                   "Desired ai backend (0: python (keras), 1: naa, 2: cuda)")
      ->check(CLI::PositiveNumber)
      ->default_val(RuntimeParameters::AI_BACKEND_DEFAULT);

  app.add_flag("-c,--copy-non-reactive", params.copy_non_reactive_regions,
               "Copy non-reactive regions instead of computing them");

  app.add_option("-f,--fn", params.function_code, "Function code for the NAA")
      ->check(CLI::PositiveNumber);

  app.add_flag("--rds", params.as_rds,
               "Save output as .rds file instead of default .qs2");

  app.add_flag("--qs", params.as_qs,
               "Save output as .qs file instead of default .qs2");

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
  params.out_ext = "qs2";
  if (params.as_rds)
    params.out_ext = "rds";
  if (params.as_qs)
    params.out_ext = "qs";

  if (MY_RANK == 0) {
    // MSG("Complete results storage is " + BOOL_PRINT(simparams.store_result));
    MSG("Output format/extension is " + params.out_ext);
    MSG("Work Package Size: " + std::to_string(params.work_package_size));
    MSG("DHT is " + BOOL_PRINT(params.use_dht));
    MSG("AI Surrogate is " + BOOL_PRINT(params.ai));

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
  R["TMP_PROPS"] = Rcpp::wrap(chem.getField().GetProps());

  std::unique_ptr<AIContext> ai_ctx = nullptr;
  size_t retrain_counter = 0;
  size_t field_size = 0;

  if (params.ai) {

    if (params.function_code == 1) {
      ai_ctx = std::make_unique<AIContext>(
          "./bench/barite/barite_trained.weights.h5");
    } else if (params.function_code == 2) {
      ai_ctx = std::make_unique<AIContext>(
          "./bench/dolo/dolomite_trained.weights.h5");
    }

    R.parseEval(
        "mean <- as.numeric(standard$mean[ai_surrogate_species_output])");
    R.parseEval(
        "scale <- as.numeric(standard$scale[ai_surrogate_species_output])");

    std::vector<float> mean = R["mean"];
    std::vector<float> scale = R["scale"];

    field_size = chem.getField().GetRequestedVecSize();
    std::cout << field_size << std::endl;

    ai_ctx->scaler.set_scaler(mean, scale);

    // initialzie training backens only if retraining is desired
    if (params.ai_backend == PYTHON_BACKEND) {
      MSG("AI Surrogate with Python/keras backend enabled.")
      std::cerr << "Not implemented" << std::endl;
    } else if (params.ai_backend == NAA_BACKEND) {
      MSG("AI Surrogate with NAA backend enabled.")
      ai_ctx->training_backend =
          std::make_unique<NAABackend<ai_type_t>>(4 * field_size);
    }

    if (!params.disable_retraining) {
      ai_ctx->training_backend->training_thread(
          ai_ctx->design_buffer, ai_ctx->results_buffer, ai_ctx->model,
          ai_ctx->meta_params, ai_ctx->scaler, ai_ctx->data_semaphore_write,
          ai_ctx->data_semaphore_read, ai_ctx->model_semaphore,
          ai_ctx->training_is_running, params.function_code);
    }
  }

  /* SIMULATION LOOP */
  double dSimTime{0};
  for (uint32_t iter = 1; iter < maxiter + 1; iter++) {
    double start_t = MPI_Wtime();

    const double &dt = params.timesteps[iter - 1];

    std::cout << std::endl;

    /* displaying iteration number, with C++ and R iterator */
    MSG("Going through iteration " + std::to_string(iter) + "/" +
        std::to_string(maxiter));

    MSG("Current time step is " + std::to_string(dt));

    /* run transport */
    diffusion.simulate(dt);

    // validity_vector:
    //    vector with length = number of grid cells (1 indicates
    //    that values are copied into the next time step (no chemical reaction)
    //    or the ai surrogate prediction was good enoug. In both cases the
    //    workers skip the exact simulation with PHREEQC) predictor_idx
    //    ai_validity_vector
    // predictors:
    //    data frame with elements that are used as input for the ai surrogate
    //    model. The data frame can be smaller than the grid size if
    //    copy_non_reactive_regions option is enabled. In this case only the
    //    reactive data are used for ai prediction.
    // targets:
    //    data frame with elements that are used as ai outputs for the
    //    retraining step.

    if (params.ai || params.copy_non_reactive_regions) {

      chem.getField().update(diffusion.getField());

      R["TMP"] = Rcpp::wrap(chem.getField().AsVector());
      R.parseEval(
          std::string("field <- setNames(data.frame(matrix(TMP, nrow=" +
                      std::to_string(chem.getField().GetRequestedVecSize()) +
                      ")), TMP_PROPS)"));

      R.parseEval("validity_vector <- rep(FALSE, nrow(field))");

      R.parseEval("length(validity_vector)");

      if (params.copy_non_reactive_regions) {
        R.parseEval(
            "validity_vector <- field[[threshold$species]] < threshold$value");
      }
    }

    MSG("Chemistry start");
    if (params.ai) {
      double ai_start_t = MPI_Wtime();

      // deep copy field
      R.parseEval("predictors <- field");
      // get only ai related species
      R.parseEval("predictors <- predictors[ai_surrogate_species_input]");

      // remove already copied values
      R.parseEval("predictors <- predictors[!validity_vector,]");

      // store row names of predictors
      R.parseEval("predictor_idx <- row.names(predictors)");

      R.parseEval("predictors_scaled <- preprocess(predictors)");

      std::vector<std::vector<float>> predictors_scaled =
          R["predictors_scaled"];

      std::vector<float> predictions_scaled =
          ai_ctx->model.predict(predictors_scaled, params.batch_size,
                                ai_ctx->model_semaphore); // features per cell
      int n_samples = R.parseEval("nrow(predictors)");
      int n_output_features = ai_ctx->model.weight_matrices.back().cols();
      std::vector<double> predictions_scaled_double(predictions_scaled.begin(),
                                                    predictions_scaled.end());

      R["TMP"] = predictions_scaled_double;
      R["n_samples"] = n_samples;
      R["n_output"] = n_output_features;

      R.parseEval("predictions_scaled <- setNames(data.frame(matrix(TMP, "
                  "nrow=n_samples, ncol=n_output, byrow=TRUE)), "
                  "ai_surrogate_species_output)");

      R.parseEval("predictions <- postprocess(predictions_scaled)");

      MSG("AI Validation");

      R.parseEval("ai_validity_vector <- validate_predictions(predictors, "
                  "predictions) ");

      // get only indices where prediction was valid
      R.parseEval("predictor_idx <- predictor_idx[ai_validity_vector]");

      // set in global validity vector all elements to true, where prediction
      // was possible
      R.parseEval("validity_vector[as.numeric(predictor_idx)] <- TRUE");

      MSG("AI TempField");
      // maybe row.names was overwritten by function calls ??
      R.parseEval("row.names(predictions) <- row.names(predictors)");
      // subset predictions to ai_validity_vector == TRUE
      R.parseEval("predictions <- predictions[ai_validity_vector,]");
      // merge predicted values into field stored in R
      R.parseEval("field[as.numeric(row.names(predictions)),ai_surrogate_"
                  "species_output] "
                  "<- predictions");

      MSG("AI Set Field");
      Field predictions_field = Field(
          R.parseEval("nrow(field)"),
          Rcpp::as<std::vector<std::vector<double>>>(R.parseEval("field")),
          R.parseEval("colnames(field)"));

      chem.getField().update(predictions_field);

      double ai_end_t = MPI_Wtime();
      R["ai_prediction_time"] = ai_end_t - ai_start_t;
    }

    if (params.copy_non_reactive_regions || params.ai) {
      MSG("Set copied or predicted values for the workers");
      R.parseEval(
          "print(paste('Number of valid cells:', sum(validity_vector)))");
      chem.set_ai_surrogate_validity_vector(R.parseEval("validity_vector"));
    }

    chem.simulate(dt);

    /* AI surrogate iterative training*/
    if (params.ai == true && params.disable_retraining == false) {
      double ai_start_t = MPI_Wtime();

      R["TMP"] = Rcpp::wrap(chem.getField().AsVector());
      R.parseEval(
          std::string("targets <- setNames(data.frame(matrix(TMP, nrow=" +
                      std::to_string(chem.getField().GetRequestedVecSize()) +
                      ")), TMP_PROPS)"));

      R.parseEval("predictors_retraining <- "
                  "get_invalid_values(predictors_scaled, ai_validity_vector)");
      R.parseEval("targets <- "
                  "targets[as.numeric(row.names(predictors_retraining)), "
                  "ai_surrogate_species_output]");

      R.parseEval("targets_retraining <- preprocess(targets)");

      std::vector<std::vector<float>> predictors_retraining =
          R["predictors_retraining"];
      std::vector<std::vector<float>> targets_retraining =
          R["targets_retraining"];

      MSG("AI: add invalid data to buffer");

      ai_ctx->data_semaphore_write.acquire();

      if (predictors_retraining[0].size() > 0 &&
          targets_retraining[0].size() > 0 && retrain_counter == 0) {
        ai_ctx->design_buffer.addData(predictors_retraining);
        ai_ctx->results_buffer.addData(targets_retraining);
      }
      size_t elements_design_buffer =
          ai_ctx->design_buffer.getSize() /
          (predictors_retraining.size() * sizeof(float));
      size_t elements_results_buffer =
          ai_ctx->results_buffer.getSize() /
          (targets_retraining.size() * sizeof(float));

      std::cout << "design_buffer_size: " << elements_design_buffer
                << std::endl;
      std::cout << "results_buffer_size: " << elements_results_buffer
                << std::endl;

      if (elements_design_buffer >= 4 * field_size &&
          elements_results_buffer >= 4 * field_size &&
          ai_ctx->training_is_running == false && retrain_counter == 0) {
        ai_ctx->data_semaphore_read.release();
        retrain_counter++;
      } else {
        ai_ctx->data_semaphore_write.release();
      }

      double ai_end_t = MPI_Wtime();
      R["ai_training_time"] = ai_end_t - ai_start_t;
    }

    // MPI_Barrier(MPI_COMM_WORLD);
    double end_t = MPI_Wtime();
    dSimTime += end_t - start_t;
    R["totaltime"] = dSimTime;

    // MDL master_iteration_end just writes on disk state_T and
    // state_C after every iteration if the cmdline option
    // --ignore-results is not given (and thus the R variable
    // store_result is TRUE)
    call_master_iter_end(R, diffusion.getField(), chem.getField());

    diffusion.getField().update(chem.getField());

    MSG("End of *coupling* iteration " + std::to_string(iter) + "/" +
        std::to_string(maxiter));
    // MSG();
  } // END SIMULATION LOOP

  std::cout << std::endl;

  if (params.ai && !params.disable_retraining) {
    ai_ctx->training_backend->stop_training(ai_ctx->data_semaphore_read);
  }

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

std::array<double, 2> getBaseTotals(Field &&field, int root, MPI_Comm comm) {
  std::array<double, 2> base_totals;

  int rank;
  MPI_Comm_rank(comm, &rank);

  const bool is_master = root == rank;

  if (is_master) {
    const auto h_col = field["H"];
    const auto o_col = field["O"];

    base_totals[0] = *std::min_element(h_col.begin(), h_col.end());
    base_totals[1] = *std::min_element(o_col.begin(), o_col.end());
    MPI_Bcast(base_totals.data(), 2, MPI_DOUBLE, root, MPI_COMM_WORLD);
    return base_totals;
  }

  MPI_Bcast(base_totals.data(), 2, MPI_DOUBLE, root, comm);

  return base_totals;
}

bool getHasID(Field &&field, int root, MPI_Comm comm) {
  bool has_id;

  int rank;
  MPI_Comm_rank(comm, &rank);

  const bool is_master = root == rank;

  if (is_master) {
    const auto ID_field = field["ID"];

    std::set<double> unique_IDs(ID_field.begin(), ID_field.end());

    has_id = unique_IDs.size() > 1;

    MPI_Bcast(&has_id, 1, MPI_C_BOOL, root, MPI_COMM_WORLD);

    return has_id;
  }

  MPI_Bcast(&has_id, 1, MPI_C_BOOL, root, comm);

  return has_id;
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
        getBaseTotals(init_list.getInitialGrid(), 0, MPI_COMM_WORLD),
        getHasID(init_list.getInitialGrid(), 0, MPI_COMM_WORLD),
        run_params.use_dht,
        run_params.dht_size,
        run_params.dht_snaps,
        run_params.out_dir,
        run_params.use_interp,
        run_params.interp_bucket_entries,
        run_params.interp_size,
        run_params.interp_min_entries,
        run_params.ai,
        run_params.copy_non_reactive_regions};

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

      if (run_params.ai || run_params.copy_non_reactive_regions) {
        /* Incorporate ai surrogate from R */
        R.parseEvalQ(ai_surrogate_r_library);
        /* Use dht species for model input and output */
        const auto &names = init_list.getChemistryInit().dht_species.getNames();
        for (const auto &name : names) {
          std::cout << name << " ";
        }
        std::cout << "\n"; //
        const std::string ai_surrogate_input_script =
            init_list.getChemistryInit().ai_surrogate_input_script;

        MSG("AI: sourcing user-provided script");
        R.parseEvalQ(ai_surrogate_input_script);

        MSG("AI: initialize AI model");
      }

      MSG("Init done on process with rank " + std::to_string(MY_RANK));

      // MPI_Barrier(MPI_COMM_WORLD);

      DiffusionModule diffusion(init_list.getDiffusionInit(),
                                init_list.getInitialGrid());

      chemistry.masterSetField(init_list.getInitialGrid());

      Rcpp::List profiling = RunMasterLoop(R, run_params, diffusion, chemistry);

      MSG("finished simulation loop");

      R["profiling"] = profiling;
      R["setup"] = *global_rt_setup;
      R["setup$out_ext"] = run_params.out_ext;

      std::string r_vis_code;
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
