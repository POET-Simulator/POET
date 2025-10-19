#include "ControlModule.hpp"
#include "IO/Datatypes.hpp"
#include "IO/HDF5Functions.hpp"
#include "IO/StatsIO.hpp"
#include <cmath>

bool poet::ControlModule::isControlIteration(uint32_t iter) {
  control_interval_enabled = (iter % control_interval == 0);
  if (control_interval_enabled) {
    MSG("[Control] Control interval triggered at iteration " +
        std::to_string(iter));
  }
  return control_interval_enabled;
}

void poet::ControlModule::beginIteration() {
  if (rollback_enabled) {
    if (sur_disabled_counter > 0) {
      sur_disabled_counter--;
      MSG("Rollback counter: " + std::to_string(sur_disabled_counter));
    } else {
      rollback_enabled = false;
    }
  }
}

void poet::ControlModule::endIteration(uint32_t iter) {
  /* Writing a checkpointing */
  if (checkpoint_interval > 0 && iter % checkpoint_interval == 0) {
    MSG("Writing checkpoint of iteration " + std::to_string(iter));
    write_checkpoint(out_dir, "checkpoint" + std::to_string(iter) + ".hdf5",
                     {.field = chem->getField(), .iteration = iter});
  }

  /* Control Logic*/
  if (control_interval_enabled && !rollback_enabled) {
    writeStatsToCSV(error_history, species_names, out_dir,
                    "stats_overview");

    if (triggerRollbackIfExceeded(*chem, *params, iter)) {
      rollback_enabled = true;
      rollback_counter++;
      sur_disabled_counter = control_interval;
      MSG("Interpolation disabled for the next " +
          std::to_string(control_interval) + ".");
    }
  }
}

void poet::ControlModule::BCastControlFlags() {
  int interp_flag = rollback_enabled ? 0 : 1;
  int dht_fill_flag = rollback_enabled ? 1 : 0;
  chem->ChemBCast(&interp_flag, 1, MPI_INT);
  chem->ChemBCast(&dht_fill_flag, 1, MPI_INT);
}

bool poet::ControlModule::triggerRollbackIfExceeded(ChemistryModule &chem,
                                                    RuntimeParameters &params,
                                                    uint32_t &iter) {

  if (error_history.empty()) {
    MSG("No error history yet; skipping rollback check.");
    return false;
  }

  const auto &mape = chem.error_history.back().mape;
  const auto &props = chem.getField().GetProps();

  for (uint32_t i = 0; i < params.mape_threshold.size(); ++i) {
    // Skip invalid entries
    if (mape[i] == 0) {
      continue;
    }
    bool mape_exceeded = mape[i] > params.mape_threshold[i];

    if (mape_exceeded) {
      uint32_t rollback_iter = ((iter - 1) / params.checkpoint_interval) *
                               params.checkpoint_interval;

      MSG("[THRESHOLD EXCEEDED] " + props[i] +
          " has MAPE = " + std::to_string(mape[i]) +
          " exceeding threshold = " + std::to_string(params.mape_threshold[i]) +
          " → rolling back to iteration " + std::to_string(rollback_iter));

      Checkpoint_s checkpoint_read{.field = chem.getField()};
      read_checkpoint(params.out_dir,
                      "checkpoint" + std::to_string(rollback_iter) + ".hdf5",
                      checkpoint_read);
      iter = checkpoint_read.iteration;
      return true;
    }
  }
  MSG("All species are within their MAPE and RRMSE thresholds.");
  return false;
}

void poet::ControlModule::computeSpeciesErrors(
    const std::vector<double> &reference_values,
    const std::vector<double> &surrogate_values, uint32_t size_per_prop) {

  SimulationErrorStats species_error_stats(species_count, params->global_iter,
                                           rollback_counter);

  for (uint32_t i = 0; i < species_count; ++i) {
    double err_sum = 0.0;
    double sqr_err_sum = 0.0;
    uint32_t base_idx = i * size_per_prop;

    for (uint32_t j = 0; j < size_per_prop; ++j) {
      const double ref_value = reference_values[base_idx + j];
      const double sur_value = surrogate_values[base_idx + j];

      if (ref_value == 0.0) {
        if (sur_value != 0.0) {
          err_sum += 1.0;
          sqr_err_sum += 1.0;
        }
        // Both zero: skip
      } else {
        double alpha = 1.0 - (sur_value / ref_value);
        err_sum += std::abs(alpha);
        sqr_err_sum += alpha * alpha;
      }
    }

    species_error_stats.mape[i] = 100.0 * (err_sum / size_per_prop);
    species_error_stats.rrmse[i] =
        (size_per_prop > 0) ? std::sqrt(sqr_err_sum / size_per_prop) : 0.0;
  }
  error_history.push_back(species_error_stats);
}
