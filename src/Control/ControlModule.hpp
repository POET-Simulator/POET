#ifndef CONTROLMODULE_H_
#define CONTROLMODULE_H_

#include "Base/Macros.hpp"
#include "Chemistry/ChemistryModule.hpp"
#include "poet.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace poet {

class ChemistryModule;

class ControlModule {

public:
  /* Control configuration*/

  // std::uint32_t global_iter = 0;
  // std::uint32_t sur_disabled_counter = 0;
  // std::uint32_t rollback_counter = 0;

  void updateControlIteration(const uint32_t iter);

  auto GetGlobalIteration() const noexcept { return global_iteration; }

  // void beginIteration();

  void endIteration(const uint32_t iter);

 // void BCastControlFlags();

  //bool triggerRollbackIfExceeded(ChemistryModule &chem,
    //                             RuntimeParameters &params, uint32_t &iter);

  struct SimulationErrorStats {
    std::vector<double> mape;
    std::vector<double> rrmse;
    uint32_t iteration; // iterations in simulation after rollbacks
    uint32_t rollback_count;

    SimulationErrorStats(uint32_t species_count, uint32_t iter, uint32_t counter)
        : mape(species_count, 0.0), rrmse(species_count, 0.0), iteration(iter),
          rollback_count(counter) {}
  };

  void computeSpeciesErrors(const std::vector<double> &reference_values,
                                   const std::vector<double> &surrogate_values,
                                   const uint32_t size_per_prop);

  std::vector<SimulationErrorStats> error_history;

  struct ControlSetup {
    std::string out_dir;
    std::uint32_t checkpoint_interval;
    std::uint32_t control_interval;
    std::vector<std::string> species_names;
    std::vector<double> mape_threshold;
  };

  void enableControlLogic(const ControlSetup &setup) {
    this->out_dir = setup.out_dir;
    this->checkpoint_interval = setup.checkpoint_interval;
    this->control_interval = setup.control_interval;
    this->species_names = setup.species_names;
    this->mape_threshold = setup.mape_threshold;
  }

  bool GetControlIntervalEnabled() const {
    return this->control_interval_enabled;
  }

  auto GetControlInterval() const { return this->control_interval; }

  std::vector<double> GetMapeThreshold() const { return this->mape_threshold; }

  /* Profiling getters */
  auto GetMasterCtrlLogicTime() const { return this->ctrl_time; }

  auto GetMasterCtrlBcastTime() const { return this->bcast_ctrl_time; }

  auto GetMasterRecvCtrlLogicTime() const { return this->recv_ctrl_time; }

private:
  bool rollback_enabled = false;
  bool control_interval_enabled = false;

  poet::ChemistryModule *chem = nullptr;

  std::uint32_t checkpoint_interval = 0;
  std::uint32_t control_interval = 0;
  std::uint32_t global_iteration = 0;
  std::vector<double> mape_threshold;

  std::vector<std::string> species_names;
  std::string out_dir;

  double ctrl_time = 0.0;
  double bcast_ctrl_time = 0.0;
  double recv_ctrl_time = 0.0;

  /* Buffer for shuffled surrogate data */
  std::vector<double> sur_shuffled;
};

} // namespace poet

#endif // CONTROLMODULE_H_