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
  ControlModule(RuntimeParameters *run_params, ChemistryModule *chem_module)
      : params(run_params), chem(chem_module) {};

  /* Control configuration*/
  std::vector<std::string> species_names;
  uint32_t species_count = 0;
  std::string out_dir;

  bool rollback_enabled = false;
  bool control_interval_enabled = false;

  std::uint32_t global_iter = 0;
  std::uint32_t sur_disabled_counter = 0;
  std::uint32_t rollback_counter = 0;
  std::uint32_t checkpoint_interval = 0;
  std::uint32_t control_interval = 0;

  std::vector<double> mape_threshold;
  std::vector<double> rrmse_threshold;

  double ctrl_t = 0.;
  double bcast_ctrl_t = 0.;
  double recv_ctrl_t = 0.;

  /* Buffer for shuffled surrogate data */
  std::vector<double> sur_shuffled;

  bool isControlIteration(uint32_t iter);

  void beginIteration();

  void endIteration(uint32_t iter);

  void BCastControlFlags();

  bool triggerRollbackIfExceeded(ChemistryModule &chem,
                                 RuntimeParameters &params, uint32_t &iter);

  struct SimulationErrorStats {
    std::vector<double> mape;
    std::vector<double> rrmse;
    uint32_t iteration; // iterations in simulation after rollbacks
    uint32_t rollback_count;

    SimulationErrorStats(size_t species_count, uint32_t iter, uint32_t counter)
        : mape(species_count, 0.0), rrmse(species_count, 0.0), iteration(iter),
          rollback_count(counter) {}
  };

  static void computeSpeciesErrors(const std::vector<double> &reference_values,
                                   const std::vector<double> &surrogate_values,
                                   uint32_t size_per_prop);

  std::vector<SimulationErrorStats> error_history;

  struct ControlSetup {
    std::string out_dir;
    std::uint32_t checkpoint_interval;
    std::uint32_t control_interval;
    std::uint32_t species_count;

    std::vector<std::string> species_names;
    std::vector<double> mape_threshold;
    std::vector<double> rrmse_threshold;
  };

  void enableControlLogic(const ControlSetup &setup) {
    out_dir = setup.out_dir;
    checkpoint_interval = setup.checkpoint_interval;
    control_interval = setup.control_interval;
    species_count = setup.species_count;

    species_names = setup.species_names;
    mape_threshold = setup.mape_threshold;
    rrmse_threshold = setup.rrmse_threshold;
  }

  /* Profiling getters */
  auto GetMasterCtrlLogicTime() const { return this->ctrl_t; }

  auto GetMasterCtrlBcastTime() const { return this->bcast_ctrl_t; }

  auto GetMasterRecvCtrlLogicTime() const { return this->recv_ctrl_t; }

private:
  RuntimeParameters *params;
  ChemistryModule *chem;
};

} // namespace poet

#endif // CONTROLMODULE_H_