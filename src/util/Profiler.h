#ifndef PROFILER_H
#define PROFILER_H

#include "../model/ChemSim.h"
#include "../model/TransportSim.h"
#include "RRuntime.h"
#include "SimParams.h"

namespace poet {
class Profiler {
 public:
  static int startProfiling(t_simparams &params, ChemMaster &chem,
                            TransportSim &trans, RRuntime &R, double simtime);
};

}  // namespace poet

#endif  // PROFILER_H