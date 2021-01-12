#ifndef TRANSPORT_SIM_H
#define TRANSPORT_SIM_H

#include "../util/RRuntime.h"

namespace poet {
class TransportSim {
 public:
  TransportSim(RRuntime &R);

  void runIteration();

  double getTransportTime();

 private:
  RRuntime &R;

  double transport_t = 0.f;
};
}  // namespace poet

#endif  // TRANSPORT_SIM_H