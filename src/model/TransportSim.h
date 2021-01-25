#ifndef TRANSPORT_SIM_H
#define TRANSPORT_SIM_H

#include <RRuntime.h>

namespace poet {
class TransportSim {
 public:
  TransportSim(RRuntime &R);

  void run();
  void end();

  double getTransportTime();

 private:
  RRuntime &R;

  double transport_t = 0.f;
};
}  // namespace poet

#endif  // TRANSPORT_SIM_H