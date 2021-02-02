#ifndef TRANSPORT_SIM_H
#define TRANSPORT_SIM_H

#include <RRuntime.h>

namespace poet {
/**
 * @brief Class describing transport simulation
 *
 * Offers simple methods to run an iteration and end the simulation.
 *
 */
class TransportSim {
 public:
  /**
   * @brief Construct a new TransportSim object
   *
   * The instance will only be initialized with given R object.
   *
   * @param R RRuntime object
   */
  TransportSim(RRuntime &R);

  /**
   * @brief Run simulation for one iteration
   *
   * This will simply call the R function 'master_advection'
   *
   */
  void run();

  /**
   * @brief End simulation
   *
   * All measured timings are distributed to the R runtime
   *
   */
  void end();

  /**
   * @brief Get the transport time
   *
   * @return double time spent in transport
   */
  double getTransportTime();

 private:
  /**
   * @brief Instance of RRuntime
   *
   */
  RRuntime &R;

  /**
   * @brief time spent for transport
   *
   */
  double transport_t = 0.f;
};
}  // namespace poet

#endif  // TRANSPORT_SIM_H