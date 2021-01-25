#include "TransportSim.h"

#include <mpi.h>

using namespace poet;

TransportSim::TransportSim(RRuntime &R_) : R(R_) {}

void TransportSim::run() {
  double sim_a_transport, sim_b_transport;

  sim_b_transport = MPI_Wtime();
  R.parseEvalQ("mysetup <- master_advection(setup=mysetup)");
  sim_a_transport = MPI_Wtime();

  transport_t += sim_a_transport - sim_b_transport;
}

void TransportSim::end() {
  R["simtime_transport"] = transport_t;
  R.parseEvalQ("profiling$simtime_transport <- simtime_transport");
}

double TransportSim::getTransportTime() { return this->transport_t; }
