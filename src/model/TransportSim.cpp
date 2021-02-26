/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2021 Marco De Lucia (GFZ Potsdam)
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
