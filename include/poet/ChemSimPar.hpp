/*
** Copyright (C) 2018-2021 Alexander Lindemann, Max Luebke (University of
** Potsdam)
**
** Copyright (C) 2018-2022 Marco De Lucia, Max Luebke (GFZ Potsdam)
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

#ifndef CHEMSIMPAR_H
#define CHEMSIMPAR_H

#include "ChemSimSeq.hpp"
#include "DHT_Wrapper.hpp"
#include "Grid.hpp"
#include "RInside.h"
#include "SimParams.hpp"
#include <array>
#include <bits/stdint-uintn.h>
#include <cstdint>
#include <mpi.h>

#include <string>
#include <vector>

/** Number of data elements that are kept free at each work package */
#define BUFFER_OFFSET 5

/** Message tag indicating work */
#define TAG_WORK 42
/** Message tag indicating to finish loop */
#define TAG_FINISH 43
/** Message tag indicating timing profiling */
#define TAG_TIMING 44
/** Message tag indicating collecting DHT performance */
#define TAG_DHT_PERF 45
/** Message tag indicating simulation reached the end of an itertation */
#define TAG_DHT_ITER 47

namespace poet {

/**
 * @brief Class providing execution of master chemistry
 *
 * Providing member functions to run an iteration and to end a simulation. Also
 * a loop to send and recv pkgs from workers is implemented.
 *
 */
class ChemMaster : public ChemSim {
public:
  /**
   * @brief Construct a new ChemMaster object
   *
   * The following steps are executed to create a new object of ChemMaster:
   * -# all needed simulation parameters are extracted
   * -# memory is allocated
   * -# distribution of work packages is calculated
   *
   * @param params Simulation parameters as SimParams object
   * @param R_ R runtime
   * @param grid_ Grid object
   */
  ChemMaster(SimParams &params, RInside &R_, Grid &grid_);

  /**
   * @brief Destroy the ChemMaster object
   *
   * By freeing ChemMaster all buffers allocated in the Constructor are freed.
   *
   */
  ~ChemMaster();

  /**
   * @brief Run iteration of simulation in parallel mode
   *
   * To run the chemistry simulation parallel following steps are done:
   *
   * -# 'Shuffle' the grid by previously calculated distribution of work
   * packages. Convert R grid to C memory area.
   * -# Start the send/recv loop.
   * Detailed description in sendPkgs respectively in recvPkgs.
   * -# 'Unshuffle'
   * the grid and convert C memory area to R grid.
   * -# Run 'master_chemistry'
   *
   * The main tasks are instrumented with time measurements.
   *
   */
  void masterLoop(double dt);

  /**
   * @brief End chemistry simulation.
   *
   * Notify the worker to finish their 'work'-loop. This is done by sending
   * every worker an empty message with the tag TAG_FINISH. Now the master will
   * receive measured times and DHT metrics from all worker one after another.
   * Finally he will write all data to the R runtime and return this function.
   *
   */
  void endLoop();

  /**
   * @brief Get the send time
   *
   * Time spent in send loop.
   *
   * @return double sent time in seconds
   */
  double getSendTime();

  /**
   * @brief Get the recv time
   *
   * Time spent in recv loop.
   *
   * @return double recv time in seconds
   */
  double getRecvTime();

  /**
   * @brief Get the idle time
   *
   * Time master was idling in MPI_Probe of recv loop.
   *
   * @return double idle time in seconds
   */
  double getIdleTime();

  /**
   * @brief Get the Worker time
   *
   * Time spent in whole send/recv loop.
   *
   * @return double worker time in seconds
   */
  double getWorkerTime();

  /**
   * @brief Get the ChemMaster time
   *
   * Time spent in 'master_chemistry' R function.
   *
   * @return double ChemMaster time in seconds
   */
  double getChemMasterTime();

  /**
   * @brief Get the sequential time
   *
   * Time master executed code which must be run sequential.
   *
   * @return double seqntial time in seconds.
   */
  double getSeqTime();

private:
  void printProgressbar(int count_pkgs, int n_wp, int barWidth = 70);
  inline void sendPkgs(int &pkg_to_send, int &count_pkgs, int &free_workers,
                       double *work_pointer, const double &dt);
  inline void recvPkgs(int &pkg_to_recv, bool to_send, int &free_workers);
  void shuffleField(const std::vector<double> &in_field, uint32_t size_per_prop,
                    uint32_t prop_count, double *out_buffer);
  void unshuffleField(const double *in_buffer, uint32_t size_per_prop,
                      uint32_t prop_count, std::vector<double> &out_field);

  bool dht_enabled;
  unsigned int wp_size;

  double send_t = 0.f;
  double recv_t = 0.f;
  double master_idle = 0.f;
  double worker_t = 0.f;
  double chem_master = 0.f;
  double seq_t = 0.f;

  typedef struct {
    char has_work;
    double *send_addr;
  } worker_struct;

  worker_struct *workerlist;
  double *send_buffer;
  double *mpi_buffer;
};

/**
 * @brief Class providing execution of worker chemistry
 *
 * Providing mainly a function to loop and wait for messages from the master.
 *
 */
class ChemWorker : public ChemSim {
public:
  /**
   * @brief Construct a new ChemWorker object
   *
   * The following steps are executed to create a new object of ChemWorker:
   * -# all needed simulation parameters are extracted
   * -# memory is allocated
   * -# Preparetion to create a DHT
   * -# and finally create a new DHT_Wrapper
   *
   * @param params Simulation parameters as SimParams object
   * @param R_ R runtime
   * @param grid_ Grid object
   * @param dht_comm Communicator addressing all processes marked as worker
   */
  ChemWorker(SimParams &params, RInside &R_, Grid &grid_, MPI_Comm dht_comm);

  /**
   * @brief Destroy the ChemWorker object
   *
   * Therefore all buffers are freed and the DHT_Wrapper object is destroyed.
   *
   */
  ~ChemWorker();

  /**
   * @brief Start the 'work' loop
   *
   * Loop in an endless loop. At the beginning probe for a message from the
   * master process. If there is a receivable message evaluate the message tag.
   *
   */
  void loop();

private:
  void doWork(MPI_Status &probe_status);
  void postIter();
  void finishWork();
  void writeFile();
  void readFile();

  bool dht_enabled;
  bool dt_differ;
  int dht_snaps;
  std::string dht_file;
  unsigned int dht_size_per_process;
  std::vector<bool> dht_flags;
  double *mpi_buffer_results;
  DHT_Wrapper *dht;

  std::array<double, 3> timing;
  double idle_t = 0.f;
  int phreeqc_count = 0;

  double *mpi_buffer;
};
} // namespace poet
#endif // CHEMSIMPAR_H
