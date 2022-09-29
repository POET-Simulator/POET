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

#ifndef CHEMSIM_H
#define CHEMSIM_H

#include "DHT_Wrapper.hpp"
#include "RRuntime.hpp"
#include "SimParams.hpp"
#include "Grid.hpp"
#include <mpi.h>

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
 * @brief Base class of the chemical simulation
 *
 * Providing member functions to run an iteration and to end a simulation. Also
 * containing basic parameters for simulation.
 *
 */
class ChemSim {
public:
  /**
   * @brief Construct a new ChemSim object
   *
   * Creating a new instance of class ChemSim will just extract simulation
   * parameters from SimParams object.
   *
   * @param params SimParams object
   * @param R_ R runtime
   * @param grid_ Initialized grid
   */
  ChemSim(SimParams &params, RRuntime &R_, Grid &grid_);

  /**
   * @brief Run iteration of simulation in sequential mode
   *
   * This will call the correspondending R function slave_chemistry, followed by
   * the execution of master_chemistry.
   *
   * @todo change function name. Maybe 'slave' to 'seq'.
   *
   */
  virtual void run();

  /**
   * @brief End simulation
   *
   * End the simulation by distribute the measured runtime of simulation to the
   * R runtime.
   *
   */
  virtual void end();

  /**
   * @brief Get the Chemistry Time
   *
   * @return double Runtime of sequential chemistry simulation in seconds
   */
  double getChemistryTime();

protected:
  /**
   * @brief Current simulation time or 'age' of simulation
   *
   */
  double current_sim_time = 0;

  /**
   * @brief Current iteration
   *
   */
  int iteration = 0;

  /**
   * @brief Current simulation timestep
   *
   */
  int dt = 0;

  /**
   * @brief Rank of process in MPI_COMM_WORLD
   *
   */
  int world_rank;

  /**
   * @brief Size of communicator MPI_COMM_WORLD
   *
   */
  int world_size;

  /**
   * @brief Number of grid cells in each work package
   *
   */
  unsigned int wp_size;

  /**
   * @brief Instance of RRuntime object
   *
   */
  RRuntime &R;

  /**
   * @brief Initialized grid object
   *
   */
  Grid &grid;

  /**
   * @brief Stores information about size of the current work package
   *
   */
  std::vector<int> wp_sizes_vector;

  /**
   * @brief Absolute path to output path
   *
   */
  std::string out_dir;

  /**
   * @brief Pointer to sending buffer
   *
   */
  double *send_buffer;

  /**
   * @brief Worker struct
   *
   * This struct contains information which worker as work and which work
   * package he is working on.
   *
   */
  typedef struct {
    char has_work;
    double *send_addr;
  } worker_struct;

  /**
   * @brief Pointer to worker_struct
   *
   */
  worker_struct *workerlist;

  /**
   * @brief Pointer to mpi_buffer
   *
   * Typically for the master this is a continous C memory area containing the
   * grid. For worker the memory area will just have the size of a work package.
   *
   */
  double *mpi_buffer;

  /**
   * @brief Total chemistry runtime
   *
   */
  double chem_t = 0.f;
};

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
  ChemMaster(SimParams &params, RRuntime &R_, Grid &grid_);

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
  void run() override;

  /**
   * @brief End chemistry simulation.
   *
   * Notify the worker to finish their 'work'-loop. This is done by sending
   * every worker an empty message with the tag TAG_FINISH. Now the master will
   * receive measured times and DHT metrics from all worker one after another.
   * Finally he will write all data to the R runtime and return this function.
   *
   */
  void end() override;

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
  /**
   * @brief Print a progressbar
   *
   * Prints a progressbar to stdout according to count of processed work
   * packages in this iteration.
   *
   * @param count_pkgs Last processed index of work package
   * @param n_wp Number of work packages
   * @param barWidth Width of the progressbar/Count of characters to display the
   * bar
   */
  void printProgressbar(int count_pkgs, int n_wp, int barWidth = 70);

  /**
   * @brief Start send loop
   *
   * Send a work package to every free worker, which are noted in a worker
   * struct. After a work package was sent move pointer on work grid to the next
   * work package. Use MPI_Send to transfer work package to worker.
   *
   * @param pkg_to_send Pointer to variable containing how much work packages
   * are still to send
   * @param count_pkgs Pointer to variable indexing the current work package
   * @param free_workers Pointer to variable with the count of free workers
   */
  void sendPkgs(int &pkg_to_send, int &count_pkgs, int &free_workers);

  /**
   * @brief Start recv loop
   *
   * Receive processed work packages by worker. This is done by first probing
   * for a message. If a message is receivable, receive it and put result into
   * respective memory area. Continue, but now with a non blocking MPI_Probe. If
   * a message is receivable or if no work packages are left to send, receive
   * it. Otherwise or if all remaining work packages are received exit loop.
   *
   * @param pkg_to_recv Pointer to variable counting the to receiving work
   * packages
   * @param to_send Bool indicating if there are still work packages to send
   * @param free_workers Pointer to worker to variable holding the number of
   * free workers
   */
  void recvPkgs(int &pkg_to_recv, bool to_send, int &free_workers);

  /**
   * @brief Indicating usage of DHT
   *
   */
  bool dht_enabled;

  /**
   * @brief Default number of grid cells in each work package
   *
   */
  unsigned int wp_size;

  /**
   * @brief Pointer to current to be processed work package
   *
   */
  double *work_pointer;

  /**
   * @brief Time spent in send loop
   *
   */
  double send_t = 0.f;

  /**
   * @brief Time spent in recv loop
   *
   */
  double recv_t = 0.f;

  /**
   * @brief Time master is idling in MPI_Probe
   *
   */
  double master_idle = 0.f;

  /**
   * @brief Time spent in send/recv loop
   *
   */
  double worker_t = 0.f;

  /**
   * @brief Time spent in sequential chemistry part
   *
   */
  double chem_master = 0.f;

  /**
   * @brief Time spent in sequential instructions
   *
   */
  double seq_t = 0.f;
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
  ChemWorker(SimParams &params, RRuntime &R_, Grid &grid_, MPI_Comm dht_comm);

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
  /**
   * @brief Evaluating message to receive as work package
   *
   * These steps are done:
   *
   * -# Receive message
   * -# Evaluate message header containing information about work package size,
   * current iteration and timestep, simulation age
   * -# if DHT is enabled check DHT for previously simulated results
   * -# run simulation of work package
   * -# send results back to master
   * -# if DHT is enabled write simulated grid points to DHT
   *
   * @param probe_status message status of produced by MPI_Probe in loop
   */
  void doWork(MPI_Status &probe_status);

  /**
   * @brief Action to do after iteration
   *
   * If DHT is enabled print statistics and if dht_snaps is set to 2 write DHT
   * snapshots.
   *
   */
  void postIter();

  /**
   * @brief Message tag evaluates to TAG_FINISH
   *
   * Send all the collected timings and (possbile) DHT metrics to the master.
   *
   */
  void finishWork();

  /**
   * @brief Write DHT snapshot
   *
   */
  void writeFile();

  /**
   * @brief Read DHT snapshot
   *
   */
  void readFile();

  /**
   * @brief Indicates usage of DHT
   *
   */
  bool dht_enabled;

  /**
   * @brief Boolean if dt differs between iterations
   *
   */
  bool dt_differ;

  /**
   * @brief Value between 0 and 2, indicating when to write DHT snapshots
   *
   */
  int dht_snaps;

  /**
   * @brief Absolute path to DHT snapshot file
   *
   */
  std::string dht_file;

  /**
   * @brief Count of bytes each process should allocate for the DHT
   *
   */
  unsigned int dht_size_per_process;

  /**
   * @brief Indicates which grid cells were previously simulated and don't need
   * to be simulated now
   *
   */
  std::vector<bool> dht_flags;

  /**
   * @brief simulated results are stored here
   *
   */
  double *mpi_buffer_results;

  /**
   * @brief Instance of DHT_Wrapper
   *
   */
  DHT_Wrapper *dht;

  /**
   * @brief Array to store timings
   *
   * The values are stored in following order
   *
   * -# PHREEQC time
   * -# DHT_get time
   * -# DHT_fill time
   *
   */
  double timing[3];

  /**
   * @brief Time worker is idling in MPI_Probe
   *
   */
  double idle_t = 0.f;

  /**
   * @brief Count of PHREEQC calls
   *
   */
  int phreeqc_count = 0;
};
} // namespace poet
#endif // CHEMSIM_H
