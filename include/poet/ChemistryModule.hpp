//  Time-stamp: "Last modified 2023-04-24 14:30:06 mluebke"

#ifndef CHEMISTRYMODULE_H_
#define CHEMISTRYMODULE_H_

#include "Field.hpp"
#include "IrmResult.h"
#include "PhreeqcRM.h"
#include "poet/DHT_Wrapper.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <map>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>

namespace poet {
/**
 * \brief Wrapper around PhreeqcRM to provide POET specific parallelization with
 * easy access.
 */
class ChemistryModule : public PhreeqcRM {
public:
#ifdef POET_USE_PRM
  /**
   * Creates a new instance of Chemistry module with given grid cell count and
   * using MPI communicator.
   *
   * Constructor is just a wrapper around PhreeqcRM's constructor, so just calls
   * the base class constructor.
   *
   * \param nxyz Count of grid cells.
   * \param communicator MPI communicator where work will be distributed.
   */
  ChemistryModule(uint32_t nxyz, MPI_Comm communicator)
      : PhreeqcRM(nxyz, communicator), n_cells(nxyz){};
#else
  /**
   * Creates a new instance of Chemistry module with given grid cell count, work
   * package size and communicator.
   *
   * This constructor shall only be called by the master. To create workers, see
   * ChemistryModule::createWorker .
   *
   * When the use of parallelization is intended, the nxyz value shall be set to
   * 1 to save memory and only one node is needed for initialization.
   *
   * \param nxyz Count of grid cells to allocate and initialize for each
   * process. For parellel use set to 1 at the master.
   * \param wp_size Count of grid cells to fill each work package at maximum.
   * \param communicator MPI communicator to distribute work in.
   */
  ChemistryModule(uint32_t nxyz, uint32_t wp_size, MPI_Comm communicator);

  /**
   * Deconstructor, which frees DHT data structure if used.
   */
  ~ChemistryModule();
#endif

  /**
   * Parses the input script and extract information needed during runtime.
   *
   * **Only run by master**.
   *
   * Database must be loaded beforehand.
   *
   * \param input_script_path Path to input script to parse.
   */
  void RunInitFile(const std::string &input_script_path);

  /**
   * Run the chemical simulation with parameters set.
   */
  void RunCells();

  /**
   * Returns the chemical field.
   */
  auto GetField() const { return this->chem_field; }
  /**
   * Returns all known species names, including not only aqueous species, but
   * also equilibrium, exchange, surface and kinetic reactants.
   */
  auto GetPropNames() const { return this->prop_names; }

  /**
   * Return the accumulated runtime in seconds for chemical simulation.
   */
  auto GetChemistryTime() const { return this->chem_t; }

#ifndef POET_USE_PRM

  /**
   * Create a new worker instance inside given MPI communicator.
   *
   * Wraps communication needed before instanciation can take place.
   *
   * \param communicator MPI communicator to distribute work in.
   *
   * \returns A worker instance with fixed work package size.
   */
  static ChemistryModule createWorker(MPI_Comm communicator);

  /**
   * Default work package size.
   */
  static constexpr uint32_t CHEM_DEFAULT_WP_SIZE = 5;

  /**
   * Intended to alias input parameters for grid initialization with a single
   * value per species.
   */
  using SingleCMap = std::unordered_map<std::string, double>;

  /**
   * Intended to alias input parameters for grid initialization with mutlitple
   * values per species.
   */
  using VectorCMap = std::unordered_map<std::string, std::vector<double>>;

  /**
   * Enumerating DHT file options
   */
  enum {
      DHT_FILES_DISABLED = 0, //!< disabled file output
      DHT_FILES_SIMEND,       //!< only output of snapshot after simulation
      DHT_FILES_ITEREND       //!< output snapshots after each iteration
  };

  /**
   * **This function has to be run!**
   *
   * Merge initial values from existing module with the chemistry module and set
   * according internal variables.
   *
   * \param other Field to merge chemistry with. Most likely it is something
   * like the diffusion field.
   */
  void initializeField(const Field &other);

  /**
   * **Only called by workers!** Start the worker listening loop.
   */
  void WorkerLoop();

  /**
   * **Called by master** Advise the workers to break the loop.
   */
  void MasterLoopBreak();

  /**
   * **Master only** Enables DHT usage for the workers.
   *
   * If called multiple times enabling the DHT, the current state of DHT will be
   * overwritten with a new instance of the DHT.
   *
   * \param enable Enables or disables the usage of the DHT.
   * \param size_mb Size in megabyte to allocate for the DHT if enabled.
   * \param key_species Name of species to be used for key creation.
   */
  void SetDHTEnabled(bool enable, uint32_t size_mb,
                     const std::vector<std::string> &key_species);
  /**
   * **Master only** Set DHT snapshots to specific mode.
   *
   * \param type DHT snapshot mode.
   * \param out_dir Path to output DHT snapshots.
   */
  void SetDHTSnaps(int type, const std::string &out_dir);
  /**
   * **Master only** Set the vector with significant digits to round before
   * inserting into DHT.
   *
   * \param signif_vec Vector defining significant digit for each species. Order
   * is defined by prop_type vector (ChemistryModule::GetPropNames).
   */
  void SetDHTSignifVector(std::vector<uint32_t> &signif_vec);
  /**
   * **Master only** Set the DHT rounding type of each species. See
   * DHT_PROP_TYPES enumemartion for explanation.
   *
   * \param proptype_vec Vector defining DHT prop type for each species.
   */
  void SetDHTPropTypeVector(std::vector<uint32_t> proptype_vec);
  /**
   * **Master only** Load the state of the DHT from given file.
   *
   * \param input_file File to load data from.
   */
  void ReadDHTFile(const std::string &input_file);

  /**
   * **Master only** Return count of grid cells.
   */
  auto GetNCells() const { return this->n_cells; }
  /**
   * **Master only** Return work package size.
   */
  auto GetWPSize() const { return this->wp_size; }
  /**
   * **Master only** Return the time in seconds the master spent waiting for any
   * free worker.
   */
  auto GetMasterIdleTime() const { return this->idle_t; }
  /**
   * **Master only** Return the time in seconds the master spent in sequential
   * part of the simulation, including times for shuffling/unshuffling field
   * etc.
   */
  auto GetMasterSequentialTime() const { return this->seq_t; }
  /**
   * **Master only** Return the time in seconds the master spent in the
   * send/receive loop.
   */
  auto GetMasterLoopTime() const { return this->send_recv_t; }

  /**
   * **Master only** Collect and return all accumulated timings recorded by
   * workers to run Phreeqc simulation.
   *
   * \return Vector of all accumulated Phreeqc timings.
   */
  std::vector<double> GetWorkerPhreeqcTimings() const;
  /**
   * **Master only** Collect and return all accumulated timings recorded by
   * workers to get values from the DHT.
   *
   * \return Vector of all accumulated DHT get times.
   */
  std::vector<double> GetWorkerDHTGetTimings() const;
  /**
   * **Master only** Collect and return all accumulated timings recorded by
   * workers to write values to the DHT.
   *
   * \return Vector of all accumulated DHT fill times.
   */
  std::vector<double> GetWorkerDHTFillTimings() const;
  /**
   * **Master only** Collect and return all accumulated timings recorded by
   * workers waiting for work packages from the master.
   *
   * \return Vector of all accumulated waiting times.
   */
  std::vector<double> GetWorkerIdleTimings() const;

  /**
   * **Master only** Collect and return DHT hits of all workers.
   *
   * \return Vector of all count of DHT hits.
   */
  std::vector<uint32_t> GetWorkerDHTHits() const;

  /**
   * **Master only** Collect and return DHT misses of all workers.
   *
   * \return Vector of all count of DHT misses.
   */
  std::vector<uint32_t> GetWorkerDHTMiss() const;

  /**
   * **Master only** Collect and return DHT evictions of all workers.
   *
   * \return Vector of all count of DHT evictions.
   */
  std::vector<uint32_t> GetWorkerDHTEvictions() const;

  /**
   * **Master only** Returns the current state of the chemical field.
   *
   * \return Reference to the chemical field.
   */
  Field &getField() { return this->chem_field; }
#endif
protected:
#ifdef POET_USE_PRM
  void PrmToPoetField(std::vector<double> &field);
#else
  enum {
    CHEM_INIT,
    CHEM_WP_SIZE,
    CHEM_INIT_SPECIES,
    CHEM_DHT_ENABLE,
    CHEM_DHT_SIGNIF_VEC,
    CHEM_DHT_PROP_TYPE_VEC,
    CHEM_DHT_SNAPS,
    CHEM_DHT_READ_FILE,
    CHEM_WORK_LOOP,
    CHEM_PERF,
    CHEM_BREAK_MAIN_LOOP
  };

  enum { LOOP_WORK, LOOP_END };

  enum {
    WORKER_PHREEQC,
    WORKER_DHT_GET,
    WORKER_DHT_FILL,
    WORKER_IDLE,
    WORKER_DHT_HITS,
    WORKER_DHT_MISS,
    WORKER_DHT_EVICTIONS
  };

  struct worker_s {
    double phreeqc_t = 0.;
    double dht_get = 0.;
    double dht_fill = 0.;
    double idle_t = 0.;
  };

  struct worker_info_s {
    char has_work = 0;
    double *send_addr;
  };

  using worker_list_t = std::vector<struct worker_info_s>;
  using workpointer_t = std::vector<double>::iterator;

  void MasterRunParallel();
  void MasterRunSequential();

  void MasterSendPkgs(worker_list_t &w_list, workpointer_t &work_pointer,
                      int &pkg_to_send, int &count_pkgs, int &free_workers,
                      double dt, uint32_t iteration,
                      const std::vector<uint32_t> &wp_sizes_vector);
  void MasterRecvPkgs(worker_list_t &w_list, int &pkg_to_recv, bool to_send,
                      int &free_workers);

  std::vector<double> MasterGatherWorkerTimings(int type) const;
  std::vector<uint32_t> MasterGatherWorkerMetrics(int type) const;

  void WorkerProcessPkgs(struct worker_s &timings, uint32_t &iteration);

  void WorkerDoWork(MPI_Status &probe_status, int double_count,
                    struct worker_s &timings);
  void WorkerPostIter(MPI_Status &prope_status, uint32_t iteration);
  void WorkerPostSim(uint32_t iteration);

  void WorkerWriteDHTDump(uint32_t iteration);
  void WorkerReadDHTDump(const std::string &dht_input_file);

  void WorkerPerfToMaster(int type, const struct worker_s &timings);
  void WorkerMetricsToMaster(int type);

  IRM_RESULT WorkerRunWorkPackage(std::vector<double> &vecWP,
                                  std::vector<std::uint32_t> &vecMapping,
                                  double dSimTime, double dTimestep);

  std::vector<uint32_t> CalculateWPSizesVector(uint32_t n_cells,
                                               uint32_t wp_size) const;

  std::vector<double> shuffleField(const std::vector<double> &in_field,
                                   uint32_t size_per_prop, uint32_t prop_count,
                                   uint32_t wp_count);
  void unshuffleField(const std::vector<double> &in_buffer,
                      uint32_t size_per_prop, uint32_t prop_count,
                      uint32_t wp_count, std::vector<double> &out_field);
  std::vector<std::uint32_t>
  parseDHTSpeciesVec(const std::vector<std::string> &species_vec) const;

  int comm_size, comm_rank;
  MPI_Comm group_comm;

  bool is_sequential;
  bool is_master;

  uint32_t wp_size{CHEM_DEFAULT_WP_SIZE};
  bool dht_enabled{false};
  int dht_snaps_type{DHT_FILES_DISABLED};
  std::string dht_file_out_dir;

  poet::DHT_Wrapper *dht = nullptr;

  static constexpr uint32_t BUFFER_OFFSET = 5;

  inline void ChemBCast(void *buf, int count, MPI_Datatype datatype) const {
    MPI_Bcast(buf, count, datatype, 0, this->group_comm);
  }

  inline void PropagateFunctionType(int &type) const {
    ChemBCast(&type, 1, MPI_INT);
  }
  double simtime = 0.;
  double idle_t = 0.;
  double seq_t = 0.;
  double send_recv_t = 0.;

  std::array<double, 2> base_totals{0};

#endif

  double chem_t = 0.;

  uint32_t n_cells = 0;
  uint32_t prop_count = 0;
  std::vector<std::string> prop_names;

  Field chem_field{0};

  static constexpr int MODULE_COUNT = 5;

  std::array<std::uint32_t, MODULE_COUNT> speciesPerModule{};
};
} // namespace poet

#endif // CHEMISTRYMODULE_H_
