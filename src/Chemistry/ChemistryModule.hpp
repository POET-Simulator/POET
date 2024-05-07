
#ifndef CHEMISTRYMODULE_H_
#define CHEMISTRYMODULE_H_

#include "DataStructures/Field.hpp"
#include "DataStructures/NamedVector.hpp"

#include "ChemistryDefs.hpp"

#include "Init/InitialList.hpp"
#include "SurrogateModels/DHT_Wrapper.hpp"
#include "SurrogateModels/Interpolation.hpp"

#include "PhreeqcEngine.hpp"
#include <array>
#include <cstdint>
#include <map>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

namespace poet {
/**
 * \brief Wrapper around PhreeqcRM to provide POET specific parallelization with
 * easy access.
 */
class ChemistryModule {
public:
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
  ChemistryModule(uint32_t wp_size,
                  const InitialList::ChemistryInit chem_params,
                  MPI_Comm communicator);

  /**
   * Deconstructor, which frees DHT data structure if used.
   */
  ~ChemistryModule();

  void masterSetField(Field<TugType> &&field);
  /**
   * Run the chemical simulation with parameters set.
   */
  void simulate(double dt);

  /**
   * Returns all known species names, including not only aqueous species, but
   * also equilibrium, exchange, surface and kinetic reactants.
   */
  // auto GetPropNames() const { return this->prop_names; }

  /**
   * Return the accumulated runtime in seconds for chemical simulation.
   */
  auto GetChemistryTime() const { return this->chem_t; }

  void setFilePadding(std::uint32_t maxiter) {
    this->file_pad =
        static_cast<std::uint8_t>(std::ceil(std::log10(maxiter + 1)));
  }

  struct SurrogateSetup {
    std::vector<std::string> prop_names;

    bool dht_enabled;
    std::uint32_t dht_size_mb;

    bool interp_enabled;
    std::uint32_t interp_bucket_size;
    std::uint32_t interp_size_mb;
    std::uint32_t interp_min_entries;
  };

  void masterEnableSurrogates(const SurrogateSetup &setup) {
    // FIXME: This is a hack to get the prop_names and prop_count from the setup
    this->prop_names = setup.prop_names;
    this->prop_count = setup.prop_names.size();

    this->dht_enabled = setup.dht_enabled;
    this->interp_enabled = setup.interp_enabled;

    if (this->dht_enabled || this->interp_enabled) {
      this->initializeDHT(setup.dht_size_mb, this->params.dht_species);
    }

    if (this->interp_enabled) {
      this->initializeInterp(setup.interp_bucket_size, setup.interp_size_mb,
                             setup.interp_min_entries,
                             this->params.interp_species);
    }
  }

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
    DHT_SNAPS_DISABLED = 0, //!< disabled file output
    DHT_SNAPS_SIMEND,       //!< only output of snapshot after simulation
    DHT_SNAPS_ITEREND       //!< output snapshots after each iteration
  };

  /**
   * **Only called by workers!** Start the worker listening loop.
   */
  void WorkerLoop();

  /**
   * **Called by master** Advise the workers to break the loop.
   */
  void MasterLoopBreak();

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
  Field<TugType> &getField() { return *this->chem_field; }

  /**
   * **Master only** Enable/disable progress bar.
   *
   * \param enabled True if print progressbar, false if not.
   */
  void setProgressBarPrintout(bool enabled) {
    this->print_progessbar = enabled;
  };

  std::vector<uint32_t> GetWorkerInterpolationCalls() const;

  std::vector<double> GetWorkerInterpolationWriteTimings() const;
  std::vector<double> GetWorkerInterpolationReadTimings() const;
  std::vector<double> GetWorkerInterpolationGatherTimings() const;
  std::vector<double> GetWorkerInterpolationFunctionCallTimings() const;

  std::vector<uint32_t> GetWorkerPHTCacheHits() const;

protected:
  void initializeDHT(uint32_t size_mb,
                     const NamedVector<std::uint32_t> &key_species);
  void setDHTSnapshots(int type, const std::string &out_dir);
  void setDHTReadFile(const std::string &input_file);

  void initializeInterp(std::uint32_t bucket_size, std::uint32_t size_mb,
                        std::uint32_t min_entries,
                        const NamedVector<std::uint32_t> &key_species);

  enum {
    CHEM_FIELD_INIT,
    CHEM_DHT_ENABLE,
    CHEM_DHT_SIGNIF_VEC,
    CHEM_DHT_SNAPS,
    CHEM_DHT_READ_FILE,
    CHEM_IP_ENABLE,
    CHEM_IP_MIN_ENTRIES,
    CHEM_IP_SIGNIF_VEC,
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
    WORKER_IP_WRITE,
    WORKER_IP_READ,
    WORKER_IP_GATHER,
    WORKER_IP_FC,
    WORKER_DHT_HITS,
    WORKER_DHT_EVICTIONS,
    WORKER_PHT_CACHE_HITS,
    WORKER_IP_CALLS
  };

  std::vector<uint32_t> interp_calls;
  std::vector<uint32_t> dht_hits;
  std::vector<uint32_t> dht_evictions;

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

  void MasterRunParallel(double dt);
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

  void WorkerRunWorkPackage(WorkPackage &work_package, double dSimTime,
                            double dTimestep);

  std::vector<uint32_t> CalculateWPSizesVector(uint32_t n_cells,
                                               uint32_t wp_size) const;

  std::vector<std::int32_t>
  parseDHTSpeciesVec(const NamedVector<std::uint32_t> &key_species,
                     const std::vector<std::string> &to_compare) const;

  void BCastStringVec(std::vector<std::string> &io);

  int comm_size, comm_rank;
  MPI_Comm group_comm;

  bool is_sequential;
  bool is_master;

  uint32_t wp_size;
  bool dht_enabled{false};
  int dht_snaps_type{DHT_SNAPS_DISABLED};
  std::string dht_file_out_dir;

  poet::DHT_Wrapper *dht = nullptr;

  bool interp_enabled{false};
  std::unique_ptr<poet::InterpolationModule> interp;

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

  bool print_progessbar{false};

  std::uint8_t file_pad{1};

  double chem_t{0.};

  uint32_t n_cells = 0;
  uint32_t prop_count = 0;
  std::vector<std::string> prop_names;

  std::unique_ptr<Field<TugType>> chem_field;

  const InitialList::ChemistryInit params;

  std::map<int, std::unique_ptr<PhreeqcEngine>> phreeqc_instances;
};
} // namespace poet

#endif // CHEMISTRYMODULE_H_
