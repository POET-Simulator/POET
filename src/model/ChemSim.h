#ifndef CHEMSIM_H
#define CHEMSIM_H

#include <DHT_Wrapper.h>
#include <RRuntime.h>
#include <SimParams.h>
#include <mpi.h>

#include <vector>

#include "Grid.h"

#define BUFFER_OFFSET 5

#define TAG_WORK 42
#define TAG_FINISH 43
#define TAG_TIMING 44
#define TAG_DHT_PERF 45
#define TAG_DHT_STATS 46
#define TAG_DHT_STORE 47
#define TAG_DHT_ITER 48

namespace poet {
class ChemSim {
 public:
  ChemSim(t_simparams *params, RRuntime &R_, Grid &grid_);

  virtual void run();
  virtual void end();

  double getChemistryTime();

 protected:
  double current_sim_time = 0;
  int iteration = 0;
  int dt = 0;

  int world_rank;
  int world_size;
  unsigned int wp_size;

  RRuntime &R;
  Grid &grid;

  std::vector<int> wp_sizes_vector;
  std::string out_dir;

  double *send_buffer;
  typedef struct {
    char has_work;
    double *send_addr;
  } worker_struct;
  worker_struct *workerlist;

  double *mpi_buffer;

  double chem_t = 0.f;
};

class ChemMaster : public ChemSim {
 public:
  ChemMaster(t_simparams *params, RRuntime &R_, Grid &grid_);
  ~ChemMaster();

  void run() override;
  void end() override;

  double getSendTime();
  double getRecvTime();
  double getIdleTime();
  double getWorkerTime();
  double getChemMasterTime();
  double getSeqTime();

 private:
  void printProgressbar(int count_pkgs, int n_wp, int barWidth = 70);
  void sendPkgs(int &pkg_to_send, int &count_pkgs, int &free_workers);
  void recvPkgs(int &pkg_to_recv, bool to_send, int &free_workers);

  bool dht_enabled;
  unsigned int wp_size;
  double *work_pointer;

  double send_t = 0.f;
  double recv_t = 0.f;
  double master_idle = 0.f;
  double worker_t = 0.f;
  double chem_master = 0.f;
  double seq_t = 0.f;
};

class ChemWorker : public ChemSim {
 public:
  ChemWorker(t_simparams *params_, RRuntime &R_, Grid &grid_,
             MPI_Comm dht_comm);
  ~ChemWorker();

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
  t_simparams *params;

  double *mpi_buffer_results;

  DHT_Wrapper *dht;

  double timing[3];
  double idle_t = 0.f;
  int phreeqc_count = 0;
};
}  // namespace poet
#endif  // CHEMSIM_H
