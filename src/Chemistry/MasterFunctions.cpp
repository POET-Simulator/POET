#include "ChemistryModule.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <mpi.h>
#include <vector>

std::vector<uint32_t>
poet::ChemistryModule::MasterGatherWorkerMetrics(int type) const
{
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);

  uint32_t dummy;
  std::vector<uint32_t> metrics(this->comm_size);

  MPI_Gather(&dummy, 1, MPI_UINT32_T, metrics.data(), 1, MPI_UINT32_T, 0,
             this->group_comm);

  metrics.erase(metrics.begin());
  return metrics;
}

std::vector<double>
poet::ChemistryModule::MasterGatherWorkerTimings(int type) const
{
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);

  double dummy;
  std::vector<double> timings(this->comm_size);

  MPI_Gather(&dummy, 1, MPI_DOUBLE, timings.data(), 1, MPI_DOUBLE, 0,
             this->group_comm);

  timings.erase(timings.begin());
  return timings;
}

std::vector<double> poet::ChemistryModule::GetWorkerPhreeqcTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_PHREEQC);
}

std::vector<double> poet::ChemistryModule::GetWorkerDHTGetTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_DHT_GET);
}

std::vector<double> poet::ChemistryModule::GetWorkerDHTFillTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_DHT_FILL);
}

std::vector<double> poet::ChemistryModule::GetWorkerIdleTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_IDLE);
}

std::vector<uint32_t> poet::ChemistryModule::GetWorkerDHTHits() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  type = WORKER_DHT_HITS;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);

  MPI_Status probe;
  MPI_Probe(MPI_ANY_SOURCE, WORKER_DHT_HITS, this->group_comm, &probe);
  int count;
  MPI_Get_count(&probe, MPI_UINT32_T, &count);

  std::vector<uint32_t> ret(count);
  MPI_Recv(ret.data(), count, MPI_UINT32_T, probe.MPI_SOURCE, WORKER_DHT_HITS,
           this->group_comm, NULL);

  return ret;
}

std::vector<uint32_t> poet::ChemistryModule::GetWorkerDHTEvictions() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  type = WORKER_DHT_EVICTIONS;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);

  MPI_Status probe;
  MPI_Probe(MPI_ANY_SOURCE, WORKER_DHT_EVICTIONS, this->group_comm, &probe);
  int count;
  MPI_Get_count(&probe, MPI_UINT32_T, &count);

  std::vector<uint32_t> ret(count);
  MPI_Recv(ret.data(), count, MPI_UINT32_T, probe.MPI_SOURCE,
           WORKER_DHT_EVICTIONS, this->group_comm, NULL);

  return ret;
}

std::vector<double>
poet::ChemistryModule::GetWorkerInterpolationWriteTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_IP_WRITE);
}

std::vector<double>
poet::ChemistryModule::GetWorkerInterpolationReadTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_IP_READ);
}

std::vector<double>
poet::ChemistryModule::GetWorkerInterpolationGatherTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_IP_GATHER);
}

std::vector<double>
poet::ChemistryModule::GetWorkerInterpolationFunctionCallTimings() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  return MasterGatherWorkerTimings(WORKER_IP_FC);
}

std::vector<uint32_t>
poet::ChemistryModule::GetWorkerInterpolationCalls() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  type = WORKER_IP_CALLS;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);

  MPI_Status probe;
  MPI_Probe(MPI_ANY_SOURCE, WORKER_IP_CALLS, this->group_comm, &probe);
  int count;
  MPI_Get_count(&probe, MPI_UINT32_T, &count);

  std::vector<uint32_t> ret(count);
  MPI_Recv(ret.data(), count, MPI_UINT32_T, probe.MPI_SOURCE, WORKER_IP_CALLS,
           this->group_comm, NULL);

  return ret;
}

std::vector<uint32_t> poet::ChemistryModule::GetWorkerPHTCacheHits() const
{
  int type = CHEM_PERF;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
  type = WORKER_PHT_CACHE_HITS;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);

  MPI_Status probe;
  MPI_Probe(MPI_ANY_SOURCE, type, this->group_comm, &probe);
  int count;
  MPI_Get_count(&probe, MPI_UINT32_T, &count);

  std::vector<uint32_t> ret(count);
  MPI_Recv(ret.data(), count, MPI_UINT32_T, probe.MPI_SOURCE, type,
           this->group_comm, NULL);

  return ret;
}
void poet::ChemistryModule::computeStats(const std::vector<double> &pqc_vector,
                                         const std::vector<double> &sur_vector,
                                         uint32_t size_per_prop, uint32_t species_count,
                                         error_stats &stats)
{
  for (uint32_t i = 0; i < species_count; i++)
  {
    double err_sum = 0.0;
    double sqr_err_sum = 0.0;

    for (uint32_t j = 0; j < size_per_prop; j++)
    {

      const double &pqc_value = pqc_vector[i * size_per_prop + j];
      const double &sur_value = sur_vector[i * size_per_prop + j];

      if (pqc_value == 0 && sur_value == 0) {
        //
      } else if (pqc_value == 0 && sur_value != 0) {
        std::cout << "NOOOO! pqc = " << pqc_value << ", sur = " << sur_value << "\n";
        err_sum += 1.;
        sqr_err_sum += 1.;
      }  else {
        const double alpha = 1 - (sur_value/pqc_value);
        err_sum += std::abs(alpha);
        sqr_err_sum += alpha * alpha;
      }

      // if (pqc_value != 0)
      // {
      //   double rel_err = (pqc_value - sur_value) / pqc_value;
      //   err_sum += std::abs(rel_err);
      //   sqr_err_sum += rel_err * rel_err;
      // }
      // if (pqc_value == 0 && sur_value != 0)
      // {
      //   err_sum += 1.0;
      //   sqr_err_sum += 1.0;
      // }
      // else: both cases are zero, skip (no error)
      if (i == 6 && (j % 1000 == 0)) {
        std::cout << "pqc = " << pqc_value << ", sur = " << sur_value << "\n";
      }
    }
    if (i == 0)
    {
      std::cout << "computeStats, i==0, err_sum: " << err_sum << std::endl;
      std::cout << "computeStats, i==0, sqr_err_sum: " << sqr_err_sum << std::endl;
    }
    stats.mape[i] = 100.0  * (err_sum / size_per_prop);
    stats.rrsme[i] = (size_per_prop > 0) ? std::sqrt(sqr_err_sum / size_per_prop) : 0.0;
  }

}

inline std::vector<int> shuffleVector(const std::vector<int> &in_vector,
                                      uint32_t size_per_prop,
                                      uint32_t wp_count)
{
  std::vector<int> out_buffer(in_vector.size());
  uint32_t write_i = 0;
  for (uint32_t i = 0; i < wp_count; i++)
  {
    for (uint32_t j = i; j < size_per_prop; j += wp_count)
    {
      out_buffer[write_i] = in_vector[j];
      write_i++;
    }
  }
  return out_buffer;
}

inline std::vector<double> shuffleField(const std::vector<double> &in_field,
                                        uint32_t size_per_prop,
                                        uint32_t species_count,
                                        uint32_t wp_count)
{
  std::vector<double> out_buffer(in_field.size());
  uint32_t write_i = 0;
  for (uint32_t i = 0; i < wp_count; i++)
  {
    for (uint32_t j = i; j < size_per_prop; j += wp_count)
    {
      for (uint32_t k = 0; k < species_count; k++)
      {
        out_buffer[(write_i * species_count) + k] =
            in_field[(k * size_per_prop) + j];
      }
      write_i++;
    }
  }
  return out_buffer;
}

inline void unshuffleField(const std::vector<double> &in_buffer,
                           uint32_t size_per_prop, uint32_t species_count,
                           uint32_t wp_count, std::vector<double> &out_field)
{
  uint32_t read_i = 0;

  for (uint32_t i = 0; i < wp_count; i++)
  {
    for (uint32_t j = i; j < size_per_prop; j += wp_count)
    {
      for (uint32_t k = 0; k < species_count; k++)
      {
        out_field[(k * size_per_prop) + j] =
            in_buffer[(read_i * species_count) + k];
      }
      read_i++;
    }
  }
}

inline void printProgressbar(int count_pkgs, int n_wp, int barWidth = 70)
{
  /* visual progress */
  double progress = (float)(count_pkgs + 1) / n_wp;

  std::cout << "[";
  int pos = barWidth * progress;
  for (int iprog = 0; iprog < barWidth; ++iprog)
  {
    if (iprog < pos)
      std::cout << "=";
    else if (iprog == pos)
      std::cout << ">";
    else
      std::cout << " ";
  }
  std::cout << "] " << int(progress * 100.0) << " %\r";
  std::cout.flush();
  /* end visual progress */
}

inline void poet::ChemistryModule::MasterSendPkgs(
    worker_list_t &w_list, workpointer_t &work_pointer, int &pkg_to_send,
    int &count_pkgs, int &free_workers, double dt, uint32_t iteration, uint32_t control_iteration,
    const std::vector<uint32_t> &wp_sizes_vector)
{
  /* declare variables */
  int local_work_package_size;

  /* search for free workers and send work */
  for (int p = 0; p < this->comm_size - 1; p++)
  {
    if (w_list[p].has_work == 0 && pkg_to_send > 0) /* worker is free */
    {
      /* to enable different work_package_size, set local copy of
       * work_package_size to pre-calculated work package size vector */

      local_work_package_size = (int)wp_sizes_vector[count_pkgs];
      count_pkgs++;

      /* note current processed work package in workerlist */
      w_list[p].send_addr = work_pointer.base();

      /* push work pointer to next work package */
      const uint32_t end_of_wp = local_work_package_size * this->prop_count;
      std::vector<double> send_buffer(end_of_wp + this->BUFFER_OFFSET);
      std::copy(work_pointer, work_pointer + end_of_wp, send_buffer.begin());

      work_pointer += end_of_wp;

      // fill send buffer starting with work_package ...
      // followed by: work_package_size
      send_buffer[end_of_wp] = (double)local_work_package_size;
      // current iteration of simulation
      send_buffer[end_of_wp + 1] = (double)iteration;
      // size of timestep in seconds
      send_buffer[end_of_wp + 2] = dt;
      // current time of simulation (age) in seconds
      send_buffer[end_of_wp + 3] = this->simtime;
      // current work package start location in field
      uint32_t wp_start_index = std::accumulate(wp_sizes_vector.begin(), std::next(wp_sizes_vector.begin(), count_pkgs), 0);
      send_buffer[end_of_wp + 4] = wp_start_index;
      // whether this iteration is a control iteration
      send_buffer[end_of_wp + 5] = control_iteration;

      /* ATTENTION Worker p has rank p+1 */
      // MPI_Send(send_buffer, end_of_wp + BUFFER_OFFSET, MPI_DOUBLE, p + 1,
      //          LOOP_WORK, this->group_comm);
      MPI_Send(send_buffer.data(), send_buffer.size(), MPI_DOUBLE, p + 1,
               LOOP_WORK, this->group_comm);

      /* Mark that worker has work to do */
      w_list[p].has_work = 1;
      free_workers--;
      pkg_to_send -= 1;
    }
  }
}

inline void poet::ChemistryModule::MasterRecvPkgs(worker_list_t &w_list,
                                                  int &pkg_to_recv,
                                                  bool to_send,
                                                  int &free_workers)
{
  /* declare most of the variables here */
  int need_to_receive = 1;
  double idle_a, idle_b;
  int p, size;

  MPI_Status probe_status;
  // master_recv_a = MPI_Wtime();
  /* start to loop as long there are packages to recv and the need to receive
   */
  while (need_to_receive && pkg_to_recv > 0)
  {
    // only of there are still packages to send and free workers are available
    if (to_send && free_workers > 0)
      // non blocking probing
      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &need_to_receive,
                 &probe_status);
    else
    {
      idle_a = MPI_Wtime();
      // blocking probing
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &probe_status);
      idle_b = MPI_Wtime();
      this->idle_t += idle_b - idle_a;
    }

    /* if need_to_receive was set to true above, so there is a message to
     * receive */
    if (need_to_receive)
    {
      p = probe_status.MPI_SOURCE;
      if (probe_status.MPI_TAG == LOOP_WORK)
      {
        MPI_Get_count(&probe_status, MPI_DOUBLE, &size);
        MPI_Recv(w_list[p - 1].send_addr, size, MPI_DOUBLE, p, LOOP_WORK,
                 this->group_comm, MPI_STATUS_IGNORE);
        w_list[p - 1].has_work = 0;
        pkg_to_recv -= 1;
        free_workers++;
      }
      if (probe_status.MPI_TAG == LOOP_CTRL)
      {
        MPI_Get_count(&probe_status, MPI_DOUBLE, &size);

        // layout of buffer is [phreeqc][surrogate]
        std::vector<double> recv_buffer(size);

        MPI_Recv(recv_buffer.data(), size, MPI_DOUBLE, p, LOOP_CTRL,
                 this->group_comm, MPI_STATUS_IGNORE);

        std::copy(recv_buffer.begin(), recv_buffer.begin() + (size / 2),
                  w_list[p - 1].send_addr);

        sur_shuffled.insert(sur_shuffled.end(), recv_buffer.begin() + (size / 2),
                            recv_buffer.begin() + size);

        w_list[p - 1].has_work = 0;
        pkg_to_recv -= 1;
        free_workers++;
      }
    }
  }
}

void poet::ChemistryModule::simulate(double dt)
{
  double start_t{MPI_Wtime()};
  if (this->is_sequential)
  {
    MasterRunSequential();
    return;
  }

  MasterRunParallel(dt);
  double end_t{MPI_Wtime()};
  this->chem_t += end_t - start_t;
}

void poet::ChemistryModule::MasterRunSequential()
{
  // std::vector<double> shuffled_field =
  //     shuffleField(chem_field.AsVector(), n_cells, prop_count, 1);

  // std::vector<std::vector<double>> input;
  // for (std::size_t i = 0; i < n_cells; i++) {
  //   input.push_back(
  //       std::vector<double>(shuffled_field.begin() + (i * prop_count),
  //                           shuffled_field.begin() + ((i + 1) *
  //                           prop_count)));
  // }

  // this->setDumpedField(input);
  // PhreeqcRM::RunCells();
  // this->getDumpedField(input);

  // shuffled_field.clear();
  // for (std::size_t i = 0; i < n_cells; i++) {
  //   shuffled_field.insert(shuffled_field.end(), input[i].begin(),
  //                         input[i].end());
  // }

  // std::vector<double> out_vec{shuffled_field};
  // unshuffleField(shuffled_field, n_cells, prop_count, 1, out_vec);
  // chem_field = out_vec;
}

void poet::ChemistryModule::MasterRunParallel(double dt)
{
  /* declare most of the needed variables here */
  double seq_a, seq_b, seq_c, seq_d;
  double worker_chemistry_a, worker_chemistry_b;
  double sim_e_chemistry, sim_f_chemistry;
  int pkg_to_send, pkg_to_recv;
  int free_workers;
  int i_pkgs;
  int ftype;

  const std::vector<uint32_t> wp_sizes_vector =
      CalculateWPSizesVector(this->n_cells, this->wp_size);

  if (this->ai_surrogate_enabled)
  {
    ftype = CHEM_AI_BCAST_VALIDITY;
    PropagateFunctionType(ftype);
    this->ai_surrogate_validity_vector = shuffleVector(this->ai_surrogate_validity_vector,
                                                       this->n_cells,
                                                       wp_sizes_vector.size());
    ChemBCast(&this->ai_surrogate_validity_vector.front(), this->n_cells, MPI_INT);
  }

  ftype = CHEM_WORK_LOOP;
  PropagateFunctionType(ftype);

  MPI_Barrier(this->group_comm);

  static uint32_t iteration = 0;
  uint32_t control_iteration = static_cast<uint32_t>(this->runtime_params->control_iteration_active ? 1 : 0);
  if (control_iteration)
  {
    sur_shuffled.clear();
    sur_shuffled.reserve(this->n_cells * this->prop_count);
  }

  /* start time measurement of sequential part */
  seq_a = MPI_Wtime();

  /* shuffle grid */
  // grid.shuffleAndExport(mpi_buffer);

  std::vector<double> mpi_buffer =
      shuffleField(chem_field.AsVector(), this->n_cells, this->prop_count,
                   wp_sizes_vector.size());

  /* setup local variables */
  pkg_to_send = wp_sizes_vector.size();
  pkg_to_recv = wp_sizes_vector.size();

  workpointer_t work_pointer = mpi_buffer.begin();
  worker_list_t worker_list(this->comm_size - 1);

  free_workers = this->comm_size - 1;
  i_pkgs = 0;

  /* end time measurement of sequential part */
  seq_b = MPI_Wtime();
  seq_t += seq_b - seq_a;

  /* start time measurement of chemistry time needed for send/recv loop */
  worker_chemistry_a = MPI_Wtime();

  /* start send/recv loop */
  // while there are still packages to recv
  while (pkg_to_recv > 0)
  {
    // print a progressbar to stdout
    if (print_progessbar)
    {
      printProgressbar((int)i_pkgs, (int)wp_sizes_vector.size());
    }
    // while there are still packages to send
    if (pkg_to_send > 0)
    {
      // send packages to all free workers ...
      MasterSendPkgs(worker_list, work_pointer, pkg_to_send, i_pkgs,
                     free_workers, dt, iteration, control_iteration, wp_sizes_vector);
    }
    // ... and try to receive them from workers who has finished their work
    MasterRecvPkgs(worker_list, pkg_to_recv, pkg_to_send > 0, free_workers);
  }

  // Just to complete the progressbar
  std::cout << std::endl;

  /* stop time measurement of chemistry time needed for send/recv loop */
  worker_chemistry_b = MPI_Wtime();
  this->send_recv_t += worker_chemistry_b - worker_chemistry_a;

  /* start time measurement of sequential part */
  seq_c = MPI_Wtime();

  /* unshuffle grid */
  // grid.importAndUnshuffle(mpi_buffer);
  std::vector<double> out_vec{mpi_buffer};
  unshuffleField(mpi_buffer, this->n_cells, this->prop_count,
                 wp_sizes_vector.size(), out_vec);
  chem_field = out_vec;

  /* do master stuff */

  if (control_iteration)
  {
    control_iteration_counter++;

    std::vector<double> sur_unshuffled{sur_shuffled};

    unshuffleField(sur_shuffled, this->n_cells, this->prop_count,
                   wp_sizes_vector.size(), sur_unshuffled);

    error_stats stats(this->prop_count, control_iteration_counter * runtime_params->control_iteration);

    computeStats(out_vec, sur_unshuffled, this->n_cells, this->prop_count, stats);
    error_stats_history.push_back(stats);

    // to do: control values to epsilon
  }

  /* start time measurement of master chemistry */
  sim_e_chemistry = MPI_Wtime();

  /* end time measurement of sequential part */
  seq_d = MPI_Wtime();
  this->seq_t += seq_d - seq_c;

  /* end time measurement of whole chemistry simulation */

  /* advise workers to end chemistry iteration */
  for (int i = 1; i < this->comm_size; i++)
  {
    MPI_Send(NULL, 0, MPI_DOUBLE, i, LOOP_END, this->group_comm);
  }

  this->simtime += dt;
  iteration++;
}

void poet::ChemistryModule::MasterLoopBreak()
{
  int type = CHEM_BREAK_MAIN_LOOP;
  MPI_Bcast(&type, 1, MPI_INT, 0, this->group_comm);
}

std::vector<uint32_t>
poet::ChemistryModule::CalculateWPSizesVector(uint32_t n_cells,
                                              uint32_t wp_size) const
{
  bool mod_pkgs = (n_cells % wp_size) != 0;
  uint32_t n_packages =
      (uint32_t)(n_cells / wp_size) + static_cast<int>(mod_pkgs);

  std::vector<uint32_t> wp_sizes_vector(n_packages, 0);

  for (int i = 0; i < n_cells; i++)
  {
    wp_sizes_vector[i % n_packages] += 1;
  }

  return wp_sizes_vector;
}

void poet::ChemistryModule::masterSetField(Field field)
{
  this->chem_field = field;
  this->prop_count = field.GetProps().size();

  int ftype = CHEM_FIELD_INIT;
  PropagateFunctionType(ftype);

  ChemBCast(&this->prop_count, 1, MPI_UINT32_T);
}
