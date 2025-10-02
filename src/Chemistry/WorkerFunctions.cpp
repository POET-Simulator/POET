#include "ChemistryModule.hpp"
#include "SurrogateModels/DHT_Wrapper.hpp"
#include "SurrogateModels/Interpolation.hpp"

#include "Chemistry/ChemistryDefs.hpp"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <mpi.h>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace poet
{

  inline std::string get_string(int root, MPI_Comm communicator)
  {
    int count;
    MPI_Bcast(&count, 1, MPI_INT, root, communicator);

    char *buffer = new char[count + 1];
    MPI_Bcast(buffer, count, MPI_CHAR, root, communicator);

    buffer[count] = '\0';

    std::string ret_str(buffer);
    delete[] buffer;

    return ret_str;
  }

    void poet::ChemistryModule::WorkerLoop()
    {
      struct worker_s timings;

      // HACK: defining the worker iteration count here, which will increment after
      // each CHEM_ITER_END message
      uint32_t iteration = 1;
      bool loop = true;

      while (loop)
      {
        int func_type;
        PropagateFunctionType(func_type);

        switch (func_type)
        {
        case CHEM_FIELD_INIT:
        {
          ChemBCast(&this->prop_count, 1, MPI_UINT32_T);
          if (this->ai_surrogate_enabled)
          {
            this->ai_surrogate_validity_vector.resize(
                this->n_cells); // resize statt reserve?
          }
          break;
        }
        case CHEM_AI_BCAST_VALIDITY:
        {
          // Receive the index vector of valid ai surrogate predictions
          MPI_Bcast(&this->ai_surrogate_validity_vector.front(), this->n_cells,
                    MPI_INT, 0, this->group_comm);
          break;
        }
        case CHEM_INTERP:
        {
          int interp_flag;
          ChemBCast(&interp_flag, 1, MPI_INT);
          this->interp_enabled = (interp_flag == 1);
          break;
        }
        case CHEM_WORK_LOOP:
        {
          WorkerProcessPkgs(timings, iteration);
          break;
        }
        case CHEM_PERF:
        {
          int type;
          ChemBCast(&type, 1, MPI_INT);
          if (type < WORKER_DHT_HITS)
          {
            WorkerPerfToMaster(type, timings);
            break;
          }
          WorkerMetricsToMaster(type);
          break;
        }
        case CHEM_BREAK_MAIN_LOOP:
        {
          WorkerPostSim(iteration);
          loop = false;
          break;
        }
        default:
        {
          throw std::runtime_error("Worker received unknown tag from master.");
        }
        }
      }
    }

    void poet::ChemistryModule::WorkerProcessPkgs(struct worker_s &timings,
                                                  uint32_t &iteration)
    {
      MPI_Status probe_status;
      bool loop = true;

      MPI_Barrier(this->group_comm);

      while (loop)
      {
        double idle_a = MPI_Wtime();
        MPI_Probe(0, MPI_ANY_TAG, this->group_comm, &probe_status);
        double idle_b = MPI_Wtime();

        switch (probe_status.MPI_TAG)
        {
        case LOOP_WORK:
        {
          timings.idle_t += idle_b - idle_a;
          int count;
          MPI_Get_count(&probe_status, MPI_DOUBLE, &count);

          WorkerDoWork(probe_status, count, timings);
          break;
        }
        case LOOP_END:
        {
          WorkerPostIter(probe_status, iteration);
          iteration++;
          loop = false;
          break;
        }
        }
      }
    }

  void poet::ChemistryModule::WorkerDoWork(MPI_Status &probe_status,
                                           int double_count,
                                           struct worker_s &timings)
  {
    static int counter = 1;

    double dht_get_start, dht_get_end;
    double phreeqc_time_start, phreeqc_time_end;
    double dht_fill_start, dht_fill_end;

    uint32_t iteration;
    double dt;
    double current_sim_time;
    uint32_t wp_start_index;
    int count = double_count;
    bool control_iteration_active = false;
    std::vector<double> mpi_buffer(count);

    /* receive */
    MPI_Recv(mpi_buffer.data(), count, MPI_DOUBLE, 0, LOOP_WORK, this->group_comm,
             MPI_STATUS_IGNORE);

    /* decrement count of work_package by BUFFER_OFFSET */
    count -= BUFFER_OFFSET;
    /* check for changes on all additional variables given by the 'header' of
     * mpi_buffer */

    // work_package_size
    poet::WorkPackage s_curr_wp(mpi_buffer[count]);

    // current iteration of simulation
    iteration = mpi_buffer[count + 1];

    // current timestep size
    dt = mpi_buffer[count + 2];

    // current simulation time ('age' of simulation)
    current_sim_time = mpi_buffer[count + 3];

    // current work package start location in field
    wp_start_index = mpi_buffer[count + 4];

    control_iteration_active = (mpi_buffer[count + 5] == 1);

    for (std::size_t wp_i = 0; wp_i < s_curr_wp.size; wp_i++)
    {
      s_curr_wp.input[wp_i] =
          std::vector<double>(mpi_buffer.begin() + this->prop_count * wp_i,
                              mpi_buffer.begin() + this->prop_count * (wp_i + 1));
    }

    // std::cout << this->comm_rank << ":" << counter++ << std::endl;
    if (dht_enabled || interp_enabled)
    {
      dht->prepareKeys(s_curr_wp.input, dt);
    }

    if (dht_enabled)
    {
      /* check for values in DHT */
      dht_get_start = MPI_Wtime();
      dht->checkDHT(s_curr_wp);
      dht_get_end = MPI_Wtime();
      timings.dht_get += dht_get_end - dht_get_start;
    }

    if (interp_enabled)
    {
      interp->tryInterpolation(s_curr_wp);
    }

    if (this->ai_surrogate_enabled)
    {
      // Map valid predictions from the ai surrogate in the workpackage
      for (int i = 0; i < s_curr_wp.size; i++)
      {
        if (this->ai_surrogate_validity_vector[wp_start_index + i] == 1)
        {
          s_curr_wp.mapping[i] = CHEM_AISURR;
        }
      }
    }

    /* if control iteration: create copy surrogate results (output and mappings) and then set them to zero,
      give this to phreeqc */

    poet::WorkPackage s_curr_wp_control = s_curr_wp;

    if (control_iteration_active)
    {
      for (std::size_t wp_i = 0; wp_i < s_curr_wp_control.size; wp_i++)
      {
        s_curr_wp_control.output[wp_i] = std::vector<double>(this->prop_count, 0.0);
        s_curr_wp_control.mapping[wp_i] = 0;
      }
    }

    phreeqc_time_start = MPI_Wtime();

    WorkerRunWorkPackage(control_iteration_active ? s_curr_wp_control : s_curr_wp, current_sim_time, dt);

    phreeqc_time_end = MPI_Wtime();

    if (control_iteration_active)
    {     
      std::size_t sur_wp_offset = s_curr_wp.size * this->prop_count;
 
      mpi_buffer.resize(count + sur_wp_offset);
      
      for (std::size_t wp_i = 0; wp_i < s_curr_wp_control.size; wp_i++)
      {
        std::copy(s_curr_wp_control.output[wp_i].begin(), s_curr_wp_control.output[wp_i].end(),
                  mpi_buffer.begin() + this->prop_count * wp_i);
      }

      // s_curr_wp only contains the interpolated data
      // copy surrogate output after the the pqc output, mpi_buffer[pqc][interp]

      for (std::size_t wp_i = 0; wp_i < s_curr_wp.size; wp_i++)
      {
      if (s_curr_wp.mapping[wp_i] != CHEM_PQC) // only copy if surrogate was used
      {
        std::copy(s_curr_wp.output[wp_i].begin(), s_curr_wp.output[wp_i].end(),
                  mpi_buffer.begin() + sur_wp_offset + this->prop_count * wp_i);
      } else 
      {
        // if pqc was used, copy pqc results again
        std::copy(s_curr_wp_control.output[wp_i].begin(), s_curr_wp_control.output[wp_i].end(),
                  mpi_buffer.begin() + sur_wp_offset + this->prop_count * wp_i);
      }

      }

      count += sur_wp_offset;
    }
    else
    {
      for (std::size_t wp_i = 0; wp_i < s_curr_wp.size; wp_i++)
      {
        std::copy(s_curr_wp.output[wp_i].begin(), s_curr_wp.output[wp_i].end(),
                  mpi_buffer.begin() + this->prop_count * wp_i);
      }
    }

    /* send results to master */
    MPI_Request send_req;

    int mpi_tag = control_iteration_active ? LOOP_CTRL : LOOP_WORK;
    MPI_Isend(mpi_buffer.data(), count, MPI_DOUBLE, 0, mpi_tag, MPI_COMM_WORLD, &send_req);

    if (dht_enabled || interp_enabled)
    {
      /* write results to DHT */
      dht_fill_start = MPI_Wtime();
      dht->fillDHT(control_iteration_active ? s_curr_wp_control : s_curr_wp);
      dht_fill_end = MPI_Wtime();

      if (interp_enabled)
      {
        interp->writePairs();
      }
      timings.dht_fill += dht_fill_end - dht_fill_start;
    }

    timings.phreeqc_t += phreeqc_time_end - phreeqc_time_start;

    MPI_Wait(&send_req, MPI_STATUS_IGNORE);
  }

  void poet::ChemistryModule::WorkerPostIter(MPI_Status &prope_status,
                                             uint32_t iteration)
  {
    MPI_Recv(NULL, 0, MPI_DOUBLE, 0, LOOP_END, this->group_comm,
             MPI_STATUS_IGNORE);

    if (this->dht_enabled)
    {
      dht_hits.push_back(dht->getHits());
      dht_evictions.push_back(dht->getEvictions());
      dht->resetCounter();

      if (this->dht_snaps_type == DHT_SNAPS_ITEREND)
      {
        WorkerWriteDHTDump(iteration);
      }
    }

    if (this->interp_enabled)
    {
      std::stringstream out;
      interp_calls.push_back(interp->getInterpolationCount());
      interp->resetCounter();
      interp->writePHTStats();
      if (this->dht_snaps_type == DHT_SNAPS_ITEREND)
      {
        out << this->dht_file_out_dir << "/iter_" << std::setfill('0')
            << std::setw(this->file_pad) << iteration << ".pht";
        interp->dumpPHTState(out.str());
      }

      const auto max_mean_idx =
          DHT_get_used_idx_factor(this->interp->getDHTObject(), 1);

      if (max_mean_idx >= 2)
      {
        DHT_flush(this->interp->getDHTObject());
        DHT_flush(this->dht->getDHT());
        if (this->comm_rank == 2)
        {
          std::cout << "Flushed both DHT and PHT!\n\n";
        }
      }
    }

    RInsidePOET::getInstance().parseEvalQ("gc()");
  }

  void poet::ChemistryModule::WorkerPostSim(uint32_t iteration)
  {
    if (this->dht_enabled && this->dht_snaps_type >= DHT_SNAPS_ITEREND)
    {
      WorkerWriteDHTDump(iteration);
    }
    if (this->interp_enabled && this->dht_snaps_type >= DHT_SNAPS_ITEREND)
    {
      std::stringstream out;
      out << this->dht_file_out_dir << "/iter_" << std::setfill('0')
          << std::setw(this->file_pad) << iteration << ".pht";
      interp->dumpPHTState(out.str());
    }
  }

  void poet::ChemistryModule::WorkerWriteDHTDump(uint32_t iteration)
  {
    std::stringstream out;
    out << this->dht_file_out_dir << "/iter_" << std::setfill('0')
        << std::setw(this->file_pad) << iteration << ".dht";
    int res = dht->tableToFile(out.str().c_str());
    if (res != DHT_SUCCESS && this->comm_rank == 2)
      std::cerr
          << "CPP: Worker: Error in writing current state of DHT to file.\n";
    else if (this->comm_rank == 2)
      std::cout << "CPP: Worker: Successfully written DHT to file " << out.str()
                << "\n";
  }

  void poet::ChemistryModule::WorkerReadDHTDump(
      const std::string &dht_input_file)
  {
    int res = dht->fileToTable((char *)dht_input_file.c_str());
    if (res != DHT_SUCCESS)
    {
      if (res == DHT_WRONG_FILE)
      {
        if (this->comm_rank == 1)
          std::cerr
              << "CPP: Worker: Wrong file layout! Continue with empty DHT ...\n";
      }
      else
      {
        if (this->comm_rank == 1)
          std::cerr << "CPP: Worker: Error in loading current state of DHT from "
                       "file. Continue with empty DHT ...\n";
      }
    }
    else
    {
      if (this->comm_rank == 2)
        std::cout << "CPP: Worker: Successfully loaded state of DHT from file "
                  << dht_input_file << "\n";
    }
  }

  void poet::ChemistryModule::WorkerRunWorkPackage(WorkPackage &work_package,
                                                   double dSimTime,
                                                   double dTimestep)
  {

    std::vector<std::vector<double>> inout_chem = work_package.input;
    std::vector<std::size_t> to_ignore;

    for (std::size_t wp_id = 0; wp_id < work_package.size; wp_id++)
    {
      if (work_package.mapping[wp_id] != CHEM_PQC)
      {
        to_ignore.push_back(wp_id);
      }
    }
    this->pqc_runner->run(inout_chem, dTimestep, to_ignore);

    for (std::size_t wp_id = 0; wp_id < work_package.size; wp_id++)
    {
      if (work_package.mapping[wp_id] == CHEM_PQC)
      {
        work_package.output[wp_id] = inout_chem[wp_id];
      }
    }
  }

  void poet::ChemistryModule::WorkerPerfToMaster(int type,
                                                 const struct worker_s &timings)
  {
    switch (type)
    {
    case WORKER_PHREEQC:
    {
      MPI_Gather(&timings.phreeqc_t, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
                 this->group_comm);
      break;
    }
    case WORKER_DHT_GET:
    {
      MPI_Gather(&timings.dht_get, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
                 this->group_comm);
      break;
    }
    case WORKER_DHT_FILL:
    {
      MPI_Gather(&timings.dht_fill, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
                 this->group_comm);
      break;
    }
    case WORKER_IDLE:
    {
      MPI_Gather(&timings.idle_t, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0,
                 this->group_comm);
      break;
    }
    case WORKER_IP_WRITE:
    {
      double val = interp->getPHTWriteTime();
      MPI_Gather(&val, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0, this->group_comm);
      break;
    }
    case WORKER_IP_READ:
    {
      double val = interp->getPHTReadTime();
      MPI_Gather(&val, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0, this->group_comm);
      break;
    }
    case WORKER_IP_GATHER:
    {
      double val = interp->getDHTGatherTime();
      MPI_Gather(&val, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0, this->group_comm);
      break;
    }
    case WORKER_IP_FC:
    {
      double val = interp->getInterpolationTime();
      MPI_Gather(&val, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0, this->group_comm);
      break;
    }
    default:
    {
      throw std::runtime_error("Unknown perf type in master's message.");
    }
    }
  }

  void poet::ChemistryModule::WorkerMetricsToMaster(int type)
  {
    MPI_Comm worker_comm = dht->getCommunicator();
    int worker_rank;
    MPI_Comm_rank(worker_comm, &worker_rank);

    MPI_Comm &group_comm = this->group_comm;

    auto reduce_and_send = [&worker_rank, &worker_comm, &group_comm](
                               std::vector<std::uint32_t> &send_buffer, int tag)
    {
      std::vector<uint32_t> to_master(send_buffer.size());
      MPI_Reduce(send_buffer.data(), to_master.data(), send_buffer.size(),
                 MPI_UINT32_T, MPI_SUM, 0, worker_comm);

      if (worker_rank == 0)
      {
        MPI_Send(to_master.data(), to_master.size(), MPI_UINT32_T, 0, tag,
                 group_comm);
      }
    };

    switch (type)
    {
    case WORKER_DHT_HITS:
    {
      reduce_and_send(dht_hits, WORKER_DHT_HITS);
      break;
    }
    case WORKER_DHT_EVICTIONS:
    {
      reduce_and_send(dht_evictions, WORKER_DHT_EVICTIONS);
      break;
    }
    case WORKER_IP_CALLS:
    {
      reduce_and_send(interp_calls, WORKER_IP_CALLS);
      return;
    }
    case WORKER_PHT_CACHE_HITS:
    {
      std::vector<std::uint32_t> input = this->interp->getPHTLocalCacheHits();
      reduce_and_send(input, WORKER_PHT_CACHE_HITS);
      return;
    }
    default:
    {
      throw std::runtime_error("Unknown perf type in master's message.");
    }
    }
  }

} // namespace poet
