//  Time-stamp: "Last modified 2023-08-01 23:18:45 mluebke"

#include "poet/DHT_Wrapper.hpp"
#include "poet/HashFunctions.hpp"
#include "poet/Interpolation.hpp"
#include "poet/LookupKey.hpp"
#include "poet/Rounding.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <iterator>
#include <memory>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>
extern "C" {
#include "poet/DHT.h"
}

namespace poet {

InterpolationModule::InterpolationModule(
    std::uint32_t entries_per_bucket, std::uint64_t size_per_process,
    std::uint32_t min_entries_needed, DHT_Wrapper &dht,
    const NamedVector<std::uint32_t> &interp_key_signifs,
    const std::vector<std::int32_t> &dht_key_indices)
    : dht_instance(dht), key_signifs(interp_key_signifs),
      key_indices(dht_key_indices), min_entries_needed(min_entries_needed) {

  initPHT(this->key_signifs.size(), entries_per_bucket, size_per_process,
          dht.getCommunicator());

  pht->setSourceDHT(dht.getDHT());
}

void InterpolationModule::initPHT(std::uint32_t key_count,
                                  std::uint32_t entries_per_bucket,
                                  std::uint32_t size_per_process,
                                  MPI_Comm communicator) {
  uint32_t key_size = key_count * sizeof(Lookup_Keyelement) + sizeof(double);
  uint32_t data_size = sizeof(DHT_Location);

  pht = std::make_unique<ProximityHashTable>(
      key_size, data_size, entries_per_bucket, size_per_process, communicator);
}

void InterpolationModule::writePairs(const DHT_Wrapper::DHT_ResultObject &in) {
  for (int i = 0; i < in.length; i++) {
    if (in.needPhreeqc[i]) {
      const auto coarse_key = roundKey(in.keys[i]);
      pht->writeLocationToPHT(coarse_key, in.locations[i]);
    }
  }
}

void InterpolationModule::tryInterpolation(
    DHT_Wrapper::DHT_ResultObject &dht_results,
    std::vector<std::uint32_t> &curr_mapping) {
  interp_result.status.resize(dht_results.length, NOT_NEEDED);
  interp_result.results.resize(dht_results.length, {});

  for (int i = 0; i < dht_results.length; i++) {
    if (!dht_results.needPhreeqc[i]) {
      interp_result.status[i] = NOT_NEEDED;
      continue;
    }

    auto pht_result =
        pht->query(roundKey(dht_results.keys[i]), this->key_signifs.getValues(),
                   this->min_entries_needed, dht_instance.getInputCount(),
                   dht_instance.getOutputCount());

    int pht_i = 0;

    while (pht_i < pht_result.size) {
      if (pht_result.size < this->min_entries_needed) {
        break;
      }

      auto in_it = pht_result.in_values.begin() + pht_i;
      auto out_it = pht_result.out_values.begin() + pht_i;

      bool same_sig_calcite = (pht_result.in_values[pht_i][7] == 0) ==
                              (dht_results.results[i][7] == 0);
      bool same_sig_dolomite = (pht_result.in_values[pht_i][8] == 0) ==
                               (dht_results.results[i][9] == 0);
      if (!same_sig_calcite || !same_sig_dolomite) {
        pht_result.size -= 1;
        pht_result.in_values.erase(in_it);
        pht_result.out_values.erase(out_it);
        continue;
      }

      pht_i += 1;
    }

    if (pht_result.size < this->min_entries_needed) {
      interp_result.status[i] = INSUFFICIENT_DATA;
      continue;
    }

#ifdef POET_PHT_ADD
    this->pht->incrementReadCounter(roundKey(dht_results.keys[i]));
#endif

    double start_fc = MPI_Wtime();
    // mean water
    // double mean_water = 0;
    // for (int out_i = 0; out_i < pht_result.size; out_i++) {
    //   mean_water += pht_result.out_values[out_i][0];
    // }
    // mean_water /= pht_result.size;

    interp_result.results[i] =
        f_interpolate(dht_instance.getKeyElements(), dht_results.results[i],
                      pht_result.in_values, pht_result.out_values);

    if (interp_result.results[i][7] < 0 || interp_result.results[i][9] < 0) {
      interp_result.status[i] = INSUFFICIENT_DATA;
      continue;
    }

    // interp_result.results[i][0] = mean_water;
    this->interp_t += MPI_Wtime() - start_fc;

    this->interpolations++;

    curr_mapping.erase(std::remove(curr_mapping.begin(), curr_mapping.end(), i),
                       curr_mapping.end());

    dht_results.needPhreeqc[i] = false;
    interp_result.status[i] = RES_OK;
  }
}

void InterpolationModule::resultsToWP(std::vector<double> &work_package) const {
  for (uint32_t i = 0; i < interp_result.status.size(); i++) {
    if (interp_result.status[i] == RES_OK) {
      const std::size_t length =
          interp_result.results[i].end() - interp_result.results[i].begin();
      std::copy(interp_result.results[i].begin(),
                interp_result.results[i].end(),
                work_package.begin() + (length * i));
    }
  }
}

} // namespace poet
