#include "Init/InitialList.hpp"
#include "Interpolation.hpp"

#include "DHT_Wrapper.hpp"
#include "DataStructures/NamedVector.hpp"
#include "LookupKey.hpp"

#include <Rcpp.h>
#include <Rcpp/proxy/ProtectedProxy.h>
#include <Rinternals.h>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <mpi.h>
#include <string>
#include <vector>

#include <LUCX/PHT.h>

namespace poet {

InterpolationModule::InterpolationModule(
    std::uint32_t entries_per_bucket, std::uint64_t size_per_process,
    std::uint32_t min_entries_needed, DHT_Wrapper &dht,
    const NamedVector<std::uint32_t> &interp_key_signifs,
    const std::vector<std::int32_t> &dht_key_indices,
    const std::vector<std::string> &_out_names,
    const InitialList::ChemistryHookFunctions &_hooks)
    : dht_instance(dht), key_signifs(interp_key_signifs),
      key_indices(dht_key_indices), min_entries_needed(min_entries_needed),
      dht_names(dht.getKeySpecies().getNames()), out_names(_out_names),
      hooks(_hooks) {

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

void InterpolationModule::writePairs() {
  const auto in = this->dht_instance.getDHTResults();
  for (int i = 0; i < in.filledDHT.size(); i++) {
    if (in.filledDHT[i]) {
      const auto coarse_key = roundKey(in.keys[i]);
      pht->writeLocationToPHT(coarse_key, in.locations[i]);
    }
  }
}

void InterpolationModule::tryInterpolation(WorkPackage &work_package) {
  interp_result.status.resize(work_package.size, NOT_NEEDED);

  const auto dht_results = this->dht_instance.getDHTResults();

  for (int wp_i = 0; wp_i < work_package.size; wp_i++) {
    if (work_package.input[wp_i][0] != 2) {
      interp_result.status[wp_i] = INSUFFICIENT_DATA;
      continue;
    }

    if (work_package.mapping[wp_i] != CHEM_PQC) {
      interp_result.status[wp_i] = NOT_NEEDED;
      continue;
    }

    const auto rounded_key = roundKey(dht_results.keys[wp_i]);

    auto pht_result =
        pht->query(rounded_key, this->min_entries_needed,
                   dht_instance.getInputCount(), dht_instance.getOutputCount());

    if (pht_result.size < this->min_entries_needed) {
      interp_result.status[wp_i] = INSUFFICIENT_DATA;
      continue;
    }

    if (hooks.interp_pre.isValid()) {
      NamedVector<double> nv_in(this->out_names, work_package.input[wp_i]);

      std::vector<int> rm_indices = Rcpp::as<std::vector<int>>(
          hooks.interp_pre(nv_in, pht_result.in_values));

      pht_result.size -= rm_indices.size();

      if (pht_result.size < this->min_entries_needed) {
        interp_result.status[wp_i] = INSUFFICIENT_DATA;
        continue;
      }

      for (const auto &index : rm_indices) {
        pht_result.in_values.erase(
            std::next(pht_result.in_values.begin(), index - 1));
        pht_result.out_values.erase(
            std::next(pht_result.out_values.begin(), index - 1));
      }
    }

#ifdef POET_PHT_ADD
    this->pht->incrementReadCounter(roundKey(rounded_key));
#endif

    const int cell_id = static_cast<int>(work_package.input[wp_i][0]);

    if (!to_calc_cache.contains(cell_id)) {
      const std::vector<std::int32_t> &to_calc = dht_instance.getKeyElements();
      std::vector<std::int32_t> keep_indices;

      for (std::size_t i = 0; i < to_calc.size(); i++) {
        if (!std::isnan(work_package.input[wp_i][to_calc[i]])) {
          keep_indices.push_back(to_calc[i]);
        }
      }

      to_calc_cache[cell_id] = keep_indices;
    }

    double start_fc = MPI_Wtime();

    work_package.output[wp_i] =
        f_interpolate(to_calc_cache[cell_id], work_package.input[wp_i],
                      pht_result.in_values, pht_result.out_values);

    if (hooks.interp_post.isValid()) {
      NamedVector<double> nv_result(this->out_names, work_package.output[wp_i]);
      if (Rcpp::as<bool>(hooks.interp_post(nv_result))) {
        interp_result.status[wp_i] = INSUFFICIENT_DATA;
        continue;
      }
    }

    // interp_result.results[i][0] = mean_water;
    this->interp_t += MPI_Wtime() - start_fc;

    this->interpolations++;

    work_package.mapping[wp_i] = CHEM_INTERP;
    interp_result.status[wp_i] = RES_OK;
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
