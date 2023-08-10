//  Time-stamp: "Last modified 2023-08-09 12:52:37 mluebke"

#ifndef INTERPOLATION_H_
#define INTERPOLATION_H_

#include "DHT.h"
#include "DHT_Wrapper.hpp"
#include "DataStructures.hpp"
#include "LookupKey.hpp"
#include "poet/DHT_Wrapper.hpp"
#include "poet/Rounding.hpp"
#include <cassert>
#include <iostream>
#include <list>
#include <memory>
#include <mpi.h>
#include <string>
#include <utility>
extern "C" {
#include "poet/DHT.h"
}

#include "poet/LookupKey.hpp"
#include <cstdint>
#include <functional>
#include <unordered_map>
#include <vector>

namespace poet {

class ProximityHashTable {
public:
  using bucket_indicator = std::uint64_t;

  ProximityHashTable(uint32_t key_size, uint32_t data_size,
                     uint32_t entry_count, uint32_t size_per_process,
                     MPI_Comm communicator);
  ~ProximityHashTable();

  // delete assign and copy operator
  ProximityHashTable &operator=(const ProximityHashTable *) = delete;
  ProximityHashTable(const ProximityHashTable &) = delete;

  struct PHT_Result {
    std::uint32_t size;
    std::vector<std::vector<double>> in_values;
    std::vector<std::vector<double>> out_values;

    std::uint32_t getSize() const {
      std::uint32_t sum{0};
      for (const auto &results : out_values) {
        sum += results.size() * sizeof(double);
      }
      return sum;
    }
  };

  void setSourceDHT(DHT *src) {
    this->source_dht = src;
    this->dht_key_count = src->key_size / sizeof(Lookup_Keyelement);
    this->dht_data_count = src->data_size / sizeof(double);
    this->dht_buffer.resize(src->data_size + src->key_size);
  }

  void writeLocationToPHT(LookupKey key, DHT_Location location);

  const PHT_Result &query(const LookupKey &key,
                          const std::uint32_t min_entries_needed,
                          const std::uint32_t input_count,
                          const std::uint32_t output_count);

  std::uint64_t getLocations(const LookupKey &key);

  void getEntriesFromLocation(const PHT_Result &locations,
                              const std::vector<uint32_t> &signif);

  void writeStats() { DHT_print_statistics(this->prox_ht); }

  DHT *getDHTObject() { return this->prox_ht; }

  auto getPHTWriteTime() const -> double { return this->pht_write_t; };
  auto getPHTReadTime() const -> double { return this->pht_read_t; };
  auto getDHTGatherTime() const -> double { return this->pht_gather_dht_t; };

  auto getLocalCacheHits() const -> std::vector<std::uint32_t> {
    return this->all_cache_hits;
  }

  void storeAndResetCounter() {
    all_cache_hits.push_back(cache_hits);
    cache_hits = 0;
  }

#if POET_PHT_ADD
  void incrementReadCounter(const LookupKey &key);
#endif

private:
  enum { INTERP_CB_OK, INTERP_CB_FULL, INTERP_CB_ALREADY_IN };

  static int PHT_callback_function(int in_data_size, void *in_data,
                                   int out_data_size, void *out_data);

  static std::vector<double> convertKeysFromDHT(Lookup_Keyelement *keys_in,
                                                std::uint32_t key_size);

  static bool similarityCheck(const LookupKey &fine, const LookupKey &coarse,
                              const std::vector<uint32_t> &signif);

  char *bucket_store;

  class Cache
      : private std::unordered_map<LookupKey, PHT_Result, LookupKeyHasher> {
  public:
    void operator()(const LookupKey &key, const PHT_Result val);

    std::pair<bool, PHT_Result> operator[](const LookupKey &key);
    void flush() { this->clear(); }

  protected:
  private:
    static constexpr std::int64_t MAX_CACHE_SIZE = 100E6;

    std::int64_t free_mem{MAX_CACHE_SIZE};

    std::list<LookupKey> lru_queue;
    using lru_iterator = typename std::list<LookupKey>::iterator;

    std::unordered_map<LookupKey, lru_iterator, LookupKeyHasher> keyfinder;
  };

  Cache localCache;
  DHT *prox_ht;

  std::uint32_t dht_evictions = 0;

  DHT *source_dht = nullptr;

  PHT_Result lookup_results;
  std::vector<char> dht_buffer;

  std::uint32_t dht_key_count;
  std::uint32_t dht_data_count;

  MPI_Comm communicator;

  double pht_write_t = 0.;
  double pht_read_t = 0.;
  double pht_gather_dht_t = 0.;

  std::uint32_t cache_hits{0};
  std::vector<std::uint32_t> all_cache_hits{};
};

class InterpolationModule {
public:
  using InterpFunction = std::vector<double> (*)(
      const std::vector<std::int32_t> &, const std::vector<double> &,
      const std::vector<std::vector<double>> &,
      const std::vector<std::vector<double>> &);

  InterpolationModule(std::uint32_t entries_per_bucket,
                      std::uint64_t size_per_process,
                      std::uint32_t min_entries_needed, DHT_Wrapper &dht,
                      const NamedVector<std::uint32_t> &interp_key_signifs,
                      const std::vector<std::int32_t> &dht_key_indices);

  enum result_status { RES_OK, INSUFFICIENT_DATA, NOT_NEEDED };

  struct InterpolationResult {
    std::vector<std::vector<double>> results;
    std::vector<result_status> status;

    void ResultsToWP(std::vector<double> &currWP);
  };

  void setInterpolationFunction(InterpFunction func) {
    this->f_interpolate = func;
  }

  void setMinEntriesNeeded(std::uint32_t entries) {
    this->min_entries_needed = entries;
  }

  auto getMinEntriesNeeded() { return this->min_entries_needed; }

  void writePairs();

  void tryInterpolation(WorkPackage &work_package);

  void resultsToWP(std::vector<double> &work_package) const;

  auto getPHTWriteTime() const { return pht->getPHTWriteTime(); };
  auto getPHTReadTime() const { return pht->getPHTReadTime(); };
  auto getDHTGatherTime() const { return pht->getDHTGatherTime(); };
  auto getInterpolationTime() const { return this->interp_t; };

  auto getInterpolationCount() const -> std::uint32_t {
    return this->interpolations;
  }

  auto getPHTLocalCacheHits() const -> std::vector<std::uint32_t> {
    return this->pht->getLocalCacheHits();
  }

  void resetCounter() {
    this->interpolations = 0;
    this->pht->storeAndResetCounter();
  }

  void writePHTStats() { this->pht->writeStats(); }
  void dumpPHTState(const std::string &filename) {
    DHT_to_file(this->pht->getDHTObject(), filename.c_str());
  }

  static constexpr std::uint32_t COARSE_DIFF = 2;
  static constexpr std::uint32_t COARSE_SIGNIF_DEFAULT =
      DHT_Wrapper::DHT_KEY_SIGNIF_DEFAULT - COARSE_DIFF;

private:
  void initPHT(std::uint32_t key_count, std::uint32_t entries_per_bucket,
               std::uint32_t size_per_process, MPI_Comm communicator);

  static std::vector<double> dummy(const std::vector<std::int32_t> &,
                                   const std::vector<double> &,
                                   const std::vector<std::vector<double>> &,
                                   const std::vector<std::vector<double>> &) {
    return {};
  }

  double interp_t = 0.;

  std::uint32_t interpolations{0};

  InterpFunction f_interpolate = dummy;

  std::uint32_t min_entries_needed = 5;

  std::unique_ptr<ProximityHashTable> pht;

  DHT_Wrapper &dht_instance;

  NamedVector<std::uint32_t> key_signifs;
  std::vector<std::int32_t> key_indices;

  InterpolationResult interp_result;
  PHT_Rounder rounder;

  LookupKey roundKey(const LookupKey &in_key) {
    LookupKey out_key;

    for (std::uint32_t i = 0; i < key_indices.size(); i++) {
      out_key.push_back(rounder.round(in_key[key_indices[i]], key_signifs[i]));
    }

    // timestep
    out_key.push_back(in_key.back());

    return out_key;
  }
};
} // namespace poet

#endif // INTERPOLATION_H_
