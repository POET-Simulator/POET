#include "Interpolation.hpp"

#include "DHT_Wrapper.hpp"
#include "HashFunctions.hpp"
#include "LUCX/DHT.h"
#include "LookupKey.hpp"
#include "Rounding.hpp"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <vector>

#include <LUCX/PHT.h>

namespace poet {

ProximityHashTable::ProximityHashTable(uint32_t key_size, uint32_t data_size,
                                       uint32_t entry_count,
                                       uint32_t size_per_process,
                                       MPI_Comm communicator_)
    : communicator(communicator_), prox_ht(new PHT) {

  bucket_store = std::make_unique<char[]>(data_size * entry_count);

  uint32_t buckets_per_process =
      static_cast<std::uint32_t>(size_per_process / (data_size * entry_count));

  int status =
      PHT_create(key_size, data_size, buckets_per_process, entry_count,
                 &poet::Murmur2_64A, communicator, this->prox_ht.get());

  if (status != PHT_SUCCESS) {
    throw std::runtime_error("Failed to create PHT.");
  }
}

ProximityHashTable::~ProximityHashTable() {
  PHT_free(this->prox_ht.get(), NULL);
}

void ProximityHashTable::writeLocationToPHT(LookupKey key,
                                            DHT_Location location) {
  double start = MPI_Wtime();

  int status =
      PHT_write(this->prox_ht.get(), key.data(), &location, NULL, NULL);

  if (status == PHT_WRITE_SUCCESS_WITH_EVICTION) {
    this->dht_evictions++;
  }

  this->pht_write_t += MPI_Wtime() - start;
}

const ProximityHashTable::PHT_Result &ProximityHashTable::query(
    const LookupKey &key, const std::uint32_t min_entries_needed,
    const std::uint32_t input_count, const std::uint32_t output_count) {

  double start_r = MPI_Wtime();
  const auto cache_ret = localCache[key];
  if (cache_ret.first) {
    cache_hits++;
    return (lookup_results = cache_ret.second);
  }

  std::uint32_t slots_read;
  int status =
      PHT_read(prox_ht.get(), key.data(), bucket_store.get(), &slots_read);

  this->pht_read_t += MPI_Wtime() - start_r;

  if (status == PHT_READ_MISS || slots_read < min_entries_needed) {
    this->lookup_results.size = 0;
    return lookup_results;
  }

  DHT_Location *bucket_elements =
      reinterpret_cast<DHT_Location *>(bucket_store.get());

  lookup_results.size = slots_read;

  lookup_results.in_values.clear();
  lookup_results.in_values.resize(slots_read);

  lookup_results.out_values.resize(slots_read);

  for (std::uint32_t i = 0; i < slots_read; i++) {
    DHT_Location &loc = bucket_elements[i];
    int status =
        DHT_read_location(source_dht, loc.first, loc.second, dht_buffer.data());

    if (status == DHT_READ_MISS) {
      continue;
    }

    auto *buffer = reinterpret_cast<double *>(dht_buffer.data());

    lookup_results.in_values[i] =
        std::vector<double>(buffer, buffer + input_count);

    buffer += input_count;
    lookup_results.out_values[i] =
        std::vector<double>(buffer, buffer + output_count);
  }

  if (lookup_results.size != 0) {
    localCache(key, lookup_results);
  }

  return lookup_results;
}

inline bool
ProximityHashTable::similarityCheck(const LookupKey &fine,
                                    const LookupKey &coarse,
                                    const std::vector<uint32_t> &signif) {

  PHT_Rounder rounder;

  for (int i = 0; i < signif.size(); i++) {
    if (!(rounder.round(fine[i], signif[i]) == coarse[i])) {
      return false;
    }
  }

  return true;
}

inline std::vector<double>
ProximityHashTable::convertKeysFromDHT(Lookup_Keyelement *keys_in,
                                       std::uint32_t key_size) {
  std::vector<double> output(key_size);
  DHT_Rounder rounder;
  for (int i = 0; i < key_size; i++) {
    output[i] = rounder.convert(keys_in[i]);
  }

  return output;
}

void ProximityHashTable::Cache::operator()(const LookupKey &key,
                                           const PHT_Result val) {
  const auto elemIt = this->find(key);

  if (elemIt == this->end()) {

    if (this->free_mem < 0) {
      const LookupKey &to_del = this->lru_queue.back();
      const auto elem_d = this->find(to_del);
      this->free_mem += elem_d->second.getMemSize();
      this->erase(to_del);
      this->keyfinder.erase(to_del);
      this->lru_queue.pop_back();
    }

    this->insert({key, val});
    this->lru_queue.emplace_front(key);
    this->keyfinder[key] = lru_queue.begin();
    this->free_mem -= val.getMemSize();
    return;
  }

  elemIt->second = val;
}

std::pair<bool, ProximityHashTable::PHT_Result>
ProximityHashTable::Cache::operator[](const LookupKey &key) {
  const auto elemIt = this->find(key);

  if (elemIt == this->end()) {
    return {false, {}};
  }

  this->lru_queue.splice(lru_queue.begin(), lru_queue, this->keyfinder[key]);
  return {true, elemIt->second};
}

#ifdef POET_PHT_ADD
static int PHT_increment_counter(int in_data_size, void *in_data,
                                 int out_data_size, void *out_data) {
  char *start = reinterpret_cast<char *>(out_data);
  std::uint64_t *counter = reinterpret_cast<std::uint64_t *>(
      start + out_data_size - sizeof(std::uint64_t));
  *counter += 1;

  return 0;
}

void ProximityHashTable::incrementReadCounter(const LookupKey &key) {
  auto *old_func_ptr = this->prox_ht->accumulate_callback;
  DHT_set_accumulate_callback(prox_ht, PHT_increment_counter);
  int ret, dummy;
  DHT_write_accumulate(prox_ht, key.data(), 0, NULL, NULL, NULL, &ret);
  DHT_set_accumulate_callback(prox_ht, old_func_ptr);
}
#endif
} // namespace poet
