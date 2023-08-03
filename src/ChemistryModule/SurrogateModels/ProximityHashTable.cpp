//  Time-stamp: "Last modified 2023-08-01 17:11:42 mluebke"

#include "poet/DHT_Wrapper.hpp"
#include "poet/HashFunctions.hpp"
#include "poet/Interpolation.hpp"
#include "poet/LookupKey.hpp"
#include "poet/Rounding.hpp"
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <memory>
#include <unordered_set>
#include <vector>
extern "C" {
#include "poet/DHT.h"
}

namespace poet {

ProximityHashTable::ProximityHashTable(uint32_t key_size, uint32_t data_size,
                                       uint32_t entry_count,
                                       uint32_t size_per_process,
                                       MPI_Comm communicator_)
    : communicator(communicator_) {

  data_size *= entry_count;
  data_size += sizeof(bucket_indicator);

#ifdef POET_PHT_ADD
  data_size += sizeof(std::uint64_t);
#endif

  bucket_store = new char[data_size];

  uint32_t buckets_per_process =
      static_cast<std::uint32_t>(size_per_process / (data_size + key_size));

  this->prox_ht = DHT_create(communicator, buckets_per_process, data_size,
                             key_size, &poet::Murmur2_64A);

  DHT_set_accumulate_callback(this->prox_ht, PHT_callback_function);
}

ProximityHashTable::~ProximityHashTable() {
  delete[] bucket_store;
  if (prox_ht) {
    DHT_free(prox_ht, NULL, NULL);
  }
}

int ProximityHashTable::PHT_callback_function(int in_data_size, void *in_data,
                                              int out_data_size,
                                              void *out_data) {
  const int max_elements_per_bucket =
      static_cast<int>((out_data_size - sizeof(bucket_indicator)
#ifdef POET_PHT_ADD
                        - sizeof(std::uint64_t)
#endif
                            ) /
                       in_data_size);
  DHT_Location *input = reinterpret_cast<DHT_Location *>(in_data);

  bucket_indicator *occupied_buckets =
      reinterpret_cast<bucket_indicator *>(out_data);
  DHT_Location *pairs = reinterpret_cast<DHT_Location *>(occupied_buckets + 1);

  if (*occupied_buckets == max_elements_per_bucket) {
    return INTERP_CB_FULL;
  }

  for (bucket_indicator i = 0; i < *occupied_buckets; i++) {
    if (pairs[i] == *input) {
      return INTERP_CB_ALREADY_IN;
    }
  }

  pairs[(*occupied_buckets)++] = *input;

  return INTERP_CB_OK;
}

void ProximityHashTable::writeLocationToPHT(LookupKey key,
                                            DHT_Location location) {

  double start = MPI_Wtime();

  // if (localCache[key].first) {
  //   return;
  // }

  int ret_val;

  int status = DHT_write_accumulate(prox_ht, key.data(), sizeof(location),
                                    &location, NULL, NULL, &ret_val);

  if (status == DHT_WRITE_SUCCESS_WITH_EVICTION) {
    this->dht_evictions++;
  }

  // if (ret_val == INTERP_CB_FULL) {
  //   localCache(key, {});
  // }

  this->pht_write_t += MPI_Wtime() - start;
}

const ProximityHashTable::PHT_Result &ProximityHashTable::query(
    const LookupKey &key, const std::vector<std::uint32_t> &signif,
    std::uint32_t min_entries_needed, std::uint32_t input_count,
    std::uint32_t output_count) {

  double start_r = MPI_Wtime();
  const auto cache_ret = localCache[key];
  if (cache_ret.first) {
    cache_hits++;
    return (lookup_results = cache_ret.second);
  }

  int res = DHT_read(prox_ht, key.data(), bucket_store);
  this->pht_read_t += MPI_Wtime() - start_r;

  if (res != DHT_SUCCESS) {
    this->lookup_results.size = 0;
    return lookup_results;
  }

  auto *bucket_element_count =
      reinterpret_cast<bucket_indicator *>(bucket_store);
  auto *bucket_elements =
      reinterpret_cast<DHT_Location *>(bucket_element_count + 1);

  if (*bucket_element_count < min_entries_needed) {
    this->lookup_results.size = 0;
    return lookup_results;
  }

  lookup_results.size = *bucket_element_count;
  auto locations = std::vector<DHT_Location>(
      bucket_elements, bucket_elements + *(bucket_element_count));

  lookup_results.in_values.clear();
  lookup_results.in_values.reserve(*bucket_element_count);

  lookup_results.out_values.clear();
  lookup_results.out_values.reserve(*bucket_element_count);

  for (const auto &loc : locations) {
    double start_g = MPI_Wtime();
    DHT_read_location(source_dht, loc.first, loc.second, dht_buffer.data());
    this->pht_gather_dht_t += MPI_Wtime() - start_g;

    auto *buffer = reinterpret_cast<double *>(dht_buffer.data());

    lookup_results.in_values.push_back(
        std::vector<double>(buffer, buffer + input_count));

    buffer += input_count;
    lookup_results.out_values.push_back(
        std::vector<double>(buffer, buffer + output_count));

    // if (!similarityCheck(check_key, key, signif)) {
    //   // TODO: original stored location in PHT was overwritten in DHT.
    //   Need to
    //   // handle this!
    //   lookup_results.size--;
    //   if (lookup_results.size < min_entries_needed) {
    //     lookup_results.size = 0;
    //     break;
    //   }
    //   continue;
    // }

    // auto input = convertKeysFromDHT(buffer_start, dht_key_count);
    // // remove timestep from the key
    // input.pop_back();
    // lookup_results.in_keys.push_back(input);

    // auto *data = reinterpret_cast<double *>(buffer + dht_key_count);
    // lookup_results.out_values.push_back(
    //     std::vector<double>(data, data + dht_data_count));
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
      this->free_mem += elem_d->second.getSize();
      this->erase(to_del);
      this->keyfinder.erase(to_del);
      this->lru_queue.pop_back();
    }

    this->insert({key, val});
    this->lru_queue.emplace_front(key);
    this->keyfinder[key] = lru_queue.begin();
    this->free_mem -= val.getSize();
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
