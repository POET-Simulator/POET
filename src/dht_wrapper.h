#pragma once
#include <RInside.h> 
#include <string>
#include <vector>
#include <math.h>
#include "DHT.h"

using namespace std;
using namespace Rcpp;

/*Functions*/
uint64_t get_md5(int key_size, void* key);
void fuzz_for_dht(RInside &R, int var_count, void *key, double dt);
void check_dht(RInside &R, int length, std::vector<bool> &out_result_index, double *work_package);
void fill_dht(RInside &R, int length, std::vector<bool> &result_index, double *work_package, double *results);
void print_statistics();
int table_to_file(char* filename);
int file_to_table(char* filename);

/*globals*/
extern bool dht_enabled;
extern int dht_snaps;
extern std::string dht_file;
extern bool dt_differ;

//Default DHT Size per process in Byte (defaults to 1 GiB)
#define DHT_SIZE_PER_PROCESS 1073741824

//sets default dht access and distribution strategy
#define DHT_STRATEGY 0
// 0 -> DHT is on workers, access from workers only
// 1 -> DHT is on workers + master, access from master only !NOT IMPLEMENTED YET!

#define ROUND(value,signif) (((int) (pow(10.0, (double) signif) * value)) * pow(10.0, (double) -signif))

extern int dht_strategy;
extern int dht_significant_digits;
extern std::vector<int> dht_significant_digits_vector;
extern std::vector<string> prop_type_vector;
extern bool dht_logarithm;
extern uint64_t dht_size_per_process;

//global DHT object, can be NULL if not initialized, check strategy
extern DHT* dht_object;

//DHT Performance counter
extern uint64_t dht_hits, dht_miss, dht_collision;

extern double* fuzzing_buffer;
extern std::vector<bool> dht_flags;
