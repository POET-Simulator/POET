#include "dht_wrapper.h"
#include <openssl/md5.h>

using namespace poet;

/*init globals*/
bool dht_enabled;
int dht_snaps;
int dht_strategy;
int dht_significant_digits;
std::string dht_file;
std::vector<int> dht_significant_digits_vector;
std::vector<string> prop_type_vector;
bool dht_logarithm;
uint64_t dht_size_per_process;
uint64_t dht_hits, dht_miss, dht_collision;
RInside *R_DHT;
std::vector<bool> dht_flags;
DHT *dht_object;

double *fuzzing_buffer;

bool dt_differ;
/*functions*/

uint64_t get_md5(int key_size, void *key) {
    MD5_CTX ctx;
    unsigned char sum[MD5_DIGEST_LENGTH];
    uint64_t retval, *v1, *v2;

    MD5_Init(&ctx);
    MD5_Update(&ctx, key, key_size);
    MD5_Final(sum, &ctx);

    v1 = (uint64_t *)&sum[0];
    v2 = (uint64_t *)&sum[8];
    retval = *v1 ^ *v2;

    return retval;
}

double Round_off(RRuntime R, double N, double n) {
    double result;
    R["roundsig"] = n;
    R["roundin"] = N;

    result = R.parseEval("signif(roundin, digits=roundsig)");

    return result;
}

/*
*   Stores fuzzed version of key in fuzzing_buffer
*/
void fuzz_for_dht(RRuntime R, int var_count, void *key, double dt) {
    unsigned int i = 0;
    //introduce fuzzing to allow more hits in DHT
    for (i = 0; i < (unsigned int)var_count; i++) {
        if (prop_type_vector[i] == "act") {
            //with log10
            if (dht_logarithm) {
                if (((double *)key)[i] < 0)
                    cerr << "dht_wrapper.cpp::fuzz_for_dht(): Warning! Negative value at key!" << endl;
                else if (((double *)key)[i] == 0)
                    fuzzing_buffer[i] = 0;
                else
                    //fuzzing_buffer[i] = Round_off(R, std::log10(((double *)key)[i]), dht_significant_digits_vector[i] - 1);
                    fuzzing_buffer[i] = ROUND(-(std::log10(((double *)key)[i])), dht_significant_digits_vector[i]);
            } else {
                //without log10
                //fuzzing_buffer[i] = Round_off(R, ((double *)key)[i], dht_significant_digits_vector[i]);
                fuzzing_buffer[i] = ROUND((((double *)key)[i]), dht_significant_digits_vector[i]);
            }
        } else if (prop_type_vector[i] == "logact") {
            //fuzzing_buffer[i] = Round_off(R, ((double *)key)[i], dht_significant_digits_vector[i]);
            fuzzing_buffer[i] = ROUND((((double *)key)[i]), dht_significant_digits_vector[i]);
        } else if (prop_type_vector[i] == "ignore") {
            fuzzing_buffer[i] = 0;
        } else {
            cerr << "dht_wrapper.cpp::fuzz_for_dht(): Warning! Probably wrong prop_type!" << endl;
        }
    }
    if (dt_differ)
        fuzzing_buffer[var_count] = dt;
}

void check_dht(RRuntime R, int length, std::vector<bool> &out_result_index, double *work_package) {
    void *key;
    int res;
    int var_count = prop_type_vector.size();
    double dt;

    dt = R.parseEval("mysetup$dt");

    for (int i = 0; i < length; i++) {
        key = (void *)&(work_package[i * var_count]);

        //fuzz data (round, logarithm etc.)
        fuzz_for_dht(R, var_count, key, dt);

        //overwrite input with data from DHT, IF value is found in DHT
        res = DHT_read(dht_object, fuzzing_buffer, key);

        if (res == DHT_SUCCESS) {
            //flag that this line is replaced by DHT-value, do not simulate!!
            out_result_index[i] = false;
            dht_hits++;
        } else if (res == DHT_READ_ERROR) {
            //this line is untouched, simulation is needed
            out_result_index[i] = true;
            dht_miss++;
        } else {
            //MPI ERROR ... WHAT TO DO NOW?
            //RUNNING CIRCLES WHILE SCREAMING
        }
    }
}

void fill_dht(RRuntime R, int length, std::vector<bool> &result_index, double *work_package, double *results) {
    void *key;
    void *data;
    int res;
    int var_count = prop_type_vector.size();
    double dt;

    dt = R.parseEval("mysetup$dt");

    for (int i = 0; i < length; i++) {
        key = (void *)&(work_package[i * var_count]);
        data = (void *)&(results[i * var_count]);

        if (result_index[i]) {
            //If true -> was simulated, needs to be inserted into dht

            //fuzz data (round, logarithm etc.)
            fuzz_for_dht(R, var_count, key, dt);

            res = DHT_write(dht_object, fuzzing_buffer, data);

            if (res != DHT_SUCCESS) {
                if (res == DHT_WRITE_SUCCESS_WITH_COLLISION) {
                    dht_collision++;
                } else {
                    //MPI ERROR ... WHAT TO DO NOW?
                    //RUNNING CIRCLES WHILE SCREAMING
                }
            }
        }
    }
}

void print_statistics() {
    int res;

    res = DHT_print_statistics(dht_object);

    if (res != DHT_SUCCESS) {
        //MPI ERROR ... WHAT TO DO NOW?
        //RUNNING CIRCLES WHILE SCREAMING
    }
}

int table_to_file(char *filename) {
    int res = DHT_to_file(dht_object, filename);
    return res;
}

int file_to_table(char *filename) {

    int res = DHT_from_file(dht_object, filename);
    if (res != DHT_SUCCESS) 
        return res;

    DHT_print_statistics(dht_object);

    return DHT_SUCCESS;
}
