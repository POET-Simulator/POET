/*
 * File:   DHT.h
 * Author: max luebke
 *
 * Created on 16. November 2017, 09:14
 */

#ifndef DHT_H
#define DHT_H

#include <mpi.h>
#include <stdint.h>

#define DHT_MPI_ERROR -1
#define DHT_READ_ERROR -2
#define DHT_SUCCESS 0
#define DHT_WRITE_SUCCESS_WITH_COLLISION 1

#define DHT_WRONG_FILE 11
#define DHT_FILE_IO_ERROR 12
#define DHT_FILE_READ_ERROR 13
#define DHT_FILE_WRITE_ERROR 14

#define DHT_HEADER_SIZE 8

typedef struct {;
    int *writes_local, old_writes;
    int read_misses, collisions;
    int w_access, r_access;
} DHT_stats;

typedef struct {
    MPI_Win window;
    int data_size;
    int key_size;
    unsigned int table_size;
    MPI_Comm communicator;
    int comm_size;
    uint64_t(*hash_func) (int, void*);
    void* recv_entry;
    void* send_entry;
    void* mem_alloc;
    int read_misses;
    int collisions;
    unsigned int *index;
    unsigned int index_count;
    DHT_stats *stats;
} DHT;



/*
 * parameters:
 * MPI_Comm comm - communicator of processes that are holding the DHT
 * int size_per_process - number of buckets each process will create
 * int data_size - size of data in bytes
 * int key_size - size of key in bytes
 * *hash_func - pointer to hashfunction
 *
 * return:
 * NULL if error during initialization
 * DHT* if success
 */
extern DHT* DHT_create(MPI_Comm comm, unsigned int size_per_process, int data_size, int key_size, uint64_t(*hash_func)(int, void*));

/*
 * parameters:
 * DHT *table - DHT_object created by DHT_create
 * void* data - pointer to data
 * void* - pointer to key
 *
 * return:
 * error value (see above)
 */
extern int DHT_write(DHT *table, void* key, void* data);

/*
 * parameters:
 * DHT *table - DHT_object created by DHT_create
 * void* key - pointer to key
 * void* destination - pointer which will hold the resulting data from DHT
 *
 * return:
 * error value (see above)
 */
extern int DHT_read(DHT *table, void* key, void* destination);

extern int DHT_to_file(DHT *table, char* filename);

extern int DHT_from_file(DHT *table, char* filename);

/*
 * parameters:
 * DHT *table - DHT_object created by DHT_create
 * int* collision_counter - pointer which will hold the total count of collisions
 * int* readerror_counter - pointer which will hold the total count of read errors
 *
 * return:
 * error value (see above)
 */
extern int DHT_free(DHT *table, int* collision_counter, int* readerror_counter);

/*
 * parameters:
 * DHT *table - DHT_object created by DHT_create
 * 
 * return:
 * error value (DHT_SUCCESS or DHT_MPI_ERROR)
 */
extern int DHT_print_statistics(DHT *table);

#endif /* DHT_H */