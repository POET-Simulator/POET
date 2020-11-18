#include "DHT.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>

/*
 * determining destination rank and index by hash of key
 *
 * return values by reference
 */
static void determine_dest(uint64_t hash, int comm_size, unsigned int table_size, unsigned int *dest_rank, unsigned int *index, unsigned int index_count) {
    uint64_t tmp;
    int char_hop = 9-index_count;
    unsigned int i;
    for (i = 0; i < index_count ; i++) {
        tmp = 0;
        memcpy(&tmp,(unsigned char *)&hash+i, char_hop);
        index[i] = (unsigned int) (tmp % table_size);
    }
    *dest_rank = (unsigned int) (hash % comm_size);
}

/**
 * set write flag to 1
 */
static void set_flag(char* flag_byte) {
    *flag_byte = 0;
    *flag_byte |= (1 << 0);
} 

/**
 * return 1 if write flag is set
 * else 0
 */
static int read_flag(char flag_byte) {
    if ((flag_byte & 0x01) == 0x01) {
        return 1;
    } else return 0;
}

/*
 * allocating memory for DHT object and buckets.
 * creating MPI window for OSC
 * filling DHT object with passed by parameters, window, 2 counters for R/W errors and 2 pointers with allocated memory for further use
 * return DHT object
 */
DHT* DHT_create(MPI_Comm comm, unsigned int size, int data_size, int key_size, uint64_t(*hash_func) (int, void*)) {
    DHT *object;
    MPI_Win window;
    void* mem_alloc;
    int comm_size, tmp;

    tmp = (int) ceil(log2(size));
    if (tmp%8 != 0) tmp = tmp + (8-(tmp%8));

    object = (DHT*) malloc(sizeof(DHT));
    if (object == NULL) return NULL;

    //every memory allocation has 1 additional byte for flags etc.
    if (MPI_Alloc_mem(size * (1 + data_size + key_size), MPI_INFO_NULL, &mem_alloc) != 0) return NULL;
    if (MPI_Comm_size(comm, &comm_size) != 0) return NULL;

    memset(mem_alloc, '\0', size * (1 + data_size + key_size));

    if (MPI_Win_create(mem_alloc, size * (1 + data_size + key_size), (1 + data_size + key_size), MPI_INFO_NULL, comm, &window) != 0) return NULL;

    object->data_size = data_size;
    object->key_size = key_size;
    object->table_size = size;
    object->window = window;
    object->hash_func = hash_func;
    object->comm_size = comm_size;
    object->communicator = comm;
    object->read_misses = 0;
    object->collisions = 0;
    object->recv_entry = malloc(1 + data_size + key_size);
    object->send_entry = malloc(1 + data_size + key_size);
    object->index_count = 9-(tmp/8);
    object->index = (unsigned int*) malloc((9-(tmp/8))*sizeof(int));
    object->mem_alloc = mem_alloc;

    DHT_stats *stats;

    stats = (DHT_stats*) malloc(sizeof(DHT_stats));
    if (stats == NULL) return NULL;

    object->stats = stats;
    object->stats->writes_local = (int*) calloc(comm_size, sizeof(int));
    object->stats->old_writes = 0;
    object->stats->read_misses = 0;
    object->stats->collisions = 0;
    object->stats->w_access = 0;
    object->stats->r_access = 0;

    return object;
}

/*
 * puts passed by data with key to DHT
 *
 * returning DHT_MPI_ERROR = -1 if MPI error occurred
 * else  DHT_SUCCESS = 0 if success
 */
int DHT_write(DHT *table, void* send_key, void* send_data) {
    unsigned int dest_rank, i;
    int result = DHT_SUCCESS;

    table->stats->w_access++;

    //determine destination rank and index by hash of key
    determine_dest(table->hash_func(table->key_size, send_key), table->comm_size, table->table_size, &dest_rank, table->index, table->index_count);

    //concatenating key with data to write entry to DHT
    set_flag((char *) table->send_entry);
    memcpy((char *) table->send_entry + 1, (char *) send_key, table->key_size);
    memcpy((char *) table->send_entry + table->key_size + 1, (char *) send_data, table->data_size);

    //locking window of target rank with exclusive lock
    if (MPI_Win_lock(MPI_LOCK_EXCLUSIVE, dest_rank, 0, table->window) != 0)
        return DHT_MPI_ERROR;
    for (i = 0; i < table->index_count; i++)
    {
        if (MPI_Get(table->recv_entry, 1 + table->data_size + table->key_size, MPI_BYTE, dest_rank, table->index[i], 1 + table->data_size + table->key_size, MPI_BYTE, table->window) != 0) return DHT_MPI_ERROR;
        if (MPI_Win_flush(dest_rank, table->window) != 0) return DHT_MPI_ERROR;

        //increment collision counter if receiving key doesn't match sending key
        //,entry has write flag + last index is reached

        if (read_flag(*(char *)table->recv_entry)) {
            if (memcmp(send_key, (char *) table->recv_entry + 1, table->key_size) != 0) {
                if (i == (table->index_count)-1) {
                    table->collisions += 1;
                    table->stats->collisions += 1;
                    result = DHT_WRITE_SUCCESS_WITH_COLLISION;
                    break;
                }
            } else break;
        } else {
            table->stats->writes_local[dest_rank]++;
            break;
        }
    }

    //put data to DHT
    if (MPI_Put(table->send_entry, 1 + table->data_size + table->key_size, MPI_BYTE, dest_rank, table->index[i], 1 + table->data_size + table->key_size, MPI_BYTE, table->window) != 0) return DHT_MPI_ERROR;
    //unlock window of target rank
    if (MPI_Win_unlock(dest_rank, table->window) != 0) return DHT_MPI_ERROR;

    return result;
}

/*
 * gets data from the DHT by key
 *
 * return  DHT_SUCCESS = 0 if success
 * DHT_MPI_ERROR = -1 if MPI error occurred
 * DHT_READ_ERROR = -2 if receiving key doesn't match sending key
 */
int DHT_read(DHT *table, void* send_key, void* destination) {
    unsigned int dest_rank, i;

    table->stats->r_access++;

    //determine destination rank and index by hash of key
    determine_dest(table->hash_func(table->key_size, send_key), table->comm_size, table->table_size, &dest_rank, table->index, table->index_count);

    //locking window of target rank with shared lock
    if (MPI_Win_lock(MPI_LOCK_SHARED, dest_rank, 0, table->window) != 0) return DHT_MPI_ERROR;
    //receive data
    for (i = 0; i < table->index_count; i++) {
        if (MPI_Get(table->recv_entry, 1 + table->data_size + table->key_size, MPI_BYTE, dest_rank, table->index[i], 1 + table->data_size + table->key_size, MPI_BYTE, table->window) != 0) return DHT_MPI_ERROR;
        if (MPI_Win_flush(dest_rank, table->window) != 0) return DHT_MPI_ERROR;
        
        //increment read error counter if write flag isn't set or key doesn't match passed by key + last index reached
        //else copy data to dereference of passed by destination pointer

        if ((read_flag(*(char *) table->recv_entry)) == 0) {
            table->read_misses += 1;
            table->stats->read_misses += 1;
            if (MPI_Win_unlock(dest_rank, table->window) != 0) return DHT_MPI_ERROR;
            return DHT_READ_ERROR;
        }

        if (memcmp(((char*)table->recv_entry) + 1, send_key, table->key_size) != 0) {
            if (i == (table->index_count)-1) {
                table->read_misses += 1;
                table->stats->read_misses += 1;
                if (MPI_Win_unlock(dest_rank, table->window) != 0) return DHT_MPI_ERROR;
                return DHT_READ_ERROR;
            }
        } else break;
    }

    //unlock window of target rank
    if (MPI_Win_unlock(dest_rank, table->window) != 0) return DHT_MPI_ERROR;

    memcpy((char *) destination, (char *) table->recv_entry + table->key_size + 1, table->data_size);

    return DHT_SUCCESS;
}

int DHT_to_file(DHT* table, char* filename) {
    //open file
    MPI_File file;
    if (MPI_File_open(table->communicator, filename, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &file) != 0) return DHT_FILE_IO_ERROR;

    int rank;
    MPI_Comm_rank(table->communicator, &rank);

    //write header (key_size and data_size)
    if (rank == 0) {
        if (MPI_File_write(file, &table->key_size, 1, MPI_INT, MPI_STATUS_IGNORE) != 0) return DHT_FILE_WRITE_ERROR;
        if (MPI_File_write(file, &table->data_size, 1, MPI_INT, MPI_STATUS_IGNORE) != 0) return DHT_FILE_WRITE_ERROR;
    }

    if (MPI_File_seek_shared(file, DHT_HEADER_SIZE, MPI_SEEK_SET) != 0) return DHT_FILE_IO_ERROR;

    char* ptr;
    int bucket_size = table->key_size + table->data_size + 1;

    //iterate over local memory
    for (unsigned int i = 0; i < table->table_size; i++) {
        ptr = (char *) table->mem_alloc + (i * bucket_size);
        //if bucket has been written to (checked by written_flag)...
        if (read_flag(*ptr)) {
            //write key and data to file
            if (MPI_File_write_shared(file, ptr + 1, bucket_size - 1, MPI_BYTE, MPI_STATUS_IGNORE) != 0) return DHT_FILE_WRITE_ERROR;
        }
    }
    //close file
    if (MPI_File_close(&file) != 0) return DHT_FILE_IO_ERROR;

    return DHT_SUCCESS;  
}

int DHT_from_file(DHT* table, char* filename) {
    MPI_File file;
    MPI_Offset f_size;
    int e_size, m_size, cur_pos, rank, offset;
    char* buffer;
    void* key;
    void* data;

    if (MPI_File_open(table->communicator, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &file) != 0) return DHT_FILE_IO_ERROR;

    if (MPI_File_get_size(file, &f_size) != 0) return DHT_FILE_IO_ERROR;

    MPI_Comm_rank(table->communicator, &rank);

    e_size = table->key_size + table->data_size;
    m_size = e_size > DHT_HEADER_SIZE ? e_size : DHT_HEADER_SIZE;
    buffer = (char *) malloc(m_size);

    if (MPI_File_read(file, buffer, DHT_HEADER_SIZE, MPI_BYTE, MPI_STATUS_IGNORE) != 0) return DHT_FILE_READ_ERROR;

    if (*(int *) buffer != table->key_size) return DHT_WRONG_FILE;
    if (*(int *) (buffer + 4) != table->data_size) return DHT_WRONG_FILE;

    offset = e_size*table->comm_size;

    if (MPI_File_seek(file, DHT_HEADER_SIZE, MPI_SEEK_SET) != 0) return DHT_FILE_IO_ERROR;
    cur_pos = DHT_HEADER_SIZE + (rank * e_size);

    while(cur_pos < f_size) {
        if (MPI_File_seek(file, cur_pos, MPI_SEEK_SET) != 0) return DHT_FILE_IO_ERROR;
        MPI_Offset tmp;
        MPI_File_get_position(file, &tmp);
        if (MPI_File_read(file, buffer, e_size, MPI_BYTE, MPI_STATUS_IGNORE) != 0) return DHT_FILE_READ_ERROR;
        key = buffer;
        data = (buffer+table->key_size);
        if (DHT_write(table, key, data) == DHT_MPI_ERROR) return DHT_MPI_ERROR;

        cur_pos += offset;
    }

    free (buffer);
    if (MPI_File_close(&file) != 0) return DHT_FILE_IO_ERROR; 

    return DHT_SUCCESS;
}

/*
 * frees up memory and accumulate counter
 */
int DHT_free(DHT* table, int* collision_counter, int* readerror_counter) {
    int buf;

    if (collision_counter != NULL) {
        buf = 0;
        if (MPI_Reduce(&table->collisions, &buf, 1, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;
        *collision_counter = buf;
    }
    if (readerror_counter != NULL) {
        buf = 0;
        if (MPI_Reduce(&table->read_misses, &buf, 1, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;
        *readerror_counter = buf;
    }
    if (MPI_Win_free(&(table->window)) != 0) return DHT_MPI_ERROR;
    if (MPI_Free_mem(table->mem_alloc) != 0) return DHT_MPI_ERROR;
    free(table->recv_entry);
    free(table->send_entry);
    free(table->index);

    free(table->stats->writes_local);
    free(table->stats);

    free(table);

    return DHT_SUCCESS;
}

/*
 * prints a table with statistics about current use of DHT 
 * for each participating process and summed up results containing:
    *  1. occupied buckets (in respect to the memory of this process)
    *  2. free buckets (in respect to the memory of this process)
    *  3. calls of DHT_write (w_access)
    *  4. calls of DHT_read (r_access)
    *  5. read misses (see DHT_READ_ERROR)
    *  6. collisions (see DHT_WRITE_SUCCESS_WITH_COLLISION)
 * 3-6 will reset with every call of this function
 * finally the amount of new written entries is printed out (in relation to last call of this funtion)
 */
int DHT_print_statistics(DHT *table) {
    int *written_buckets;
    int *read_misses, sum_read_misses;
    int *collisions, sum_collisions;
    int sum_w_access, sum_r_access, *w_access, *r_access;
    int rank;

    MPI_Comm_rank(table->communicator, &rank);

    //disable possible warning of unitialized variable, which is not the case
    //since we're using MPI_Gather to obtain all values only on rank 0
    #pragma GCC diagnostic push
    #pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

    //obtaining all values from all processes in the communicator
    if (rank == 0) read_misses = (int*) malloc(table->comm_size*sizeof(int));
    if (MPI_Gather(&table->stats->read_misses, 1, MPI_INT, read_misses, 1, MPI_INT, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    if (MPI_Reduce(&table->stats->read_misses, &sum_read_misses, 1, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    table->stats->read_misses = 0;
    
    if (rank == 0) collisions = (int*) malloc(table->comm_size*sizeof(int));
    if (MPI_Gather(&table->stats->collisions, 1, MPI_INT, collisions, 1, MPI_INT, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    if (MPI_Reduce(&table->stats->collisions, &sum_collisions, 1, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    table->stats->collisions = 0;

    if (rank == 0) w_access = (int*) malloc(table->comm_size*sizeof(int));
    if (MPI_Gather(&table->stats->w_access, 1, MPI_INT, w_access, 1, MPI_INT, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    if (MPI_Reduce(&table->stats->w_access, &sum_w_access, 1, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    table->stats->w_access = 0;

    if (rank == 0) r_access = (int*) malloc(table->comm_size*sizeof(int));
    if (MPI_Gather(&table->stats->r_access, 1, MPI_INT, r_access, 1, MPI_INT, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    if (MPI_Reduce(&table->stats->r_access, &sum_r_access, 1, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;
    table->stats->r_access = 0;
    
    if (rank == 0) written_buckets = (int*) calloc(table->comm_size, sizeof(int));
    if (MPI_Reduce(table->stats->writes_local, written_buckets, table->comm_size, MPI_INT, MPI_SUM, 0, table->communicator) != 0) return DHT_MPI_ERROR;

    if (rank == 0) { //only process with rank 0 will print out results as a table
        int sum_written_buckets = 0;

        for (int i=0; i < table->comm_size; i++) {
            sum_written_buckets += written_buckets[i];
        }

        int members = 7;
        int padsize = (members*12)+1;
        char pad[padsize+1];

        memset(pad, '-', padsize*sizeof(char));
        pad[padsize]= '\0';
        printf("\n");
        printf("%-35s||resets with every call of this function\n", " ");
        printf("%-11s|%-11s|%-11s||%-11s|%-11s|%-11s|%-11s\n", 
            "rank", 
            "occupied", 
            "free", 
            "w_access",
            "r_access",
            "read misses",
            "collisions");
        printf("%s\n", pad);
        for (int i = 0; i < table->comm_size; i++) {
            printf("%-11d|%-11d|%-11d||%-11d|%-11d|%-11d|%-11d\n", 
                i, 
                written_buckets[i], 
                table->table_size-written_buckets[i], 
                w_access[i],
                r_access[i],
                read_misses[i],
                collisions[i]);
        }
        printf("%s\n", pad);
        printf("%-11s|%-11d|%-11d||%-11d|%-11d|%-11d|%-11d\n", 
            "sum", 
            sum_written_buckets, 
            (table->table_size*table->comm_size)-sum_written_buckets, 
            sum_w_access,
            sum_r_access,
            sum_read_misses,
            sum_collisions);

        printf("%s\n", pad);
        printf("%s %d\n", 
            "new entries:", 
            sum_written_buckets - table->stats->old_writes);

        printf("\n");
        fflush(stdout);

        table->stats->old_writes = sum_written_buckets;
    }

    //enable warning again
    #pragma GCC diagnostic pop

    MPI_Barrier(table->communicator);
    return DHT_SUCCESS;
}