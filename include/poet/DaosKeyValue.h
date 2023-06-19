
/**
 * @file DaosKeyValue.h
 * @author Nico Sauerbrei (nico.sauerbrei@uni-potsdam.de)
 * @brief API to interact with DAOS
 * @version 0.1
 * @date 01 Jun 2023
 *
 * This file implements the communication between POET and the DAOS
 * Key-Value Store
 */

#ifndef DAOS_KEY_VALUE_H
#define DAOS_KEY_VALUE_H

#include <mpi.h>
#include <stdint.h>

#include <daos.h>

#define DAOS_SUCCESS 0
#define DAOS_ERROR -1
#define DAOS_MPI_ERROR -2
#define DAOS_READ_MISS -3

/**
 * Internal struct to store statistics about read and write accesses and also
 * read misses and evictions.
 * <b>All values will be resetted to zero after a call of
 * DHT_print_statistics().</b>
 * Internal use only!
 *
 */
typedef struct
{
  /** Count of writes to specific process this process did. */
  int *writes_local;
  /** Writes after last call of DHT_print_statistics. */
  int old_writes;
  /** How many read misses occur? */
  int read_misses;
  /** How many buckets where evicted? */
  int evictions;
  /** How many calls of DHT_write() did this process? */
  int w_access;
  /** How many calls of DHT_read() did this process? */
  int r_access;
} DAOSKV_stats;

/**
 * Struct which serves as a handler or so called \a DHT-object. Will
 * be created by DHT_create and must be passed as a parameter to every following
 * function. Stores all relevant data.
 * Do not touch outside DHT functions!
 */
typedef struct
{

  /** MPI communicator of all participating processes. */
  MPI_Comm communicator;
  /** Size of the MPI communicator respectively all participating processes. */
  int comm_size;
  /** Rank of the process in the MPI communicator. */
  int rank;
  /** Count of read misses over all time. */
  int read_misses;
  /** Count of evictions over all time. */
  int evictions;

  /**  Label of the DAOS container.*/
  char *cont_label;
  /**  DAOS pool handle.*/
  daos_handle_t poh;
  /**  DAOS container handle.*/
  daos_handle_t coh;
  /**  DAOS object handle.*/
  daos_handle_t oh;

#ifdef DHT_STATISTICS
  /** Detailed statistics of the usage of the DHT. */
  DAOS_stats *stats;
#endif
} DAOSKV;
#ifdef __cplusplus
extern "C"
{
#endif

  extern DAOSKV *DAOSKV_create(MPI_Comm comm);
  extern int DAOSKV_free(DAOSKV *object);
  extern int DAOSKV_write(DAOSKV *object, void *key, int key_size, void *send_data, int send_size);
  extern int DAOSKV_read(DAOSKV *object, void *key, int key_size, void *recv_data, int recv_size);
  extern int DAOSKV_remove(DAOSKV *object, void *key, int key_size);
#ifdef __cplusplus
}
#endif

#endif /* DAOS_KEY_VALUE_H */
