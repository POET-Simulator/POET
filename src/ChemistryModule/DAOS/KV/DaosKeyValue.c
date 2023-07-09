/*
** Copyright (C) 2017-2021 Max Luebke (University of Potsdam)
**
** POET is free software; you can redistribute it and/or modify it under the
** terms of the GNU General Public License as published by the Free Software
** Foundation; either version 2 of the License, or (at your option) any later
** version.
**
** POET is distributed in the hope that it will be useful, but WITHOUT ANY
** WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
** A PARTICULAR PURPOSE. See the GNU General Public License for more details.
**
** You should have received a copy of the GNU General Public License along with
** this program; if not, write to the Free Software Foundation, Inc., 51
** Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

#include "poet/DaosKeyValue.h"

#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <daos.h>

enum handleType
{
  HANDLE_POOL,
  HANDLE_CO,
};

#define ENUM_DESC_BUF 512
#define ENUM_DESC_NR 5

enum
{
  OBJ_DKEY,
  OBJ_AKEY
};

static inline void
handle_share(DAOSKV *object, int type)
{
  d_iov_t ghdl = {NULL, 0, 0};
  int rc;

  if (object->rank == 0)
  {
    /** fetch size of global handle */
    if (type == HANDLE_POOL)
      rc = daos_pool_local2global(object->poh, &ghdl);
    else
      rc = daos_cont_local2global(object->coh, &ghdl);
  }

  /** broadcast size of global handle to all peers */
  MPI_Bcast(&ghdl.iov_buf_len, 1, MPI_UINT64_T, 0, object->communicator);

  /** allocate buffer for global pool handle */
  ghdl.iov_buf = malloc(ghdl.iov_buf_len);
  ghdl.iov_len = ghdl.iov_buf_len;

  if (object->rank == 0)
  {
    /** generate actual global handle to share with peer tasks */
    if (type == HANDLE_POOL)
      rc = daos_pool_local2global(object->poh, &ghdl);
    else
      rc = daos_cont_local2global(object->coh, &ghdl);
  }

  /** broadcast global handle to all peers */
  MPI_Bcast(ghdl.iov_buf, ghdl.iov_len, MPI_BYTE, 0, object->communicator);

  if (object->rank != 0)
  {
    /** unpack global handle */
    if (type == HANDLE_POOL)
    {
      /* NB: Only pool_global2local are different */
      rc = daos_pool_global2local(ghdl, &object->poh);
    }
    else
    {
      rc = daos_cont_global2local(object->poh, ghdl, &object->coh);
    }
  }

  free(ghdl.iov_buf);

  MPI_Barrier(object->communicator);
}

DAOSKV *DAOSKV_create(MPI_Comm comm)
{

  DAOSKV *object;
  int rc;

  object = (DAOSKV *)malloc(sizeof(DAOSKV));
  if (object == NULL)
    return NULL;

  // fill daos-object
  object->communicator = comm;
  object->read_misses = 0;
  object->evictions = 0;

  MPI_Comm_rank(comm, &object->rank);
  MPI_Comm_size(comm, &object->comm_size);

  // if set, initialize dht_stats
#ifdef DHT_STATISTICS
  DAOSKV_stats *stats;

  stats = (DAOSKV_stats *)malloc(sizeof(DAOSKV_stats));
  if (stats == NULL)
    return NULL;

  object->stats = stats;
  object->stats->writes_local = (int *)calloc(object->comm_size, sizeof(int));
  object->stats->old_writes = 0;
  object->stats->read_misses = 0;
  object->stats->read_hits = 0;
  object->stats->evictions = 0;
  object->stats->w_access = 0;
  object->stats->r_access = 0;
#endif



  /** initialize the local DAOS stack */
  if (daos_init() != 0)
    return NULL;

  /** Call connect on rank 0 only and broadcast handle to others */
  if (object->rank == 0)
  {
    char *pool_name = getenv("DAOS_POOL");
    if(pool_name == NULL)
      printf("Pool label invalid \n");
    if (daos_pool_connect(pool_name, NULL, DAOS_PC_RW, &object->poh,
                          NULL, NULL) != 0)
      return NULL;
  }

  /** share pool handle with peer tasks */
  handle_share(object, HANDLE_POOL);

  /*
   * Create and open container on rank 0 and share the handle.
   *
   * Alternatively, one could create the container outside of this program
   * using the daos utility: daos cont create --pool=puuid
   * and pass the uuid to the app.
   */
  if (object->rank == 0)
  {
    /** create container */
    if (getenv("DAOS_CONT") != NULL)
    {
      object->cont_label= getenv("DAOS_CONT");
    }
    else
    {
      object->cont_label = "Poet_Pool2";
    }

    /** check & open container if it already exist */
    if (0 != daos_cont_open(object->poh, object->cont_label, DAOS_COO_RW, &object->coh, NULL, NULL))
    {
      /** create & open container*/
      daos_cont_create_with_label(object->poh, object->cont_label, NULL, NULL, NULL);

      /** open container */
      if (daos_cont_open(object->poh, object->cont_label, DAOS_COO_RW, &object->coh, NULL, NULL) != 0)
        return NULL;
    }
  }
  /** share container handle with peer tasks */
  handle_share(object, HANDLE_CO);

  /** open object */

  daos_obj_id_t oid;

  oid.hi = 0;
  oid.lo = 4;

  daos_obj_generate_oid(object->coh, &oid, DAOS_OT_KV_HASHED, OC_SX, 0, 0);

  daos_kv_open(object->coh, oid, DAOS_OO_RW, &object->oh, NULL);

  return object;
}

int DAOSKV_free(DAOSKV *object)
{
  MPI_Barrier(object->communicator);

  if (daos_kv_close(object->oh, NULL) != 0)
    return DAOS_ERROR;
  if (daos_cont_close(object->coh, NULL) != 0)
    return DAOS_ERROR;

  if (object->rank == 0)
  {
    daos_cont_destroy(object->poh, "simple_obj", 0, NULL);
  }
  if (daos_pool_disconnect(object->poh, NULL) != 0)
    return DAOS_ERROR;

  if (daos_fini() != 0)
    return DAOS_ERROR;

  return DAOS_SUCCESS;
}

int DAOSKV_write(DAOSKV *object, void *key, int key_size, void *send_data, int send_size)
{
	#ifdef DHT_STATISTICS
		object->stats->w_access++;
	#endif
	
   //Turn key into a string
    char* keyString[(key_size*2)+1];
    
    keyToString(keyString,key,key_size);
    
  if (daos_kv_put(object->oh, DAOS_TX_NONE, 0, keyString, send_size, send_data, NULL) != 0)
    return DAOS_ERROR;

  #ifdef DHT_STATISTICS
  object->stats->writes_local[object->rank]++;
  #endif

  return DAOS_SUCCESS;
}

int DAOSKV_read(DAOSKV *object, void *key, int key_size, void *recv_data, int recv_size)
{

  #ifdef DHT_STATISTICS
    object->stats->r_access++;
  #endif

   //Turn key into a string
    char* keyString[(key_size*2)+1];
    
    keyToString(keyString,key,key_size);
    
    
    daos_size_t size = recv_size;
    int rc;

	rc = daos_kv_get(object->oh, DAOS_TX_NONE , DAOS_COND_DKEY_FETCH, keyString, &size, recv_data, NULL);



  if (rc == -DER_NONEXIST){
	#ifdef DHT_STATISTICS
      object->stats->read_misses += 1;
    #endif
    return DAOS_READ_MISS;
	}
  else if (rc != 0)
	return DAOS_ERROR;
  #ifdef DHT_STATISTICS
    object->stats->read_hits += 1;
  #endif
  return DAOS_SUCCESS;
}

int DAOSKV_remove(DAOSKV *object, void *key, int key_size)
{
	
	 //Turn key into a string
    char* keyString[(key_size*2)+1];
    
    keyToString(keyString,key,key_size);
  int rc;

  if (daos_kv_remove(object->oh, DAOS_TX_NONE, 0, keyString, NULL) != 0)
    return DAOS_ERROR;

  return DAOS_SUCCESS;
}

int keyToString(char* output, void* key, int key_size){
	
	int i;
	int offset = 0;
	for (i = 0; i < key_size; i++)
	{
		
		sprintf((char*) output + offset, "%02X", ((char*) key)[i]);
		offset += 2;
	}
		
	

	output[offset++] = '\0';
	
	return 0;
	
	
	}

int DAOSKV_print_statistics(DAOSKV *object){
  #ifdef DHT_STATISTICS
  int *written_buckets;
  int *read_misses, sum_read_misses;
  int *read_hits, sum_read_hits;
  int *evictions, sum_evictions;
  int sum_w_access, sum_r_access, *w_access, *r_access;
  int rank;
  MPI_Comm_rank(object->communicator, &rank);

// disable possible warning of unitialized variable, which is not the case
// since we're using MPI_Gather to obtain all values only on rank 0
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wmaybe-uninitialized"

  // obtaining all values from all processes in the communicator
  if (rank == 0) read_misses = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->read_misses, 1, MPI_INT, read_misses, 1,
                 MPI_INT, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->read_misses, &sum_read_misses, 1, MPI_INT,
                 MPI_SUM, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->read_misses = 0;

  if (rank == 0) read_hits = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->read_hits, 1, MPI_INT, read_hits, 1,
                 MPI_INT, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->read_hits, &sum_read_hits, 1, MPI_INT,
                 MPI_SUM, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->read_hits = 0;

  if (rank == 0) evictions = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->evictions, 1, MPI_INT, evictions, 1, MPI_INT, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->evictions, &sum_evictions, 1, MPI_INT, MPI_SUM,
                 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->evictions = 0;

  if (rank == 0) w_access = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->w_access, 1, MPI_INT, w_access, 1, MPI_INT, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->w_access, &sum_w_access, 1, MPI_INT, MPI_SUM, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->w_access = 0;

  if (rank == 0) r_access = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->r_access, 1, MPI_INT, r_access, 1, MPI_INT, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->r_access, &sum_r_access, 1, MPI_INT, MPI_SUM, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->r_access = 0;

  if (rank == 0) written_buckets = (int *)calloc(object->comm_size, sizeof(int));
  if (MPI_Reduce(object->stats->writes_local, written_buckets, object->comm_size,
                 MPI_INT, MPI_SUM, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;

  if (rank == 0) {  // only process with rank 0 will print out results as a
                    // object
    int sum_written_buckets = 0;

    for (int i = 0; i < object->comm_size; i++) {
      sum_written_buckets += written_buckets[i];
    }

    int members = 7;
    int padsize = (members * 12) + 1;
    char pad[padsize + 1];

    memset(pad, '-', padsize * sizeof(char));
    pad[padsize] = '\0';
    printf("\n");
    printf("%-35s||resets with every call of this function\n", " ");
    printf("%-11s|%-11s|%-11s||%-11s|%-11s|%-11s|%-11s|%-11s\n", "rank", "occupied",
           "free", "w_access", "r_access", "read misses","read hits", "evictions");
    printf("%s\n", pad);
    for (int i = 0; i < object->comm_size; i++) {
      printf("%-11d|%-11d|%-11d||%-11d|%-11d|%-11d|%-11d|%-11d\n", i,
             written_buckets[i], 0,
             w_access[i], r_access[i], read_misses[i],read_hits[i], evictions[i]);
    }
    printf("%s\n", pad);
    printf("%-11s|%-11d|%-11d||%-11d|%-11d|%-11d|%-11d|%-11d\n", "sum",
           sum_written_buckets,
           0,
           sum_w_access, sum_r_access, sum_read_misses,sum_read_hits, sum_evictions);

    printf("%s\n", pad);
    printf("%s %d\n",
           "new entries:", sum_written_buckets - object->stats->old_writes);

    printf("\n");
    fflush(stdout);

    object->stats->old_writes = sum_written_buckets;
  }

// enable warning again
#pragma GCC diagnostic pop

  MPI_Barrier(object->communicator);
  return DAOS_SUCCESS;
#endif

}
