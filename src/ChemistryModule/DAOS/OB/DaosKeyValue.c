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
    if (pool_name == NULL)
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
      object->cont_label = getenv("DAOS_CONT");
    }
    else
    {
      object->cont_label = "Poet_Container";
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

  daos_obj_generate_oid(object->coh, &oid, 0, OC_SX, 0, 0);

  daos_obj_open(object->coh, oid, DAOS_OO_RW, &object->oh, NULL);

  return object;
}

int DAOSKV_free(DAOSKV *object)
{
  MPI_Barrier(object->communicator);

  if (daos_obj_close(object->oh, NULL) != 0)
    return DAOS_ERROR;
  if (daos_cont_close(object->coh, NULL) != 0)
    return DAOS_ERROR;

  if (object->rank == 0)
  {
    daos_cont_destroy(object->poh, object->cont_label, 0, NULL);
  }
  if (daos_pool_disconnect(object->poh, NULL) != 0)
    return DAOS_ERROR;

  if (daos_fini() != 0)
    return DAOS_ERROR;

  return DAOS_SUCCESS;
}

int DAOSKV_write(DAOSKV *object, void *key, int key_size, void *send_data, int send_size)
{
  int rc;

#ifdef DHT_STATISTICS
  object->stats->w_access++;
#endif

  d_iov_t dkey;
  d_sg_list_t sgl;
  d_iov_t sg_iov;
  daos_iod_t iod;

  // set dkey
  d_iov_set(&dkey, key, key_size);

  d_iov_set(&sg_iov, send_data, send_size);
  sgl.sg_nr = 1;
  sgl.sg_nr_out = 0;
  sgl.sg_iovs = &sg_iov;

  // set akey
  // Here for all the same value since dkey differ
  // Maybe other way around is better? Or a mix?
  int akey = 1;
  d_iov_set(&iod.iod_name, &akey, sizeof(int));

  iod.iod_nr = 1;                 /** has to be 1 for single value */
  iod.iod_size = send_size;       /** size of the single value */
  iod.iod_recxs = NULL;           /** recx is ignored for single value */
  iod.iod_type = DAOS_IOD_SINGLE; /** value type of the akey */

  rc = daos_obj_update(object->oh, DAOS_TX_NONE, 0, &dkey, 1, &iod, &sgl,
                       NULL);

  // No space left in storage
  if (rc == -DER_NOSPACE && object->rank == 0)
  {
    trim_Space(object, 10, send_size, key_size);
  }
  if (rc != 0)
    return DAOS_ERROR;

#ifdef DHT_STATISTICS
  object->stats->writes_local[object->rank]++;
#endif

  return DAOS_SUCCESS;
}

int DAOSKV_read(DAOSKV *object, void *key, int key_size, void *recv_data, int recv_size)
{
  int rc;

#ifdef DHT_STATISTICS
  object->stats->r_access++;
#endif

  d_iov_t dkey;
  d_sg_list_t sgl;
  d_iov_t sg_iov;
  daos_iod_t iod;

  // set dkey
  d_iov_set(&dkey, key, key_size);

  d_iov_set(&sg_iov, recv_data, recv_size);
  sgl.sg_nr = 1;
  sgl.sg_nr_out = 0;
  sgl.sg_iovs = &sg_iov;

  // set akey
  // Here for all the same value since dkey differ
  // Maybe other way around is better? Or a mix?
  int akey = 1;
  d_iov_set(&iod.iod_name, &akey, sizeof(int));

  iod.iod_nr = 1;                 /** 1 for single value */
  iod.iod_size = DAOS_REC_ANY;    /** size of the single value, set to DAOS_REC_ANY to check if key was written or not */
  iod.iod_recxs = NULL;           /** recx is ignored for single value */
  iod.iod_type = DAOS_IOD_SINGLE; /** value type of the akey */

  /** fetch a dkey */
  rc = daos_obj_fetch(object->oh, DAOS_TX_NONE, 0, &dkey, 1, &iod, &sgl,
                      NULL, NULL);
  if (rc != 0)
  {
    return DAOS_ERROR;
  }

  if (iod.iod_size == 0)
  {
#ifdef DHT_STATISTICS
    object->stats->read_misses += 1;
#endif
    return DAOS_READ_MISS;
  }
#ifdef DHT_STATISTICS
  object->stats->read_hits += 1;
#endif
  return DAOS_SUCCESS;
}

int DAOSKV_remove(DAOSKV *object, void *key, int key_size)
{
  d_iov_t dkey;

  // set dkey
  d_iov_set(&dkey, key, key_size);
  if (daos_obj_punch_dkeys(object->oh, DAOS_TX_NONE, 0, 1, &dkey, NULL) != 0)
    return DAOS_ERROR;

  return DAOS_SUCCESS;
}

int enumerate_key(DAOSKV *object, int *total_nr, int key_size)
{

  char *buf;
  daos_key_desc_t kds[5];
  daos_anchor_t anchor = {0};
  d_sg_list_t sgl;
  d_iov_t sg_iov;
  int key_nr = 0;
  int rc;

  buf = malloc(key_size);
  d_iov_set(&sg_iov, buf, key_size);
  sgl.sg_nr = 1;
  sgl.sg_nr_out = 0;
  sgl.sg_iovs = &sg_iov;

  while (!daos_anchor_is_eof(&anchor))
  {
    uint32_t nr = 5;

    memset(buf, 0, key_size);

    rc = daos_obj_list_dkey(object->oh, DAOS_TX_NONE, &nr, kds,
                            &sgl, &anchor, NULL);

        //If there is no key, break the loop
    if(buf[0] == '\0'){
      break;
    }
    printf("Enumareted over dkey: %s\n", buf);

    if (rc != 0)
    {
      printf("Error retrieving Key %d \n", rc);
      return DAOS_ERROR;
    }
    if (nr == 0)
      continue;
    key_nr += nr;
  }

  *total_nr = key_nr;

  return DAOS_SUCCESS;
}

int delete_n_entries(DAOSKV *object, int toDelete, int key_size)
{
  daos_handle_t th = DAOS_TX_NONE;

  char *buf;
  daos_key_desc_t kds[5];
  daos_anchor_t anchor = {0};
  d_sg_list_t sgl;
  d_iov_t sg_iov;

  int rc;

  buf = malloc(key_size);
  d_iov_set(&sg_iov, buf, key_size);
  sgl.sg_nr = 1;
  sgl.sg_nr_out = 0;
  sgl.sg_iovs = &sg_iov;
  memset(buf, 0, key_size);

  /* allocate transaction */
  rc = daos_tx_open(object->coh, &th, 0, NULL);

  int key_nr = 0;
  while (!daos_anchor_is_eof(&anchor) && key_nr < toDelete)
  {
    uint32_t nr = 5;

    rc = daos_obj_list_dkey(object->oh, DAOS_TX_NONE, &nr, kds,
                            &sgl, &anchor, NULL);
    
        //If there is no key, break the loop
    if(buf[0] == '\0'){
      break;
    }

    // Add delete of key to transaction th
    printf("Delete dkey: %s\n", buf);
    d_iov_t dkey;
    // set dkey
    d_iov_set(&dkey, buf, key_size);
    if (daos_obj_punch_dkeys(object->oh, th, 0, 1, &dkey, NULL) != 0)
      printf("Delete n Key Error");

    if (rc != 0)
    {
      printf("Error retrieving Key %d \n", rc);
      return DAOS_ERROR;
    }
    if (nr == 0)
      continue;
    key_nr += nr;
  }

  // commit transaction, retry if failure
  rc = daos_tx_commit(th, NULL);
  if (rc)
  {
    printf("Commit error: %d\n", rc);
    if (rc == -DER_TX_RESTART)
    {
      /* conflict with another transaction, try again */
      rc = daos_tx_restart(th, NULL);
    }
  }

  // free transaction resources
  rc = daos_tx_close(th, NULL);
  return DAOS_SUCCESS;
}

struct daos_space get_pool_size(DAOSKV *object)
{
  int rc;
  daos_pool_info_t pinfo = {0};
  struct daos_pool_space *ps = &pinfo.pi_space;

  // query only the space, replace with DPI_ALL for all infos
  pinfo.pi_bits = DPI_SPACE;
  rc = daos_pool_query(object->poh, NULL, &pinfo, NULL, NULL);
  // size of storage
  // printf("Total Size:%d\n", ps->ps_space.s_total[DAOS_MEDIA_SCM]+ps->ps_space.s_total[DAOS_MEDIA_NVME]);
  // printf("Free Size:%d\n", ps->ps_space.s_free[DAOS_MEDIA_SCM]+ ps->ps_space.s_free[DAOS_MEDIA_NVME]);

  return ps->ps_space;
}

int trim_Space(DAOSKV *object, float deletePercentage, int dataSize, int keySize)
{
  // Get current usage of the storage space
  struct daos_space space = get_pool_size(object);

  long int total_size = space.s_total[DAOS_MEDIA_SCM] + space.s_total[DAOS_MEDIA_NVME];

  // Estimate, total number of entries
  int totalNumberOfEntries = total_size / (dataSize + keySize);
  // Calculate how many keys to delete
  int toDeleteEntries = totalNumberOfEntries * deletePercentage / 100;

  delete_n_entries(object, toDeleteEntries, keySize);

  return DAOS_SUCCESS;
}

int DAOSKV_print_statistics(DAOSKV *object)
{
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
  if (rank == 0)
    read_misses = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->read_misses, 1, MPI_INT, read_misses, 1,
                 MPI_INT, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->read_misses, &sum_read_misses, 1, MPI_INT,
                 MPI_SUM, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->read_misses = 0;

  if (rank == 0)
    read_hits = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->read_hits, 1, MPI_INT, read_hits, 1,
                 MPI_INT, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->read_hits, &sum_read_hits, 1, MPI_INT,
                 MPI_SUM, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->read_hits = 0;

  if (rank == 0)
    evictions = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->evictions, 1, MPI_INT, evictions, 1, MPI_INT, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->evictions, &sum_evictions, 1, MPI_INT, MPI_SUM,
                 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->evictions = 0;

  if (rank == 0)
    w_access = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->w_access, 1, MPI_INT, w_access, 1, MPI_INT, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->w_access, &sum_w_access, 1, MPI_INT, MPI_SUM, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->w_access = 0;

  if (rank == 0)
    r_access = (int *)malloc(object->comm_size * sizeof(int));
  if (MPI_Gather(&object->stats->r_access, 1, MPI_INT, r_access, 1, MPI_INT, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  if (MPI_Reduce(&object->stats->r_access, &sum_r_access, 1, MPI_INT, MPI_SUM, 0,
                 object->communicator) != 0)
    return DAOS_MPI_ERROR;
  object->stats->r_access = 0;

  if (rank == 0)
    written_buckets = (int *)calloc(object->comm_size, sizeof(int));
  if (MPI_Reduce(object->stats->writes_local, written_buckets, object->comm_size,
                 MPI_INT, MPI_SUM, 0, object->communicator) != 0)
    return DAOS_MPI_ERROR;

  if (rank == 0)
  { // only process with rank 0 will print out results as a
    // object
    int sum_written_buckets = 0;

    for (int i = 0; i < object->comm_size; i++)
    {
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
           "free", "w_access", "r_access", "read misses", "read hits", "evictions");
    printf("%s\n", pad);
    for (int i = 0; i < object->comm_size; i++)
    {
      printf("%-11d|%-11d|%-11d||%-11d|%-11d|%-11d|%-11d|%-11d\n", i,
             written_buckets[i], 0,
             w_access[i], r_access[i], read_misses[i], read_hits[i], evictions[i]);
    }
    printf("%s\n", pad);
    printf("%-11s|%-11d|%-11d||%-11d|%-11d|%-11d|%-11d|%-11d\n", "sum",
           sum_written_buckets,
           0,
           sum_w_access, sum_r_access, sum_read_misses, sum_read_hits, sum_evictions);

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
