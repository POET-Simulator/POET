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

  // if set, initialize dht_stats
#ifdef DHT_STATISTICS
  DHT_stats *stats;

  stats = (DHT_stats *)malloc(sizeof(DHT_stats));
  if (stats == NULL)
    return NULL;

  object->stats = stats;
  object->stats->writes_local = (int *)calloc(comm_size, sizeof(int));
  object->stats->old_writes = 0;
  object->stats->read_misses = 0;
  object->stats->evictions = 0;
  object->stats->w_access = 0;
  object->stats->r_access = 0;
#endif

  MPI_Comm_rank(comm, &object->rank);
  MPI_Comm_size(comm, &object->comm_size);

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

  if (rc != 0)
    return DAOS_ERROR;

  return DAOS_SUCCESS;
}

int DAOSKV_read(DAOSKV *object, void *key, int key_size, void *recv_data, int recv_size)
{
  int rc;

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
    return DAOS_ERROR;

  if (iod.iod_size == 0)
    return DAOS_READ_MISS;

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
