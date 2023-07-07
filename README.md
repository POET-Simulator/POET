<!--
    Time-stamp: "Last modified 2023-01-19 12:06:10 delucia"
-->

# POET Daos

## Branch Information

This branch replaces POET's MPI-based key-value store with the one provided by
the DAOS-API, which is meant to be run in an environment which is connected with
a DAOS server.

## Installation

To install this version of POET please follow the instructions in the main branche.
The DAOS library must be build on the system.

To build POET with DAOS as the key-value store either the POET_DAOSKV or POET_DAOSOB flag needs to be true,
where POET_DAOSKV uses DAOS KV API, and POET_DAOSOB uses DAOS object API. For example:

```sh
cmake -DPOET_DAOSOB=True -DCMAKE_INSTALL_PREFIX=/home/dtest/poet ..
```

## Running

To run POET with DAOS the --dht (WIP) parameter needs to be set, furthermore,
-x DAOS_POOL needs to be given as well, with DAOS_POOL being an environment variable with a running pool name,  
for example:

```sh
export DAOS_POOL="<pool_name>"

mpirun -n 4 -x DAOS_POOL ./poet --dht ../bench/dolo_diffu_inner/dolo_diffu_inner.R output
```

## Differences to the main branch

- src/ChemistryModule/CMakeLists.txt was modified to compile DAOS code
- src/ChemistryModule/DaosKeyValue.c implements the underlying framework to connect and disconnect from the DAOS server, as well as call the read and write operations
- src/ChemistryModule/DHT_Wrapper.cpp was slightly modified to use DaosKeyValue.c instead of DHT.c

- include/poet/ was with the header file of DaosKeyValue extend and DHT_Wrapper.hpp was modified as well