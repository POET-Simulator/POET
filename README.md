
<!--
    Time-stamp: "Last modified 2021-02-08 13:46:00 mluebke"
-->

# POET

POET is a coupled reactive transport simulator implementing a parallel
architecture and a fast distributed hash table.

## External Libraries

We use external libraries:

- **argh** - https://github.com/adishavit/argh (BSD license)

## Installation

### Requirements

To compile POET you need several software to be installed:

- C/C++ compiler (tested with GCC)
- MPI-Implementation (tested with OpenMPI and MVAPICH)
- R language and environment 
- CMake 3.9+

The following R libraries must then be installed:

- [devtools](https://www.r-project.org/nosvn/pandoc/devtools.html)
- [Rcpp](https://cran.r-project.org/web/packages/Rcpp/index.html)
- [RInside](https://cran.r-project.org/web/packages/RInside/index.html)
- [RedModRphree](https://gitext.gfz-potsdam.de/delucia/RedModRphree)
- [Rmufits](https://gitext.gfz-potsdam.de/delucia/Rmufits)

### Compiling source code

The generation of makefiles is done with CMake. So, running

```
cmake . -B build
```

will create the directory `build`. `cd` into it and run `make` to start build
process.

If everything went well you'll find the executable at `build/src/poet`. 

During the generation of Makefiles, various options can be specified via `cmake
-D <option>=<value> [...]`. Currently there are the following available options:

- **DHT_Debug**=*boolean* - toggles the output of detailed statistics about DHT
  usage (`cmake -D DHT_Debug=ON`). Defaults to *OFF*.
  
### Example: Build from scratch

Assuming that only the C/C++ compiler, MPI libraries, R runtime environment and
CMake have been installed, POET can be installed as follows:

``` sh
# start R environment
$ R

# install R dependencies
> install.packages(c("devtools", "Rcpp", "RInside"))
> devtools::install_git("https://gitext.gfz-potsdam.de/delucia/RedModRphree.git")
> devtools::install_git("https://gitext.gfz-potsdam.de/delucia/Rmufits.git")
> q(save="no")

# cd into POET project root
$ cd <POET_dir>

# Build process
$ mkdir build && cd build
$ cmake -DCMAKE_INSTALL_PREFIX=/home/<user>/poet ..
$ make -j<max_numprocs>
$ make install
```

This will install a POET project structure into `/home/<user>/poet` which is
called hereinafter `<POET_INSTALL_DIR>`. With this version of POET I would not
recommend to install to hierarchies like `/usr/local/` etc.

The correspondending directory tree would be look like this:

``` sh
.
└── poet/
    ├── bin/
    │   └── poet
    ├── data/
    │   └── SimDol2D.R
    └── R_lib/
        ├── kin_r_library.R
        └── parallel_r_library.R
```

The R libraries will be read in during runtime and the paths are hardcoded
absolute paths inside `poet.cpp`. So, if you consider to move `bin/poet` either
change paths of the R source files and recompile POET or also move `R_lib/*`
according to the binary.

## Running

Run POET by `mpirun ./poet <OPTIONS> <SIMFILE> <OUTPUT_DIRECTORY>` where:

- **OPTIONS** - runtime parameters (explained below)
- **SIMFILE** - simulation described as R script (currently supported:
  `<POET_INSTALL_DIR>/data/SimDol2D.R`)
- **OUTPUT_DIRECTORY** - path, where all output of POET should be stored

### Runtime options

The following parameters can be set:

| Option                  | Value        | Description                                                     |
|-------------------------|--------------|-----------------------------------------------------------------|
| **-work-package-size=** | *1..n*       | size of work packages (defaults to *5*)                         |
| **-ignore-result**      |              | disables store of simulation resuls                             |
| **-dht**                |              | enabling DHT usage (defaults to *OFF*)                          |
| **-dht-nolog**          |              | disabling applying of logarithm before rounding                 |
| **-dht-signif=**        | *1..n*       | set rounding to number of significant digits (defaults to  *5*) |
| **-dht-strategy=**      | *0-1*        | change DHT strategy. **NOT IMPLEMENTED YET** (Defaults to *0*)  |
| **-dht-size=**          | *1-n*        | size of DHT per process involved in byte (defaults to *1 GiB*)  |
| **-dht-snaps=**         | *0-2*        | disable or enable storage of DHT snapshots                      |
| **-dht-file=**          | `<SNAPSHOT>` | initializes DHT with the given snapshot file                    |

#### Additions to `dht-signif`

Only used if no vector is given in setup file. For individual valuies per column
use R vector `signif_vector` in `SIMFILE`.

#### Additions to `dht-snaps`
Following values can be set:
    - *0* = snapshots are disabled
    - *1* = only stores snapshot at the end of the simulation with name
      `<OUTPUT_DIRECTORY>.dht`
    - *2* = stores snapshot at the end and after each iteration iteration
      snapshot files are stored in <DIRECTORY>/iter<n>.dht
    

### Example: Running from scratch


We will continue the above example and start a simulation with `SimDol2D.R`,
which is the only simulation supported at this moment. To start the simulation
with 4 processes `cd` into your previously installed POET-dir
`<POET_INSTALL_DIR>/bin` and run:

``` sh
mpirun -n 4 ./poet ../data/SimDol2D.R output
```

After a finished simulation all data generated by POET will be found in the
directory `output`.

You might want to use the DHT to cache previously simulated data points and
reuse them in further time-steps. Just append `-dht` to the options of POET to
activate the usage of the DHT. The resulting call would look like this:

``` sh
mpirun -n 4 ./poet -dht SimDol2D.R output
```

## Examples included (more to come)

- **SimDol2D.R** - simple chemistry (Calcite/Dolomite) on a 50x50 2D grid, 20
time steps 2)
- ~~ **SimDolKtz.R** - simple chemistry (Calcite/Dolomite) on Ketzin grid (~650k
elements), 20 time steps The flow snapshots are **NOT INCLUDED** in project
directory but must be provided separately.~~

## About the usage of MPI_Wtime()
Implemented time measurement functions use `MPI_Wtime()`. Some important
informations from the OpenMPI Man Page:

For example, on platforms that support it, the clock_gettime() function will be
used to obtain a monotonic clock value with whatever precision is supported on
that platform (e.g., nanoseconds).
