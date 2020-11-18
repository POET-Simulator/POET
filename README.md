
<!--
    Time-stamp: "Last modified 2020-02-01 18:14:13 delucia"
-->

# install libraries from MDL

    library(devtools)
    devtools::install_git("https://gitext.gfz-potsdam.de/delucia/RedModRphree.git")
    devtools::install_git("https://gitext.gfz-potsdam.de/delucia/Rmufits.git")
  
# USAGE

mpirun ./kin <OPTIONS> <simfile.R> <DIRECTORY>

OPTIONS:

--work-package-size=<1-n>     ... size of work packages (default 5)

--ignore-result         ... disables store of simulation resuls

DHT:

--dht                   ... enable dht (default is off)

--dht-log               ... enable logarithm application before rounding (default is off)

--dht-signif=<1-n>      ... set rounding to number of significant digits (default 5) 
                            (only used if no vector is given in setup file)
                            (for individual values per column use R vector "signif_vector" in setup file)

--dht-strategy=<0-1>    ... change dht strategy, not implemented yet (default 0, dht on workers)

--dht-size=<1-n>        ... size of dht per process involved (see dht-strategy) in byte (default 1GiB)

--dht-snaps=<0-2>	... enable or disable storage of DHT snapshots
			    0 = snapshots are disabled
			    1 = only stores snapshot at the end of the simulation with name <DIRECTORY>.dht
			    2 = stores snapshot at the end and after each iteration
			        iteration snapshot files are stored in <DIRECTORY>/iter<n>.dht

--dht-file=<snapshot> 	... initializes DHT with the given snapshot file

###############################################################################


# about the usage of MPI_Wtime()
From the OpenMPI Man Page:

For example, on platforms that support it, the clock_gettime() function will be used 
to obtain a monotonic clock value with whatever precision is supported on that platform (e.g., nanoseconds). 

# External Libraries
Cmdline Parsing -> https://github.com/adishavit/argh


# Examples included (more to come)
1) SimDol2D.R     ... simple chemistry (Calcite/Dolomite) on a 50x50 2D grid, 20 time steps
2) SimDolKtz.R    ... simple chemistry (Calcite/Dolomite) on Ketzin grid (~650k elements), 20 time steps
                      The flow snapshots are NOT INCLUDED in svn but must be provided separately



