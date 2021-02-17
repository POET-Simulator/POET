# Output Files

POET will place all simulation data and other files inside the given
`<OUTPUT_DIRECTORY>`. The directory will look like this:

``` {.example}
.
└── <OUTPUT_DIRECTORY>/
    ├── iter000.rds
    ├── iter000.dht
    ├── ...
    ├── iter<n>.rds
    ├── iter<n>.dht
    ├── setup.rds
    └── timings.rds
```

## Description

All `.rds` file can be read inside a R runtime using e.g.
`readRDS("<FILE>")`. The following description can be given to the
files:

| File          | Description                                                                                                       |
|---------------|-------------------------------------------------------------------------------------------------------------------|
| iter<*n*>.rds | Defines the state of the grid after *n* iteration, especially the state after transport (`T`) and chemistry (`C`) |
| iter<*n*>.dht | DHT-snapshot of the *n* th iteration                                                                              |
| setup.rds     | Summary of all simulation parameters given at startup                                                             |
| timings.rds   | Various measured timings by POET                                                                                  |

## Timings

POET provides built-in time measurements of (sub) routines. The
corresponding values can be found in `<OUTPUT_DIRECTORY>/timings.rds`
and possible to read out within a R runtime with
`readRDS("timings.rds")`. There you will find the following values:

| Value              | Description                                                                |
|--------------------|----------------------------------------------------------------------------|
| simtime            | time spent in whole simulation loop without any initialization and cleanup |
| simtime\_transport | measured time in *transport* subroutine                                    |
| simtime\_chemistry | measured time in *chemistry* subroutine (actual parallelized part)         |

### chemistry subsetting

If running parallel there are also measured timings which are subsets of
*simtime\_​chemistry*.

| Value                      | Description                                        |
|----------------------------|----------------------------------------------------|
| simtime\_workers           | time spent in send/​recv loop of master             |
| simtime\_chemistry\_master | sequential part of master chemistry                |
| phreeqc                    | measured time of each worker in PHREEQC subroutine |

### DHT usage {#DHT-usage}

If running parallel and with DHT activated two more timings and also
some profiling about the DHT usage are given:

| Value           | Description                                             |
|-----------------|---------------------------------------------------------|
| dht\_fill\_time | time to write data to DHT                               |
| dht\_get\_time  | time to retreive data from DHT                          |
| dh\_hits        | count of data points retrieved from DHT                 |
| dht\_miss       | count of misses/count of data points written to DHT     |
| dht\_evictions  | count of data points evicted by another write operation |
