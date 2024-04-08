set(bench_files
    surfex.R
    ex.R
)

set(runtime_files
    surfex_rt.R 
    ex_rt.R
)

ADD_BENCH_TARGET(
    surfex_bench
    bench_files 
    runtime_files
    "surfex"
)

add_dependencies(${BENCHTARGET} surfex_bench)