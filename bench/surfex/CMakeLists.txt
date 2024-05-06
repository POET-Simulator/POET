set(bench_files
    # surfex.R
    # ex.R
    PoetEGU_surfex_500.R
)

set(runtime_files
    # surfex_rt.R 
    # ex_rt.R
    PoetEGU_surfex_500_rt.R
)

ADD_BENCH_TARGET(
    surfex_bench
    bench_files 
    runtime_files
    "surfex"
)

add_dependencies(${BENCHTARGET} surfex_bench)