set(bench_files
    barite_fgcs_4.R
    dolo_fgcs_3.R
)

set(runtime_files
    barite_fgcs_4_rt.R
    dolo_fgcs_3_rt.R
)

ADD_BENCH_TARGET(
    fgcs_bench
    bench_files 
    runtime_files
    "fgcs"
)

add_dependencies(${BENCHTARGET} fgcs_bench)