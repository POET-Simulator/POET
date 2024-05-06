set(bench_files
    dolo_inner_large.R
    dolo_interp.R
)

set(runtime_files
    dolo_inner_large_rt.R
    dolo_interp_rt.R
)

ADD_BENCH_TARGET(
    dolo_bench
    bench_files 
    runtime_files
    "dolo"
)

add_dependencies(${BENCHTARGET} dolo_bench)