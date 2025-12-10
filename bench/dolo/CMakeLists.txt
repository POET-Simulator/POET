set(bench_files
    dolo_interp.R
)

set(runtime_files
    dolo_interp_rt_dt2000.R
)

set(weight_files
    dolomite_trained.weights.h5)

ADD_BENCH_TARGET(
    dolo_bench
    bench_files 
    runtime_files
    weight_files
    "dolo"
)

add_dependencies(${BENCHTARGET} dolo_bench)