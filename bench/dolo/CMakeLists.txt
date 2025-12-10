set(bench_files
    dolo_interp.R
)

set(runtime_files
    dolo_interp_rt_dt2000.R
)

set(model_weights
    dolomite_trained.weights.h5)

ADD_BENCH_TARGET(
    dolo_bench
    bench_files 
    runtime_files
    model_weights
    "dolo"
)

add_dependencies(${BENCHTARGET} dolo_bench)