# Create a list of files 
set(bench_files
    barite_200.R
)

set(runtime_files
    barite_200_rt.R
)

set(weight_files
    barite_trained.weights.h5)

# add_custom_target(barite_bench DEPENDS ${bench_files} ${runtime_files})

ADD_BENCH_TARGET(barite_bench 
    bench_files 
    runtime_files
    weight_files
    "barite"
)

add_dependencies(${BENCHTARGET} barite_bench)