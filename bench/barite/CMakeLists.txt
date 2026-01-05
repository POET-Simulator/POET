# Create a list of files
set(bench_files
    ${PROJECT_SOURCE_DIR}/ext/ai-benchmarks-utils/data/barite/abschluss/barite_200_abschluss.R
)

set(runtime_files
    ${PROJECT_SOURCE_DIR}/ext/ai-benchmarks-utils/data/barite/abschluss/runtime_barite_200_abschluss.R
    ${PROJECT_SOURCE_DIR}/ext/ai-benchmarks-utils/data/barite/abschluss/barite_trained.weights.h5
    ${PROJECT_SOURCE_DIR}/ext/ai-benchmarks-utils/data/barite/abschluss/barite_200ai_surrogate_input_script.R
)

message("Project source: ${PROJECT_SOURCE_DIR}")

# add_custom_target(barite_bench DEPENDS ${bench_files} ${runtime_files})

ADD_BENCH_TARGET(barite_bench
    bench_files
    runtime_files
    "barite"
)

add_dependencies(${BENCHTARGET} barite_bench)
