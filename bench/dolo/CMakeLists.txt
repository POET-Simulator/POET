set(bench_files
    ${PROJECT_SOURCE_DIR}/ext/ai-benchmarks-utils/data/dolomite/abschluss/
)

set(runtime_files
    ${PROJECT_SOURCE_DIR}/ext/ai-benchmarks-utils/data/dolomite/abschluss/
)

ADD_BENCH_TARGET(
    dolo_bench
    bench_files
    runtime_files
    "dolo"
)

add_dependencies(${BENCHTARGET} dolo_bench)
