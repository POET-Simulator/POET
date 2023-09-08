




taskset --cpu-list 0-23:2 mpirun --allow-run-as-root -x DAOS_POOL=test_pool -n 12 ./poet ../share/poet/bench/daos/dolo_diffu_edge.R ../../POETR >> ../../POET.out
#    taskset --cpu-list 0-23:2 mpirun --allow-run-as-root -x DAOS_POOL=test_pool -n $i ./build/src/kivibench-DAOSKV -x 10000 -y 3 -m 1000 -k 10 -v 12000 --csv >> benchmarks/clientbig$i.csv

