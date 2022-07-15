set -e

cd work
setup_IC io_bench
mpirun -np 8 -hostfile ../hosts -map-by node io_bench_argot io_bench-init/io_bench-init io_bench.pars
