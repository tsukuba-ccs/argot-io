set -e
LANG=C

cd work
setup_IC io_bench
mpirun -np 8 -hostfile $HOME/hosts -map-by node io_bench_argot io_bench-init/io_bench-init io_bench.pars

echo see work/io_bench-out and remove work/\{dmp,io_bench-*,io_bench.diag}
echo OK
