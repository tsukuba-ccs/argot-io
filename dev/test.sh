set -e
LANG=C

export CHFS_CHUNK_SIZE=$((512*1024))
export CHFS_BUF_SIZE=$((1024*1024))

eval $(chfsctl -h $HOME/hosts -L work/log start)
chlist

cd work
rm -rf dmp io_bench-* io_bench.diag

PROB="-x 64 -y 64 -z 32 -X 2 -Y 2 -Z 1"

#setup_IC $PROB io_bench
#mpirun -np 4 -hostfile $HOME/hosts -map-by node io_bench_argot $PROB io_bench-init/io_bench-init io_bench.pars

echo setup_IC
setup_IC -a chfs $PROB io_bench
echo io_bench_io
mpirun -x CHFS_SERVER -x CHFS_CHUNK_SIZE -x CHFS_BUF_SIZE \
	-np 4 -hostfile $HOME/hosts -map-by node io_bench_argot -a chfs $PROB \
	io_bench-init/io_bench-init io_bench.pars

echo chfsctl stop
chfsctl -h $HOME/hosts stop
echo see work/io_bench-out and remove work/\{dmp,io_bench-*,io_bench.diag,log}
echo OK
