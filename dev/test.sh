set -e
LANG=C

export CHFS_CHUNK_SIZE=$((512*1024))
export CHFS_BUF_SIZE=$((1024*1024))

eval $(chfsctl -h $HOME/hosts -L work/log start)
chlist

cd work
rm -rf dmp io_bench-* io_bench.diag

#setup_IC io_bench
#mpirun -x LD_LIBRARY_PATH -np 8 -hostfile $HOME/hosts -map-by node io_bench_argot io_bench-init/io_bench-init io_bench.pars

setup_IC -a chfs io_bench
mpirun -x LD_LIBRARY_PATH -x CHFS_SERVER -x CHFS_CHUNK_SIZE -x CHFS_BUF_SIZE \
	-np 8 -hostfile $HOME/hosts -map-by node io_bench_argot -a chfs \
	io_bench-init/io_bench-init io_bench.pars

chfsctl -h $HOME/hosts stop
echo see work/io_bench-out and remove work/\{dmp,io_bench-*,io_bench.diag,log}
echo OK
