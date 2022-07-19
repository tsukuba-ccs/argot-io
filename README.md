# ARGOT-IO Benchmark

## About the benchmark

This benchmark is an I/O kernel of ARGOT space radiate transfer
application.

## How to install

1. Edit the problem size and the number of processes in
   src/run_param.h

* Problem size
```
#define NMESH_X_TOTAL (128)
#define NMESH_Y_TOTAL (128)
#define NMESH_Z_TOTAL (128)
```
* The number of processes
```
#define NNODE_X (2)
#define NNODE_Y (2)
#define NNODE_Z (2)
```

2. Compile the program by configure and make.
```
% autoreconf -i
% ./configure
% make
```

## How to execute

1. Create initial data
```
% setup_IC io_bench
```
The initial data will be generated in ./io_bench-init directory.

2. Execute io_bench_argot with io_bench-init/io_bench-init io_bench.pars
```
% mpiexec -n 8 io_bench_argot io_bench-init/io_bench-init io_bench.pars
```
The input file io_bench.pars is included in the work directory.

3. The performance number is the last line of io_bench-out/out_000_000_000
   file.
