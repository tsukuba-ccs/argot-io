# ARGOT-IO Benchmark

## About the benchmark

This benchmark is an I/O kernel of ARGOT space radiate transfer
application.

## How to install

1. Compile the program by configure and make.
```
% autoreconf -i
% ./configure
% make
```

## How to execute

1. Create initial data
```
% setup_IC -x 128 -y 128 -z 128 -X 2 -Y 2 -Z 2 io_bench
```
-x, -y, and -z options specify the problem size, and -X, -Y, and, -Z
options specify the number of processes to execute io_bench_argot.

The initial data will be generated in ./io_bench-init directory.

1. Execute io_bench_argot with io_bench-init/io_bench-init io_bench.pars
```
% mpiexec -n 8 io_bench_argot -x 128 -y 128 -z 128 -X 2 -Y 2 -Z 2 io_bench-init/io_bench-init io_bench.pars
```
The input file io_bench.pars is included in the work directory.

1. The performance number is the last line of io_bench-out/out_000_000_000
   file.
