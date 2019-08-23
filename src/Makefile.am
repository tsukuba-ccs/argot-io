bin_PROGRAMS = io_bench_argot setup_IC
io_bench_argot_SOURCES = init_mpi.c init_run.c input_data.c \
	input_params.c io_bench_argot.c cosmology.c \
	mpi_rank.c make_directory.c output_data.c timing.c
setup_IC_SOURCES = uniform_medium_IC.c output_data.c \
	mpi_rank.c make_directory.c source.c timing.c
AM_CFLAGS = -O3 -funroll-loops -std=c99 -mcmodel=medium
ACLOCAL_AMFLAGS = -I m4

CC = $(MPICC)
