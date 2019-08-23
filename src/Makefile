CC=mpicc
CFLAGS=-O3 -funroll-loops -std=c99 -mcmodel=medium 

IO_BENCH_ARGOT_OBJ = init_mpi.o init_run.o input_data.o input_params.o io_bench_argot.o cosmology.o mpi_rank.o make_directory.o output_data.o timing.o

SETUP_IC_OBJ = uniform_medium_IC.o output_data.o mpi_rank.o make_directory.o source.o timing.o

IO_BENCH_ARGOT_DEP := $(IO_BENCH_ARGOT_OBJ)
SETUP_IC_DEP := $(SETUP_IC_OBJ)

all: io_bench_argot setup_IC

io_bench_argot : $(IO_BENCH_ARGOT_DEP)
	$(CC) $(CFLAGS) -o $@ $(IO_BENCH_ARGOT_OBJ) $(CLIBS) -lm

setup_IC : $(SETUP_IC_OBJ)
	$(CC) $(CFLAGS) -o $@ $(SETUP_IC_OBJ) $(CLIBS) -lm

clean: 
	-rm -rf io_bench_argot setup_IC
	-rm -rf *.o

distclean: clean
	-rm -rf *~