bin_PROGRAMS = io_bench_argot setup_IC argot-io
io_bench_argot_SOURCES = init_mpi.c init_run.c input_data.c \
	input_params.c io_bench_argot.c cosmology.c \
	mpi_rank.c make_directory.c output_data.c timing.c \
	pad_region.c set_fluid_pad.c fluid_hllc.c calc_timestep.c \
	io_posix.c verify_output.c \
	run_param.h cosmology.h source.h fluid.h chemistry.h cross_section.h \
	mpi_param.h prototype.h constants.h radiation.h \
	io_ops.h
setup_IC_SOURCES = uniform_medium_IC.c output_data.c input_data.c \
	mpi_rank.c make_directory.c source.c timing.c \
	io_posix.c verify_output.c \
	run_param.h cosmology.h source.h fluid.h chemistry.h cross_section.h \
	mpi_param.h prototype.h constants.h radiation.h \
	io_ops.h
argot_io_SOURCES = argot-io.c io_posix.c io_ops.h
AM_CFLAGS = -O3 -funroll-loops -std=c99 -mcmodel=medium -D__SECOND_ORDER_RUNGE_KUTTA__
if GFARM
io_bench_argot_SOURCES += io_gfarm.c
setup_IC_SOURCES += io_gfarm.c
argot_io_SOURCES += io_gfarm.c
endif
if CHFS
io_bench_argot_SOURCES += io_chfs.c
setup_IC_SOURCES += io_chfs.c
argot_io_SOURCES += io_chfs.c
AM_CFLAGS += @CHFS_CFLAGS@
AM_LDFLAGS = @CHFS_LIBS@
endif
ACLOCAL_AMFLAGS = -I m4

CC = $(MPICC)
