#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mpi.h>
#include <unistd.h>
#include <getopt.h>

#include "constants.h"
#include "run_param.h"
#include "mpi_param.h"
#include "radiation.h"
#include "source.h"

#include "prototype.h"
#include "io_ops.h"

int main(int argc, char **argv)
{
  
  static struct run_param this_run;
  static struct mpi_param this_mpi;
  
  struct radiation_src *src;
  struct fluid_mesh *mesh;
  char *api = "posix";
  int c;

  while ((c = getopt(argc, argv, "a:")) != -1) {
    switch (c) {
    case 'a':
      api = optarg;
      break;
    }
  }
  argc -= optind;
  argv += optind;

  init_io_ops(api);
  io_ops->init();

  src = (struct radiation_src *)malloc(sizeof(struct radiation_src)*NSOURCE_MAX);
  mesh = (struct fluid_mesh *)malloc(sizeof(struct fluid_mesh)*NMESH_LOCAL);
  
  MPI_Init(&argc, &argv);

  init_mpi(&this_run, &this_mpi);
  input_data(mesh, src, &this_run, argv[0]);
  input_params(&this_run, argv[1]);
  init_run(&this_run);

  fprintf(this_run.proc_file, 
	  "# reading mesh files : %14.6e [sec] : %14.6e [GB/sec])\n",
	  this_run.input_mesh_wt, this_run.input_mesh_tp);
  //  fprintf(this_run.proc_file, 
  //	  "# reading a source file : %14.6e [sec] : %14.6e [GB/sec])\n",
  //	  this_run.input_src_wt, this_run.input_src_tp);
	  

  float tnow, dtime;

  this_run.tnow = 0.0;

  this_run.output_wt = 0.0;
  this_run.output_file_size = 0;
  
  dtime = 0.01;

  while(this_run.tnow < 0.5) {

    sleep(1);

    this_run.tnow += dtime;
    MPI_Barrier(MPI_COMM_WORLD);
    output_data_in_run(mesh, src, &this_run, this_run.model_name);
  }

  fprintf(this_run.proc_file,
	  "# output throughput : %14.6e [GB/sec]\n",
	  (float)this_run.output_file_size/this_run.output_wt/1.0e9);

  MPI_Finalize();
  io_ops->term();
}
