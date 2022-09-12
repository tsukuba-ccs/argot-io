#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "run_param.h"
#include "cross_section.h"

#include "prototype.h"
#include "io_ops.h"

#define FILENAME_LENGTH (256)
void make_directory(char*);

void init_run(struct run_param *this_run)
{
  static char dirname[FILENAME_LENGTH];
  static char proc_filename[FILENAME_LENGTH];
  struct io_ops *saved_io_ops = io_ops;

  sprintf(dirname,"%s-out",this_run->model_name);
  /* XXX - hack to use POSIX */
  init_io_ops("posix");
  make_directory(dirname);
  io_ops = saved_io_ops;
  
  sprintf(proc_filename, "%s-out/out_%03d_%03d_%03d", this_run->model_name,
	  this_run->rank_x, this_run->rank_y, this_run->rank_z);

  this_run->proc_file = fopen(proc_filename,"a");

  /* checking parameter configuration of this run */
  fprintf(this_run->proc_file,"\n\n\n\n#==================================================\n");
  fprintf(this_run->proc_file,"# model name : %s\n", this_run->model_name);
  fprintf(this_run->proc_file,"# total number of mesh along X-axis : %d\n", 
          this_run->nmesh_x_total);
  fprintf(this_run->proc_file,"# total number of mesh along Y-axis : %d\n", 
          this_run->nmesh_y_total);
  fprintf(this_run->proc_file,"# total number of mesh along Z-axis : %d\n", 
          this_run->nmesh_z_total);
  fprintf(this_run->proc_file,"# local number of mesh along X-axis : %d\n", 
          this_run->nmesh_x_local);
  fprintf(this_run->proc_file,"# local number of mesh along Y-axis : %d\n", 
          this_run->nmesh_y_local);
  fprintf(this_run->proc_file,"# local number of mesh along Z-axis : %d\n", 
          this_run->nmesh_z_local);
  fprintf(this_run->proc_file,"# number of chemical species : %d\n", 
          this_run->nspecies);


}
