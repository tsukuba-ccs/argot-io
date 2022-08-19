#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/times.h>

#include "run_param.h"
#include "fluid.h"
#include "source.h"
#include "prototype.h"
#include "io_ops.h"

#define FILENAME_LENGTH (256)

void make_directory(char*);

void output_mesh_single(struct fluid_mesh *mesh, struct run_param *this_run,
			char *filename)
{
  void *output_fp;

  io_ops->write_open(filename, &output_fp);
  io_ops->write(this_run, sizeof(*this_run), output_fp);

  int imesh;
  for(imesh=0;imesh<NMESH_X_LOCAL*NMESH_Y_LOCAL*NMESH_Z_LOCAL;imesh++) {
    struct fluid_mesh_io tmp_mesh;

    tmp_mesh.dens = mesh[imesh].dens;
    tmp_mesh.eneg = mesh[imesh].eneg;
    tmp_mesh.momx = mesh[imesh].momx;
    tmp_mesh.momy = mesh[imesh].momy;
    tmp_mesh.momz = mesh[imesh].momz;
    tmp_mesh.chem = mesh[imesh].chem;

    io_ops->write(&tmp_mesh, sizeof(tmp_mesh), output_fp);
  }
  io_ops->close(output_fp);
}

void output_mesh(struct fluid_mesh *mesh, struct run_param *this_run,
		 char *prefix)
{

  static char filename[FILENAME_LENGTH];

  sprintf(filename,"%s_%03d_%03d_%03d",prefix, 
	  this_run->rank_x, this_run->rank_y, this_run->rank_z);

  output_mesh_single(mesh, this_run, filename);
}


void output_src(struct radiation_src *src, struct run_param *this_run,
		char *prefix)
{
  void *output_fp;
  static char filename[FILENAME_LENGTH];

  sprintf(filename,"%s_src.dat",prefix);

  io_ops->write_open(filename, &output_fp);
  io_ops->write(&(this_run->nsrc), sizeof(uint64_t), output_fp);
  io_ops->write(&this_run->freq, sizeof(struct freq_param), output_fp);
  io_ops->write(src, sizeof(struct radiation_src) * this_run->nsrc, output_fp);
  io_ops->close(output_fp);
}

void output_data(struct fluid_mesh *mesh, 
		 struct radiation_src *src,
		 struct run_param *this_run, 
		 char *prefix)
{
  output_mesh(mesh, this_run, prefix);
  //  if(this_run->mpi_rank == 0) output_src(src, this_run, prefix);
}

void output_data_in_run(struct fluid_mesh *mesh,
			struct radiation_src *src,
			struct run_param *this_run,
			char *prefix)
{
  static char prefix_stamp[256];
  static char dirname[128];

  struct timeval start_tv, stop_tv;
  struct tms start_tms, stop_tms;

  int64_t mesh_data_size;

  int64_t data_size;
  float walltime;

  mesh_data_size = sizeof(struct fluid_mesh_io)*NMESH_LOCAL + sizeof(struct run_param);

  if(this_run->output_indx >= this_run->noutput) {
    return;
  }

  if(this_run->step % 20 == 0) {

    times(&start_tms);
    gettimeofday(&start_tv, NULL);

    make_directory("dmp");
    sprintf(prefix_stamp, "dmp/%s-dmp",prefix);
    output_data(mesh, src, this_run, prefix_stamp);

    times(&stop_tms);
    gettimeofday(&stop_tv, NULL);
    
    walltime = wallclock_timing(start_tv, stop_tv);
    MPI_Allreduce(MPI_IN_PLACE, &walltime, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    this_run->output_wt += walltime;

    data_size = mesh_data_size;
    MPI_Allreduce(MPI_IN_PLACE, &data_size, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    this_run->output_file_size += data_size;

    int io_check = verify_output(mesh, this_run, prefix_stamp);
    if(io_check != 0) {
      fprintf(this_run->proc_file,
	      "# Invalid IO detected.\n");
      fflush(this_run->proc_file);
    }    

    fprintf(this_run->proc_file, 
	    "# dumping files at t = %14.6e : %14.6e [sec] : %14.6e [GB/sec] \n",
	    this_run->tnow, walltime, mesh_data_size/walltime/1.0e9);
    fflush(this_run->proc_file);

  }

  if(this_run->tnow > this_run->output_timing[this_run->output_indx]){

    times(&start_tms);
    gettimeofday(&start_tv, NULL);

    sprintf(dirname,"%s-%02d",prefix, this_run->output_indx);
    make_directory(dirname);

    sprintf(prefix_stamp, "%s/%s-%02d",
            dirname,prefix,this_run->output_indx);

    output_data(mesh, src, this_run, prefix_stamp);
    this_run->output_indx++;

    times(&stop_tms);
    gettimeofday(&stop_tv, NULL);
    
    walltime = wallclock_timing(start_tv, stop_tv);
    MPI_Allreduce(MPI_IN_PLACE, &walltime, 1, MPI_FLOAT, MPI_MAX, MPI_COMM_WORLD);
    this_run->output_wt += walltime;

    data_size = mesh_data_size;
    MPI_Allreduce(MPI_IN_PLACE, &data_size, 1, MPI_LONG_LONG_INT, MPI_SUM, MPI_COMM_WORLD);
    this_run->output_file_size += data_size;

    int io_check = verify_output(mesh, this_run, prefix_stamp);
    if(io_check != 0) {
      fprintf(this_run->proc_file,
	      "# Invalid IO detected.\n");
    }else{
      fprintf(this_run->proc_file, "# IO verification passed\n");
    }

    fprintf(this_run->proc_file, 
	    "# dumping files at t = %14.6e : %14.6e [sec] : %14.6e [GB/sec] \n",
	    this_run->tnow, walltime, mesh_data_size/walltime/1.0e9);
    fflush(this_run->proc_file);
  }
}
