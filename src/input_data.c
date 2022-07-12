#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <sys/times.h>
#include <inttypes.h>

#include "run_param.h"
#include "fluid.h"
#include "source.h"
#include "prototype.h"
#include "io_ops.h"

#define FILENAME_LENGTH (256)

void setup_data(struct run_param *this_run)
{
  float dx, dy, dz;

  this_run->xmin = this_run->ymin = this_run->zmin = 0.0;
  this_run->xmax = this_run->ymax = this_run->zmax = 1.0;

  dx = (this_run->xmax-this_run->xmin)/(float)NNODE_X;
  dy = (this_run->ymax-this_run->ymin)/(float)NNODE_Y;
  dz = (this_run->zmax-this_run->zmin)/(float)NNODE_Z;

  this_run->xmin_local = this_run->rank_x*dx;
  this_run->xmax_local = this_run->xmin_local+dx;
  
  this_run->ymin_local = this_run->rank_y*dy;
  this_run->ymax_local = this_run->ymin_local+dy;

  this_run->xmin_local = this_run->rank_x*dx;
  this_run->xmax_local = this_run->xmin_local+dx;
  
  this_run->ymin_local = this_run->rank_y*dy;
  this_run->ymax_local = this_run->ymin_local+dy;

  this_run->nmesh_x_total = this_run->nmesh_y_total = this_run->nmesh_z_total = 128;
  this_run->nmesh_x_local = this_run->nmesh_y_local = this_run->nmesh_z_local = 32;

  this_run->delta_x = (this_run->xmax-this_run->xmin)/(float)this_run->nmesh_x_total;
  this_run->delta_y = (this_run->ymax-this_run->ymin)/(float)this_run->nmesh_y_total;
  this_run->delta_z = (this_run->zmax-this_run->zmin)/(float)this_run->nmesh_z_total;
}

void input_mesh_header_prefix(struct run_param *this_run, char *prefix)
{
  static char filename[FILENAME_LENGTH];

  sprintf(filename, "%s_%03d_%03d_%03d", prefix, 
	  this_run->rank_x, this_run->rank_y, this_run->rank_z);

  input_mesh_header(this_run, filename);
}

void input_mesh_header(struct run_param *this_run, char *filename)
{
  struct run_param input_run;
  void *input_fp;
  
  io_ops->read_open(filename, &input_fp);
  io_ops->read(&input_run, sizeof(struct run_param), input_fp);
  io_ops->close(input_fp);

  if(input_run.nmesh_x_total != NMESH_X_TOTAL ||
     input_run.nmesh_y_total != NMESH_Y_TOTAL ||
     input_run.nmesh_z_total != NMESH_Z_TOTAL) {
    fprintf(stderr, "# Inconsistent size of the global mesh\n");
    fprintf(stderr, "# input_run.nmesh_x_total = %d\n", input_run.nmesh_x_total);
    fprintf(stderr, "# input_run.nmesh_y_total = %d\n", input_run.nmesh_y_total);
    fprintf(stderr, "# input_run.nmesh_z_total = %d\n", input_run.nmesh_z_total);
    exit(EXIT_FAILURE);
  }

  if(input_run.nmesh_x_local != NMESH_X_LOCAL ||
     input_run.nmesh_y_local != NMESH_Y_LOCAL ||
     input_run.nmesh_z_local != NMESH_Z_LOCAL) {
    fprintf(stderr, "# Inconsistent size of the local mesh\n");
    fprintf(stderr, "# input_run.nmesh_x_local = %d\n", input_run.nmesh_x_local);
    fprintf(stderr, "# input_run.nmesh_y_local = %d\n", input_run.nmesh_y_local);
    fprintf(stderr, "# input_run.nmesh_z_local = %d\n", input_run.nmesh_z_local);
    exit(EXIT_FAILURE);
  }

  if(input_run.mpi_nproc != NNODE_X*NNODE_Y*NNODE_Z) {
    fprintf(stderr, "# Inconsistent number of MPI processes\n");
    fprintf(stderr, "# input_run.mpi_nproc = %d\n", input_run.mpi_nproc);
    fprintf(stderr, "# this_run.mpi_nproc = %d\n", NNODE_X*NNODE_Y*NNODE_Z);
    exit(EXIT_FAILURE);
  }

  if(input_run.nnode_x != NNODE_X ||
     input_run.nnode_y != NNODE_Y ||
     input_run.nnode_z != NNODE_Z) {
    fprintf(stderr, "# Inconsistent domain decomposition\n");
    fprintf(stderr, "# input_run.nnode_x = %d\n",input_run.nnode_x);
    fprintf(stderr, "# input_run.nnode_y = %d\n",input_run.nnode_y);
    fprintf(stderr, "# input_run.nnode_z = %d\n",input_run.nnode_z);
    exit(EXIT_FAILURE);
  }

  if(input_run.nspecies != NSPECIES) {
    fprintf(stderr, "# Inconsistent number of chemical species \n");
    fprintf(stderr, "# input_run.nspecies = %d\n", input_run.nspecies);
    exit(EXIT_FAILURE);
  }

  if(input_run.nchannel != NCHANNEL) {
    fprintf(stderr, "# Inconsistent channel number of radiation transfer \n");
    fprintf(stderr, "# input_run.nchannel = %d \n", input_run.nchannel);
    exit(EXIT_FAILURE);
  }

  io_ops->read_open(filename, &input_fp);
  io_ops->read(this_run, sizeof(struct run_param), input_fp);
  io_ops->close(input_fp);
}

void input_mesh_single(struct fluid_mesh *mesh, struct run_param *this_run,
		       char *filename)
{
  void *input_fp;

  struct timeval start_tv, stop_tv;
  struct tms start_tms, stop_tms;

  uint64_t mesh_data_size;
  float walltime;

  mesh_data_size = sizeof(struct fluid_mesh_io)*NMESH_LOCAL + 
    sizeof(struct run_param);

  input_mesh_header(this_run,filename);
  
  io_ops->read_open(filename, &input_fp);

  times(&start_tms);
  gettimeofday(&start_tv, NULL);
  
  io_ops->read(this_run, sizeof(struct run_param), input_fp);
#if 0
  io_ops->read(mesh, sizeof(struct fluid_mesh) *
	NMESH_X_LOCAL*NMESH_Y_LOCAL*NMESH_Z_LOCAL,
	input_fp);
#else
  int imesh;
  for(imesh=0;imesh<NMESH_X_LOCAL*NMESH_Y_LOCAL*NMESH_Z_LOCAL;imesh++) {
    struct fluid_mesh_io temp_mesh;
    io_ops->read(&temp_mesh, sizeof(struct fluid_mesh_io), input_fp);
    mesh[imesh].dens = temp_mesh.dens;
    mesh[imesh].eneg = temp_mesh.eneg;
    mesh[imesh].momx = temp_mesh.momx;
    mesh[imesh].momy = temp_mesh.momy;
    mesh[imesh].momz = temp_mesh.momz;
    mesh[imesh].chem = temp_mesh.chem;
    mesh[imesh].uene = (mesh[imesh].eneg - 
			0.5*(SQR(mesh[imesh].momx)+
			     SQR(mesh[imesh].momy)+
			     SQR(mesh[imesh].momz))/mesh[imesh].dens)/mesh[imesh].dens;

    mesh[imesh].prev_chem = mesh[imesh].chem;
    mesh[imesh].prev_uene = mesh[imesh].uene;
  }
#endif

  times(&stop_tms);
  gettimeofday(&stop_tv, NULL);
  
  walltime = wallclock_timing(start_tv, stop_tv);
  
  this_run->input_mesh_wt = walltime;
  this_run->input_mesh_tp = mesh_data_size/walltime/1.0e9;

  io_ops->close(input_fp);
}

void input_mesh(struct fluid_mesh *mesh, struct run_param *this_run,
		char *prefix)
{

  static char filename[FILENAME_LENGTH];

  sprintf(filename, "%s_%03d_%03d_%03d", prefix, 
	  this_run->rank_x, this_run->rank_y, this_run->rank_z);

  input_mesh_single(mesh, this_run, filename);

}

void input_src(struct radiation_src *src, struct run_param *this_run, 
	       char *prefix)
{
  static char filename[FILENAME_LENGTH];
  void *input_fp;

  struct timeval start_tv, stop_tv;
  struct tms start_tms, stop_tms;

  uint64_t source_data_size;
  float walltime;

  uint64_t input_nsrc;

  source_data_size = sizeof(struct radiation_src)*this_run->nsrc
    +sizeof(struct freq_param);

  sprintf(filename,"%s_src.dat", prefix);

  io_ops->read_open(filename, &input_fp);

  times(&start_tms);
  gettimeofday(&start_tv, NULL);
  
  if(input_fp == NULL) {
    fprintf(stderr, "# File %s not found\n", filename);
    exit(EXIT_FAILURE);
  }
  
  io_ops->read(&(this_run->nsrc), sizeof(uint64_t), input_fp);

  if(this_run->nsrc > NSOURCE_MAX) {
    fprintf(stderr, "# Exceeds the max. number of the sources\n");
    fprintf(stderr, "# input_nsrc = %" PRIu64 "\n", this_run->nsrc);
    exit(EXIT_FAILURE);
  }

  io_ops->read(&this_run->freq, sizeof(struct freq_param), input_fp);
  io_ops->read(src, sizeof(struct radiation_src) * this_run->nsrc, input_fp);

  times(&stop_tms);
  gettimeofday(&stop_tv, NULL);
    
  walltime = wallclock_timing(start_tv, stop_tv);

  this_run->input_src_wt = walltime;
  this_run->input_src_tp = source_data_size/walltime/1.0e9;

  io_ops->close(input_fp);
}

void input_src_file(struct radiation_src *src, struct run_param *this_run, 
                    char *src_file)
{
  void *input_fp;

  uint64_t input_nsrc;

  io_ops->read_open(src_file, &input_fp);
  
  if(input_fp == NULL) {
    fprintf(stderr, "# File %s not found\n", src_file);
    exit(EXIT_FAILURE);
  }
  
  io_ops->read(&(this_run->nsrc), sizeof(uint64_t), input_fp);

  if(this_run->nsrc > NSOURCE_MAX) {
    fprintf(stderr, "# Exceeds the max. number of the sources\n");
    fprintf(stderr, "# input_nsrc = %" PRIu64 "\n", this_run->nsrc);
    exit(EXIT_FAILURE);
  }

  io_ops->read(&this_run->freq, sizeof(struct freq_param), input_fp);
  io_ops->read(src, sizeof(struct radiation_src) * this_run->nsrc, input_fp);

  io_ops->close(input_fp);
}


void input_data(struct fluid_mesh *mesh, struct radiation_src *src, 
		struct run_param *this_run, char *prefix)
{
  input_mesh(mesh, this_run, prefix);
  //  input_src(src, this_run, prefix);
}
