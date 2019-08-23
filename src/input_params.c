#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <mpi.h>

#define MAX_MODEL_NAME_LEN (128)

#include "run_param.h"

void input_params(struct run_param *this_run, char *param_filename)
{
  FILE *param_file;

  param_file = fopen(param_filename, "r");

  if(param_file == NULL) {
    if(this_run->mpi_rank == 0) 
      fprintf(stderr, "File %s not found in input_params.\n", param_filename);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* set the model name */
  this_run->model_name = (char *) malloc(sizeof(char)*MAX_MODEL_NAME_LEN);
  fscanf(param_file,"%s",this_run->model_name);

  static char diag_filename[256];
  sprintf(diag_filename,"%s.diag", this_run->model_name);

  /* Open the diag_file */
  if(this_run->mpi_rank == 0) this_run->diag_file = fopen(diag_filename, "a");

  /* set nmesh_per_loop */
  int ngrp;
  fscanf(param_file, "%d", &ngrp);
  this_run->nmesh_per_loop = NMESH_LOCAL/ngrp;
  if((NMESH_LOCAL % this_run->nmesh_per_loop) != 0) {
    fprintf(this_run->proc_file, 
	    "NMESH_LOCAL should be divisible by this_run->nmesh_per_loop\n.");
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* # of output timing */
  fscanf(param_file,"%d",&(this_run->noutput));

  this_run->output_timing = (float *) malloc(sizeof(float)*this_run->noutput);

  this_run->output_indx = 0;

  int iout;
  float prev_output_timing;

  for(iout=0;iout<this_run->noutput;iout++) {

    fscanf(param_file,"%f",&(this_run->output_timing[iout]));

#ifdef __COSMOLOGICAL__
    if(iout>0 && this_run->output_timing[iout] > prev_output_timing) {
      fprintf(this_run->err_file,
              "Output timing must be in reducing order.\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if(this_run->output_timing[iout] > this_run->znow) {
#else
    if(iout>0 && this_run->output_timing[iout] < prev_output_timing) {
      fprintf(this_run->err_file,
	      "Output timing must be in increasing order.\n");
      MPI_Finalize();
      exit(EXIT_FAILURE);
    }
    if(this_run->output_timing[iout] < this_run->tnow) {
#endif
      this_run->output_indx = iout+1;
    }
    prev_output_timing = this_run->output_timing[iout];
  }

#ifdef __COSMOLOGICAL__
  this_run->tend = ztotime(this_run->output_timing[this_run->noutput-1], 
			   this_run->cosm);
#else
  this_run->tend = this_run->output_timing[this_run->noutput-1];
#endif
      
  if(this_run->output_indx == this_run->noutput || 
     this_run->tnow > this_run->tend) {
    if(this_run->mpi_rank == 0){
      fprintf(stderr,"This run already finished.\n");
      fprintf(stderr,"No outputs requested.\n");
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

}
