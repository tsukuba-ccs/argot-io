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

int verify_mesh_data(struct fluid_mesh *mesh1, struct fluid_mesh *mesh2)
{
  int result = 0;

#pragma omp parallel for reduce(result:+)
  for(int im=0;im<NMESH_LOCAL;im++) {
    if(mesh1[im].dens != mesh2[im].dens) result += 1;
    if(mesh1[im].eneg != mesh2[im].eneg) result += 1;
    if(mesh1[im].momx != mesh2[im].momx) result += 1;
    if(mesh1[im].momy != mesh2[im].momy) result += 1;
    if(mesh1[im].momz != mesh2[im].momz) result += 1;
    if(mesh1[im].uene != mesh2[im].uene) result += 1;
#if 1
    if(mesh1[im].chem.wmol  != mesh2[im].chem.wmol ) result += 1;
    if(mesh1[im].chem.felec != mesh2[im].chem.felec) result += 1;
    if(mesh1[im].chem.fHI   != mesh2[im].chem.fHI  ) result += 1;
    if(mesh1[im].chem.fHII  != mesh2[im].chem.fHII ) result += 1;
#ifdef __HYDROGEN_MOL__
    if(mesh1[im].chem.fH2I   != mesh2[im].chem.fH2I ) result += 1;
    if(mesh1[im].chem.fH2II  != mesh2[im].chem.fH2II) result += 1;
    if(mesh1[im].chem.fHM    != mesh2[im].chem.fHM  ) result += 1;
#endif
#ifdef __HELIUM__
    if(mesh1[im].chem.fHeI   != mesh2[im].chem.fHeI  ) result += 1;
    if(mesh1[im].chem.fHeII  != mesh2[im].chem.fHeII ) result += 1;
    if(mesh1[im].chem.fHeIII != mesh2[im].chem.fHeIII) result += 1;
#endif
#endif
  }

  return result;
}

int verify_output(struct fluid_mesh *mesh, struct run_param *this_run,
		  char *prefix)
{
  struct fluid_mesh *mesh_io;
  struct radiation_src *src_io;
  struct run_param run_io;

  mesh_io = (struct fluid_mesh *)malloc(sizeof(struct fluid_mesh)*NMESH_LOCAL);
  src_io = (struct radiation_src *)malloc(sizeof(struct radiation_src)*NSOURCE_MAX);

  run_io.rank_x = this_run->rank_x;
  run_io.rank_y = this_run->rank_y;
  run_io.rank_z = this_run->rank_z;
  run_io.mpi_rank = this_run->mpi_rank;
  //run_io.nchannel = this_run->nchannel;
  
  input_data(mesh_io, src_io, &run_io, prefix);

  int result = verify_mesh_data(mesh, mesh_io);

  return result;
}
