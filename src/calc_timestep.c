#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <mpi.h>

#include "constants.h"
#include "run_param.h"
#include "fluid.h"

#include "prototype.h"

#define MESH(ix,iy,iz) (mesh[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))])

float calc_timestep_fluid(struct fluid_mesh *mesh, struct run_param *this_run)
{
  float dtime_min;

  dtime_min = FLT_MAX;

#pragma omp parallel for schedule(auto) reduction(min:dtime_min)
  for(int ix=0;ix<NMESH_X_LOCAL;ix++) {
    for(int iy=0;iy<NMESH_Y_LOCAL;iy++){
      for(int iz=0;iz<NMESH_Z_LOCAL;iz++){
        struct fluid_mesh *target_mesh;
        float cs, gamma, velx, vely, velz;
        float dtime_x, dtime_y, dtime_z, dtime;

        target_mesh = &MESH(ix,iy,iz);

        gamma = gamma_total(target_mesh, this_run);
        cs = sqrtf((gamma-1.0)*gamma*target_mesh->uene);
        velx = fabsf(target_mesh->momx/target_mesh->dens);
        vely = fabsf(target_mesh->momy/target_mesh->dens);
        velz = fabsf(target_mesh->momz/target_mesh->dens);

        dtime_x = this_run->delta_x/(velx+cs);
        dtime_y = this_run->delta_y/(vely+cs);
        dtime_z = this_run->delta_z/(velz+cs);

        dtime = fminf(dtime_x, fminf(dtime_y, dtime_z));
        dtime_min = fminf(dtime, dtime_min);
      }
    }
  }

  #ifndef __SERIAL__
  float dtime_min_local;
  dtime_min_local = dtime_min;

  MPI_Allreduce(&dtime_min_local, &dtime_min, 1, MPI_FLOAT, MPI_MIN, MPI_COMM_WORLD);
#endif

  return (COURANT_FACT*dtime_min);

}
