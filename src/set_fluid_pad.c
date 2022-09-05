#include <stdio.h>
#include <stdlib.h>

#include "run_param.h"
#include "fluid.h"

#define MESH(ix,iy,iz)   (mesh[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))])

void set_pad_x(struct pad_region *pad, struct fluid_mesh *mesh,
	       struct fluid_mesh *mesh_im2, struct fluid_mesh *mesh_im1,
	       struct fluid_mesh *mesh_ip1, struct fluid_mesh *mesh_ip2,
	       int ix, int iy, int iz,
	       struct run_param *this_run)
{
  struct fluid_mesh reflect_lo1, reflect_lo2;
  struct fluid_mesh reflect_hi1, reflect_hi2;
  int ixm2, ixm1, ixp1, ixp2;
  
  ixm2 = ix-2;  ixm1 = ix-1;
  ixp1 = ix+1;  ixp2 = ix+2;

  if(ixm2 < 0) {
#ifdef __ISOLATED__
    if(this_run->rank_x == 0) { /* at the boundary of the simulation box */
#ifdef REFLECT_BOUNDARY_X_LO
      ixm2 = -1-ixm2;
      reflect_lo2 = MESH(ixm2, iy, iz);
      reflect_lo2.momx *= -1.0;
      *mesh_im2 = reflect_lo2;
#elif defined(OUTFLOW_BOUNDARY_X_LO)
      ixm2 = 0;
      *mesh_im2 = MESH(ixm2, iy, iz);
#endif
    }else{ /* at the boundary of the local volume */
      *mesh_im2 = pad->PAD_X_LO(ixm2+2, iy, iz);
    }
#else /* ! __ISOLATED__  */
    *mesh_im2 = pad->PAD_X_LO(ixm2+2, iy, iz);
#endif
  }else{
    *mesh_im2 = MESH(ixm2, iy, iz);
  }

  if(ixm1 < 0) {
#ifdef __ISOLATED__
    if(this_run->rank_x == 0) {
#ifdef REFLECT_BOUNDARY_X_LO
      ixm1 = -1-ixm1;
      reflect_lo1 = MESH(ixm1, iy, iz);
      reflect_lo1.momx *= -1.0;
      *mesh_im1 = reflect_lo1;
#elif defined(OUTFLOW_BOUNDARY_X_LO)
      ixm1 = 0;
      *mesh_im1 = MESH(ixm1, iy, iz);
#endif
    }else{
      *mesh_im1 = pad->PAD_X_LO(ixm1+2, iy, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_im1 = pad->PAD_X_LO(ixm1+2, iy, iz);
#endif
  }else{
    *mesh_im1 = MESH(ixm1, iy, iz);
  }

  if(ixp1>NMESH_X_LOCAL-1) {
#ifdef __ISOLATED__
    if(this_run->rank_x == NNODE_X-1) {
#ifdef REFLECT_BOUNDARY_X_HI
      ixp1 = 2*NMESH_X_LOCAL-ixp1-1;
      reflect_hi1 = MESH(ixp1, iy, iz);
      reflect_hi1.momx *= -1.0;
      *mesh_ip1 = reflect_hi1;
#elif defined(OUTFLOW_BOUNDARY_X_HI)
      ixp1 = NMESH_X_LOCAL-1;
      *mesh_ip1 = MESH(ixp1, iy, iz);
#endif
    }else{
      *mesh_ip1 = pad->PAD_X_HI(ixp1-NMESH_X_LOCAL, iy, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_ip1 = pad->PAD_X_HI(ixp1-NMESH_X_LOCAL, iy, iz);
#endif
  }else{
    *mesh_ip1 = MESH(ixp1, iy, iz);
  }

  if(ixp2>NMESH_X_LOCAL-1) {
#ifdef __ISOLATED__
    if(this_run->rank_x == NNODE_X-1) {
#ifdef REFLECT_BOUNDARY_X_HI
      ixp2 = 2*NMESH_X_LOCAL-ixp2-1;
      reflect_hi2 = MESH(ixp2, iy, iz);
      reflect_hi2.momx *= -1.0;
      *mesh_ip2 = reflect_hi2;
#elif defined(OUTFLOW_BOUNDARY_X_HI)
      ixp2 = NMESH_X_LOCAL-1;
      *mesh_ip2 = MESH(ixp2, iy, iz);
#endif
    }else{
      *mesh_ip2 = pad->PAD_X_HI(ixp2-NMESH_X_LOCAL, iy, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_ip2 = pad->PAD_X_HI(ixp2-NMESH_X_LOCAL, iy, iz);
#endif
  }else{
    *mesh_ip2 = MESH(ixp2, iy, iz);
  }
}

void set_pad_y(struct pad_region *pad, struct fluid_mesh *mesh,
	       struct fluid_mesh *mesh_im2, struct fluid_mesh *mesh_im1,
	       struct fluid_mesh *mesh_ip1, struct fluid_mesh *mesh_ip2,
	       int ix, int iy, int iz,
	       struct run_param *this_run)
{
  struct fluid_mesh reflect_lo1, reflect_lo2;
  struct fluid_mesh reflect_hi1, reflect_hi2;
  int iym2, iym1, iyp1, iyp2;
  
  iym2 = iy-2;  iym1 = iy-1;
  iyp1 = iy+1;  iyp2 = iy+2;
  
  if(iym2 < 0) {
#ifdef __ISOLATED__
    if(this_run->rank_y == 0) { /* at the boundary of the simulation box */
#ifdef REFLECT_BOUNDARY_Y_LO
      iym2 = -1-iym2;
      reflect_lo2 = MESH(ix, iym2, iz);
      reflect_lo2.momy *= -1.0;
      *mesh_im2 = reflect_lo2;
#elif defined(OUTFLOW_BOUNDARY_Y_LO)
      iym2 = 0;
      *mesh_im2 = MESH(ix, iym2, iz);
#endif
    }else{ /* at the boundary of the local volume */
      *mesh_im2 = pad->PAD_Y_LO(ix, iym2+2, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_im2 = pad->PAD_Y_LO(ix, iym2+2, iz);
#endif
  }else{
    *mesh_im2 = MESH(ix, iym2, iz);
  }

  if(iym1 < 0) {
#ifdef __ISOLATED__
    if(this_run->rank_y == 0) {
#ifdef REFLECT_BOUNDARY_Y_LO
      iym1 = -1-iym1;
      reflect_lo1 = MESH(ix, iym1, iz);
      reflect_lo1.momy *= -1.0;
      *mesh_im1 = reflect_lo1;
#elif defined(OUTFLOW_BOUNDARY_Y_LO)
      iym1 = 0;
      *mesh_im1 = MESH(ix, iym1, iz);
#endif
    }else{
      *mesh_im1 = pad->PAD_Y_LO(ix, iym1+2, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_im1 = pad->PAD_Y_LO(ix, iym1+2, iz);
#endif
  }else{
    *mesh_im1 = MESH(ix, iym1, iz);
  }

  if(iyp1>NMESH_Y_LOCAL-1) {
#ifdef __ISOLATED__
    if(this_run->rank_y == NNODE_Y-1) {
#ifdef REFLECT_BOUNDARY_Y_HI
      iyp1 = 2*NMESH_Y_LOCAL-iyp1-1;
      reflect_hi1 = MESH(ix, iyp1, iz);
      reflect_hi1.momy *= -1.0;
      *mesh_ip1 = reflect_hi1;
#elif defined(OUTFLOW_BOUNDARY_Y_HI)
      iyp1 = NMESH_Y_LOCAL-1;
      *mesh_ip1 = MESH(ix, iyp1, iz);
#endif
    }else{
      *mesh_ip1 = pad->PAD_Y_HI(ix, iyp1-NMESH_Y_LOCAL, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_ip1 = pad->PAD_Y_HI(ix, iyp1-NMESH_Y_LOCAL, iz);
#endif
  }else{
    *mesh_ip1 = MESH(ix, iyp1, iz);
  }

  if(iyp2>NMESH_Y_LOCAL-1) {
#ifdef __ISOLATED__
    if(this_run->rank_y == NNODE_Y-1) {
#ifdef REFLECT_BOUNDARY_Y_HI
      iyp2 = 2*NMESH_Y_LOCAL-iyp2-1;
      reflect_hi2 = MESH(ix, iyp2, iz);
      reflect_hi2.momy *= -1.0;
      *mesh_ip2 = reflect_hi2;
#elif defined(OUTFLOW_BOUNDARY_Y_HI)
      iyp2 = NMESH_Y_LOCAL-1;
      *mesh_ip2 = MESH(ix, iyp2, iz);
#endif
    }else{
      *mesh_ip2 = pad->PAD_Y_HI(ix, iyp2-NMESH_Y_LOCAL, iz);
    }
#else /* ! __ISOLATED__ */
    *mesh_ip2 = pad->PAD_Y_HI(ix, iyp2-NMESH_Y_LOCAL, iz);
#endif
  }else{
    *mesh_ip2 = MESH(ix, iyp2, iz);
  }
}

void set_pad_z(struct pad_region *pad, struct fluid_mesh *mesh,
	       struct fluid_mesh *mesh_im2, struct fluid_mesh *mesh_im1,
	       struct fluid_mesh *mesh_ip1, struct fluid_mesh *mesh_ip2,
	       int ix, int iy, int iz,
	       struct run_param *this_run)
{
  struct fluid_mesh reflect_lo1, reflect_lo2;
  struct fluid_mesh reflect_hi1, reflect_hi2;
  int izm2, izm1, izp1, izp2;

  izm2 = iz-2;  izm1 = iz-1;
  izp1 = iz+1;  izp2 = iz+2;

  if(izm2 < 0) {
#ifdef __ISOLATED__
    if(this_run->rank_z == 0) { /* at the boundary of the simulation box */
#ifdef REFLECT_BOUNDARY_Z_LO
      izm2 = -1-izm2;
      reflect_lo2 = MESH(ix, iy, izm2);
      reflect_lo2.momz *= -1.0;
      *mesh_im2 = reflect_lo2;
#elif defined(OUTFLOW_BOUNDARY_Z_LO)
      izm2 = 0;
      *mesh_im2 = MESH(ix, iy, izm2);
#endif
    }else{ /* at the boundary of the local volume */
      *mesh_im2 = pad->PAD_Z_LO(ix, iy, izm2+2);
    }
#else /* ! __ISOLATED__ */
    *mesh_im2 = pad->PAD_Z_LO(ix, iy, izm2+2);
#endif
  }else{
    *mesh_im2 = MESH(ix, iy, izm2);
  }

  if(izm1 < 0) {
#ifdef __ISOLATED__
    if(this_run->rank_z == 0) {
#ifdef REFLECT_BOUNDARY_Z_LO
      izm1 = -1-izm1;
      reflect_lo1 = MESH(ix, iy, izm1);
      reflect_lo1.momz *= -1.0;
      *mesh_im1 = reflect_lo1;
#elif defined(OUTFLOW_BOUNDARY_Z_LO)
      izm1 = 0;
      *mesh_im1 = MESH(ix, iy, izm1);
#endif
    }else{
      *mesh_im1 = pad->PAD_Z_LO(ix, iy, izm1+2);
    }
#else /* ! __ISOLATED__ */
    *mesh_im1 = pad->PAD_Z_LO(ix, iy, izm1+2);
#endif
  }else{
    *mesh_im1 = MESH(ix, iy, izm1);
  }

  if(izp1>NMESH_Z_LOCAL-1) {
#ifdef __ISOLATED__
    if(this_run->rank_z == NNODE_Z-1) {
#ifdef REFLECT_BOUNDARY_Z_HI
      izp1 = 2*NMESH_Z_LOCAL-izp1-1;
      reflect_hi1 = MESH(ix, iy, izp1);
      reflect_hi1.momz *= -1.0;
      *mesh_ip1 = reflect_hi1;
#elif defined(OUTFLOW_BOUNDARY_Z_HI)
      izp1 = NMESH_Z_LOCAL-1;
      *mesh_ip1 = MESH(ix, iy, izp1);
#endif
    }else{
      *mesh_ip1 = pad->PAD_Z_HI(ix, iy, izp1-NMESH_Z_LOCAL);
    }
#else /* ! __ISOLATED__ */
    *mesh_ip1 = pad->PAD_Z_HI(ix, iy, izp1-NMESH_Z_LOCAL);
#endif
  }else{
    *mesh_ip1 = MESH(ix, iy, izp1);
  }

  if(izp2>NMESH_Z_LOCAL-1) {
#ifdef __ISOLATED__
    if(this_run->rank_z == NNODE_Z-1) {
#ifdef REFLECT_BOUNDARY_Z_HI
      izp2 = 2*NMESH_Z_LOCAL-izp2-1;
      reflect_hi2 = MESH(ix, iy, izp2);
      reflect_hi2.momz *= -1.0;
      *mesh_ip2 = reflect_hi2;
#elif defined(OUTFLOW_BOUNDARY_Z_HI)
      izp2 = NMESH_Z_LOCAL-1;
      *mesh_ip2 = MESH(ix, iy, izp2);
#endif
    }else{
      *mesh_ip2 = pad->PAD_Z_HI(ix, iy, izp2-NMESH_Z_LOCAL);
    }
#else /* ! __ISOLATED__ */
    *mesh_ip2 = pad->PAD_Z_HI(ix, iy, izp2-NMESH_Z_LOCAL);
#endif
  }else{
    *mesh_ip2 = MESH(ix, iy, izp2);
  }
}
