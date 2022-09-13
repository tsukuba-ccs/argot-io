#ifndef __ARGOT_FLUID__
#define __ARGOT_FLUID__

#include "chemistry.h"

#define __X_BOUNDARY_OUTFLOW__
#define __Y_BOUNDARY_OUTFLOW__
#define __Z_BOUNDARY_OUTFLOW__

#ifdef __X_BOUNDARY_REFLECT__
#define REFLECT_BOUNDARY_X_LO
#define REFLECT_BOUNDARY_X_HI
#endif

#ifdef __X_BOUNDARY_OUTFLOW__
#define OUTFLOW_BOUNDARY_X_LO
#define OUTFLOW_BOUNDARY_X_HI
#endif

#ifdef __Y_BOUNDARY_REFLECT__
#define REFLECT_BOUNDARY_Y_LO
#define REFLECT_BOUNDARY_Y_HI
#endif

#ifdef __Y_BOUNDARY_OUTFLOW__
#define OUTFLOW_BOUNDARY_Y_LO
#define OUTFLOW_BOUNDARY_Y_HI
#endif

#ifdef __Z_BOUNDARY_REFLECT__
#define REFLECT_BOUNDARY_Z_LO
#define REFLECT_BOUNDARY_Z_HI
#endif

#ifdef __Z_BOUNDARY_OUTFLOW__
#define OUTFLOW_BOUNDARY_Z_LO
#define OUTFLOW_BOUNDARY_Z_HI
#endif

#define NPAD (2)

#define PAD_Z_LO(ix, iy, iz) \
  pad_z_lo[(iz)+NPAD*((iy)+NMESH_Y_LOCAL*(ix))]
#define PAD_Z_HI(ix, iy, iz) \
  pad_z_hi[(iz)+NPAD*((iy)+NMESH_Y_LOCAL*(ix))]
#define PAD_X_LO(ix, iy, iz) \
  pad_x_lo[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))]
#define PAD_X_HI(ix, iy, iz) \
  pad_x_hi[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))]
#define PAD_Y_LO(ix, iy, iz) \
  pad_y_lo[(iz)+NMESH_Z_LOCAL*((iy)+NPAD*(ix))]
#define PAD_Y_HI(ix, iy, iz) \
  pad_y_hi[(iz)+NMESH_Z_LOCAL*((iy)+NPAD*(ix))]

#define HIGH_MACH_THRESHOLD (0.1)
#define SHOCK_PRESS_THRESHOLD (0.3)
#define COURANT_FACT (0.1)

// HLLC flag
//#define __ACOUSTIC_TYPE__
#define __TWO_RAREFACTION__
//#define __TWO_SCHOCK__

/* parameters for limiter */
#define KAPPA    (-2.5)
#define KPLUS    (1.0+KAPPA)
#define KMINUS   (1.0-KAPPA)
#define BPARAM   ((3.0-KAPPA)/(1.0-KAPPA))

/* paramters for AUSM+ scheme */
#define alpha (0.1875) /* 3.0/16.0 */
#define beta  (0.125) /* 1.0/8.0 */

#define ABS(A) ( ((A)>=0.0) ? (A):(-A) )
#define SIGN(A) ( ((A)>0.0) ? (1.0):((A)<0.0) ? (-1.0):(0.0) )

#define Mp_beta(M) ( ( 0.25)*SQR((M)+1.0)+(beta)*SQR((M)*(M)-1.0) )
#define Mm_beta(M) ( (-0.25)*SQR((M)-1.0)-(beta)*SQR((M)*(M)-1.0) )

#define Pp_alpha(M) (( 0.25)*SQR((M)+1.0)*(2.0-(M)) + (alpha)*(M)*SQR((M)*(M)-1.0) )
#define Pm_alpha(M) (( 0.25)*SQR((M)-1.0)*(2.0+(M)) - (alpha)*(M)*SQR((M)*(M)-1.0) )

#define M_p(M) ( (fabsf(M)>=1.0) ? ((0.5)*((M) + fabsf(M))):(Mp_beta(M)) )
#define M_m(M) ( (fabsf(M)>=1.0) ? ((0.5)*((M) - fabsf(M))):(Mm_beta(M)) )
#define P_p(M) ( (fabsf(M)> 1.0) ? ((0.5)*(1.0 + SIGN(M))):(Pp_alpha(M)) )
#define P_m(M) ( (fabsf(M)> 1.0) ? ((0.5)*(1.0 - SIGN(M))):(Pm_alpha(M)) )

#define MUSCL_L(Ujm1,Uj,Ujp1) ( (Uj) + 0.25*( (KMINUS)*minmod(((Uj)-(Ujm1)) ,(BPARAM)*((Ujp1)-(Uj))) )+ 0.25*( (KPLUS)*minmod(((Ujp1)-(Uj)) ,(BPARAM)*((Uj)-(Ujm1))) ) )
#define MUSCL_R(Uj,Ujp1,Ujp2) ( (Ujp1) - 0.25*( (KMINUS)*minmod(((Ujp2)-(Ujp1)) ,(BPARAM)*((Ujp1)-(Uj))) )- 0.25*( (KPLUS)*minmod(((Ujp1)-(Uj)) ,(BPARAM)*((Ujp2)-(Ujp1))) ) )

#define minmod(x,y) ( (((x)*(y))>0.0)&&(ABS(x)<=ABS(y)) ? (x):((((x)*(y))>0.0)&&(ABS(x)>ABS(y))) ? (y):(0.0) )

#define FLUID_TINY (1.0e-32)

struct fluid_mesh_raw {
  float dens;
  float eneg;
  float momx, momy, momz;
  float uene, duene, durad;
  float etrp;
  float pot;
  float dens_prev; /* density in the previous step */

  struct prim_chem chem;
  struct prim_chem prev_chem; /* chemical composition in the previous iteration */
  float  prev_uene; /* specific energy in the previous iteration */

  short under_shock;
  short high_mach;
};

struct fluid_mesh {
  float dens;
  float eneg;
  float momx, momy, momz;
  float uene, duene, durad;
  float etrp;
  float pot;
  float dens_prev; /* density in the previous step */

  struct prim_chem chem;
  struct prim_chem prev_chem; /* chemical composition in the previous iteration */
  float  prev_uene; /* specific energy in the previous iteration */

  short under_shock;
  short high_mach;
  char pad[128 - sizeof(struct fluid_mesh_raw)];
};

struct fluid_mesh_io_raw {
  float dens;
  float eneg;
  float momx, momy, momz;
  float uene, pot;
  struct prim_chem chem;
};

struct fluid_mesh_io {
  float dens;
  float eneg;
  float momx, momy, momz;
  float uene, pot;
  struct prim_chem chem;
  char pad[64 - sizeof(struct fluid_mesh_io_raw)];
};

struct pad_region {
  struct fluid_mesh *pad_x_lo, *pad_x_hi;
  struct fluid_mesh *pad_y_lo, *pad_y_hi;
  struct fluid_mesh *pad_z_lo, *pad_z_hi;
};


#endif /* __ARGOT_FLUID__ */
