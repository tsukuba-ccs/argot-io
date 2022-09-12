#ifdef __cplusplus
extern "C" {
#endif

#ifndef __ARGOT_RUN_PARAM__
#define __ARGOT_RUN_PARAM__

#include <stdint.h>

#define NMESH_X_TOTAL (128)
#define NMESH_Y_TOTAL (128)
#define NMESH_Z_TOTAL (128)

#define NNODE_X (2)
#define NNODE_Y (2)
#define NNODE_Z (2)

#define NMESH_Z_TOTAL_P2 (NMESH_Z_TOTAL+2)

#ifdef __ISOLATED_GRAV__
#define NMESH_X_POTEN    (NMESH_X_TOTAL*2)
#define NMESH_Y_POTEN    (NMESH_Y_TOTAL*2)
#define NMESH_Z_POTEN    (NMESH_Z_TOTAL*2)
#define NMESH_Z_POTEN_P2 (NMESH_Z_POTEN+2)
#define NMESH_X_GREEN    (NMESH_X_POTEN)
#define NMESH_Y_GREEN    (NMESH_Y_POTEN)
#define NMESH_Z_GREEN    (NMESH_Z_POTEN)
#define NMESH_Z_GREEN_P2 (NMESH_Z_GREEN+2)
#else /* !__ISOLATED_GRAV__ or periodic boundary */
#define NMESH_X_POTEN    (NMESH_X_TOTAL)
#define NMESH_Y_POTEN    (NMESH_Y_TOTAL)
#define NMESH_Z_POTEN    (NMESH_Z_TOTAL)
#define NMESH_Z_POTEN_P2 (NMESH_Z_POTEN+2)
#define NMESH_X_GREEN    (NMESH_X_TOTAL/2+1)
#define NMESH_Y_GREEN    (NMESH_Y_TOTAL/2+1)
#define NMESH_Z_GREEN    (NMESH_Z_TOTAL/2+1)
#endif

#define NNODE (NNODE_X*NNODE_Y*NNODE_Z)
#define NSEG_PER_RAY (5) /* max. number of segments for a single light ray */

#define NMESH_X_LOCAL (NMESH_X_TOTAL/NNODE_X)
#define NMESH_Y_LOCAL (NMESH_Y_TOTAL/NNODE_Y)
#define NMESH_Z_LOCAL (NMESH_Z_TOTAL/NNODE_Z)
#define NMESH_LOCAL (NMESH_X_LOCAL*NMESH_Y_LOCAL*NMESH_Z_LOCAL)

#define NSOURCE_MAX (32768)

#define NGRID_NU (32)

#define ARGOT_THETA_CRIT (0.7)

// number of chemical species 

#define NSPECIES_HYDROGEN_ATOM (2)
#define NSPECIES_HELIUM (3)
#define NSPECIES_HYDROGEN_MOL (3)

#ifdef __HYDROGEN_MOL__
#define NSPECIES_HYDROGEN (NSPECIES_HYDROGEN_ATOM+NSPECIES_HYDROGEN_MOL)
#else /* ! __HYDROGEN_MOL__ */
#define NSPECIES_HYDROGEN (NSPECIES_HYDROGEN_ATOM)
#endif 

#ifdef __HELIUM__
#define NSPECIES (NSPECIES_HYDROGEN+NSPECIES_HELIUM)
#else /* ! __HELIUM__ */
#define NSPECIES (NSPECIES_HYDROGEN)
#endif


// number of radiation transfer channel 
#define NCHANNEL_HYDROGEN (1)
#define NCHANNEL_HELIUM (2)
#define NCHANNEL_HYDROGEN_MOL (5)

#ifdef __HELIUM__

#ifdef __HYDROGEN_MOL__
#define NCHANNEL (NCHANNEL_HYDROGEN+NCHANNEL_HELIUM+NCHANNEL_HYDROGEN_MOL)
#else  /* ! __HYDROGEN_MOL__ */
#define NCHANNEL (NCHANNEL_HYDROGEN+NCHANNEL_HELIUM)
#endif /* __HYDROGEN_MOL__ */

#else /* ! __HELIUM__ */

#ifdef __HYDROGEN_MOL__
#define NCHANNEL (NCHANNEL_HYDROGEN+NCHANNEL_HYDROGEN_MOL)
#else  /* ! __HYDROGEN_MOL__ */
#define NCHANNEL (NCHANNEL_HYDROGEN)
#endif /* __HYDROGEN_MOL__ */

#endif

//#define MAX_NRAY_PER_TARGET (8192)
#define MAX_NRAY_PER_TARGET (NMESH_LOCAL)
//#define MAX_NRAY_PER_TARGET (NSOURCE_MAX)
#define NMESH_PER_LOOP (NMESH_LOCAL)
#define ALLOWED_MEM_SIZE_IN_GB (20.0)

#define DTFACT_RAD (5.0) /* ratio of the radiation dt relative to the chemical one */

#ifdef __USE_GPU__
#define NMAX_CUDA_DEV (4) /* number of CUDA devices per node */
#define NMESH_PER_BLOCK (256)
#define NSEG_MAX_PER_BLOCK (512)
#define NSEG_MAX_PER_DEV (16384)
#else /* !__USE_GPU__ */
#define NMAX_CUDA_DEV (0) /* dummy */
#endif /* __USE_GPU__ */

#define N_SIDE (8)

#define TORR_DIST_FACT (1.0e-2)

#include <stdio.h>
#include "cosmology.h"
#include "source.h"
#include "fluid.h"
#include "cross_section.h"

struct run_param {

  /* model name */
  char *model_name;

  FILE *proc_file, *diag_file, *err_file;

  int mpi_rank, mpi_nproc;
  int rank_x, rank_y, rank_z;
  int nnode_x, nnode_y, nnode_z;

  double lunit, munit, tunit,eunit;
  double denstonh, uenetok, masstonh;

  float anow, tnow, znow, hnow;
  struct cosmology cosm;
  float tend;
  float dtime;

  struct freq_param freq;
  struct cross_section csect[NGRID_NU];

  int nspecies; /* number of chemical species */
  int nchannel; /* number of radiation channels */
  int ngrid_nu;

  int step;

  int output_indx;
  int noutput;
  float *output_timing;

  uint64_t  nsrc;  /* number of the total radiating sources */
  uint64_t  nray;  /* number of the rays targeted to meshes in this domain */
  uint64_t  nseg;  /* number of the ray segments that get across this domain */

  /* global min and max of the x-, y- and z-coordinate */
  float  xmin, ymin, zmin;
  float  xmax, ymax, zmax;

  /* local min and max of the x-, y- and z-coordinate */
  float  xmin_local, ymin_local, zmin_local;
  float  xmax_local, ymax_local, zmax_local;

  /* number of the global mesh */
  int    nmesh_x_total, nmesh_y_total, nmesh_z_total;

  /* number of the local mesh */
  int    nmesh_x_local, nmesh_y_local, nmesh_z_local;

  /* dx dy dz */
  float  delta_x, delta_y, delta_z;
  
  /* cirital theta parameter */
  float theta_crit;

  /* number of target meshes to simultaneously compute radiation field */
  int nmesh_per_loop;

  /* pad region for hydro code */
  struct pad_region pad;

  float input_mesh_tp, input_src_tp;
  float input_mesh_wt, input_src_wt;

  uint64_t output_file_size;
  float output_wt;
};

#ifndef MIN
#define MIN(a,b) ((a) > (b) ? (b) : (a))
#endif

#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif

#ifndef CUBE
#define CUBE(x) ((x)*(x)*(x))
#endif

#ifndef QUAD
#define QUAD(x) ((x)*(x)*(x)*(x))
#endif

#ifndef NORML2
#define NORML2(x,y,z) (SQR(x)+SQR(y)+SQR(z))
#endif

#define PRINT_SEG(s) \
  fprintf(this_run->proc_file,						\
	  "%llu%d %d %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %d %d %d %d %d %d\n",	\
	  s->ray_indx, s->target_rank, s->local_rank,			\
	  s->xpos_start, s->ypos_start, s->zpos_start,			\
	  s->xpos_end, s->ypos_end, s->zpos_end,                        \
	  this_run->xmin_local+x_cur, \
	  this_run->ymin_local+y_cur, \
	  this_run->zmin_local+z_cur, \
	  abs_dist, ix_cur, iy_cur, iz_cur, ix_end, iy_end, iz_end)

#endif  /* __ARGOT_RUN_PARAM__ */

#ifdef __cplusplus
}
#endif
