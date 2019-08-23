#ifdef __cplusplus
extern "C" {
#endif

#include "source.h"

#ifndef __ARGOT_RADIATION__
#define __ARGOT_RADIATION__

struct photoion_rate {
  float GammaHI;
  float HeatHI;
#ifdef __HELIUM__
  float GammaHeI, GammaHeII;
  float HeatHeI, HeatHeII;
#endif
#ifdef __HYDROGEN_MOL__
  float GammaHM, GammaH2I_I, GammaH2I_II, GammaH2II_I, GammaH2II_II;
  float HeatHM, HeatH2I_I, HeatH2I_II, HeatH2II_I, HeatH2II_II;
#endif
};

struct ray_segment {
  uint64_t ray_indx;  /* index of the light ray */
  int   local_rank;
  int   target_rank;
  float xpos_start,ypos_start,zpos_start;
  float xpos_end,  ypos_end,  zpos_end;

  float optical_depth_HI;
#ifdef __HELIUM__
  float optical_depth_HeI, optical_depth_HeII;
#endif /* __HELIUM__ */
#ifdef __HYDROGEN_MOL__
  float optical_depth_HM;
  float optical_depth_H2I_I, optical_depth_H2I_II;
  float optical_depth_H2II_I, optical_depth_H2II_II;
#endif /* __HYDROGEN_MOL__ */

};

struct light_ray {
  /* radiation source */
  struct radiation_src src;

  /* ray segment */
  int    num_segment;
  struct ray_segment segment[NSEG_PER_RAY];

  /* global indices of the target mesh */
  int ix_target, iy_target, iz_target;

  float optical_depth_HI;
#ifdef __HELIUM__
  float optical_depth_HeI, optical_depth_HeII;
#endif /* __HELIUM__ */
#ifdef __HYDROGEN_MOL__
  float optical_depth_HM;
  float optical_depth_H2I_I, optical_depth_H2I_II;
  float optical_depth_H2II_I, optical_depth_H2II_II;
#endif /* __HYDROGEN_MOL__ */
};

struct light_ray_IO {
  /* radiation source */
  struct radiation_src src;

  /* global indices of the target mesh */
  int ix_target, iy_target, iz_target;

  float optical_depth_HI;
#ifdef __HELIUM__
  float optical_depth_HeI, optical_depth_HeII;
#endif /* __HELIUM__ */
#ifdef __HYDROGEN_MOL__
  float optical_depth_HM;
  float optical_depth_H2I_I, optical_depth_H2I_II;
  float optical_depth_H2II_I, optical_depth_H2II_II;
#endif /* __HYDROGEN_MOL__ */
};

#endif   /* __ARGOT_RADIATION__ */

#ifdef __cplusplus
}
#endif
