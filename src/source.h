#ifndef __ARGOT_SOURCE__
#define __ARGOT_SOURCE__

#include "run_param.h"

struct freq_param {
  double log_nu_min, log_nu_max, dlog_nu;
  double nu[NGRID_NU];
};

struct radiation_src {
  int type; /* 0:blackbody / 1:power law / >=2:other */
  float param; 
  float xpos, ypos, zpos;
  double photon_rate[NGRID_NU];
};

/* prototypes */
void setup_freq_param(struct freq_param*);
void setup_photon_rate(struct freq_param*, struct radiation_src*, double);

#endif /* __ARGOT_SOURCE__ */
