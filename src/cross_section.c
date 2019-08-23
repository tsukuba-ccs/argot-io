#include "run_param.h"
#include "cross_section.h"
#include "chemistry.h"

void setup_cross_section(struct cross_section *csect, struct freq_param *freq)
{
  for(int inu=0;inu<NGRID_NU;inu++) {
    csect[inu].csect_HI   = csectHI(freq->nu[inu]);
#ifdef __HELIUM__
    csect[inu].csect_HeI  = csectHeI(freq->nu[inu]);
    csect[inu].csect_HeII = csectHeII(freq->nu[inu]);
#endif /* __HELIUM__ */
  }
}
