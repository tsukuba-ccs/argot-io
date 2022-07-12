#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>
#include <sys/time.h>
#include <unistd.h>

#include "run_param.h"

#ifndef CLK_TCK
#define CLK_TCK sysconf(_SC_CLK_TCK)
#endif

float timing(struct tms from, struct tms to) 
{
  return ((float)(to.tms_utime+to.tms_stime)-(float)(from.tms_utime+from.tms_stime))/(float)CLK_TCK;
}

float wallclock_timing(struct timeval from, struct timeval to)
{
  return ((to.tv_sec+to.tv_usec*1.e-6)-(from.tv_sec+from.tv_usec*1.e-6));
}

int start_timing(struct timeval *start_tv, struct tms *start_tms)
{
  times(start_tms);
  gettimeofday(start_tv, NULL);

  return 0;
}

int end_timing(struct timeval *start_tv, struct timeval *end_tv,
                struct tms *start_tms, struct tms *end_tms,
                char *label, struct run_param *this_run)
{
  times(end_tms);
  gettimeofday(end_tv, NULL);

  fprintf(this_run->proc_file,
          "# %s : %12.4e [sec] (CPU) / %12.4e [sec] (Wall) \n",label,
          timing(*start_tms, *end_tms), wallclock_timing(*start_tv, *end_tv));
  fflush(this_run->proc_file);

  return 0;
}
