#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "run_param.h"
#include "source.h"
#include "constants.h"
#include "chemistry.h"

static double blackbody(double nu, double T)
// nu : frequency in nuL
// T  : temperatureu in Kelvin
// normalized at nu=nuL
{
  double I_bb, I_bb_nuL;

  I_bb     = SQR(nu*nuL)/(exp(hplanck*nu*nuL/(kboltz*T))-1.e0);
  I_bb_nuL = SQR(nuL)/(exp(hplanck*nuL/(kboltz*T))-1.e0);

  return (I_bb/I_bb_nuL);
}

static double power_law(double nu, double pl_indx)
// nu : frequency in nuL
// pl_indx : power-law index in terms of photon number spectrum;
// normalized at nu=nuL
{
  double I_pl, I_pl_nuL, I_pl_normalized;

  I_pl = pow(nu*nuL,pl_indx);
  I_pl_nuL = pow(nuL,pl_indx);
  
  I_pl_normalized = pow(nu,pl_indx);

  return I_pl_normalized;
}

void setup_freq_param(struct freq_param *freq)
{
  freq->log_nu_min = 0.0;
  freq->log_nu_max = 2.0;
  freq->dlog_nu = (freq->log_nu_max-freq->log_nu_min)/(float)(NGRID_NU);

  int inu;
  
  for(inu=0;inu<NGRID_NU;inu++) 
    freq->nu[inu] = pow(10.0,freq->log_nu_min + (inu+0.5)*freq->dlog_nu);
  
}


#define NITER (5)
#define ALLOWED_ERROR (1.0e-8)

static double count_midpnt(double (*func)(double, double), 
                           double a, double b, int n, double param)
{
  double x,tnm,sum,del,ddel;
  static double s;
  int it,j;

  if (n == 1) {
    return (s=(b-a)*(*func)(0.5*(a+b), param));
  } else {
    for(it=1,j=1;j<n-1;j++) it *= 3;
    tnm=it;
    del=(b-a)/(3.0*tnm);
    ddel=del+del;
    x=a+0.5*del;
    sum=0.0;
    for (j=1;j<=it;j++) {
      sum += (*func)(x,param);
      x += ddel;
      sum += (*func)(x,param);
      x += del;
    }
    s=(s+(b-a)*sum/tnm)/3.0;
    return s;
  }
}

double count_photon(double (*func)(double, double), double param, 
                    double nu1, double nu2)
{
  int p;
  double N, N_old;

  for(p=1;p<=NITER;p++) {
    if(p>1) N_old = N;
    N = count_midpnt(func, nu1, nu2, p, param);
    if(p>1) {
      if(fabs(N-N_old) < ALLOWED_ERROR*fabs(N_old)) return N;
    }
  }

  return N;
}

void setup_photon_rate(struct freq_param *freq, struct radiation_src *src, 
                       double total_photon_rate)
{
  int inu;
  double count_sum;

  count_sum=0.0;
  for(inu=0;inu<NGRID_NU;inu++) {
    double nu1, nu2;

    nu1 = freq->log_nu_min + freq->dlog_nu*(double)inu;
    nu2 = nu1 + freq->dlog_nu;

    nu1 = pow(10.0, nu1);
    nu2 = pow(10.0, nu2);

    if(src->type==0) { /* blackbody */
      src->photon_rate[inu] = count_photon(blackbody, 
                                           (double) src->param, nu1, nu2);
    }else if(src->type == 1) { /* power-law */
      src->photon_rate[inu] = count_photon(power_law, 
					   (double) src->param, nu1, nu2);
    }
    count_sum += src->photon_rate[inu];
  }


  for(inu=0;inu<NGRID_NU;inu++) {
    src->photon_rate[inu] *= total_photon_rate/count_sum;
  }
  
}

#if 0
int main(int argc, char **argv) 
{
  struct radiation_src src;
  static struct freq_param freq;

  src.type=0; 
  src.param=1.0e5;

  setup_freq_param(&freq);
  setup_photon_rate(&freq, &src, 5.0e48);

  printf("# freq.dlog_nu = %14.6e\n", freq.dlog_nu);

  int inu;
  for(inu=0;inu<NGRID_NU;inu++) {
    printf("%16.8e %16.8e \n",freq.nu[inu],src.photon_rate[inu]);
  }
}
#endif
