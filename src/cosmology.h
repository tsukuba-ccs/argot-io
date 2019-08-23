#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifndef __COSMOLOGY__
#define __COSMOLOGY__

struct cosmology {
  float omega0, lambda0, omegab, hubble;
  float tend;
};

float ztotime(float znow, struct cosmology cosm);
float atotime(float znow, struct cosmology cosm);
float timetoz(float tnow, struct cosmology cosm);
float timetoa(float tnow, struct cosmology cosm);
void  funcd(double x, double *f, double *df, double tau);
float rtsafe(double x1, double x2, double xacc, double tau);

#endif /* __COSMOLOGY__ */

#ifdef __cplusplus
}
#endif
