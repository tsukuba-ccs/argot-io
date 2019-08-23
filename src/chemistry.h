#ifdef __cplusplus
extern "C" {
#endif

#ifndef __CHEMISTRY__
#define __CHEMISTRY__

#include "constants.h"

#define K_to_eV  (8.61623e-5)
#define eV_to_Hz (2.41838e14)
#define nuL      (3.28e15)

#define GAMMA_MONOATOMIC (1.6666666)
#define GAMM1_MONOATOMIC (0.6666666)
#define GAMMA_DIATOMIC   (1.4)
#define GAMM1_DIATOMIC   (0.4)

#ifdef __HELIUM__
#define XHYDROGEN (0.755)
#define YHELIUM  (0.245)
#define HELIUM_FACT (0.0811258278)
#else 
#define XHYDROGEN (1.0)
#define YHELIUM  (0.0)
#define HELIUM_FACT (0.0)
#endif

#define HI_EDGE (1.0)
#define HeI_EDGE (1.8125)
#define HeII_EDGE (4.0)

struct prim_chem{
  float wmol;                    // mean molecular weight
  float felec;                   // ne/nH

  float fHI, fHII;               // number fraction relative to nH
#ifdef __HYDROGEN_MOL__
  float fH2I, fH2II, fHM;        // number fraction relative to nH
#endif
#ifdef __HELIUM__
  float fHeI, fHeII, fHeIII;     // number fraction relative to nHe
#endif

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

#ifdef __HELIUM__

#ifdef __HYDROGEN_MOL__
#define WMOL(c) (1.0/(XHYDROGEN*(1.0-(c).fH2I-(c).fH2II+(c).felec)+0.25*YHELIUM))
#else 
#define WMOL(c) (1.0/((1.0+(c).felec)*XHYDROGEN + 0.25*YHELIUM))
#endif

#else /* ! __HELIUM__ */

#ifdef __HYDRIGEN_MOL__
#define WMOL(c) (1.0/(1.0-(c).fH2I-(c).fH2II+(c).felec))
#else
#define WMOL(c) (1.0/(1.0+(c).felec))
#endif

#endif

extern double GammaHI(float);
extern double GammaHM(float);
extern double GammaH2I_I(float);
extern double GammaH2I_II(float);
extern double GammaH2II_I(float);
extern double GammaH2II_II(float);
extern double GammaHeI(float);
extern double GammaHeII(float);

extern double HeatHI(float);
extern double HeatHeI(float);
extern double HeatHeII(float);
// extern double HeatHM(float zred); 
// extern double HeatH2I_I(float zred); 
// extern double HeatH2I_II(float zred); 
// extern double HeatH2II_I(float zred); 
// extern double HeatH2II_II(float zred);  

#ifndef __CUDACC__
extern double k01(double);
extern double k02(double);
extern double k03(double);
extern double k04(double);
extern double k05(double);
extern double k06(double);
extern double k07(double);
extern double k08(double);
extern double k09(double);
extern double k10(double);
extern double k11(double);
extern double k12(double);
extern double k13(double);
extern double k14(double);
extern double k15(double);
extern double k16(double);
extern double k17(double);
extern double k18(double);
extern double k19(double);
extern double k20(double);
extern double k21(double);
#endif

extern double creh2(double);
extern double crehe2(double);
extern double crexihe2(double);
extern double crehe3(double);
extern double cioh1(double);
extern double ciohe1(double);
extern double ciohe2(double);
extern double cexh1(double);
extern double cexhe1(double);
extern double cexhe2(double);
extern double brems(double);
extern double compt(double, float);
extern double h2mol(double, double, double);

extern double calc_cooling_rate(struct prim_chem *, float, double, double);
extern double calc_heating_rate(struct prim_chem *, float, double, double);
extern double calc_heatcool_rate(struct prim_chem *, float, double, double);
extern void calc_ioneq(struct prim_chem*, double, double, float);

extern double csectHI(double);
extern double csectHeI(double);
extern double csectHeII(double);

#endif /* __CHEMISTRY__ */

#ifdef __cplusplus
}
#endif
