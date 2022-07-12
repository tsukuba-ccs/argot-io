#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <float.h>
#include <sys/time.h>
#include <sys/times.h>

#include "run_param.h"
#include "fluid.h"
#include "chemistry.h"
#include "mpi.h"
#include "prototype.h"

struct fluid_flux {
  float dens_flux;
  float momx_flux;
  float momy_flux;
  float momz_flux;
  float eneg_flux;
  float etrp_flux;

  float ndens_flux[NSPECIES];
};

#define MESH(ix,iy,iz) (mesh[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))])
#define MESH_MID(ix,iy,iz) (mesh_mid[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))])

float gamma_total(struct fluid_mesh*, struct run_param*);
void set_pad_x(struct pad_region*, struct fluid_mesh*,
	       struct fluid_mesh**, struct fluid_mesh**, struct fluid_mesh**,
	       struct fluid_mesh**, int, int, int, struct run_param*);
void set_pad_y(struct pad_region*, struct fluid_mesh*,
	       struct fluid_mesh**, struct fluid_mesh**, struct fluid_mesh**,
	       struct fluid_mesh**, int, int, int, struct run_param*);
void set_pad_z(struct pad_region*, struct fluid_mesh*, struct fluid_mesh**,
	       struct fluid_mesh**, struct fluid_mesh**,
	       struct fluid_mesh**, int, int, int, struct run_param*);

float gamma_total(struct fluid_mesh *mesh, struct run_param *this_run)
{
  return 1.666666;
}



void correct_ndens(struct fluid_mesh *mesh, float *ndens)
{
  int max_indx_H, max_indx_He; /* index of species with maximum number density */

  /* hydrogen */
  float ndens_max_H, nH;
  ndens_max_H = ndens[0];
  nH = 0.0;
  max_indx_H=0;
  for(int ispec=0;ispec<NSPECIES_HYDROGEN;ispec++) {
    if(ndens_max_H < ndens[ispec]) {
      max_indx_H = ispec;
      ndens_max_H = ndens[ispec];
    }
    if(ispec<3) {
      nH += ndens[ispec];
    }else{
      nH += 2.0*ndens[ispec];
    }
  }

  if(max_indx_H <3) {
    nH -= ndens[max_indx_H];
  }else{
    nH -= 2.0*ndens[max_indx_H];
  }

  float spec_im;
  spec_im = mesh->dens*XHYDROGEN - nH;

  if(spec_im > 0.0) {
    if(max_indx_H <3) {
      ndens[max_indx_H] = spec_im;
    }else{
      ndens[max_indx_H] = 0.5*spec_im;
    }
  } else {

    nH = 0.0;
    for(int ispec=0;ispec<NSPECIES_HYDROGEN;ispec++) {
      if(ispec<3) {
	nH += ndens[ispec];
      }else{
	nH += 2.0*ndens[ispec];
      }
    }
    
    nH = mesh->dens*XHYDROGEN/nH;
    for(int ispec=0;ispec<NSPECIES_HYDROGEN;ispec++) {
      ndens[ispec] = ndens[ispec] * nH;
    }
  }

#ifdef __HELIUM__  
  /* helium */
  int indx_He = NSPECIES_HYDROGEN;
  float ndens_max_He, nHe;
  ndens_max_He =  ndens[indx_He];
  nHe = 0.0;
  max_indx_He = indx_He;
  for(int ispec=indx_He;ispec<indx_He+NSPECIES_HELIUM;ispec++){
    if(ndens_max_He < ndens[ispec]) {
      max_indx_He = ispec;
      ndens_max_He = ndens[ispec];
    }
    nHe += ndens[ispec];
  }

  nHe -= ndens[max_indx_He];
  spec_im = mesh->dens*YHELIUM*0.25 - nHe;

  if(spec_im > 0.0) {
    ndens[max_indx_He] = spec_im;
  }else{
    
    nHe = 0.0;
    for(int ispec=indx_He;ispec<indx_He+NSPECIES_HELIUM;ispec++){
      nHe += ndens[ispec];
    }

    nHe = mesh->dens*YHELIUM*0.25/nHe;
    for(int ispec=indx_He;ispec<indx_He+NSPECIES_HELIUM;ispec++){
      ndens[ispec] = ndens[ispec] * nHe;
    }
  }

#endif

}


void set_fraction(struct fluid_mesh *mesh, float *ndens)
{
  mesh->chem.fHI  = ndens[0]/(mesh->dens*XHYDROGEN);
  mesh->chem.fHII = ndens[1]/(mesh->dens*XHYDROGEN);

  if(mesh->chem.fHI < FLT_MIN)  mesh->chem.fHI  = FLT_MIN;
  if(mesh->chem.fHII < FLT_MIN) mesh->chem.fHII = FLT_MIN;
  
  int ispec=2;
#ifdef __HYDROGEN_MOL__
  mesh->chem.fHM   = ndens[ispec++]/(mesh->dens*XHYDROGEN);
  mesh->chem.fH2I  = ndens[ispec++]/(mesh->dens*XHYDROGEN);
  mesh->chem.fH2II = ndens[ispec++]/(mesh->dens*XHYDROGEN);

  if(mesh->chem.fHM < FLT_MIN)   mesh->chem.fHM  = FLT_MIN;
  if(mesh->chem.fH2I < FLT_MIN)  mesh->chem.fH2I = FLT_MIN;
  if(mesh->chem.fH2II < FLT_MIN) mesh->chem.fH2II = FLT_MIN;
#endif /* __HYDROGEN_MOL__ */

#ifdef __HELIUM__
  mesh->chem.fHeI   = ndens[ispec++]/(mesh->dens*YHELIUM*0.25);
  mesh->chem.fHeII  = ndens[ispec++]/(mesh->dens*YHELIUM*0.25);
  mesh->chem.fHeIII = ndens[ispec++]/(mesh->dens*YHELIUM*0.25);

  if(mesh->chem.fHeI < FLT_MIN)   mesh->chem.fHeI   = FLT_MIN;
  if(mesh->chem.fHeII < FLT_MIN)  mesh->chem.fHeII  = FLT_MIN;
  if(mesh->chem.fHeIII < FLT_MIN) mesh->chem.fHeIII = FLT_MIN;
#endif /* __HELIUM__ */

  mesh->chem.felec = mesh->chem.fHII;
#ifdef __HYDROGEN_MOL__
  mesh->chem.felec -= mesh->chem.fHM;
  mesh->chem.felec += mesh->chem.fH2II;
#endif
#ifdef __HELIUM__
  mesh->chem.felec += HELIUM_FACT*(mesh->chem.fHeII+2.0*mesh->chem.fHeIII);
#endif

  if(mesh->chem.felec < FLT_MIN)   mesh->chem.felec = FLT_MIN;
}

void set_ndens(struct fluid_mesh *mesh, float *ndens) 
{
  ndens[0] = mesh->chem.fHI*mesh->dens*XHYDROGEN;
  ndens[1] = mesh->chem.fHII*mesh->dens*XHYDROGEN;

  int ispec=2;
#ifdef __HYDROGEN_MOL__
  ndens[ispec++] = mesh->chem.fHM*mesh->dens*XHYDROGEN;
  ndens[ispec++] = mesh->chem.fH2I*mesh->dens*XHYDROGEN;
  ndens[ispec++] = mesh->chem.fH2II*mesh->dens*XHYDROGEN;
#endif /* __HYDROGEN_MOL__ */

#ifdef __HELIUM__
  ndens[ispec++] = mesh->chem.fHeI*mesh->dens*YHELIUM*0.25;
  ndens[ispec++] = mesh->chem.fHeII*mesh->dens*YHELIUM*0.25;
  ndens[ispec++] = mesh->chem.fHeIII*mesh->dens*YHELIUM*0.25;
#endif /* __HELIUM__ */

}

void calc_flux_x(struct fluid_mesh *mesh, struct pad_region *pad, 
		 int ix, int iy, int iz,
		 struct fluid_flux *flux_plus, struct fluid_flux *flux_minus, 
		 struct run_param *this_run, int *shock)
{
  struct fluid_mesh *mesh_im2,*mesh_im1,*mesh_i,*mesh_ip1,*mesh_ip2;
  struct fluid_mesh reflect_lo, reflect_hi;

  float pres_im2,pres_im1,pres_i,pres_ip1,pres_ip2;
  float velx_im2,velx_im1,velx_i,velx_ip1,velx_ip2;
  float vely_im2,vely_im1,vely_i,vely_ip1,vely_ip2;
  float velz_im2,velz_im1,velz_i,velz_ip1,velz_ip2;
  float enth_im2,enth_im1,enth_i,enth_ip1,enth_ip2;
  float cs,E_thermal;

  float pres_L,dens_L,velx_L,vely_L,velz_L;
  float enth_L,etrp_L,cs_L,eneg_L;
  float pres_R,dens_R,velx_R,vely_R,velz_R;
  float enth_R,etrp_R,cs_R,eneg_R;
  float momx_L, momx_R, momy_L, momy_R, momz_L,momz_R;
  
  int ixm2, ixm1, ixp1, ixp2;

  float gamma_im2, gamma_im1, gamma_i, gamma_ip1, gamma_ip2;
  float gamm1_im2, gamm1_im1, gamm1_i, gamm1_ip1, gamm1_ip2;
  float gamma_L, gamma_R, gamm1_L, gamm1_R;
  float gamma, gamma1;
  
  float pres_star;
  float qL, qR;
  float SL, Sstar, SR;
  float dens_flux;
  
  float ndens_im2[NSPECIES]; 
  float ndens_im1[NSPECIES];
  float ndens_i[NSPECIES];
  float ndens_ip1[NSPECIES];
  float ndens_ip2[NSPECIES];
  float ndens_L[NSPECIES];
  float ndens_R[NSPECIES];

  mesh_i = &MESH(ix,iy,iz);
  
  set_pad_x(pad, mesh, &mesh_im2, &mesh_im1, &mesh_ip1, &mesh_ip2,
	    ix, iy, iz, this_run);
  
  /* adabatic index */
  gamma_im2 = gamma_total(mesh_im2, this_run); gamm1_im2 = gamma_im2-1.0;
  gamma_im1 = gamma_total(mesh_im1, this_run); gamm1_im1 = gamma_im1-1.0;
  gamma_i   = gamma_total(mesh_i  , this_run); gamm1_i   = gamma_i  -1.0;
  gamma_ip1 = gamma_total(mesh_ip1, this_run); gamm1_ip1 = gamma_ip1-1.0;
  gamma_ip2 = gamma_total(mesh_ip2, this_run); gamm1_ip2 = gamma_ip2-1.0;

  gamma  = gamma_i;
  gamma1 = gamm1_i;
  
  /*-- pressure --*/
  pres_im2 = fmaxf(gamm1_im2*mesh_im2->uene*mesh_im2->dens, FLUID_TINY);
  pres_im1 = fmaxf(gamm1_im1*mesh_im1->uene*mesh_im1->dens, FLUID_TINY);
  pres_i   = fmaxf(gamm1_i  *mesh_i->uene*mesh_i->dens    , FLUID_TINY);
  pres_ip1 = fmaxf(gamm1_ip1*mesh_ip1->uene*mesh_ip1->dens, FLUID_TINY);
  pres_ip2 = fmaxf(gamm1_ip2*mesh_ip2->uene*mesh_ip2->dens, FLUID_TINY);

  /*-- velocity --*/
  velx_im2 = mesh_im2->momx/mesh_im2->dens;
  velx_im1 = mesh_im1->momx/mesh_im1->dens;
  velx_i   = mesh_i->momx/mesh_i->dens;
  velx_ip1 = mesh_ip1->momx/mesh_ip1->dens;
  velx_ip2 = mesh_ip2->momx/mesh_ip2->dens;

  vely_im2 = mesh_im2->momy/mesh_im2->dens;
  vely_im1 = mesh_im1->momy/mesh_im1->dens;
  vely_i   = mesh_i->momy/mesh_i->dens;
  vely_ip1 = mesh_ip1->momy/mesh_ip1->dens;
  vely_ip2 = mesh_ip2->momy/mesh_ip2->dens;

  velz_im2 = mesh_im2->momz/mesh_im2->dens;
  velz_im1 = mesh_im1->momz/mesh_im1->dens;
  velz_i   = mesh_i->momz/mesh_i->dens;
  velz_ip1 = mesh_ip1->momz/mesh_ip1->dens;
  velz_ip2 = mesh_ip2->momz/mesh_ip2->dens;

  /* chemical species */
  set_ndens(mesh_im2, ndens_im2);
  set_ndens(mesh_im1, ndens_im1);
  set_ndens(mesh_i,   ndens_i);
  set_ndens(mesh_ip1, ndens_ip1);
  set_ndens(mesh_ip2, ndens_ip2);

  /* shock detection */
  *shock = 0;
  if(velx_ip1 - velx_im1 < 0.0) *shock = 1;
  if(fabsf(pres_ip1-pres_im1) > SHOCK_PRESS_THRESHOLD*pres_i) *shock = 1;

  /*-- numerical flux at i+1/2 --*/
  /* left values */
  dens_L = MUSCL_L(mesh_im1->dens, mesh_i->dens, mesh_ip1->dens);
  velx_L = MUSCL_L(velx_im1,       velx_i,       velx_ip1);
  vely_L = MUSCL_L(vely_im1,       vely_i,       vely_ip1);
  velz_L = MUSCL_L(velz_im1,       velz_i,       velz_ip1);
  pres_L = MUSCL_L(pres_im1,       pres_i,       pres_ip1);
  
  gamma_L = MUSCL_L(gamma_im1,gamma_i,gamma_ip1);
  gamm1_L = gamma_L-1.0;

  eneg_L = 0.5*dens_L*NORML2(velx_L,vely_L,velz_L) + pres_L/gamm1_L;
  enth_L = (eneg_L + pres_L)/dens_L;
  etrp_L = pres_L/pow(dens_L,gamm1_L);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_L[ispec] = MUSCL_L(ndens_im1[ispec], 
			     ndens_i[ispec],
			     ndens_ip1[ispec]);
  }
				      

  /* right values */
  dens_R = MUSCL_R(mesh_i->dens, mesh_ip1->dens, mesh_ip2->dens);
  velx_R = MUSCL_R(velx_i,       velx_ip1,       velx_ip2);
  vely_R = MUSCL_R(vely_i,       vely_ip1,       vely_ip2);
  velz_R = MUSCL_R(velz_i,       velz_ip1,       velz_ip2);
  pres_R = MUSCL_R(pres_i,       pres_ip1,       pres_ip2);

  gamma_R = MUSCL_R(gamma_i,gamma_ip1,gamma_ip2);
  gamm1_R = gamma_R-1.0;
    
  eneg_R = 0.5*dens_R*NORML2(velx_R,vely_R,velz_R) + pres_R/gamm1_R;
  enth_R = (eneg_R + pres_R)/dens_R;
  etrp_R = pres_R/pow(dens_R,gamm1_R);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_R[ispec] = MUSCL_R(ndens_i[ispec], 
			     ndens_ip1[ispec],
			     ndens_ip2[ispec]);
  }

  cs_L = sqrt(gamma_L*pres_L/dens_L);
  cs_R = sqrt(gamma_R*pres_R/dens_R);

#ifdef __ACOUSTIC_TYPE__
    /* acoustic type arpproximation */
    float pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#ifdef __TWO_RAREFACTION__ 
    /* two-rarefaction riemann solver */
    float zindx = 0.5*(gamma1)/gamma;
    pres_star = pow((cs_L+cs_R-0.5*gamma1*(velx_R-velx_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#ifdef __TWO_SCHOCK__
    /* two-shock riemann solver */
    float pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    float pres_0 = fmax(0.0, pres_pvrs);
    float gL, gR;
    gL = sqrt(2.0/((gamma_L+1.0)*dens_L)/(pres_0+(gamm1_L)/(gamma_L+1.0)*pres_L));
    gR = sqrt(2.0/((gamma_R+1.0)*dens_R)/(pres_0+(gamm1_R)/(gamma_R+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif

  if(pres_star <= pres_L) {
    qL = 1.0;
  }else{
    qL = sqrt(1.0+0.5*(gamma_L+1.0)/gamma_L*(pres_star/pres_L-1.0));
  }
  
  if(pres_star <= pres_R) {
    qR = 1.0;
  }else{
    qR = sqrt(1.0+0.5*(gamma_R+1.0)/gamma_R*(pres_star/pres_R-1.0));
  }

  SL = fmin(velx_L - cs_L*qL, velx_R - cs_R*qR);
  SR = fmax(velx_R + cs_R*qR, velx_L + cs_L*qL);
  Sstar = (pres_R-pres_L + dens_L*velx_L*(SL-velx_L) - dens_R*velx_R*(SR-velx_R))/(dens_L*(SL-velx_L)-dens_R*(SR-velx_R));

  momx_L = dens_L*velx_L;  momy_L = dens_L*vely_L;  momz_L = dens_L*velz_L;
  momx_R = dens_R*velx_R;  momy_R = dens_R*vely_R;  momz_R = dens_R*velz_R;
  
  /* numerical flux */
  if(0.0<=SL){

    dens_flux = dens_L*velx_L;

    flux_plus->dens_flux  = dens_flux;
    flux_plus->momx_flux  = pres_L + dens_L*SQR(velx_L);
    flux_plus->momy_flux  = dens_flux*vely_L;
    flux_plus->momz_flux  = dens_flux*velz_L;
    flux_plus->eneg_flux = (eneg_L+pres_L)*velx_L;
    flux_plus->etrp_flux = etrp_L*velx_L;

    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = velx_L*ndens_L[ispec];
    }
    
  }else if(SL<0.0 && 0.0<=Sstar) {
    float pres_L_star  = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
    float pres_R_star  = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_plus->dens_flux = Sstar*(SL*dens_L - dens_L*velx_L)/(SL-Sstar);
    flux_plus->momx_flux = (Sstar*(SL*momx_L - (pres_L+dens_L*SQR(velx_L))) + SL*pres_L_star)/(SL-Sstar);
    flux_plus->momy_flux = (Sstar*(SL*momy_L - (dens_L*velx_L*vely_L)))/(SL-Sstar);
    flux_plus->momz_flux = (Sstar*(SL*momz_L - (dens_L*velx_L*velz_L)))/(SL-Sstar);
    flux_plus->eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velx_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
    flux_plus->etrp_flux = Sstar*(SL*etrp_L - etrp_L*velx_L)/(SL-Sstar);
   
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = Sstar*(SL*ndens_L[ispec] - ndens_L[ispec]*velx_L)/(SL-Sstar);
    }
    
  }else if(Sstar < 0.0 && 0.0 <=SR) {
    float pres_L_star  = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
    float pres_R_star  = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_plus->dens_flux = Sstar*(SR*dens_R - dens_R*velx_R)/(SR-Sstar);
    flux_plus->momx_flux = (Sstar*(SR*momx_R - (pres_R+dens_R*SQR(velx_R))) + SR*pres_R_star)/(SR-Sstar);
    flux_plus->momy_flux = (Sstar*(SR*momy_R - (dens_R*velx_R*vely_R)))/(SR-Sstar);
    flux_plus->momz_flux = (Sstar*(SR*momz_R - (dens_R*velx_R*velz_R)))/(SR-Sstar);
    flux_plus->eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velx_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
    flux_plus->etrp_flux = Sstar*(SR*etrp_R - etrp_R*velx_R)/(SR-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = Sstar*(SR*ndens_R[ispec] - ndens_R[ispec]*velx_R)/(SR-Sstar);
    }
    
  }else if(SR<0.0) {

    dens_flux = dens_R*velx_R;
    
    flux_plus->dens_flux = dens_flux;
    flux_plus->momx_flux = pres_R + dens_R*SQR(velx_R);
    flux_plus->momy_flux = dens_flux*vely_R;
    flux_plus->momz_flux = dens_flux*velz_R;
    flux_plus->eneg_flux = (eneg_R+pres_R)*velx_R;
    flux_plus->etrp_flux = etrp_R*velx_R;
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = velx_R*ndens_R[ispec];
    }  
  }
  

  /*-- numerical flux at i-1/2 --*/
  /* left values */
  dens_L = MUSCL_L(mesh_im2->dens, mesh_im1->dens, mesh_i->dens);
  velx_L = MUSCL_L(velx_im2      , velx_im1      , velx_i      );
  vely_L = MUSCL_L(vely_im2      , vely_im1      , vely_i      );
  velz_L = MUSCL_L(velz_im2      , velz_im1      , velz_i      );
  pres_L = MUSCL_L(pres_im2      , pres_im1      , pres_i      );

  gamma_L = MUSCL_L(gamma_im2,gamma_im1,gamma_i);
  gamm1_L = gamma_L-1.0;
  
  eneg_L = 0.5*dens_L*NORML2(velx_L,vely_L,velz_L)+pres_L/gamm1_L;
  enth_L = (eneg_L + pres_L)/dens_L;
  etrp_L = pres_L/pow(dens_L,gamm1_L);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_L[ispec] = MUSCL_L(ndens_im2[ispec], 
			     ndens_im1[ispec],
			     ndens_i[ispec]);
  }

  /* right values */
  dens_R = MUSCL_R(mesh_im1->dens, mesh_i->dens, mesh_ip1->dens);
  velx_R = MUSCL_R(velx_im1      , velx_i      , velx_ip1      );
  vely_R = MUSCL_R(vely_im1      , vely_i      , vely_ip1      );
  velz_R = MUSCL_R(velz_im1      , velz_i      , velz_ip1      );
  pres_R = MUSCL_R(pres_im1      , pres_i      , pres_ip1      );

  gamma_R = MUSCL_R(gamma_im1,gamma_i,gamma_ip1);
  gamm1_R = gamma_R-1.0;
  
  eneg_R = 0.5*dens_R*NORML2(velx_R,vely_R,velz_R) + pres_R/gamm1_R;
  enth_R = (eneg_R + pres_R)/dens_R;
  etrp_R = pres_R/pow(dens_R,gamm1_R);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_R[ispec] = MUSCL_R(ndens_im1[ispec], 
			     ndens_i[ispec],
			     ndens_ip1[ispec]);
  }

  cs_L = sqrt(gamma*pres_L/dens_L);
  cs_R = sqrt(gamma*pres_R/dens_R);

#ifdef __ACOUSTIC_TYPE__
    /* acoustic type arpproximation */
    pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#ifdef __TWO_RAREFACTION__ 
    /* two-rarefaction riemann solver */
    zindx = 0.5*(gamma1)/gamma;
    pres_star = pow((cs_L+cs_R-0.5*gamma1*(velx_R-velx_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#ifdef __TWO_SCHOCK__
    /* two-shock riemann solver */
    pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(velx_R-velx_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_0 = fmax(0.0, pres_pvrs);
    gL = sqrt(2.0/((gamma_L+1.0)*dens_L)/(pres_0+(gamm1_L)/(gamma_L+1.0)*pres_L));
    gR = sqrt(2.0/((gamma_R+1.0)*dens_R)/(pres_0+(gamm1_R)/(gamma_R+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif
  
  if(pres_star <= pres_L) {
    qL = 1.0;
  }else{
    qL = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_L-1.0));
  }
    
  if(pres_star <= pres_R) {
    qR = 1.0;
  }else{
    qR = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_R-1.0));
  }

  SL = fmin(velx_L - cs_L*qL, velx_R - cs_R*qR);
  SR = fmax(velx_R + cs_R*qR, velx_L + cs_L*qL);
  Sstar = (pres_R-pres_L + dens_L*velx_L*(SL-velx_L) - dens_R*velx_R*(SR-velx_R))/(dens_L*(SL-velx_L)-dens_R*(SR-velx_R));

  momx_L = dens_L*velx_L;  momy_L = dens_L*vely_L;  momz_L = dens_L*velz_L;
  momx_R = dens_R*velx_R;  momy_R = dens_R*vely_R;  momz_R = dens_R*velz_R;
  
  /* numerical flux */
  if(0.0<=SL){

    dens_flux = dens_L*velx_L;

    flux_minus->dens_flux  = dens_flux;
    flux_minus->momx_flux  = pres_L + dens_L*SQR(velx_L);
    flux_minus->momy_flux  = dens_flux*vely_L;
    flux_minus->momz_flux  = dens_flux*velz_L;
    flux_minus->eneg_flux  = (eneg_L+pres_L)*velx_L;
    flux_minus->etrp_flux  = etrp_L*velx_L;

    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = velx_L*ndens_L[ispec];
    }
    
  }else if(SL<0.0 && 0.0<=Sstar) {
    float pres_L_star  = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
    float pres_R_star  = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_minus->dens_flux = Sstar*(SL*dens_L - dens_L*velx_L)/(SL-Sstar);
    flux_minus->momx_flux = (Sstar*(SL*momx_L - (pres_L+dens_L*SQR(velx_L))) + SL*pres_L_star)/(SL-Sstar);
    flux_minus->momy_flux = (Sstar*(SL*momy_L - (dens_L*velx_L*vely_L)))/(SL-Sstar);
    flux_minus->momz_flux = (Sstar*(SL*momz_L - (dens_L*velx_L*velz_L)))/(SL-Sstar);
    flux_minus->eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velx_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
    flux_minus->etrp_flux = Sstar*(SL*etrp_L - etrp_L*velx_L)/(SL-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = Sstar*(SL*ndens_L[ispec] - ndens_L[ispec]*velx_L)/(SL-Sstar);
    }
    
  }else if(Sstar < 0.0 && 0.0 <=SR) {
    float pres_L_star  = pres_L + dens_L*(SL-velx_L)*(Sstar-velx_L);
    float pres_R_star  = pres_R + dens_R*(SR-velx_R)*(Sstar-velx_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_minus->dens_flux = Sstar*(SR*dens_R - dens_R*velx_R)/(SR-Sstar);
    flux_minus->momx_flux = (Sstar*(SR*momx_R - (pres_R+dens_R*SQR(velx_R))) + SR*pres_R_star)/(SR-Sstar);
    flux_minus->momy_flux = (Sstar*(SR*momy_R - (dens_R*velx_R*vely_R)))/(SR-Sstar);
    flux_minus->momz_flux = (Sstar*(SR*momz_R - (dens_R*velx_R*velz_R)))/(SR-Sstar);
    flux_minus->eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velx_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
    flux_minus->etrp_flux = Sstar*(SR*etrp_R - etrp_R*velx_R)/(SR-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = Sstar*(SR*ndens_R[ispec] - ndens_R[ispec]*velx_R)/(SR-Sstar);
    }
    
  }else if(SR<0.0) {

    dens_flux = dens_R*velx_R;
    
    flux_minus->dens_flux = dens_flux;
    flux_minus->momx_flux = pres_R + dens_R*SQR(velx_R);
    flux_minus->momy_flux = dens_flux*vely_R;
    flux_minus->momz_flux = dens_flux*velz_R;
    flux_minus->eneg_flux = (eneg_R+pres_R)*velx_R;
    flux_minus->etrp_flux = etrp_R*velx_R;
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = velx_R*ndens_R[ispec];
    }  
  }

#if 0
  *dt_x = COURANT_FACT*this_run->delta_x/((vel_plus_iph - vel_minus_imh)+FLUID_TINY);
#endif

}

void calc_flux_y(struct fluid_mesh *mesh, struct pad_region *pad, 
		 int ix, int iy, int iz, 
		 struct fluid_flux *flux_plus, struct fluid_flux *flux_minus,
		 struct run_param *this_run, int *shock)
{
  struct fluid_mesh *mesh_im2,*mesh_im1,*mesh_i,*mesh_ip1,*mesh_ip2;
  struct fluid_mesh reflect_lo, reflect_hi;

  float pres_im2,pres_im1,pres_i,pres_ip1,pres_ip2;
  float velx_im2,velx_im1,velx_i,velx_ip1,velx_ip2;
  float vely_im2,vely_im1,vely_i,vely_ip1,vely_ip2;
  float velz_im2,velz_im1,velz_i,velz_ip1,velz_ip2;
  float enth_im2,enth_im1,enth_i,enth_ip1,enth_ip2;
  float cs,E_thermal;
  
  float pres_L,dens_L,velx_L,vely_L,velz_L;
  float enth_L,etrp_L,cs_L,eneg_L;
  float pres_R,dens_R,velx_R,vely_R,velz_R;
  float enth_R,etrp_R,cs_R,eneg_R;
  float momx_L, momx_R, momy_L, momy_R, momz_L,momz_R;
  
  int iym2, iym1, iyp1, iyp2;

  float gamma_im2, gamma_im1, gamma_i, gamma_ip1, gamma_ip2;
  float gamm1_im2, gamm1_im1, gamm1_i, gamm1_ip1, gamm1_ip2;
  float gamma_L, gamma_R, gamm1_L, gamm1_R;
  float gamma, gamma1;
  
  float pres_star;
  float qL, qR;
  float SL, Sstar, SR;
  float dens_flux;
  
  float ndens_im2[NSPECIES]; 
  float ndens_im1[NSPECIES];
  float ndens_i[NSPECIES];
  float ndens_ip1[NSPECIES];
  float ndens_ip2[NSPECIES];
  float ndens_L[NSPECIES];
  float ndens_R[NSPECIES];

  mesh_i = &MESH(ix,iy,iz);

  set_pad_y(pad, mesh, &mesh_im2, &mesh_im1, &mesh_ip1, &mesh_ip2,
	    ix, iy, iz, this_run);

  /* adabatic index */
  gamma_im2 = gamma_total(mesh_im2, this_run); gamm1_im2 = gamma_im2-1.0;
  gamma_im1 = gamma_total(mesh_im1, this_run); gamm1_im1 = gamma_im1-1.0;
  gamma_i   = gamma_total(mesh_i  , this_run); gamm1_i   = gamma_i  -1.0;
  gamma_ip1 = gamma_total(mesh_ip1, this_run); gamm1_ip1 = gamma_ip1-1.0;
  gamma_ip2 = gamma_total(mesh_ip2, this_run); gamm1_ip2 = gamma_ip2-1.0;

  gamma  = gamma_i;
  gamma1 = gamm1_i;
  
  /*-- pressure --*/
  pres_im2 = fmaxf(gamm1_im2*mesh_im2->uene*mesh_im2->dens, FLUID_TINY);
  pres_im1 = fmaxf(gamm1_im1*mesh_im1->uene*mesh_im1->dens, FLUID_TINY);
  pres_i   = fmaxf(gamm1_i  *mesh_i->uene*mesh_i->dens    , FLUID_TINY);
  pres_ip1 = fmaxf(gamm1_ip1*mesh_ip1->uene*mesh_ip1->dens, FLUID_TINY);
  pres_ip2 = fmaxf(gamm1_ip2*mesh_ip2->uene*mesh_ip2->dens, FLUID_TINY);

  /*-- velocity --*/
  velx_im2 = mesh_im2->momx/mesh_im2->dens;
  velx_im1 = mesh_im1->momx/mesh_im1->dens;
  velx_i   = mesh_i->momx/mesh_i->dens;
  velx_ip1 = mesh_ip1->momx/mesh_ip1->dens;
  velx_ip2 = mesh_ip2->momx/mesh_ip2->dens;

  vely_im2 = mesh_im2->momy/mesh_im2->dens;
  vely_im1 = mesh_im1->momy/mesh_im1->dens;
  vely_i   = mesh_i->momy/mesh_i->dens;
  vely_ip1 = mesh_ip1->momy/mesh_ip1->dens;
  vely_ip2 = mesh_ip2->momy/mesh_ip2->dens;

  velz_im2 = mesh_im2->momz/mesh_im2->dens;
  velz_im1 = mesh_im1->momz/mesh_im1->dens;
  velz_i   = mesh_i->momz/mesh_i->dens;
  velz_ip1 = mesh_ip1->momz/mesh_ip1->dens;
  velz_ip2 = mesh_ip2->momz/mesh_ip2->dens;

  /* chemical species */
  set_ndens(mesh_im2, ndens_im2);
  set_ndens(mesh_im1, ndens_im1);
  set_ndens(mesh_i,   ndens_i);
  set_ndens(mesh_ip1, ndens_ip1);
  set_ndens(mesh_ip2, ndens_ip2);

  /* shock detection */
  *shock = 0;
  if(vely_ip1 - vely_im1 < 0.0) *shock = 1;
  if(fabsf(pres_ip1-pres_im1) > SHOCK_PRESS_THRESHOLD*pres_i) *shock = 1;

  /*-- numerical flux at i+1/2 --*/

  /* left values */
  dens_L = MUSCL_L(mesh_im1->dens, mesh_i->dens, mesh_ip1->dens);
  velx_L = MUSCL_L(velx_im1, velx_i, velx_ip1);
  vely_L = MUSCL_L(vely_im1, vely_i, vely_ip1);
  velz_L = MUSCL_L(velz_im1, velz_i, velz_ip1);
  pres_L = MUSCL_L(pres_im1, pres_i, pres_ip1);

  gamma_L = MUSCL_L(gamma_im1,gamma_i,gamma_ip1);
  gamm1_L = gamma_L-1.0;

  eneg_L = 0.5*dens_L*NORML2(velx_L,vely_L,velz_L) + pres_L/gamm1_L;
  enth_L = (eneg_L + pres_L)/dens_L;
  etrp_L = pres_L/pow(dens_L,gamm1_L);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_L[ispec] = MUSCL_L(ndens_im1[ispec], 
			     ndens_i[ispec],
			     ndens_ip1[ispec]);
  }

  /* right values */
  dens_R = MUSCL_R(mesh_i->dens, mesh_ip1->dens, mesh_ip2->dens);
  velx_R = MUSCL_R(velx_i, velx_ip1, velx_ip2);
  vely_R = MUSCL_R(vely_i, vely_ip1, vely_ip2);
  velz_R = MUSCL_R(velz_i, velz_ip1, velz_ip2);
  pres_R = MUSCL_R(pres_i, pres_ip1, pres_ip2);

  gamma_R = MUSCL_R(gamma_i,gamma_ip1,gamma_ip2);
  gamm1_R = gamma_R-1.0;

  eneg_R = 0.5*dens_R*NORML2(velx_R,vely_R,velz_R) + pres_R/gamm1_R;
  enth_R = (eneg_R + pres_R)/dens_R;
  etrp_R = pres_R/pow(dens_R,gamm1_R);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_R[ispec] = MUSCL_R(ndens_i[ispec], 
			     ndens_ip1[ispec],
			     ndens_ip2[ispec]);
  }

  cs_L = sqrt(gamma*pres_L/dens_L);
  cs_R = sqrt(gamma*pres_R/dens_R);

#ifdef __ACOUSTIC_TYPE__
    /* acoustic type arpproximation */
    float pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(vely_R-vely_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#ifdef __TWO_RAREFACTION__ 
    /* two-rarefaction riemann solver */
    float zindx = 0.5*(gamma1)/gamma;
    pres_star = pow((cs_L+cs_R-0.5*gamma1*(vely_R-vely_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#ifdef __TWO_SCHOCK__
    /* two-shock riemann solver */
    float pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(vely_R-vely_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    float pres_0 = fmax(0.0, pres_pvrs);
    float gL, gR;
    gL = sqrt(2.0/((gamma_L+1.0)*dens_L)/(pres_0+(gamm1_L)/(gamma_L+1.0)*pres_L));
    gR = sqrt(2.0/((gamma_R+1.0)*dens_R)/(pres_0+(gamm1_R)/(gamma_R+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif

  if(pres_star <= pres_L) {
    qL = 1.0;
  }else{
    qL = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_L-1.0));
  }
  
  if(pres_star <= pres_R) {
    qR = 1.0;
  }else{
    qR = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_R-1.0));
  }

  SL = fmin(vely_L - cs_L*qL, vely_R - cs_R*qR);
  SR = fmax(vely_R + cs_R*qR, vely_L + cs_L*qL);
  Sstar = (pres_R-pres_L + dens_L*vely_L*(SL-vely_L) - dens_R*vely_R*(SR-vely_R))/(dens_L*(SL-vely_L)-dens_R*(SR-vely_R));

  momx_L = dens_L*velx_L;  momy_L = dens_L*vely_L;  momz_L = dens_L*velz_L;
  momx_R = dens_R*velx_R;  momy_R = dens_R*vely_R;  momz_R = dens_R*velz_R;
  
  /* numerical flux */
  if(0.0<=SL){

    dens_flux = dens_L*vely_L;
    
    flux_plus->dens_flux = dens_flux;
    flux_plus->momx_flux = dens_flux*velx_L;
    flux_plus->momy_flux = pres_L + dens_L*SQR(vely_L);
    flux_plus->momz_flux = dens_flux*velz_L;
    flux_plus->eneg_flux = (eneg_L+pres_L)*vely_L;
    flux_plus->etrp_flux = etrp_L*vely_L;

    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = vely_L*ndens_L[ispec];
    }
    
  }else if(SL<0.0 && 0.0<=Sstar) {
    float pres_L_star  = pres_L + dens_L*(SL-vely_L)*(Sstar-vely_L);
    float pres_R_star  = pres_R + dens_R*(SR-vely_R)*(Sstar-vely_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_plus->dens_flux = Sstar*(SL*dens_L - dens_L*vely_L)/(SL-Sstar);
    flux_plus->momx_flux = (Sstar*(SL*momx_L - (dens_L*vely_L*velx_L)))/(SL-Sstar);
    flux_plus->momy_flux = (Sstar*(SL*momy_L - (pres_L+dens_L*SQR(vely_L))) + SL*pres_L_star)/(SL-Sstar);
    flux_plus->momz_flux = (Sstar*(SL*momz_L - (dens_L*vely_L*velz_L)))/(SL-Sstar);
    flux_plus->eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*vely_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
    flux_plus->etrp_flux = Sstar*(SL*etrp_L - etrp_L*vely_L)/(SL-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = Sstar*(SL*ndens_L[ispec] - ndens_L[ispec]*vely_L)/(SL-Sstar);
    }
    
  }else if(Sstar < 0.0 && 0.0 <=SR) {
    float pres_L_star  = pres_L + dens_L*(SL-vely_L)*(Sstar-vely_L);
    float pres_R_star  = pres_R + dens_R*(SR-vely_R)*(Sstar-vely_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_plus->dens_flux = Sstar*(SR*dens_R - dens_R*vely_R)/(SR-Sstar);
    flux_plus->momx_flux = (Sstar*(SR*momx_R - (dens_R*vely_R*velx_R)))/(SR-Sstar);
    flux_plus->momy_flux = (Sstar*(SR*momy_R - (pres_R+dens_R*SQR(vely_R))) + SR*pres_R_star)/(SR-Sstar);
    flux_plus->momz_flux = (Sstar*(SR*momz_R - (dens_R*vely_R*velz_R)))/(SR-Sstar);
    flux_plus->eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*vely_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
    flux_plus->etrp_flux = Sstar*(SR*etrp_R - etrp_R*vely_R)/(SR-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = Sstar*(SR*ndens_R[ispec] - ndens_R[ispec]*vely_R)/(SR-Sstar);
    }
    
  }else if(SR<0.0) {

    dens_flux = dens_R*vely_R;
    
    flux_plus->dens_flux = dens_flux;
    flux_plus->momx_flux = dens_flux*velx_R;
    flux_plus->momy_flux = pres_R + dens_R*SQR(vely_R);
    flux_plus->momz_flux = dens_flux*velz_R;
    flux_plus->eneg_flux = (eneg_R+pres_R)*vely_R;
    flux_plus->etrp_flux = etrp_R*vely_R;
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = vely_R*ndens_R[ispec];
    }  
  }

  /*-- numerical flux at i-1/2 --*/
  /* left values */
  dens_L = MUSCL_L(mesh_im2->dens, mesh_im1->dens, mesh_i->dens);
  velx_L = MUSCL_L(velx_im2      , velx_im1      , velx_i      );
  vely_L = MUSCL_L(vely_im2      , vely_im1      , vely_i      );
  velz_L = MUSCL_L(velz_im2      , velz_im1      , velz_i      );
  pres_L = MUSCL_L(pres_im2      , pres_im1      , pres_i      );

  gamma_L = MUSCL_L(gamma_im2,gamma_im1,gamma_i);
  gamm1_L = gamma_L-1.0;

  eneg_L = 0.5*dens_L*NORML2(velx_L,vely_L,velz_L)+pres_L/gamm1_L;
  enth_L = (eneg_L + pres_L)/dens_L;
  etrp_L = pres_L/pow(dens_L,gamm1_L);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_L[ispec] = MUSCL_L(ndens_im2[ispec], 
			     ndens_im1[ispec],
			     ndens_i[ispec]);
  }

  /* right values */
  dens_R = MUSCL_R(mesh_im1->dens, mesh_i->dens, mesh_ip1->dens);
  velx_R = MUSCL_R(velx_im1      , velx_i      , velx_ip1      );
  vely_R = MUSCL_R(vely_im1      , vely_i      , vely_ip1      );
  velz_R = MUSCL_R(velz_im1      , velz_i      , velz_ip1      );
  pres_R = MUSCL_R(pres_im1      , pres_i      , pres_ip1      );

  gamma_R = MUSCL_R(gamma_im1,gamma_i,gamma_ip1);
  gamm1_R = gamma_R-1.0;

  eneg_R = 0.5*dens_R*NORML2(velx_R,vely_R,velz_R) + pres_R/gamm1_R;
  enth_R = (eneg_R + pres_R)/dens_R;
  etrp_R = pres_R/pow(dens_R,gamm1_R);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_R[ispec] = MUSCL_R(ndens_im1[ispec], 
			     ndens_i[ispec],
			     ndens_ip1[ispec]);
  }

  cs_L = sqrt(gamma*pres_L/dens_L);
  cs_R = sqrt(gamma*pres_R/dens_R);

#ifdef __ACOUSTIC_TYPE__
    /* acoustic type arpproximation */
    pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(vely_R-vely_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_star = fmax(0.0, pres_pvrs);
#endif
#ifdef __TWO_RAREFACTION__ 
    /* two-rarefaction riemann solver */
    zindx = 0.5*(gamma1)/gamma;
    pres_star = pow((cs_L+cs_R-0.5*gamma1*(vely_R-vely_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#ifdef __TWO_SCHOCK__
    /* two-shock riemann solver */
    pres_pvrs = 0.5*(pres_L+pres_R) - 
      0.5*(vely_R-vely_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
    pres_0 = fmax(0.0, pres_pvrs);
    gL = sqrt(2.0/((gamma_L+1.0)*dens_L)/(pres_0+(gamm1_L)/(gamma_L+1.0)*pres_L));
    gR = sqrt(2.0/((gamma_R+1.0)*dens_R)/(pres_0+(gamm1_R)/(gamma_R+1.0)*pres_R));

    pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif
  
  if(pres_star <= pres_L) {
    qL = 1.0;
  }else{
    qL = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_L-1.0));
  }
    
  if(pres_star <= pres_R) {
    qR = 1.0;
  }else{
    qR = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_R-1.0));
  }

  SL = fmin(vely_L - cs_L*qL, vely_R - cs_R*qR);
  SR = fmax(vely_R + cs_R*qR, vely_L + cs_L*qL);
  Sstar = (pres_R-pres_L + dens_L*vely_L*(SL-vely_L) - dens_R*vely_R*(SR-vely_R))/(dens_L*(SL-vely_L)-dens_R*(SR-vely_R));

  momx_L = dens_L*velx_L;  momy_L = dens_L*vely_L;  momz_L = dens_L*velz_L;
  momx_R = dens_R*velx_R;  momy_R = dens_R*vely_R;  momz_R = dens_R*velz_R;
  
  /* numerical flux */
  if(0.0<=SL){

    dens_flux = dens_L*vely_L;
    
    flux_minus->dens_flux = dens_flux;
    flux_minus->momx_flux = dens_flux*velx_L;
    flux_minus->momy_flux = pres_L + dens_L*SQR(vely_L);
    flux_minus->momz_flux = dens_flux*velz_L;
    flux_minus->eneg_flux = (eneg_L+pres_L)*vely_L;
    flux_minus->etrp_flux = etrp_L*vely_L;

    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = vely_L*ndens_L[ispec];
    }
    
  }else if(SL<0.0 && 0.0<=Sstar) {
    float pres_L_star  = pres_L + dens_L*(SL-vely_L)*(Sstar-vely_L);
    float pres_R_star  = pres_R + dens_R*(SR-vely_R)*(Sstar-vely_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_minus->dens_flux = Sstar*(SL*dens_L - dens_L*vely_L)/(SL-Sstar);
    flux_minus->momx_flux = (Sstar*(SL*momx_L - (dens_L*vely_L*velx_L)))/(SL-Sstar);
    flux_minus->momy_flux = (Sstar*(SL*momy_L - (pres_L+dens_L*SQR(vely_L))) + SL*pres_L_star)/(SL-Sstar);
    flux_minus->momz_flux = (Sstar*(SL*momz_L - (dens_L*vely_L*velz_L)))/(SL-Sstar);
    flux_minus->eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*vely_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
    flux_minus->etrp_flux = Sstar*(SL*etrp_L - etrp_L*vely_L)/(SL-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = Sstar*(SL*ndens_L[ispec] - ndens_L[ispec]*vely_L)/(SL-Sstar);
    }
    
  }else if(Sstar < 0.0 && 0.0 <=SR) {
    float pres_L_star  = pres_L + dens_L*(SL-vely_L)*(Sstar-vely_L);
    float pres_R_star  = pres_R + dens_R*(SR-vely_R)*(Sstar-vely_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_minus->dens_flux = Sstar*(SR*dens_R - dens_R*vely_R)/(SR-Sstar);
    flux_minus->momx_flux = (Sstar*(SR*momx_R - (dens_R*vely_R*velx_R)))/(SR-Sstar);
    flux_minus->momy_flux = (Sstar*(SR*momy_R - (pres_R+dens_R*SQR(vely_R))) + SR*pres_R_star)/(SR-Sstar);
    flux_minus->momz_flux = (Sstar*(SR*momz_R - (dens_R*vely_R*velz_R)))/(SR-Sstar);
    flux_minus->eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*vely_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
    flux_minus->etrp_flux = Sstar*(SR*etrp_R - etrp_R*vely_R)/(SR-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = Sstar*(SR*ndens_R[ispec] - ndens_R[ispec]*vely_R)/(SR-Sstar);
    }
    
  }else if(SR<0.0) {

    dens_flux = dens_R*vely_R;
    
    flux_minus->dens_flux = dens_flux;
    flux_minus->momx_flux = dens_flux*velx_R;
    flux_minus->momy_flux = pres_R + dens_R*SQR(vely_R);
    flux_minus->momz_flux = dens_flux*velz_R;
    flux_minus->eneg_flux = (eneg_R+pres_R)*vely_R;
    flux_minus->etrp_flux = etrp_R*vely_R;
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = vely_R*ndens_R[ispec];
    }  
  }
 

#if 0
  *dt_y = COURANT_FACT*this_run->delta_y/((vel_plus_iph - vel_minus_imh)+FLUID_TINY);
#endif

}

void calc_flux_z(struct fluid_mesh *mesh, struct pad_region *pad, 
		 int ix, int iy, int iz, 
		 struct fluid_flux *flux_plus, struct fluid_flux *flux_minus,
		 struct run_param *this_run, int *shock) 
{
  struct fluid_mesh *mesh_im2,*mesh_im1,*mesh_i,*mesh_ip1,*mesh_ip2;
  struct fluid_mesh reflect_lo, reflect_hi;

  float pres_im2,pres_im1,pres_i,pres_ip1,pres_ip2;
  float velx_im2,velx_im1,velx_i,velx_ip1,velx_ip2;
  float vely_im2,vely_im1,vely_i,vely_ip1,vely_ip2;
  float velz_im2,velz_im1,velz_i,velz_ip1,velz_ip2;
  float enth_im2,enth_im1,enth_i,enth_ip1,enth_ip2;
  float cs,E_thermal;

  float pres_L,dens_L,velx_L,vely_L,velz_L;
  float enth_L,etrp_L,cs_L,eneg_L;
  float pres_R,dens_R,velx_R,vely_R,velz_R;
  float enth_R,etrp_R,cs_R,eneg_R;
  float momx_L, momx_R, momy_L, momy_R, momz_L,momz_R;
  
  int izm2, izm1, izp1, izp2;

  float gamma_im2, gamma_im1, gamma_i, gamma_ip1, gamma_ip2;
  float gamm1_im2, gamm1_im1, gamm1_i, gamm1_ip1, gamm1_ip2;
  float gamma_L, gamma_R, gamm1_L, gamm1_R;
  float gamma, gamma1;
  
  float pres_star;
  float qL, qR;
  float SL, Sstar, SR;
  float dens_flux;
  
  float ndens_im2[NSPECIES]; 
  float ndens_im1[NSPECIES];
  float ndens_i[NSPECIES];
  float ndens_ip1[NSPECIES];
  float ndens_ip2[NSPECIES];
  float ndens_L[NSPECIES];
  float ndens_R[NSPECIES];

  mesh_i = &MESH(ix,iy,iz);

  set_pad_z(pad, mesh, &mesh_im2, &mesh_im1, &mesh_ip1, &mesh_ip2,
	    ix, iy, iz, this_run);

  /* adabatic index */
  gamma_im2 = gamma_total(mesh_im2, this_run); gamm1_im2 = gamma_im2-1.0;
  gamma_im1 = gamma_total(mesh_im1, this_run); gamm1_im1 = gamma_im1-1.0;
  gamma_i   = gamma_total(mesh_i  , this_run); gamm1_i   = gamma_i  -1.0;
  gamma_ip1 = gamma_total(mesh_ip1, this_run); gamm1_ip1 = gamma_ip1-1.0;
  gamma_ip2 = gamma_total(mesh_ip2, this_run); gamm1_ip2 = gamma_ip2-1.0;

  gamma  = gamma_i;
  gamma1 = gamm1_i;
  
  /*-- pressure --*/
  pres_im2 = fmaxf(gamm1_im2*mesh_im2->uene*mesh_im2->dens, FLUID_TINY);
  pres_im1 = fmaxf(gamm1_im1*mesh_im1->uene*mesh_im1->dens, FLUID_TINY);
  pres_i   = fmaxf(gamm1_i  *mesh_i->uene*mesh_i->dens    , FLUID_TINY);
  pres_ip1 = fmaxf(gamm1_ip1*mesh_ip1->uene*mesh_ip1->dens, FLUID_TINY);
  pres_ip2 = fmaxf(gamm1_ip2*mesh_ip2->uene*mesh_ip2->dens, FLUID_TINY);

  /*-- velocity --*/
  velx_im2 = mesh_im2->momx/mesh_im2->dens;
  velx_im1 = mesh_im1->momx/mesh_im1->dens;
  velx_i   = mesh_i->momx/mesh_i->dens;
  velx_ip1 = mesh_ip1->momx/mesh_ip1->dens;
  velx_ip2 = mesh_ip2->momx/mesh_ip2->dens;

  vely_im2 = mesh_im2->momy/mesh_im2->dens;
  vely_im1 = mesh_im1->momy/mesh_im1->dens;
  vely_i   = mesh_i->momy/mesh_i->dens;
  vely_ip1 = mesh_ip1->momy/mesh_ip1->dens;
  vely_ip2 = mesh_ip2->momy/mesh_ip2->dens;

  velz_im2 = mesh_im2->momz/mesh_im2->dens;
  velz_im1 = mesh_im1->momz/mesh_im1->dens;
  velz_i   = mesh_i->momz/mesh_i->dens;
  velz_ip1 = mesh_ip1->momz/mesh_ip1->dens;
  velz_ip2 = mesh_ip2->momz/mesh_ip2->dens;

  /* chemical species */
  set_ndens(mesh_im2, ndens_im2);
  set_ndens(mesh_im1, ndens_im1);
  set_ndens(mesh_i,   ndens_i);
  set_ndens(mesh_ip1, ndens_ip1);
  set_ndens(mesh_ip2, ndens_ip2);

  /* shock detection */
  *shock = 0;
  if(velz_ip1 - velz_im1 < 0.0) *shock = 1;
  if(fabsf(pres_ip1-pres_im1) > SHOCK_PRESS_THRESHOLD*pres_i) *shock = 1;

  /*-- numerical flux at i+1/2 --*/

  /* left values */
  dens_L = MUSCL_L(mesh_im1->dens, mesh_i->dens, mesh_ip1->dens);
  velx_L = MUSCL_L(velx_im1, velx_i, velx_ip1);
  vely_L = MUSCL_L(vely_im1, vely_i, vely_ip1);
  velz_L = MUSCL_L(velz_im1, velz_i, velz_ip1);
  pres_L = MUSCL_L(pres_im1, pres_i, pres_ip1);

  gamma_L = MUSCL_L(gamma_im1,gamma_i,gamma_ip1);
  gamm1_L = gamma_L-1.0;
  
  eneg_L = 0.5*dens_L*NORML2(velx_L,vely_L,velz_L) + pres_L/gamm1_L;
  enth_L = (eneg_L + pres_L)/dens_L;
  etrp_L = pres_L/pow(dens_L,gamm1_L);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_L[ispec] = MUSCL_L(ndens_im1[ispec], 
			     ndens_i[ispec],
			     ndens_ip1[ispec]);
  }

  /* right values */
  dens_R = MUSCL_R(mesh_i->dens, mesh_ip1->dens, mesh_ip2->dens);
  velx_R = MUSCL_R(velx_i, velx_ip1, velx_ip2);
  vely_R = MUSCL_R(vely_i, vely_ip1, vely_ip2);
  velz_R = MUSCL_R(velz_i, velz_ip1, velz_ip2);
  pres_R = MUSCL_R(pres_i, pres_ip1, pres_ip2);

  gamma_R = MUSCL_R(gamma_i,gamma_ip1,gamma_ip2);
  gamm1_R = gamma_R-1.0;

  eneg_R = 0.5*dens_R*NORML2(velx_R,vely_R,velz_R) + pres_R/gamm1_R;
  enth_R = (eneg_R + pres_R)/dens_R;
  etrp_R = pres_R/pow(dens_R,gamm1_R);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_R[ispec] = MUSCL_R(ndens_i[ispec], 
			     ndens_ip1[ispec],
			     ndens_ip2[ispec]);
  }

  cs_L = sqrt(gamma*pres_L/dens_L);
  cs_R = sqrt(gamma*pres_R/dens_R);

#ifdef __ACOUSTIC_TYPE__
  /* acoustic type arpproximation */
  float pres_pvrs = 0.5*(pres_L+pres_R) - 
    0.5*(velz_R-velz_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
  pres_star = fmax(0.0, pres_pvrs);
#endif
#ifdef __TWO_RAREFACTION__ 
  /* two-rarefaction riemann solver */
  float zindx = 0.5*(gamma1)/gamma;
  pres_star = pow((cs_L+cs_R-0.5*gamma1*(velz_R-velz_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#ifdef __TWO_SCHOCK__
  /* two-shock riemann solver */
  float pres_pvrs = 0.5*(pres_L+pres_R) - 
    0.5*(velz_R-velz_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
  float pres_0 = fmax(0.0, pres_pvrs);
  float gL, gR;
  gL = sqrt(2.0/((gamma_L+1.0)*dens_L)/(pres_0+(gamm1_L)/(gamma_L+1.0)*pres_L));
  gR = sqrt(2.0/((gamma_R+1.0)*dens_R)/(pres_0+(gamm1_R)/(gamma_R+1.0)*pres_R));
  
  pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif
  
  if(pres_star <= pres_L) {
    qL = 1.0;
  }else{
    qL = sqrtf(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_L-1.0));
  }
  
  if(pres_star <= pres_R) {
    qR = 1.0;
  }else{
    qR = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_R-1.0));
  }

  SL = fmin(velz_L - cs_L*qL, velz_R - cs_R*qR);
  SR = fmax(velz_R + cs_R*qR, velz_L + cs_L*qL);
  Sstar = (pres_R-pres_L + dens_L*velz_L*(SL-velz_L) - dens_R*velz_R*(SR-velz_R))/(dens_L*(SL-velz_L)-dens_R*(SR-velz_R));

  momx_L = dens_L*velx_L;  momy_L = dens_L*vely_L;  momz_L = dens_L*velz_L;
  momx_R = dens_R*velx_R;  momy_R = dens_R*vely_R;  momz_R = dens_R*velz_R;
  
  /* numerical flux */
  if(0.0<=SL){

    dens_flux = dens_L*velz_L;
    
    flux_plus->dens_flux = dens_flux;
    flux_plus->momx_flux = dens_flux*velx_L;
    flux_plus->momy_flux = dens_flux*vely_L;
    flux_plus->momz_flux = pres_L + dens_L*SQR(velz_L);
    flux_plus->eneg_flux = (eneg_L+pres_L)*velz_L;
    flux_plus->etrp_flux = etrp_L*velz_L;

    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = velz_L*ndens_L[ispec];
    }
    
  }else if(SL<0.0 && 0.0<=Sstar) {
    float pres_L_star  = pres_L + dens_L*(SL-velz_L)*(Sstar-velz_L);
    float pres_R_star  = pres_R + dens_R*(SR-velz_R)*(Sstar-velz_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_plus->dens_flux = Sstar*(SL*dens_L - dens_L*velz_L)/(SL-Sstar);
    flux_plus->momx_flux = (Sstar*(SL*momx_L - (dens_L*velz_L*velx_L)))/(SL-Sstar);
    flux_plus->momy_flux = (Sstar*(SL*momy_L - (dens_L*velz_L*vely_L)))/(SL-Sstar);
    flux_plus->momz_flux = (Sstar*(SL*momz_L - (pres_L+dens_L*SQR(velz_L))) + SL*pres_L_star)/(SL-Sstar);
    flux_plus->eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velz_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
    flux_plus->etrp_flux = Sstar*(SL*etrp_L - etrp_L*velz_L)/(SL-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = Sstar*(SL*ndens_L[ispec] - ndens_L[ispec]*velz_L)/(SL-Sstar);
    }
    
  }else if(Sstar < 0.0 && 0.0 <=SR) {
    float pres_L_star  = pres_L + dens_L*(SL-velz_L)*(Sstar-velz_L);
    float pres_R_star  = pres_R + dens_R*(SR-velz_R)*(Sstar-velz_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_plus->dens_flux = Sstar*(SR*dens_R - dens_R*velz_R)/(SR-Sstar);
    flux_plus->momx_flux = (Sstar*(SR*momx_R - (dens_R*velz_R*velx_R)))/(SR-Sstar);
    flux_plus->momy_flux = (Sstar*(SR*momy_R - (dens_R*velz_R*vely_R)))/(SR-Sstar);
    flux_plus->momz_flux = (Sstar*(SR*momz_R - (pres_R+dens_R*SQR(velz_R))) + SR*pres_R_star)/(SR-Sstar);
    flux_plus->eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velz_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
    flux_plus->etrp_flux = Sstar*(SR*etrp_R - etrp_R*velz_R)/(SR-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = Sstar*(SR*ndens_R[ispec] - ndens_R[ispec]*velz_R)/(SR-Sstar);
    }
    
  }else if(SR<0.0) {

    dens_flux = dens_R*velz_R;
    
    flux_plus->dens_flux = dens_flux;
    flux_plus->momx_flux = dens_flux*velx_R;
    flux_plus->momy_flux = dens_flux*vely_R;
    flux_plus->momz_flux = pres_R + dens_R*SQR(velz_R);
    flux_plus->eneg_flux = (eneg_R+pres_R)*velz_R;
    flux_plus->etrp_flux = etrp_R*velz_R;
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_plus->ndens_flux[ispec] = velz_R*ndens_R[ispec];
    }  
  }
  
  /*-- numerical flux at i-1/2 --*/

  /* left values */
  dens_L = MUSCL_L(mesh_im2->dens, mesh_im1->dens, mesh_i->dens);
  velx_L = MUSCL_L(velx_im2      , velx_im1      , velx_i      );
  vely_L = MUSCL_L(vely_im2      , vely_im1      , vely_i      );
  velz_L = MUSCL_L(velz_im2      , velz_im1      , velz_i      );
  pres_L = MUSCL_L(pres_im2      , pres_im1      , pres_i      );

  gamma_L = MUSCL_L(gamma_im2,gamma_im1,gamma_i);
  gamm1_L = gamma_L-1.0;

  eneg_L = 0.5*dens_L*NORML2(velx_L,vely_L,velz_L)+pres_L/gamm1_L;
  enth_L = (eneg_L + pres_L)/dens_L;
  etrp_L = pres_L/pow(dens_L,gamm1_L);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_L[ispec] = MUSCL_L(ndens_im2[ispec], 
			     ndens_im1[ispec],
			     ndens_i[ispec]);
  }

  /* right values */
  dens_R = MUSCL_R(mesh_im1->dens, mesh_i->dens, mesh_ip1->dens);
  velx_R = MUSCL_R(velx_im1      , velx_i      , velx_ip1      );
  vely_R = MUSCL_R(vely_im1      , vely_i      , vely_ip1      );
  velz_R = MUSCL_R(velz_im1      , velz_i      , velz_ip1      );
  pres_R = MUSCL_R(pres_im1      , pres_i      , pres_ip1      );

  gamma_R = MUSCL_R(gamma_im1,gamma_i,gamma_ip1);
  gamm1_R = gamma_R-1.0;

  eneg_R = 0.5*dens_R*NORML2(velx_R,vely_R,velz_R) + pres_R/gamm1_R;
  enth_R = (eneg_R + pres_R)/dens_R;
  etrp_R = pres_R/pow(dens_R,gamm1_R);

  for(int ispec=0;ispec<NSPECIES;ispec++) {
    ndens_R[ispec] = MUSCL_R(ndens_im1[ispec], 
			     ndens_i[ispec],
			     ndens_ip1[ispec]);
  }

  cs_L = sqrt(gamma*pres_L/dens_L);
  cs_R = sqrt(gamma*pres_R/dens_R);

#ifdef __ACOUSTIC_TYPE__
  /* acoustic type arpproximation */
  pres_pvrs = 0.5*(pres_L+pres_R) - 
    0.5*(velz_R-velz_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
  pres_star = fmax(0.0, pres_pvrs);
#endif
#ifdef __TWO_RAREFACTION__ 
  /* two-rarefaction riemann solver */
  zindx = 0.5*(gamma1)/gamma;
  pres_star = pow((cs_L+cs_R-0.5*gamma1*(velz_R-velz_L))/(cs_L/pow(pres_L,zindx) + cs_R/pow(pres_R,zindx)),1.0/zindx);
#endif
#ifdef __TWO_SCHOCK__
  /* two-shock riemann solver */
  pres_pvrs = 0.5*(pres_L+pres_R) - 
    0.5*(velz_R-velz_L)*0.5*(dens_L+dens_R)*0.5*(cs_L+cs_R);
  pres_0 = fmax(0.0, pres_pvrs);
  gL = sqrt(2.0/((gamma_L+1.0)*dens_L)/(pres_0+(gamm1_L)/(gamma_L+1.0)*pres_L));
  gR = sqrt(2.0/((gamma_R+1.0)*dens_R)/(pres_0+(gamm1_R)/(gamma_R+1.0)*pres_R));
  
  pres_star = (gL*pres_L + gR*pres_R - (velx_R-velx_L))/(gL+gR);
#endif
  
  if(pres_star <= pres_L) {
    qL = 1.0;
  }else{
    qL = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_L-1.0));
  }
    
  if(pres_star <= pres_R) {
    qR = 1.0;
  }else{
    qR = sqrt(1.0+0.5*(gamma+1.0)/gamma*(pres_star/pres_R-1.0));
  }

  SL = fmin(velz_L - cs_L*qL, velz_R - cs_R*qR);
  SR = fmax(velz_R + cs_R*qR, velz_L + cs_L*qL);
  Sstar = (pres_R-pres_L + dens_L*velz_L*(SL-velz_L) - dens_R*velz_R*(SR-velz_R))/(dens_L*(SL-velz_L)-dens_R*(SR-velz_R));

  momx_L = dens_L*velx_L;  momy_L = dens_L*vely_L;  momz_L = dens_L*velz_L;
  momx_R = dens_R*velx_R;  momy_R = dens_R*vely_R;  momz_R = dens_R*velz_R;
  
  /* numerical flux */
  if(0.0<=SL){

    dens_flux = dens_L*velz_L;
    
    flux_minus->dens_flux = dens_flux;
    flux_minus->momx_flux = dens_flux*velx_L;
    flux_minus->momy_flux = dens_flux*vely_L;
    flux_minus->momz_flux = pres_L + dens_L*SQR(velz_L);
    flux_minus->eneg_flux = (eneg_L+pres_L)*velz_L;
    flux_minus->etrp_flux = etrp_L*velz_L;

    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = velz_L*ndens_L[ispec];
    }
    
  }else if(SL<0.0 && 0.0<=Sstar) {
    float pres_L_star  = pres_L + dens_L*(SL-velz_L)*(Sstar-velz_L);
    float pres_R_star  = pres_R + dens_R*(SR-velz_R)*(Sstar-velz_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_minus->dens_flux = Sstar*(SL*dens_L - dens_L*velz_L)/(SL-Sstar);
    flux_minus->momx_flux = (Sstar*(SL*momx_L - (dens_L*velz_L*velx_L)))/(SL-Sstar);
    flux_minus->momy_flux = (Sstar*(SL*momy_L - (dens_L*velz_L*vely_L)))/(SL-Sstar);
    flux_minus->momz_flux = (Sstar*(SL*momz_L - (pres_L+dens_L*SQR(velz_L))) + SL*pres_L_star)/(SL-Sstar);
    flux_minus->eneg_flux = (Sstar*(SL*eneg_L - (eneg_L+pres_L)*velz_L) + SL*Sstar*pres_L_star)/(SL-Sstar);
    flux_minus->etrp_flux = Sstar*(SL*etrp_L - etrp_L*velz_L)/(SL-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = Sstar*(SL*ndens_L[ispec] - ndens_L[ispec]*velz_L)/(SL-Sstar);
    }
    
  }else if(Sstar < 0.0 && 0.0 <=SR) {
    float pres_L_star  = pres_L + dens_L*(SL-velz_L)*(Sstar-velz_L);
    float pres_R_star  = pres_R + dens_R*(SR-velz_R)*(Sstar-velz_R);
    float pres_LR_star = 0.5*(pres_L_star + pres_R_star);

    flux_minus->dens_flux = Sstar*(SR*dens_R - dens_R*velz_R)/(SR-Sstar);
    flux_minus->momx_flux = (Sstar*(SR*momx_R - (dens_R*velz_R*velx_R)))/(SR-Sstar);
    flux_minus->momy_flux = (Sstar*(SR*momy_R - (dens_R*velz_R*vely_R)))/(SR-Sstar);
    flux_minus->momz_flux = (Sstar*(SR*momz_R - (pres_R+dens_R*SQR(velz_R))) + SR*pres_R_star)/(SR-Sstar);
    flux_minus->eneg_flux = (Sstar*(SR*eneg_R - (eneg_R+pres_R)*velz_R) + SR*Sstar*pres_R_star)/(SR-Sstar);
    flux_minus->etrp_flux = Sstar*(SR*etrp_R - etrp_R*velz_R)/(SR-Sstar);
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = Sstar*(SR*ndens_R[ispec] - ndens_R[ispec]*velz_R)/(SR-Sstar);
    }
    
  }else if(SR<0.0) {

    dens_flux = dens_R*velz_R;
    
    flux_minus->dens_flux = dens_flux;
    flux_minus->momx_flux = dens_flux*velx_R;
    flux_minus->momy_flux = dens_flux*vely_R;
    flux_minus->momz_flux = pres_R + dens_R*SQR(velz_R);
    flux_minus->eneg_flux = (eneg_R+pres_R)*velz_R;
    flux_minus->etrp_flux = etrp_R*velz_R;
    
    for(int ispec=0;ispec<NSPECIES;ispec++) {
      flux_minus->ndens_flux[ispec] = velz_R*ndens_R[ispec];
    }  
  }

#if 0
  *dt_z = COURANT_FACT*this_run->delta_x/((vel_plus_iph - vel_minus_imh)+FLUID_TINY);
#endif

}


int high_mach_detector(struct fluid_mesh *mesh, struct run_param *this_run, float gamma)
{
  float therm_eneg;
  float gam1;

  gam1 = gamma-1.0;
  //  therm_eneg = (mesh->etrp)*powf(mesh->dens,gam1)/gam1;
  therm_eneg = mesh->dens*mesh->uene;
  
  if(therm_eneg < HIGH_MACH_THRESHOLD*mesh->eneg) {
    return 1;
  }else{
    return 0;
  }
}



#ifdef __GRAVITY__
void fluid_integrate(struct fluid_mesh *mesh, struct run_param *this_run,  
		     struct mpi_param *this_mpi, 
                     struct fftw_mpi_param *this_fftw_mpi,
		     float *green_func, float *dens, float dtime)
#else /* !__GRAVITY__ */
void fluid_integrate(struct fluid_mesh *mesh, struct run_param *this_run, 
		     struct mpi_param *this_mpi, float dtime)
#endif
{
  float dtdx, dtdy, dtdz;

  struct fluid_mesh *mesh_mid;
  struct pad_region pad, pad_mid;

  struct timeval start_tv, end_tv;
  struct tms start_tms, end_tms;

  struct timeval start_tv_part, end_tv_part;
  struct tms start_tms_part, end_tms_part;
  
  start_timing(&start_tv, &start_tms);

  dtdx = dtime/this_run->delta_x;
  dtdy = dtime/this_run->delta_y;
  dtdz = dtime/this_run->delta_z;

  mesh_mid = (struct fluid_mesh *) malloc(sizeof(struct fluid_mesh)*NMESH_LOCAL);
  
  init_pad_region(&pad, this_run);
  init_pad_region(&pad_mid, this_run);

  update_pad_region(mesh, &pad, this_run, this_mpi);

  start_timing(&start_tv_part, &start_tms_part);
  
#pragma omp parallel for schedule(auto) collapse(3)
  for(int ix=0;ix<NMESH_X_LOCAL;ix++) {
    for(int iy=0;iy<NMESH_Y_LOCAL;iy++) {
      for(int iz=0;iz<NMESH_Z_LOCAL;iz++) {
	
	struct fluid_mesh *target_mesh, *dest_mesh;
	struct fluid_flux flux_x_plus, flux_x_minus;
	struct fluid_flux flux_y_plus, flux_y_minus;
	struct fluid_flux flux_z_plus, flux_z_minus;
      
	float gamma, gam1;
	float therm_eneg, kinetic_eneg;
	float dt_x, dt_y, dt_z;

	int shock_x, shock_y, shock_z;
      
	target_mesh = &MESH(ix,iy,iz);
	dest_mesh = &MESH_MID(ix,iy,iz);
      
	target_mesh->under_shock = 0;
      
	calc_flux_x(mesh, &pad, ix, iy, iz, &flux_x_plus, &flux_x_minus, this_run, &shock_x);
	calc_flux_y(mesh, &pad, ix, iy, iz, &flux_y_plus, &flux_y_minus, this_run, &shock_y);
	calc_flux_z(mesh, &pad, ix, iy, iz, &flux_z_plus, &flux_z_minus, this_run, &shock_z);

	gamma = gamma_total(target_mesh, this_run);
	gam1 = gamma-1.0;
      
	target_mesh->under_shock |= shock_x;
	target_mesh->under_shock |= shock_y;
	target_mesh->under_shock |= shock_z;
	target_mesh->high_mach = high_mach_detector(target_mesh, this_run, gamma);
      
	dest_mesh->dens = target_mesh->dens 
	  - dtdx*(flux_x_plus.dens_flux - flux_x_minus.dens_flux)
	  - dtdy*(flux_y_plus.dens_flux - flux_y_minus.dens_flux)
	  - dtdz*(flux_z_plus.dens_flux - flux_z_minus.dens_flux);
	
	dest_mesh->momx = target_mesh->momx 
	  - dtdx*(flux_x_plus.momx_flux - flux_x_minus.momx_flux)
	  - dtdy*(flux_y_plus.momx_flux - flux_y_minus.momx_flux)
	  - dtdz*(flux_z_plus.momx_flux - flux_z_minus.momx_flux);

	dest_mesh->momy = target_mesh->momy 
	  - dtdx*(flux_x_plus.momy_flux - flux_x_minus.momy_flux)
	  - dtdy*(flux_y_plus.momy_flux - flux_y_minus.momy_flux)
	  - dtdz*(flux_z_plus.momy_flux - flux_z_minus.momy_flux);

	dest_mesh->momz = target_mesh->momz 
	  - dtdx*(flux_x_plus.momz_flux - flux_x_minus.momz_flux)
	  - dtdy*(flux_y_plus.momz_flux - flux_y_minus.momz_flux)
	  - dtdz*(flux_z_plus.momz_flux - flux_z_minus.momz_flux);

	dest_mesh->eneg = target_mesh->eneg 
	  - dtdx*(flux_x_plus.eneg_flux - flux_x_minus.eneg_flux)
	  - dtdy*(flux_y_plus.eneg_flux - flux_y_minus.eneg_flux)
	  - dtdz*(flux_z_plus.eneg_flux - flux_z_minus.eneg_flux);

	dest_mesh->under_shock = target_mesh->under_shock;
	dest_mesh->high_mach   = target_mesh->high_mach;

	kinetic_eneg = 0.5*NORML2(dest_mesh->momx, 
				  dest_mesh->momy, 
				  dest_mesh->momz)/dest_mesh->dens;

	float wmol = WMOL(target_mesh->chem);
	
	if( (target_mesh->high_mach) && !(target_mesh->under_shock) ) {
	  /* update the modified entropy using the master equation */

	  dest_mesh->etrp = target_mesh->etrp 
	    - dtdx*(flux_x_plus.etrp_flux - flux_x_minus.etrp_flux)
	    - dtdy*(flux_y_plus.etrp_flux - flux_y_minus.etrp_flux)
	    - dtdz*(flux_z_plus.etrp_flux - flux_z_minus.etrp_flux);

	  if(dest_mesh->etrp < FLUID_TINY) dest_mesh->etrp = FLUID_TINY;
       
	  therm_eneg = dest_mesh->etrp*powf(dest_mesh->dens,gam1)/gam1;

	  if(therm_eneg/dest_mesh->dens*this_run->uenetok*wmol < 10.0) { // 10 K
	    therm_eneg = dest_mesh->dens/(this_run->uenetok*wmol);
	  }
	  
	  dest_mesh->eneg = therm_eneg + kinetic_eneg;
	  dest_mesh->uene = therm_eneg/dest_mesh->dens;
	}else{
	  therm_eneg = dest_mesh->eneg - kinetic_eneg;

	  if(therm_eneg/dest_mesh->dens*this_run->uenetok*wmol < 10.0) { // 10 K
	    therm_eneg = dest_mesh->dens/(this_run->uenetok*wmol);
	  }

	  /*
	  if(therm_eneg < FLUID_TINY)  {
	    therm_eneg = FLUID_TINY;
	    dest_mesh->eneg = therm_eneg + kinetic_eneg;
	  }
	  */
	  dest_mesh->uene = therm_eneg/dest_mesh->dens;
	  dest_mesh->etrp = gam1*therm_eneg/powf(dest_mesh->dens,gam1);
	}

	dest_mesh->eneg = fmaxf(dest_mesh->eneg, FLUID_TINY);
	dest_mesh->uene = fmaxf(dest_mesh->uene, FLUID_TINY);

	
	float ndens[NSPECIES], ndens_mid[NSPECIES];
	set_ndens(target_mesh, ndens);
	for(int ispec=0;ispec<NSPECIES;ispec++) {

	  ndens_mid[ispec] = ndens[ispec] 
	    - dtdx*(flux_x_plus.ndens_flux[ispec] - flux_x_minus.ndens_flux[ispec])
	    - dtdy*(flux_y_plus.ndens_flux[ispec] - flux_y_minus.ndens_flux[ispec])
	    - dtdz*(flux_z_plus.ndens_flux[ispec] - flux_z_minus.ndens_flux[ispec]);

	}
	correct_ndens(dest_mesh, ndens_mid);
	set_fraction(dest_mesh, ndens_mid);
      }
    }
  }

  end_timing(&start_tv_part, &end_tv_part, &start_tms_part, &end_tms_part,
	     "fluid_integrate first_rk_loop", this_run);
  
#ifdef __GRAVITY__
  start_timing(&start_tv_part, &start_tms_part);

  /* gravitational potential at t^{n+1/2} */
  zero_out_mesh_density(dens, this_fftw_mpi);
  calc_mesh_density(dens, mesh_mid, this_run, this_fftw_mpi);
  calc_mesh_grav_pot(dens, mesh_mid, this_run, this_fftw_mpi, green_func);
#endif /* __GRAVITY__ */

  update_pad_region(mesh_mid, &pad_mid, this_run, this_mpi);

#ifdef __GRAVITY__
  add_gravity_term(mesh, mesh_mid, &pad, &pad_mid, this_run, this_mpi, dtime, 0);
  update_pad_region(mesh_mid, &pad_mid, this_run, this_mpi);

  end_timing(&start_tv_part, &end_tv_part, &start_tms_part, &end_tms_part,
	     "fluid_integrate first_calc_grav (include update_pad_region)", this_run);
#endif /* __GRAVITY__ */



#ifdef __SECOND_ORDER_RUNGE_KUTTA__

  start_timing(&start_tv_part, &start_tms_part);

#pragma omp parallel for schedule(auto) collapse(3)
  for(int ix=0;ix<NMESH_X_LOCAL;ix++) {
    for(int iy=0;iy<NMESH_Y_LOCAL;iy++) {
      for(int iz=0;iz<NMESH_Z_LOCAL;iz++) {
	struct fluid_flux flux_x_plus, flux_x_minus;
	struct fluid_flux flux_y_plus, flux_y_minus;
	struct fluid_flux flux_z_plus, flux_z_minus;

	float gamma, gam1;
	float therm_eneg, kinetic_eneg;

	int shock_x, shock_y, shock_z;

	struct fluid_mesh *target_mesh, *dest_mesh;

	target_mesh = &MESH_MID(ix,iy,iz);
	dest_mesh = &MESH(ix,iy,iz);

	calc_flux_x(mesh_mid, &pad_mid, ix, iy, iz, &flux_x_plus, &flux_x_minus, this_run, &shock_x);
	calc_flux_y(mesh_mid, &pad_mid, ix, iy, iz, &flux_y_plus, &flux_y_minus, this_run, &shock_y);
	calc_flux_z(mesh_mid, &pad_mid, ix, iy, iz, &flux_z_plus, &flux_z_minus, this_run, &shock_z);

	float ndens[NSPECIES];
	set_ndens(target_mesh, ndens);

	gamma = gamma_total(target_mesh, this_run);
	gam1 = gamma-1.0;

	dest_mesh->dens = 0.5*dest_mesh->dens + 0.5*target_mesh->dens 
	  - 0.5*dtdx*(flux_x_plus.dens_flux - flux_x_minus.dens_flux)
	  - 0.5*dtdy*(flux_y_plus.dens_flux - flux_y_minus.dens_flux)
	  - 0.5*dtdz*(flux_z_plus.dens_flux - flux_z_minus.dens_flux);
	
	dest_mesh->momx = 0.5*dest_mesh->momx + 0.5*target_mesh->momx 
	  - 0.5*dtdx*(flux_x_plus.momx_flux - flux_x_minus.momx_flux)
	  - 0.5*dtdy*(flux_y_plus.momx_flux - flux_y_minus.momx_flux)
	  - 0.5*dtdz*(flux_z_plus.momx_flux - flux_z_minus.momx_flux);

	dest_mesh->momy = 0.5*dest_mesh->momy + 0.5*target_mesh->momy 
	  - 0.5*dtdx*(flux_x_plus.momy_flux - flux_x_minus.momy_flux)
	  - 0.5*dtdy*(flux_y_plus.momy_flux - flux_y_minus.momy_flux)
	  - 0.5*dtdz*(flux_z_plus.momy_flux - flux_z_minus.momy_flux);

	dest_mesh->momz = 0.5*dest_mesh->momz + 0.5*target_mesh->momz 
	  - 0.5*dtdx*(flux_x_plus.momz_flux - flux_x_minus.momz_flux)
	  - 0.5*dtdy*(flux_y_plus.momz_flux - flux_y_minus.momz_flux)
	  - 0.5*dtdz*(flux_z_plus.momz_flux - flux_z_minus.momz_flux);

	dest_mesh->eneg = 0.5*dest_mesh->eneg + 0.5*target_mesh->eneg 
	  - 0.5*dtdx*(flux_x_plus.eneg_flux - flux_x_minus.eneg_flux)
	  - 0.5*dtdy*(flux_y_plus.eneg_flux - flux_y_minus.eneg_flux)
	  - 0.5*dtdz*(flux_z_plus.eneg_flux - flux_z_minus.eneg_flux);
	
	kinetic_eneg = 0.5*NORML2(dest_mesh->momx, 
				  dest_mesh->momy, 
				  dest_mesh->momz)/dest_mesh->dens;

	float wmol = WMOL(target_mesh->chem);
	
	if( (target_mesh->high_mach) && !(target_mesh->under_shock) ) {
	  /* update the modified entropy using the master equation */

	  dest_mesh->etrp = 0.5*dest_mesh->etrp + 0.5*target_mesh->etrp 
	    - 0.5*dtdx*(flux_x_plus.etrp_flux - flux_x_minus.etrp_flux)
	    - 0.5*dtdy*(flux_y_plus.etrp_flux - flux_y_minus.etrp_flux)
	    - 0.5*dtdz*(flux_z_plus.etrp_flux - flux_z_minus.etrp_flux);

	  if(dest_mesh->etrp < FLUID_TINY) dest_mesh->etrp = FLUID_TINY;
	
	  therm_eneg = dest_mesh->etrp*powf(dest_mesh->dens,gam1)/gam1;

	  if(therm_eneg/dest_mesh->dens*this_run->uenetok*wmol < 10.0) { // 10 K
	    therm_eneg = dest_mesh->dens/(this_run->uenetok*wmol);
	  }
	  
	  dest_mesh->eneg = therm_eneg + kinetic_eneg;
	  dest_mesh->uene = therm_eneg/dest_mesh->dens;
	}else{
	  therm_eneg = 
	    dest_mesh->eneg - kinetic_eneg;

	  if(therm_eneg/dest_mesh->dens*this_run->uenetok*wmol < 10.0) { // 10 K
	    therm_eneg = dest_mesh->dens/(this_run->uenetok*wmol);
	  }
	  /*
	  if(therm_eneg < FLUID_TINY)  {
	    therm_eneg = FLUID_TINY;
	    dest_mesh->eneg = therm_eneg + kinetic_eneg;
	  }
	  */
	  dest_mesh->uene = therm_eneg/dest_mesh->dens;
	  dest_mesh->etrp = gam1*therm_eneg/powf(dest_mesh->dens,gam1);
	}

	float ndens_mid[NSPECIES];
	set_ndens(target_mesh, ndens_mid);
	for(int ispec=0;ispec<NSPECIES;ispec++) {

	  ndens[ispec] = 0.5*ndens[ispec] + 0.5*ndens_mid[ispec]
	    - 0.5*dtdx*(flux_x_plus.ndens_flux[ispec] - flux_x_minus.ndens_flux[ispec])
	    - 0.5*dtdy*(flux_y_plus.ndens_flux[ispec] - flux_y_minus.ndens_flux[ispec])
	    - 0.5*dtdz*(flux_z_plus.ndens_flux[ispec] - flux_z_minus.ndens_flux[ispec]);

	}
	correct_ndens(dest_mesh, ndens);
	set_fraction(dest_mesh, ndens);

	/* update previous values of dest_mesh */
	dest_mesh->prev_chem = target_mesh->chem;
	dest_mesh->prev_uene = target_mesh->uene;
	
      }
    }
  }

  end_timing(&start_tv_part, &end_tv_part, &start_tms_part, &end_tms_part,
	     "fluid_integrate second_rk_loop", this_run);
  
#ifdef __GRAVITY__

  start_timing(&start_tv_part, &start_tms_part);

  /* gravitational potential at t=t^{n+1} */
  zero_out_mesh_density(dens, this_fftw_mpi);
  calc_mesh_density(dens, mesh, this_run, this_fftw_mpi);
  calc_mesh_grav_pot(dens, mesh, this_run, this_fftw_mpi, green_func);

  update_pad_region(mesh, &pad, this_run, this_mpi);

  add_gravity_term(mesh_mid, mesh, &pad_mid, &pad, this_run, this_mpi, dtime, 1);

  end_timing(&start_tv_part, &end_tv_part, &start_tms_part, &end_tms_part,
	     "fluid_integrate second_calc_grav (include update_pad_region)", this_run); 
#endif /* __GRAVITY__ */

#else /* ! __SECOND_ORDER_RUNGE_KUTTA__ */
#pragma omp parallel for schedule(auto) collapse(3)
  for(int ix=0;ix<NMESH_X_LOCAL;ix++) {
    for(int iy=0;iy<NMESH_Y_LOCAL;iy++) {
      for(int iz=0;iz<NMESH_Z_LOCAL;iz++) {
	MESH(ix,iy,iz) = MESH_MID(ix,iy,iz);

	/* update previous values of dest_mesh */
	MESH(ix,iy,iz).prev_chem = MESH_MID(ix,iy,iz).chem;
	MESH(ix,iy,iz).prev_uene = MESH_MID(ix,iy,iz).uene;
      }
    }
  }
#endif

  free_pad_region(&pad);
  free_pad_region(&pad_mid);

  free(mesh_mid);

  end_timing(&start_tv, &end_tv, &start_tms, &end_tms,
	     "fluid_integrate", this_run);
  fprintf(this_run->proc_file,"\n");

}
