#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <math.h>

#include "constants.h"
#include "run_param.h"
#include "fluid.h"
#include "source.h"
#include "prototype.h"
#include "io_ops.h"

int NMESH_X_TOTAL = 128;
int NMESH_Y_TOTAL = 128;
int NMESH_Z_TOTAL = 128;

int NNODE_X = 2;
int NNODE_Y = 2;
int NNODE_Z = 2;

#define MESH(ix,iy,iz) (mesh[(iz)+NMESH_Z_LOCAL*((iy)+NMESH_Y_LOCAL*(ix))])

void make_directory(char*);

void
usage(char *prog)
{
  fprintf(stderr, "usage: %s [-a api] "
        "[-x mesh_x_total] [-y mesh_y_total]\n\t[-z mesh_z_total] "
        "[-X nnode_x] [-Y nnode_y] [-Z nnode_z] prefix\n", prog);
  exit(EXIT_FAILURE);
}

int main(int argc, char **argv)
{
  static struct run_param this_run;

  static struct fluid_mesh *mesh;
  //  static struct radiation_src src[NSOURCE_MAX];
  struct radiation_src *src;

  static char model_name[256],label[1024],dir_name[256];

  double nH, tmpr;
  char *api = "posix", *program;
  int c;

  program = argv[0];

  while ((c = getopt(argc, argv, "a:x:X:y:Y:z:Z:")) != -1) {
    switch (c) {
    case 'a':
      api = optarg;
      break;
    case 'x':
      NMESH_X_TOTAL = atoi(optarg);
      break;
    case 'y':
      NMESH_Y_TOTAL = atoi(optarg);
      break;
    case 'z':
      NMESH_Z_TOTAL = atoi(optarg);
      break;
    case 'X':
      NNODE_X = atoi(optarg);
      break;
    case 'Y':
      NNODE_Y = atoi(optarg);
      break;
    case 'Z':
      NNODE_Z = atoi(optarg);
      break;
    default:
      usage(program);
    }
  }
  argc -= optind;
  argv += optind;

  if (argc != 1)
    usage(program);

  mesh = malloc(sizeof(*mesh) * NMESH_X_LOCAL * NMESH_Y_LOCAL * NMESH_Z_LOCAL);
  if (mesh == NULL)
    fprintf(stderr, "no memory\n"), exit(EXIT_FAILURE);

  init_io_ops(api);
  io_ops->init();

  this_run.nmesh_x_total=NMESH_X_TOTAL;
  this_run.nmesh_y_total=NMESH_Y_TOTAL;
  this_run.nmesh_z_total=NMESH_Z_TOTAL;

  this_run.nmesh_x_local=NMESH_X_LOCAL;
  this_run.nmesh_y_local=NMESH_Y_LOCAL;
  this_run.nmesh_z_local=NMESH_Z_LOCAL;

  this_run.xmin=0.0;
  this_run.ymin=0.0;
  this_run.zmin=0.0;
  
  this_run.xmax=1.0;
  this_run.ymax=1.0;
  this_run.zmax=1.0;

  this_run.delta_x=(this_run.xmax-this_run.xmin)/(float)this_run.nmesh_x_total;
  this_run.delta_y=(this_run.ymax-this_run.ymin)/(float)this_run.nmesh_y_total;
  this_run.delta_z=(this_run.zmax-this_run.zmin)/(float)this_run.nmesh_z_total;

  /* density and temperature of the medium */
  nH = 1.0e-3;
  tmpr = 1.0e4;

  this_run.lunit = 1.5e-2*mpc;
  //  this_run.lunit = 6.6e-3*mpc;
  this_run.munit = nH*CUBE(this_run.lunit)*mproton/XHYDROGEN;
  //this_run.tunit = 3.86e15; /* recombination timescale = ne*alpha_B */  
  this_run.tunit = 1.0e8*year; /* recombination timescale = ne*alpha_B */  

  this_run.denstonh = this_run.munit/CUBE(this_run.lunit)*XHYDROGEN/mproton;
  this_run.uenetok = GAMM1_MONOATOMIC*mproton/kboltz*
    SQR(this_run.lunit)/SQR(this_run.tunit);
  this_run.anow = 1.0;
  this_run.znow = 0.0;

  this_run.output_indx = -1; // just for initial conditions
  this_run.ngrid_nu = NGRID_NU;
  this_run.nspecies = NSPECIES;
  this_run.nchannel = NCHANNEL;

  /* ionization state at the initial condition */
  struct prim_chem ioneq_chem;
  //  calc_ioneq(&ioneq_chem, nH, tmpr, 0.0);
  ioneq_chem.fHI = 0.5;
  ioneq_chem.fHII = 0.5;
  ioneq_chem.felec = ioneq_chem.fHI;
  ioneq_chem.GammaHI = 0.0;
#ifdef __HELIUM__
  ioneq_chem.GammaHeI = 0.0;
  ioneq_chem.GammaHeII = 0.0;
#endif

  tmpr=1.0e2;

  //  printf("# initial fHI = %14.6e\n",ioneq_chem.fHI);

  int rank_x, rank_y, rank_z;

  this_run.nnode_x = NNODE_X;
  this_run.nnode_y = NNODE_Y;
  this_run.nnode_z = NNODE_Z;

  this_run.step = 0;
  this_run.tnow = 0.0;

  this_run.mpi_rank = 0;

  setup_freq_param(&this_run.freq);

#if 1
  this_run.nsrc = 32768;

  src = (struct radiation_src *) malloc(sizeof(struct radiation_src)*this_run.nsrc);
  src[0].xpos = 0.4999;
  src[0].ypos = 0.4999;
  src[0].zpos = 0.4999;
  src[0].type = 0; /* black body */ 
  src[0].param = 1.0e5; /* T_bb= 100000 K */
  //  src[0].photon_rate = 5.0e48;
  setup_photon_rate(&this_run.freq, &src[0], 5.0e48);
#if 0
  for(int inu=0;inu<NGRID_NU;inu++) {
    printf("%14.6e %14.6e\n", this_run.freq.nu[inu], src[0].photon_rate[inu]);
  }
#endif
#else
  this_run.nsrc = 16;
  srand(2);

  src = (struct radiation_src *) malloc(sizeof(struct radiation_src)*this_run.nsrc);
  int isrc;
  for(isrc=0;isrc<this_run.nsrc;isrc++) {
    src[isrc].xpos = (float)rand()/(float)RAND_MAX;
    src[isrc].ypos = (float)rand()/(float)RAND_MAX;
    src[isrc].zpos = (float)rand()/(float)RAND_MAX;
    src[isrc].type = 0; 
    src[isrc].param = 5.0e3;
    setup_photon_rate(&this_run.freq, &src[isrc], 5.0e48);
  }
#endif

  sprintf(model_name, "%s", argv[0]);
  sprintf(dir_name, "%s-init", argv[0]);
  make_directory(dir_name);
  sprintf(label,"%s-init/%s-init",model_name,model_name);

  output_src(src, &this_run, label);

  for(rank_x=0;rank_x<NNODE_X;rank_x++) {
    float dx_domain = (this_run.xmax-this_run.xmin)/(float)NNODE_X;
    this_run.xmin_local = this_run.xmin + (float)rank_x*dx_domain;
    this_run.xmax_local = this_run.xmin_local + dx_domain;

    this_run.rank_x = rank_x;
    for(rank_y=0;rank_y<NNODE_Y;rank_y++) {
      float dy_domain = (this_run.ymax-this_run.ymin)/(float)NNODE_Y;
      this_run.ymin_local = this_run.ymin + (float)rank_y*dy_domain;
      this_run.ymax_local = this_run.ymin_local + dy_domain;

    this_run.rank_y = rank_y;
      for(rank_z=0;rank_z<NNODE_Z;rank_z++) {
	float dz_domain = (this_run.zmax-this_run.zmin)/(float)NNODE_Z;
	this_run.zmin_local = this_run.zmin + (float)rank_z*dz_domain;
	this_run.zmax_local = this_run.zmin_local + dz_domain;

	this_run.rank_z = rank_z;

	this_run.mpi_nproc = NNODE_X*NNODE_Y*NNODE_Z;
	this_run.mpi_rank = mpi_rank(rank_x,rank_y,rank_z);

	int ix,iy,iz;
	for(ix=0;ix<NMESH_X_LOCAL;ix++) {
	  for(iy=0;iy<NMESH_Y_LOCAL;iy++) {
	    for(iz=0;iz<NMESH_Z_LOCAL;iz++) {
	      struct fluid_mesh *tgt;

	      float xpos = this_run.xmin_local + ((float)ix+0.5)*this_run.delta_x;
	      float ypos = this_run.ymin_local + ((float)iy+0.5)*this_run.delta_y;
	      float zpos = this_run.zmin_local + ((float)iz+0.5)*this_run.delta_z;
	      if( fabsf(xpos-0.5) < 0.25 &&
		  fabsf(ypos-0.5) < 0.25 &&
		  fabsf(zpos-0.5) < 0.25) {
		    

		tgt = &MESH(ix,iy,iz);

		tgt->dens = 0.1*nH/this_run.denstonh;
		tgt->momx = 0.0;
		tgt->momy = 0.0;
		tgt->momz = 0.0;
		tgt->uene = tmpr/(this_run.uenetok*WMOL(ioneq_chem));
		tgt->eneg = tgt->dens*tgt->uene;
		tgt->chem = ioneq_chem;
		tgt->prev_chem = ioneq_chem;
	      }else{
		tgt = &MESH(ix,iy,iz);

		tgt->dens = nH/this_run.denstonh;
		tgt->momx = 0.0;
		tgt->momy = 0.0;
		tgt->momz = 0.0;
		tgt->uene = tmpr/(this_run.uenetok*WMOL(ioneq_chem));
		tgt->eneg = tgt->dens*tgt->uene;
		tgt->chem = ioneq_chem;
		tgt->prev_chem = ioneq_chem;
	      }
	    }
	  }
	}

	output_mesh(mesh, &this_run, label);

      }
    }
  }

  //  printf("# initial heat capacity ratio : %14.6e\n", gamma_total(&mesh[0], &this_run));

  free(src);
  free(mesh);

  io_ops->term();
}
