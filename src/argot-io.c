#include <sys/time.h>
#include <unistd.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <getopt.h>
#include "io_ops.h"

struct cosmology {
  float omega0, lambda0, omegab, hubble;
  float tend;
};

#define NGRID_NU (32)

struct freq_param {
  double log_nu_min, log_nu_max, dlog_nu;
  double nu[NGRID_NU];
};

struct cross_section {
  float csect_HI;
#ifdef __HELIUM__
  float csect_HeI;
  float csect_HeII;
#endif /* __HELIUM__ */
};

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

struct fluid_mesh {
  float dens;
  float eneg;
  float momx, momy, momz;
  float uene, duene, durad;
  float etrp;
  float pot;
  float dens_prev; /* density in the previous step */

  struct prim_chem chem;
  struct prim_chem prev_chem; /* chemical composition in the previous iteration */
  float  prev_uene; /* specific energy in the previous iteration */

  short under_shock;
  short high_mach;
};

struct pad_region {
  struct fluid_mesh *pad_x_lo, *pad_x_hi;
  struct fluid_mesh *pad_y_lo, *pad_y_hi;
  struct fluid_mesh *pad_z_lo, *pad_z_hi;
};

struct run_param {

  /* model name */
  char *model_name;

  FILE *proc_file, *diag_file, *err_file;

  int mpi_rank, mpi_nproc;
  int rank_x, rank_y, rank_z;
  int nnode_x, nnode_y, nnode_z;

  double lunit, munit, tunit, eunit;
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
} run_param;

struct fluid_mesh_io {
  float dens;
  float eneg;
  float momx, momy, momz;
  float pot;
  struct prim_chem chem;
};

void
output_performance(char *diag, uint64_t size, float time)
{
	printf("%s: %llu %g %g\n", diag, (unsigned long long)size,
	       time, size / time / 1.0e9);
}

#define timerval_sub(t1, t2) \
	(((t1)->tv_sec - (t2)->tv_sec) \
	 + .000001f * ((t1)->tv_usec - (t2)->tv_usec))
#define SQR(x) ((x)*(x))

void
input_data(int size, struct fluid_mesh *mesh, struct run_param *param,
	   char *fn)
{
	void *fp;
	int i;
	struct timeval tv1, tv2;
	uint64_t input_size = 0;
	float t;

	gettimeofday(&tv1, NULL);
	io_ops->read_open(fn, &fp);
	io_ops->read(param, sizeof(*param), fp);
	input_size += sizeof(*param);
	for (i = 0; i < size; ++i) {
		struct fluid_mesh_io temp_mesh;

		io_ops->read(&temp_mesh, sizeof(struct fluid_mesh_io), fp);
		mesh[i].dens = temp_mesh.dens;
		mesh[i].eneg = temp_mesh.eneg;
		mesh[i].momx = temp_mesh.momx;
		mesh[i].momy = temp_mesh.momy;
		mesh[i].momz = temp_mesh.momz;
		mesh[i].chem = temp_mesh.chem;
		mesh[i].uene = (mesh[i].eneg -
		    0.5 * (SQR(mesh[i].momx) + SQR(mesh[i].momy) +
			   SQR(mesh[i].momz)) / mesh[i].dens) / mesh[i].dens;
		mesh[i].prev_chem = mesh[i].chem;
		mesh[i].prev_uene = mesh[i].uene;
	}
	input_size += sizeof(struct fluid_mesh_io) * size;
	gettimeofday(&tv2, NULL);
	t = timerval_sub(&tv2, &tv1);
	output_performance(fn, input_size, t);
}

void
output_mesh(int size, struct fluid_mesh *mesh, struct run_param *param,
	    char *fn)
{
	void *fp;
	int i;
	struct timeval tv1, tv2;
	uint64_t output_size = 0;
	float t;

	gettimeofday(&tv1, NULL);
	io_ops->write_open(fn, &fp);
	io_ops->write(param, sizeof(struct run_param), fp);
	output_size += sizeof(struct run_param);

	for (i = 0; i < size; ++i) {
		struct fluid_mesh_io tmp_mesh;

		tmp_mesh.dens = mesh[i].dens;
		tmp_mesh.eneg = mesh[i].eneg;
		tmp_mesh.momx = mesh[i].momx;
		tmp_mesh.momy = mesh[i].momy;
		tmp_mesh.momz = mesh[i].momz;
		tmp_mesh.chem = mesh[i].chem;
		io_ops->write(&tmp_mesh, sizeof(struct fluid_mesh_io), fp);
	}
	output_size += sizeof(struct fluid_mesh_io) * size;
	param->output_file_size += output_size;
	io_ops->close(fp);
	gettimeofday(&tv2, NULL);
	t = timerval_sub(&tv2, &tv1);
	param->output_wt += t;
	output_performance(fn, output_size, t);
}

void
usage(void)
{
	fprintf(stderr, "Usage: argot-io-light [-a api] [-i input_fn] "
	    "[-n iter] [-x size_x] [-y size_y] [-z size_z]\n");
	exit(EXIT_FAILURE);
}

int
main(int argc, char *argv[])
{
	struct fluid_mesh *mesh;
	uint64_t size, size_x = 128, size_y = 128, size_z = 128;
	int c, i, j, iter = 10;
	char *input_fn = NULL, fn[20], *api = "posix";

	while ((c = getopt(argc, argv, "a:i:n:x:y:z:")) != -1) {
		switch (c) {
		case 'a':
			api = optarg;
			break;
		case 'i':
			input_fn = optarg;
			break;
		case 'n':
			iter = atoi(optarg);
			break;
		case 'x':
			size_x = atoi(optarg);
			break;
		case 'y':
			size_y = atoi(optarg);
			break;
		case 'z':
			size_z = atoi(optarg);
			break;
		default:
			usage();
		}
	}
	argc -= optind;
	argv += optind;
	if (argc != 0)
		usage();

	init_io_ops(api);
	io_ops->init();
	size = size_x * size_y * size_z;
	mesh = malloc(sizeof(struct fluid_mesh) * size);
	if (input_fn != NULL)
		input_data(size, mesh, &run_param, input_fn);
	run_param.output_file_size = run_param.output_wt = 0;
	for (i = 0; i < iter; ++i) {
		sprintf(fn, "dmp0");
		for (j = 0; j < 5; ++j) {
			sleep(1);
			output_mesh(size, mesh, &run_param, fn);
		}
		sprintf(fn, "output%02d", i);
		sleep(1);
		output_mesh(size, mesh, &run_param, fn);
	}
	output_performance("total", run_param.output_file_size,
	    run_param.output_wt);
	sprintf(fn, "dmp0");
	io_ops->unlink(fn);
	for (i = 0; i < iter; ++i) {
		sprintf(fn, "output%02d", i);
		io_ops->unlink(fn);
	}
	io_ops->term();

	return (0);
}
