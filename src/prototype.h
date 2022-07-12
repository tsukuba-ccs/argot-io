#ifndef __ARGOT_PROTOTYPE__
#define __ARGOT_PROTOTYPE__

#include "run_param.h"
#include "mpi_param.h"
#include "fluid.h"
#include "radiation.h"
#include "source.h"
#include "cross_section.h"

#ifdef __DIFFUSE_RADIATION__
#include "diffuse_photon.h"
#endif /* __DIFFUSE_RADIATION__ */

#include <sys/time.h>
#include <sys/times.h>

int mpi_rank(int,int,int);
void init_mpi(struct run_param*, struct mpi_param*);
void init_run(struct run_param*);
void setup_data(struct run_param*);
void input_mesh(struct fluid_mesh*, struct run_param*, char*);
void input_mesh_single(struct fluid_mesh*, struct run_param*, char*);
void input_mesh_header(struct run_param*, char*);
void input_mesh_header_prefix(struct run_param*, char*);
void input_src(struct radiation_src*, struct run_param*, char*);
void input_src_file(struct radiation_src*, struct run_param*, char*);
void input_data(struct fluid_mesh*, struct radiation_src*, struct run_param*, char*);
void input_params(struct run_param*, char*);
void output_mesh(struct fluid_mesh*, struct run_param*, char*);
void output_mesh_single(struct fluid_mesh*, struct run_param*, char*);
void output_src(struct radiation_src*, struct run_param*, char*);
void output_data(struct fluid_mesh*, struct radiation_src*, struct run_param*, char*);
void output_data_in_run(struct fluid_mesh*, struct radiation_src*, struct run_param*, char*);
void output_diagnostics(struct fluid_mesh*, struct run_param*, float);
void calc_ray_segment(struct light_ray*, struct run_param*);
void setup_cross_section(struct cross_section*, struct freq_param*);
void setup_light_ray_long(struct light_ray**, struct radiation_src*, struct run_param*);
void setup_light_ray_tree(struct light_ray**, struct radiation_src*, struct run_param*);
void assign_ray_segment(struct light_ray*, struct ray_segment*, struct run_param*, struct mpi_param*);
void count_ray_segment(struct light_ray*, struct run_param*);
void accum_optical_depth(struct light_ray*, struct ray_segment*, struct run_param*);
void advance_reaction(struct fluid_mesh*, struct prim_chem*,struct run_param*, float);
void advance_reaction_and_heatcool(struct fluid_mesh*, float *, struct prim_chem*, struct run_param*, float, int*, int*);
void advance_heatcool(struct fluid_mesh*, float*, struct prim_chem*,struct run_param*, float, int*, int*);
float gamma_total(struct fluid_mesh*, struct run_param*);
double calc_dtime(struct fluid_mesh*, struct prim_chem*, struct run_param*);
void update_chemistry(struct fluid_mesh*, struct run_param*);
float timing(struct tms, struct tms);
float wallclock_timing(struct timeval, struct timeval);
int start_timing(struct timeval*, struct tms*);
int end_timing(struct timeval*, struct timeval*, struct tms*, struct tms*, const char*, struct run_param*);
double calc_timestep_fluid(struct fluid_mesh*, struct run_param*);
void init_pad_region(struct pad_region*, struct run_param*);
void update_pad_x(struct fluid_mesh*, struct pad_region*, struct run_param*, struct mpi_param*);
void update_pad_y(struct fluid_mesh*, struct pad_region*, struct run_param*, struct mpi_param*);
void update_pad_z(struct fluid_mesh*, struct pad_region*, struct run_param*, struct mpi_param*);
void init_pad_region(struct pad_region*, struct run_param*);
void free_pad_region(struct pad_region*);
void update_pad_region(struct fluid_mesh*, struct pad_region*, struct run_param*, struct mpi_param*);
void fluid_integrate(struct fluid_mesh*, struct run_param*, struct mpi_param*, float);
float calc_mem_size_for_radiation(struct radiation_src*, struct run_param*, int);
int get_optimal_nmesh_per_loop(struct radiation_src*, struct run_param*);
void set_optimal_nmesh_per_loop(struct radiation_src*, struct run_param*);



#ifdef __USE_GPU__
void calc_optical_depth(struct ray_segment*, struct cuda_mem_space*, struct cuda_param*, struct run_param*);
void calc_photoion_rate_at_first(struct fluid_mesh*, struct radiation_src*, struct cuda_mem_space*, struct cuda_param*, struct run_param*, struct mpi_param*);
void calc_photoion_rate(struct cuda_mem_space*, struct light_ray*, struct light_ray_IO*, struct cuda_param*, struct run_param*);
void step_radiation(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*, struct cuda_mem_space*, struct cuda_param*, float);
void smooth_photoion_rate(struct fluid_mesh*, struct run_param*, struct cuda_mem_space*, struct cuda_param*, struct mpi_param*);
double calc_timestep_chem(struct cuda_mem_space*, struct cuda_param*, struct run_param*);
double calc_timestep(struct fluid_mesh*, struct cuda_mem_space*, struct cuda_param*, struct run_param*);
void step_chemistry(struct cuda_mem_space*, struct cuda_param*, struct run_param*, float, float*, float*);
void zero_out_photoion_rate(struct cuda_mem_space*, struct cuda_param*, struct run_param*);
#ifdef __DIFFUSE_RADIATION__
void step_radiation_tree(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*, struct cuda_mem_space*, struct cuda_param*, struct host_diffuse_param*, struct cuda_diffuse_param*, float);
#else /* !__DIFFUSE_RADIATION__ */
void step_radiation_tree(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*, struct cuda_mem_space*, struct cuda_param*, float);
#endif /* __DIFFUSE_RADIATION__ */
#else /* !__USE_GPU__ */
void calc_optical_depth(struct ray_segment*, struct fluid_mesh*, struct run_param*);
void calc_photoion_rate_at_first(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*);
void calc_photoion_rate(struct fluid_mesh*, struct light_ray*, struct run_param*);
void step_radiation(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*, float);
void smooth_photoion_rate(struct fluid_mesh*, struct run_param*, struct mpi_param*);
double calc_timestep_chem(struct fluid_mesh*, struct run_param*);
double calc_timestep(struct fluid_mesh*, struct run_param*);
void step_chemistry(struct fluid_mesh*, struct run_param*, float, float*, float*);
void zero_out_photoion_rate(struct fluid_mesh*, struct run_param*);
#ifdef __DIFFUSE_RADIATION__
void step_radiation_tree(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*, struct host_diffuse_param*, float);
#else /* !__DIFFUSE_RADIATION__ */
void step_radiation_tree(struct fluid_mesh*, struct radiation_src*, struct run_param*, struct mpi_param*, float);
#endif /* __DIFFUSE_RADIATION__ */
#endif /* __USE_GPU__ */

#ifdef __USE_GPU__
void init_gpu(struct fluid_mesh*, struct cuda_mem_space*, struct cuda_param*, struct run_param*);
void send_mesh_data(struct fluid_mesh*, struct cuda_mem_space*, struct cuda_param*, struct run_param*);
void recv_mesh_data(struct fluid_mesh*, struct cuda_mem_space*, struct cuda_param*, struct run_param*);
void update_chemistry_gpu(struct cuda_mem_space*, struct cuda_param*, struct run_param*);
void allocate_pinned_segment(struct ray_segment**, uint64_t);
void deallocate_pinned_segment(struct ray_segment*);
void allocate_pinned_light_ray_IO(struct light_ray_IO**, uint64_t);
void deallocate_pinned_light_ray_IO(struct light_ray_IO*);
#endif

#ifdef __DIFFUSE_RADIATION__
#include "diffuse_photon/diffuse_prototype.h"
#endif

#endif /* __ARGOT_PROTOTYPE__ */


