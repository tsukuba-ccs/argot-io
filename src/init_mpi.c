#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <assert.h>

#include "run_param.h"
#include "mpi_param.h"

int mpi_rank(int, int,int);

void init_mpi(struct run_param *this_run, struct mpi_param *this_mpi)
{
  MPI_Comm_size(MPI_COMM_WORLD, &(this_run->mpi_nproc));
  MPI_Comm_rank(MPI_COMM_WORLD, &(this_run->mpi_rank));

  if(NNODE_X*NNODE_Y*NNODE_Z != this_run->mpi_nproc) {
    fprintf(stderr, "Number of nodes are inconsistet.\n");
    fprintf(stderr, "NNODE : %d\n", NNODE_X*NNODE_Y*NNODE_Z);
    fprintf(stderr, 
            "Number of nodes specified by MPI : %d\n", this_run->mpi_nproc);
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  /* Determine the 3-dimensional rank */
  this_run->rank_x = this_run->mpi_rank / (NNODE_Y*NNODE_Z);
  this_run->rank_y = 
    (this_run->mpi_rank-this_run->rank_x*NNODE_Y*NNODE_Z) / NNODE_Z;
  this_run->rank_z = this_run->mpi_rank - 
    this_run->rank_x*NNODE_Y*NNODE_Z - this_run->rank_y*NNODE_Z;

  assert(mpi_rank(this_run->rank_x, this_run->rank_y, this_run->rank_z) == 
         this_run->mpi_rank);

  /* compute the number of chemical channel */
  int nchannel;

  nchannel = 1; /* HI */
#ifdef __HELIUM__
  nchannel += 2; /* HeI, HeII */
#endif
#ifdef __HYDROGEN_MOL__
  nchannel += 5; /* HM, H2I_I, H2I_II, H2II_I, H2II_II */
#endif /* __HYDROGEN_MOL__ */

  this_run->nchannel = nchannel;
  assert(nchannel == NCHANNEL);

  /* define the MPI data type for the ray_segment structure */
  int blockcount[3];
  MPI_Datatype type[3];
  MPI_Aint adr[3];

  blockcount[0]=1;
  blockcount[1]=2;
  blockcount[2]=6+nchannel;
  type[0]=MPI_UNSIGNED_LONG_LONG;
  type[1]=MPI_INT;
  type[2]=MPI_FLOAT;
  adr[0]=0;
  adr[1]=8;
  adr[2]=16;

  MPI_Type_create_struct(3,blockcount, adr, type, &(this_mpi->segment_type));
  MPI_Type_commit(&(this_mpi->segment_type));

  /* define the MPI data type for the prim_chem structure */
  blockcount[0] = 6;
#ifdef __HYDROGEN_MOL__
  blockcount[0] += 13;
#endif
#ifdef __HELIUM__
  blockcount[0] += 7;
#endif
  type[0]=MPI_FLOAT;
  adr[0]=0;

  MPI_Type_create_struct(1,blockcount, adr, type, &(this_mpi->prim_chem_type));
  MPI_Type_commit(&(this_mpi->prim_chem_type));

  /* define the MPI data type for the fluid_mesh_io structure */
  blockcount[0]=7;
  blockcount[1]=1;
  type[0]=MPI_FLOAT;
  type[1]=this_mpi->prim_chem_type;
  adr[0]=0;
  adr[1]=24;

  MPI_Type_create_struct(2,blockcount, adr, type, &(this_mpi->fluid_mesh_io_type));
  MPI_Type_commit(&(this_mpi->fluid_mesh_io_type));

  /* define the MPI data type fo ther photoion_rate structure */
  blockcount[0] = 2;
#ifdef __HELIUM__
  blockcount[0] += 4;
#endif
#ifdef __HYDROGEN_MOL__
  blockcount[0] += 10;
#endif
  type[0]=MPI_FLOAT;
  adr[0]=0;

  MPI_Type_create_struct(1,blockcount, adr, type, &(this_mpi->photoion_rate_type));
  MPI_Type_commit(&(this_mpi->photoion_rate_type));
  
}

