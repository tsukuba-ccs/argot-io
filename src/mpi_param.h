#ifndef __ARGOT_MPI_PARAM__
#define __ARGOT_MPI_PARAM__

#include <mpi.h>

struct mpi_param {
  MPI_Datatype segment_type;
  MPI_Datatype prim_chem_type;
  MPI_Datatype fluid_mesh_io_type;
  MPI_Datatype photoion_rate_type;
};

#endif /* __ARGOT_MPI_PARAM__ */
