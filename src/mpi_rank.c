#include <stdio.h>
#include <math.h>

#include "run_param.h"

int mpi_rank(int rank_x, int rank_y, int rank_z)
{
  int ix, iy, iz;

#ifdef __ISOLATED__
  if(rank_x>=NNODE_X || rank_x < 0 ||
     rank_y>=NNODE_Y || rank_y < 0 ||
     rank_z>=NNODE_Z || rank_z < 0){
    return (-1);
  }else{
    return ((rank_z)+NNODE_Z*((rank_y)+NNODE_Y*(rank_x)));
  }
#else /* __PERIODIC__ */
  ix = (rank_x + NNODE_X) % NNODE_X;
  iy = (rank_y + NNODE_Y) % NNODE_Y;
  iz = (rank_z + NNODE_Z) % NNODE_Z;

  return ((iz)+NNODE_Z*((iy)+NNODE_Y*(ix)));
#endif
}
