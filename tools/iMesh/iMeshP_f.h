#ifndef IMESHP_F_H
#define IMESHP_F_H

#define iMeshP_PartitionHandle integer
#define iMeshP_PartHandle integer

#endif 

#include "iMesh_f.h"

      integer iMeshP_LOCAL
      integer iMeshP_REMOTE
      integer iMeshP_INVALID

      integer iMeshP_INTERNAL
      integer iMeshP_BOUNDARY
      integer iMeshP_GHOST

      parameter (iMeshP_LOCAL = 0)
      parameter (iMeshP_REMOTE = 1)
      parameter (iMeshP_INVALID = 2)

      parameter (iMeshP_INTERNAL = 0)
      parameter (iMeshP_BOUNDARY = 1)
      parameter (iMeshP_GHOST = 2)
