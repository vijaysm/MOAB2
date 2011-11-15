! MOAB structured mesh extension test
! 
! This test also tests fortran free-source format
!

#define ERROR(rval) if (0 .ne. rval) call exit(1)

program ScdMesh
implicit none
integer comm1, mysize,myproc,ier
#include "iMesh_f.h"
iMesh_Instance ::  mesh
iBase_EntitySetHandle :: handle
iBase_EntityHandle :: root_set
integer :: local_dims(6),global_dims(6)
integer :: geom_dim,num_regions, num_verts

! declarations

! create the Mesh instance

local_dims(1)=0
local_dims(2)=0
local_dims(3)=-1
local_dims(4)=64
local_dims(5)=64
local_dims(6)=-1

global_dims(1)=0
global_dims(2)=0
global_dims(3)=-1
global_dims(4)=64
global_dims(5)=64
global_dims(6)=-1

call iMesh_newMesh('MOAB', mesh, ier); ERROR(ier);

call iMesh_createStructuredMesh(%VAL(mesh), local_dims, global_dims, %VAL(0),%VAL(0),%VAL(0), %VAL(1), %VAL(-1), &
  %VAL(-1), %VAL(-1), %VAL(0), handle, ier); ERROR(ier);

call iMesh_getRootSet(%VAL(mesh), root_set, ier); ERROR(ier);

call iMesh_getGeometricDimension(%VAL(mesh), geom_dim, ier); ERROR(ier);

call iMesh_getNumOfType(%VAL(mesh), %VAL(root_set), %VAL(iBase_FACE), num_regions, ier); ERROR(ier);

call iMesh_getNumOfType(%VAL(mesh), %VAL(root_set), %VAL(iBase_VERTEX), num_verts, ier); ERROR(ier);

call exit(0)
end
