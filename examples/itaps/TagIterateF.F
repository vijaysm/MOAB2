c TagIterate: Get a pointer to tag memory and interact directly with that 
c 
c This program shows how to get a pointer to tag memory and assign it to
c a fortran array, allowing an application to bypass the API and/or share
c memory with MOAB.

c Usage: TagIterate <mesh_filename>

      program tag_iterate
#include "iMesh_f.h"

c declarations
      iMesh_Instance mesh
      iBase_EntityHandle root_set
      iBase_EntityArrIterator iter
      iBase_TagHandle tagh
      IBASE_HANDLE_T rpents, ents, rpdata
      real*8 tag_data
      pointer (rpdata, tag_data(0:*))
      pointer (rpents, ents(0:*))
      pointer (ipoffsets, ioffsets(0,*))
      integer ierr, ents_alloc, ents_size, num_regions, atend, count

c create the Mesh instance
      call iMesh_newMesh("MOAB", mesh, ierr)

c get the root set
      call iMesh_getRootSet(%VAL(mesh), root_set, ierr)

c load the mesh
      call iMesh_load(%VAL(mesh), %VAL(root_set), "125hex.vtk", "",ierr)

c get the number of regions
      call iMesh_getNumOfType(%VAL(mesh), %VAL(root_set), 
     1 %VAL(iBase_REGION), num_regions, ierr)

c get an iterator over regions
      call iMesh_initEntArrIter(%VAL(mesh), %VAL(root_set), 
     1 %VAL(iBase_REGION), %VAL(iMesh_ALL_TOPOLOGIES),%VAL(num_regions),
     1 %VAL(0), iter, ierr) 

c create a tag to put on the regions
      call iMesh_createTagWithOptions(%VAL(mesh), "dumtag", 
     1 "moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 
     1 %VAL(1), %VAL(iBase_DOUBLE), tagh, ierr)

c iterate over tag memory
 10   call iMesh_tagIterate(%VAL(mesh), %VAL(tagh), %VAL(iter), 
     1     rpdata, count, ierr)
    
c step the iterator over count entities
      call iMesh_stepEntArrIter(%VAL(mesh), %VAL(iter), %VAL(count), 
     1     atend, ierr)

      if (atend .eq. 0) go to 10


      call iMesh_dtor(%VAL(mesh), ierr)

      end
