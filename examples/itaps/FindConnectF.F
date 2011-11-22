c FindConnect: Interacting with iMesh
c 
c This program shows how to get more information about a mesh, by
c getting connectivity two different ways (as connectivity and as
c adjacent 0-dimensional entities).

c Usage: FindConnect

      program findconnect
#include "iMesh_f.h"

c declarations
      iMesh_Instance mesh
      iBase_EntityHandle  ents
      IBASE_HANDLE_T rpverts, rpallverts, ipoffsets
      pointer (rpents, ents(0:*))
      pointer (rpverts, verts(0:*))
      pointer (rpallverts, allverts(0:*))
      pointer (ipoffsets, ioffsets(0,*))
      integer ierr, ents_alloc, ents_size
      integer iverts_alloc, iverts_size
      integer allverts_alloc, allverts_size
      integer offsets_alloc, offsets_size

c create the Mesh instance
      call iMesh_newMesh("MOAB", mesh, ierr)

c load the mesh
      call iMesh_load(%VAL(mesh), %VAL(0), "125hex.vtk", "", ierr)

c get all 3d elements
      ents_alloc = 0
      call iMesh_getEntities(%VAL(mesh), %VAL(0), %VAL(iBase_REGION), 
     1     %VAL(iMesh_ALL_TOPOLOGIES), rpents, ents_alloc, ents_size,
     1     ierr)

      ivert_uses = 0

c iterate through them; 
      do i = 0, ents_size-1
c get connectivity
         iverts_alloc = 0
         call iMesh_getEntAdj(%VAL(mesh), %VAL(ents(i)), 
     1        %VAL(iBase_VERTEX), rpverts, iverts_alloc, iverts_size,
     1        ierr)
c sum number of vertex uses

         vert_uses = vert_uses + iverts_size

         call free(rpverts)
      end do

c now get adjacencies in one big block
      allverts_alloc = 0
      offsets_alloc = 0
      call iMesh_getEntArrAdj(%VAL(mesh), %VAL(rpents), 
     1     %VAL(ents_size), %VAL(iBase_VERTEX), rpallverts, 
     1     allverts_alloc, allverts_size, ipoffsets, offsets_alloc, 
     1     offsets_size, ierr)

      call free(rpallverts);
      call free(ipoffsets);
      call free(rpents);

c compare results of two calling methods
      if (allverts_size .ne. vert_uses) then
         write(*,'("Sizes didn''t agree!")')
      else 
         write(*,'("Sizes did agree!")')
      endif

      call iMesh_dtor(%VAL(mesh), ierr)

      end
