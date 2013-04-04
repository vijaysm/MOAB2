  ! ListSetsNTags: list sets & tags from a mesh
  ! 
  ! This program shows how to read and list sets and tags from a mesh
  !
  ! Usage: SetsNTags <mesh_file_name>
  
  
  
#define ERRORR(a) if (0 .ne. err) print *, a
  
  program ListSetsNTags

#include "iMesh_f.h"

    iMesh_Instance mesh
    iBase_EntitySetHandle root_set
    integer err

    IBASE_HANDLE_T rpsets, rptags
    pointer (rpsets, sets(0:*))
    pointer (rptags, tags(0:*))
    iBase_EntitySetHandle sets
    iBase_TagHandle tags
    integer sets_alloc, sets_size, tags_alloc, tags_size

    real*8 dbl_val
    integer int_val, tag_type
    character*128 tname, fname
    character*1024 tname2

    integer i, j, num_hops, num_commands, tname_len
    logical read_par
    data read_par/.false./

    num_commands = command_argument_count()
    if (num_commands .eq. 0) then
       fname = "../MeshFiles/125hex.vtk"
    else
       call get_command_argument(1, tname, tname_len, err)
       if (err .ne. 0) then
          ERRORR("Problem getting filename argument.")
          call exit
       endif
       fname = tname
       if (num_commands .eq. 2) then
          call get_command_argument(2, tname, tname_len, err)
          if (err .ne. 0) then
             ERRORR("Problem getting filename argument.")
             call exit
          endif
          if (tname(1:1) .eq. 'p' .or. tname(1:1) .eq. 'P') then
             read_par = .true.
          endif
       endif
    endif

    ! create the Mesh instance
    call iMesh_newMesh("", mesh, err)
    ERRORR("Error creating new mesh.")


    call iMesh_getRootSet(%VAL(mesh), root_set, err)
    ERRORR("Couldn't get root set.")

    ! load the mesh
    if (read_par) then
       call iMesh_load(%VAL(mesh), %VAL(root_set), fname, &
  " moab:PARALLEL=READ_PART moab:PARTITION=PARALLEL_PARTITION moab:PARTITION_DISTRIBUTE moab:PARALLEL_RESOLVE_SHARED_ENTS " &
            , err)
    else
       call iMesh_load(%VAL(mesh), %VAL(root_set), fname, "", err)
    endif
    ERRORR("Couldn't load mesh.")

    ! get all sets
    sets_alloc = 0
    num_hops = 1
    call iMesh_getEntSets(%VAL(mesh), %VAL(root_set), %VAL(num_hops), &
         rpsets, sets_alloc, sets_size, err)
    ERRORR("Couldn't get all sets.")

    ! iterate through them, checking whether they have tags
    do i = 0, sets_size-1
       ! get connectivity
       tags_alloc = 0
       call iMesh_getAllEntSetTags(%VAL(mesh), %VAL(sets(i)), &
            rptags, tags_alloc, tags_size, err)
       ERRORR("Failed to get ent set tags.")

       if (0 .ne. tags_size) then
          print *, "Set ", sets(i), " Tags:"
       end if

       ! list tag names on this set
       do j = 0, tags_size-1
          call iMesh_getTagName(%VAL(mesh), %VAL(tags(j)), tname, err)
          call iMesh_getTagType(%VAL(mesh), %VAL(tags(j)), tag_type, err)
          ERRORR("Failed to get tag type.")
          if (iBase_INTEGER .eq. tag_type) then
             call iMesh_getEntSetIntData(%VAL(mesh), %VAL(sets(i)), &
                  %VAL(tags(j)), int_val, err)
             ERRORR("Failed to get int data type.")
             print *, tname, int_val
          else if (iBase_DOUBLE .eq. tag_type) then
             call iMesh_getEntSetDblData(%VAL(mesh), %VAL(sets(i)), &
                  %VAL(tags(j)), dbl_val, err)
             print *, tname, dbl_val
          else
             print *, tname
          end if

       end do

       if (tags_size .ne. 0) call free(rptags)
       tags_alloc = 0
    end do

    if (sets_size .ne. 0) call free(rpsets)

    call iMesh_dtor(%VAL(mesh), err)
    ERRORR("Failed to destruct interface.")

!    return
  end program ListSetsNTags

