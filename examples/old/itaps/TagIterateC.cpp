/*
 * This program shows how to get a pointer to tag memory, allowing an application to work
 * directly with tag memory instead of calling through the api.
 */

#include "iMesh.h"
#include "iMesh_extensions.h"
#include "stdio.h"
#include "string.h"

#define CHKERR(e, m) if (iBase_SUCCESS != e) {printf(m); return 1;}
    
int main( int argc, char *argv[] )
{
  iMesh_Instance mesh;
  iBase_EntitySetHandle root_set;
  int err;
  const char *filename;
  iBase_EntityArrIterator iter;
  iBase_TagHandle tagh;
  int count, atend;
  double *tag_data;
  int num_regions;

  if (argc == 2) {
    filename = argv[1];
  }
  else {
    printf("Usage: %s <mesh_filename>\n", argv[0]);
    return 0;
  }

    /* initialize the Mesh */
  iMesh_newMesh(NULL, &mesh, &err, 0);
  CHKERR(err, "Failed to create a mesh instance.\n");
  iMesh_getRootSet(mesh, &root_set, &err);
  CHKERR(err, "Failed to return a root set.\n");

  iMesh_load(mesh, root_set, filename, NULL, &err, strlen(filename), 0);

    /* get the number of regions in the mesh */
  iMesh_getNumOfType(mesh, root_set, iBase_REGION, &num_regions, &err);
  CHKERR(err, "Failed to get number of regions.");
  
    /* get an iterator to all regions in the model */
  iMesh_initEntArrIter( mesh, root_set, iBase_REGION, iMesh_ALL_TOPOLOGIES, num_regions, 0, &iter, &err );
  CHKERR(err, "Failed to create iterator over regions.");

    /* create a tag to put on the regions */
  iMesh_createTagWithOptions(mesh, "dumtag", "moab:TAG_STORAGE_TYPE=DENSE moab:TAG_DEFAULT_VALUE=0.0", 1, iBase_DOUBLE, &tagh, &err, 6, 54);
  CHKERR(err, "Failed to create a tag.\n");
  
  atend = 0;
  
  while (!atend) {
    
      /* get a pointer to that tag data; this will allocate tag storage if it isn't allocated yet */
    iMesh_tagIterate(mesh, tagh, iter, &tag_data, &count, &err);
    CHKERR(err, "Failed to create a tag.\n");
    
      /* step the iterator over count entities */
    iMesh_stepEntArrIter(mesh, iter, count, &atend, &err);
    CHKERR(err, "Failed to step iterator.\n");

      /* operate on tag data, or store it for later */
  }
  
  iMesh_dtor(mesh, &err);
  CHKERR(err, "Failed to destroy iMesh instance.\n");

  return 0;
}
