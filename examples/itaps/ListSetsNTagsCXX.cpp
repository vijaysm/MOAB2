/** ListSetsNTags: list sets & tags from a mesh
 * 
 * This program shows how to read and list sets and tags from a mesh
 *
 * Usage: SetsNTags <mesh_file_name>
 *
 */

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <vector>
#include "iMesh.h"

#define ERRORR(a) {if (iBase_SUCCESS != err) {std::cout << a << std::endl; return err;}}

int main( int argc, char *argv[] )
{
  const char *filename;
  bool read_par = false;
  if (argc < 2 || !strcmp(argv[1], "-h")) {
    std::cout << "Usage: " << argv[0] << " <filename> [-p]" << std::endl;
    return 0;
  }
  else {
      // Check command line arg
    filename = argv[1];

    if (argc > 2 && !strcmp(argv[2], "-p")) read_par = true;
  }
  
    // create the Mesh instance
  iMesh_Instance mesh;
  int err;
  iMesh_newMesh(NULL, &mesh, &err, 0);
  ERRORR("Error creating new mesh.\n");
  
  
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &err);
  ERRORR("Couldn't get root set.");
  
    // load the mesh
  if (read_par) {
    const char *opt = " moab:PARALLEL=READ_PART moab:PARTITION=PARALLEL_PARTITION moab:PARTITION_DISTRIBUTE moab:PARALLEL_RESOLVE_SHARED_ENTS moab:DEBUG_PIO=2 moab:DEBUG_IO=2 moab:SETS=SETS ";
    iMesh_load(mesh, root_set, filename, opt, &err, strlen(filename), strlen(opt));
  }
  else
    iMesh_load(mesh, root_set, filename, NULL, &err, strlen(filename), 0);
  ERRORR("Couldn't load mesh.");

    // get all sets
  iBase_EntitySetHandle *sets = NULL;
  int sets_alloc = 0, sets_size;
  iMesh_getEntSets(mesh, root_set, 1, &sets, &sets_alloc, &sets_size, &err);
  ERRORR("Couldn't get all sets.");

    // iterate through them, checking whether they have tags
  iBase_TagHandle *tags = NULL;
  int tags_alloc = 0, tags_size;
  int i, j;
  for (i = 0; i < sets_size; i++) {
      // get connectivity
    iMesh_getAllEntSetTags(mesh, sets[i], &tags, &tags_alloc, &tags_size, &err);
    ERRORR("Failed to get ent set tags.");

    if (0 != tags_size) {
      std::cout << "Set " << sets[i] << "; Tags:" << std::endl;
      
        // list tag names on this set
      for (j = 0; j < tags_size; j++) {
        char tname[128];
        std::vector<int> int_val;
        std::vector<double> dbl_val;
        iMesh_getTagName(mesh, tags[j], tname, &err, sizeof(tname));
        tname[sizeof(tname)-1] = '\0';
        int tag_type, tag_size;
        iMesh_getTagType(mesh, tags[j], &tag_type, &err);
        ERRORR("Failed to get tag type.");
        iMesh_getTagType(mesh, tags[j], &tag_size, &err);
        ERRORR("Failed to get tag size.");
        if (iBase_INTEGER == tag_type) {
          int_val.resize(tag_size);
          iMesh_getEntSetIntData(mesh, sets[i], tags[j], &int_val[0], &err);
          ERRORR("Failed to get int data type.");
          std::cout << tname << " = " << int_val[0] ;
          if (tag_size > 1) 
            for (int k = 1; k < tag_size; i++) std::cout << ", " << int_val[k];
          std::cout << std::endl;
        }
        else if (iBase_DOUBLE == tag_type) {
          dbl_val.resize(tag_size);
          iMesh_getEntSetDblData(mesh, sets[i], tags[j], &dbl_val[0], &err);
          std::cout << tname << " = " << dbl_val[0] ;
          if (tag_size > 1) 
            for (int k = 1; k < tag_size; i++) std::cout << ", " << dbl_val[k];
          std::cout << std::endl;
        }
        else std::cout << tname << " (opaque) " << std::endl;
      }
    }
    std::cout << std::endl;
    
    free(tags);
    tags = NULL;
    tags_alloc = 0;
  }
  
  free(sets);

  iMesh_dtor(mesh, &err);
  ERRORR("Failed to destruct interface.");
  
  return 0;
}
