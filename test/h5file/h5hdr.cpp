#include "moab/H5MInterface.hpp"
#include <stdlib.h>
#include <stdio.h>

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

using namespace moab;
int main( int argc, char* argv[] )
{

  const char * filename =STRINGIFY(MESHDIR) "/64bricks_12ktet.h5m";
  if (argc != 2) {
    fprintf( stderr,"Usage: %s <filename>\n", argv[0] );
    fprintf(stderr, "use default file: %s\n", filename);
  }
  else
  {
    filename = argv[1];
  }
  H5MInterface h5intf(filename);
  int nv, ne, nf, nr, dim, npart;
  if (h5intf.get_mesh_info(&nv, &ne, &nf, &nr, &dim, &npart))
    return 1;

  printf("number of vertices: %d\n", nv);
  printf("number of edges: %d\n", ne);
  printf("number of faces: %d\n", nf);
  printf("number of regions: %d\n", nr);
  printf("number of dimensions: %d\n", dim);
  printf("number of partitions: %d\n", npart);
  return 0;
}
