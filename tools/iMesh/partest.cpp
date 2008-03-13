#include <iostream>
#include <stdio.h>
#include "mpi.h"
#include "iMesh.h"


#define IMESH_ASSERT(ierr) if (ierr!=0) printf("imesh assert\n");
#define IMESH_NULL 0

int main(int argc, char* argv[]){
  MPI_Init(&argc, &argv);
  printf("Hello\n");

  iMesh_Instance imesh;
  int ierr, num_sets;


  imesh = IMESH_NULL;
  iMesh_newMesh("PARALLEL", &imesh, &ierr, 8);
  IMESH_ASSERT(ierr);


  const char options[] = "PARALLEL=BCAST_DELETE;PARTITION=MATERIAL_SET;PARTITION_DISTRIBUTE;";
  const char filename[] = "64bricks_1mhex.h5m";

  iMesh_load(imesh,
             IMESH_NULL,
             filename,
             options,
             &ierr,
             strlen(filename),
             strlen(options));
  IMESH_ASSERT(ierr);

  
  iMesh_getNumEntSets(imesh,
                      IMESH_NULL,
                      1,
                      &num_sets,
                      &ierr);
  IMESH_ASSERT(ierr);
  printf("There's %d entity sets here\n", num_sets);


  printf("Done\n");
  MPI_Finalize(); //probably the 4th time this is called.. no big deal

}
