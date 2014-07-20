/*
 *  This program updates a manufactured tracer field from time T0 to time T1, in parallel.
 *  Input: arrival mesh, already distributed on processors, and a departure position for
 *  each vertex, saved in a tag DP
 */
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "moab_mpi.h"
#include "iMeshP.h"

#define IMESH_ASSERT(ierr) if (ierr!=0) printf("imesh assert\n");
#define IMESH_NULL 0

#ifdef __cplusplus
extern "C" {
#endif
void create_mesh(iMesh_Instance instance,
    iBase_EntitySetHandle * imesh_euler_set, double * icoords, int * icorners,
    int nc2, int nelem, int * ierr);

#ifdef __cplusplus
} // extern "C"
#endif

int main(int argc, char* argv[]) {
  MPI_Init(&argc, &argv);

  iMesh_Instance imesh;
  iMeshP_PartitionHandle partn;
  int ierr;

  iBase_EntitySetHandle root;
  imesh = IMESH_NULL;
  iMesh_newMesh(0, &imesh, &ierr, 0);
  IMESH_ASSERT(ierr);
  iMesh_getRootSet(imesh, &root, &ierr);
  IMESH_ASSERT(ierr);

  iMeshP_createPartitionAll(imesh, MPI_COMM_WORLD, &partn, &ierr);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  IMESH_ASSERT(ierr);

  iBase_EntitySetHandle euler_set;

  /*iMesh_createEntSet(imesh, 0, &euler_set, &ierr);
   IMESH_ASSERT(ierr);*/

  // read connectivity and fvm coords from file
  FILE * finp;
  if (argc > 1)
    finp = fopen(argv[1], "r");
  else
    finp = fopen("meshinp.txt", "r");
  if (finp == NULL) {
    fprintf(stderr, "Can't open input file %s\n", argv[1]);
    exit(1);
  }
  int nelem = 0, nc = 0;
  int readvals = fscanf(finp, "%d %d ", &nelem, &nc);

  if (readvals != 2)
    exit(1);

  double * coords = (double*) malloc(
      nelem * (nc + 1) * (nc + 1) * 3 * sizeof(double));
  int * corners = (int*) malloc(nelem * 4 * sizeof(int));
  int ixc = 0, ico = 0;
  for (int i = 0; i < nelem; i++) {
    fscanf(finp, "%d %d %d %d", &corners[ixc], &corners[ixc + 1],
        &corners[ixc + 2], &corners[ixc + 3]);
    ixc += 4;
    for (int i1 = 0; i1 <= nc; i1++) {
      for (int j = 0; j <= nc; j++) {
        fscanf(finp, "%lf %lf %lf", &coords[ico], &coords[ico + 1],
            &coords[ico + 2]);
        ico += 3;
      }
    }
  }

  create_mesh(imesh, &euler_set, coords, corners, nc, nelem, &ierr);

  // write everything
  const char * out_name = "out.h5m";
  const char optionswrite[] = " moab:PARALLEL=WRITE_PART ";
  iMeshP_saveAll(imesh, partn, euler_set, out_name, optionswrite, &ierr,
      strlen(out_name), strlen(optionswrite));
  IMESH_ASSERT(ierr);

  if (0 == rank)
    printf("Done\n");
  MPI_Finalize(); //probably the 4th time this is called.. no big deal

}
