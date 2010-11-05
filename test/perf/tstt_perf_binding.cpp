/** TSTT Mesh Interface brick mesh performance test
 * 
 * This program tests TSTT mesh interface functions used to create a square structured
 * mesh.  Boilerplate taken from tstt mesh interface test in MOAB and performance test in MOAB
 *
 */

// Different platforms follow different conventions for usage
#ifndef NT
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

#include <iostream>
#include <assert.h>
#include "iMesh.h"

// needed to get the proper size for handles

using namespace std;

double LENGTH = 1.0;

// forward declare some functions
void query_elem_to_vert(iMesh_Instance mesh);
void query_vert_to_elem(iMesh_Instance mesh);
void print_time(const bool print_em, double &tot_time, double &utime, double &stime,
                long &imem, long &rmem);
void build_connect(const int nelem, const int vstart, int *&connect);
void testB(iMesh_Instance mesh, const int nelem, const double *coords, int *connect);
void testC(iMesh_Instance mesh, const int nelem, const double *coords);
void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride);
void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2);
void build_coords(const int nelem, double *&coords);
void build_connect(const int nelem, const int vstart, int *&connect);

int main( int argc, char *argv[] )
{
  int nelem = 20;
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <ints_per_side> <A|B|C>" << std::endl;
    return 1;
  }
  
  char which_test = '\0';
  
  sscanf(argv[1], "%d", &nelem);
  sscanf(argv[2], "%c", &which_test);

  if (which_test != 'B' && which_test != 'C') {
      std::cout << "Must indicate B or C for test." << std::endl;
      return 1;
  }
  
  std::cout << "number of elements: " << nelem << "; test " 
            << which_test << std::endl;

    // initialize the data in native format

    // pre-build the coords array
  double *coords;
  build_coords(nelem, coords);
  assert(NULL != coords);

    // test B: create mesh using bulk interface

    // create an implementation
  iMesh_Instance mesh;
  int result;
  iMesh_newMesh(NULL, &mesh, &result, 0);
  int *connect = NULL;

  if (iBase_SUCCESS != result) {
    cerr << "Couldn't create mesh instance." << endl;
    iMesh_dtor(mesh, &result);
    return 1;
  }

  iMesh_setGeometricDimension(mesh, 3, &result);
  if (iBase_SUCCESS != result) {
    cerr << "Couldn't set geometric dimension." << endl;
    iMesh_dtor(mesh, &result);
    return 1;
  }
  
  
  switch (which_test) {
    case 'B':
        build_connect(nelem, 1, connect);

        // test B: create mesh using bulk interface
        testB(mesh, nelem, coords, connect);
        break;
      
    case 'C':
    // test C: create mesh using individual interface
      testC(mesh, nelem, coords);
      break;
  }

  free(coords);
  
  return 0;
}

void testB(iMesh_Instance mesh, 
           const int nelem, const double *coords,
           int *connect) 
{
  double utime, stime, ttime0, ttime1, ttime2, ttime3, ttime4;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4;
  
  print_time(false, ttime0, utime, stime, imem0, rmem0);
  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  int num_elems = nelem*nelem*nelem;
  
    // create vertices as a block; initialize to NULL so allocation is done in interface
  iBase_EntityHandle *vertices = NULL;
  int vertices_allocated = 0, vertices_size;
  int result;
  iMesh_createVtxArr(mesh, num_verts, iBase_BLOCKED, 
                     coords, 3*num_verts, 
                     &vertices, &vertices_allocated, &vertices_size, &result);
  if (iBase_SUCCESS != result) {
    cerr << "Couldn't create vertices in bulk call" << endl;
    free(vertices);
    return;
  }

    // need to explicitly fill connectivity array, since we don't know
    // the format of entity handles
  int nconnect = 8 * num_elems;
  iBase_EntityHandle *sidl_connect = (iBase_EntityHandle*) malloc(nconnect*sizeof(iBase_EntityHandle));
  
  for (int i = 0; i < nconnect; i++) {
      // use connect[i]-1 because we used starting vertex index (vstart) of 1
    assert(connect[i]-1 < num_verts);
    sidl_connect[i] = vertices[connect[i]-1];
  }

    // no longer need vertices and connect arrays, free here to reduce overall peak memory usage
  free(vertices);
  free(connect);
  
    // create the entities
  iBase_EntityHandle *new_hexes = NULL;
  int new_hexes_allocated = 0, new_hexes_size;
  int *status = NULL;
  int status_allocated = 0, status_size;
  
  iMesh_createEntArr(mesh, iMesh_HEXAHEDRON, sidl_connect, nconnect, 
                     &new_hexes, &new_hexes_allocated, &new_hexes_size,
                     &status, &status_allocated, &status_size, &result);
  if (iBase_SUCCESS != result) {
    cerr << "Couldn't create hex elements in bulk call" << endl;
    free(new_hexes);
    free(status);
    free(sidl_connect);
    return;
  }

  print_time(false, ttime1, utime, stime, imem1, rmem1);

  free(status);
  free(new_hexes);
  free(sidl_connect);

    // query the mesh 2 ways
  query_elem_to_vert(mesh);

  print_time(false, ttime2, utime, stime, imem2, rmem2);

  query_vert_to_elem(mesh);
  
  print_time(false, ttime3, utime, stime, imem3, rmem3);

  iMesh_dtor(mesh, &result);

  print_time(false, ttime4, utime, stime, imem4, rmem4);

  std::cout << "TSTTb/MOAB ucd blocked: nelem, construct, e_to_v query, v_to_e query, after dtor = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << ", " 
            << ttime4-ttime3 << " seconds" 
            << std::endl;
  std::cout << "TSTTb/MOAB ucd blocked memory (rss): initial, after v/e construction, e-v query, v-e query, after dtor:" 
            << rmem0 << ", " << rmem1 << ", " << rmem2 << ", " << rmem3 <<  ", " << rmem4 << " kb" << std::endl;

}

void testC(iMesh_Instance mesh, const int nelem, const double *coords) 
{
  double utime, stime, ttime0, ttime1, ttime2, ttime3, ttime4;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4;
  print_time(false, ttime0, utime, stime, imem0, rmem0);

    // need some dimensions
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int num_verts = numv*numv*numv;
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))

    // array to hold vertices created individually
  iBase_EntityHandle *sidl_vertices = (iBase_EntityHandle*) malloc(num_verts*sizeof(iBase_EntityHandle));
  int result;

  for (int i = 0; i < num_verts; i++) {
      // create the vertex
    iMesh_createVtx(mesh, coords[i], coords[i+num_verts], coords[i+2*num_verts],
                    sidl_vertices+i, &result);
    if (iBase_SUCCESS != result) {
      cerr << "Couldn't create vertex in individual call" << endl;
      return;
    }

  }
  
  iBase_EntityHandle tmp_conn[8], new_hex;

  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        int vijk = VINDEX(i,j,k);
        tmp_conn[0] = sidl_vertices[vijk];
        tmp_conn[1] = sidl_vertices[vijk+1];
        tmp_conn[2] = sidl_vertices[vijk+1+numv];
        tmp_conn[3] = sidl_vertices[vijk+numv];
        tmp_conn[4] = sidl_vertices[vijk+numv*numv];
        tmp_conn[5] = sidl_vertices[vijk+1+numv*numv];
        tmp_conn[6] = sidl_vertices[vijk+1+numv+numv*numv];
        tmp_conn[7] = sidl_vertices[vijk+numv+numv*numv];
        
          // create the entity
  
        int status;
        iMesh_createEnt(mesh, iMesh_HEXAHEDRON, tmp_conn, 8, 
                        &new_hex, &status, &result);
        if (iBase_SUCCESS != result) {
          cerr << "Couldn't create hex element in individual call" << endl;
          return;
        }
      }
    }
  }

  print_time(false, ttime1, utime, stime, imem1, rmem1);

  free(sidl_vertices);
  
    // query the mesh 2 ways
  query_elem_to_vert(mesh);

  print_time(false, ttime2, utime, stime, imem2, rmem2);

  query_vert_to_elem(mesh);
  
  print_time(false, ttime3, utime, stime, imem3, rmem3);

  iMesh_dtor(mesh, &result);

  print_time(false, ttime4, utime, stime, imem4, rmem4);

  std::cout << "TSTTb/MOAB ucd indiv: nelem, construct, e_to_v query, v_to_e query, after dtor = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << ", " 
            << ttime4-ttime3 << " seconds" 
            << std::endl;
  std::cout << "TSTTb/MOAB ucd indiv memory (rss): initial, after v/e construction, e-v query, v-e query, after dtor:" 
            << rmem0 << ", " << rmem1 << ", " << rmem2 << ", " << rmem3 <<  ", " << rmem4 << " kb" << std::endl;

}

void query_elem_to_vert(iMesh_Instance mesh)
{
  iBase_EntityHandle *all_hexes = NULL;
  int all_hexes_size, all_hexes_allocated = 0;

    // get all the hex elements
  int success;
  iBase_EntitySetHandle root_set;
  iMesh_getRootSet(mesh, &root_set, &success);
  if (iBase_SUCCESS != success) {
    cerr << "Couldn't get root set." << endl;
    return;
  }
  
  iMesh_getEntities(mesh, root_set, iBase_REGION, 
                    iMesh_HEXAHEDRON, 
                    &all_hexes, &all_hexes_allocated, 
                    &all_hexes_size, &success);
  if (iBase_SUCCESS != success) {
    cerr << "Couldn't get all hex elements in query_mesh" << endl;
    return;
  }

    // now loop over elements
  iBase_EntityHandle *dum_connect = NULL;
  int dum_connect_allocated = 0, dum_connect_size;
  double *dum_coords = NULL;
  int dum_coords_size, dum_coords_allocated = 0;
  int order;
  iMesh_getDfltStorage(mesh, &order, &success);
  if (iBase_SUCCESS != success) return;

  for (int i = 0; i < all_hexes_size; i++) {
      // get the connectivity of this element; will allocate space on 1st iteration,
      // but will have correct size on subsequent ones
    iMesh_getEntAdj(mesh, all_hexes[i], iBase_VERTEX,
                    &dum_connect, &dum_connect_allocated, &dum_connect_size, 
                    &success);

    if (iBase_SUCCESS == success) {
        // get vertex coordinates; ; will allocate space on 1st iteration,
        // but will have correct size on subsequent ones
      iMesh_getVtxArrCoords(mesh, dum_connect, dum_connect_size, 
                            order,
                            &dum_coords, &dum_coords_allocated, 
                            &dum_coords_size, &success);

      double centroid[3] = {0.0, 0.0, 0.0};
      if (order == iBase_BLOCKED) {
        for (int j = 0; j < 8; j++) {
          centroid[0] += dum_coords[j];
          centroid[1] += dum_coords[8+j];
          centroid[2] += dum_coords[16+j];
        }
      }
      else {
        for (int j = 0; j < 8; j++) {
          centroid[0] += dum_coords[3*j];
          centroid[1] += dum_coords[3*j+1];
          centroid[2] += dum_coords[3*j+2];
        }
      }
    }
      
    if (iBase_SUCCESS != success) {
      cerr << "Problem getting connectivity or vertex coords." << endl;
      return;
    }
  }
  
  free(all_hexes);
  free(dum_connect);
  free(dum_coords);
}

void query_vert_to_elem(iMesh_Instance mesh)
{
  iBase_EntityHandle *all_verts = NULL;
  int all_verts_allocated = 0, all_verts_size;

  iBase_EntitySetHandle root_set;
  int success;
  iMesh_getRootSet(mesh, &root_set, &success);
  if (iBase_SUCCESS != success) {
    cerr << "Couldn't get root set." << endl;
    return;
  }

    // get all the vertices elements
  iMesh_getEntities(mesh, root_set, iBase_VERTEX, 
                    iMesh_POINT, &all_verts, 
                    &all_verts_allocated, &all_verts_size, &success);
  if (iBase_SUCCESS != success) {
    cerr << "Couldn't get all vertices in query_vert_to_elem" << endl;
    return;
  }

    // for this mesh, should never be more than 8 hexes connected to a vertex
  iBase_EntityHandle *dum_hexes = (iBase_EntityHandle*) calloc(8, sizeof(iBase_EntityHandle));
  
  int dum_hexes_allocated = 8, dum_hexes_size;

    // now loop over vertices
  for (int i = 0; i < all_verts_size; i++) {

        // get the connectivity of this element; will have to allocate space on every
        // iteration, since size can vary
    iMesh_getEntAdj(mesh, all_verts[i], iBase_REGION,
                    &dum_hexes, &dum_hexes_allocated, &dum_hexes_size,
                    &success);
    if (iBase_SUCCESS != success) {
      cerr << "Problem getting connectivity or vertex coords." << endl;
      return;
    }
  }

  free(dum_hexes);
  free(all_verts);
}

void print_time(const bool print_em, double &tot_time, double &utime, double &stime,
                long &imem, long &rmem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  utime = (double)r_usage.ru_utime.tv_sec +
     ((double)r_usage.ru_utime.tv_usec/1.e6);
  stime = (double)r_usage.ru_stime.tv_sec +
     ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = utime + stime;
  if (print_em)
    std::cout << "User, system, total time = " << utime << ", " << stime 
              << ", " << tot_time << std::endl;

  imem = rmem = 0;
#ifndef LINUX
  imem = r_usage.ru_idrss;
  rmem = r_usage.ru_maxrss;
  std::cout << "Max resident set size = " << r_usage.ru_maxrss << " kbytes" << std::endl;
  std::cout << "Int resident set size = " << r_usage.ru_idrss << " kbytes" << std::endl;
#else
  system("ps o args,drs,rss | grep perf | grep -v grep");  // RedHat 9.0 doesnt fill in actual memory data 
#endif
    //delete [] hex_array;
}

void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride) 
{
  for (int i = 1; i < nelem; i++) {
    start[i*stride] = start[0]+i*xint;
    start[nelem+1+i*stride] = start[nelem+1]+i*xint;
    start[2*(nelem+1)+i*stride] = start[2*(nelem+1)]+i*xint;
  }
}

void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2) 
{
    // 2D TFI on a face starting at a, with strides stride1 in ada and stride2 in tse
  for (int j = 1; j < nelem; j++) {
    double tse = j * xint;
    for (int i = 1; i < nelem; i++) {
      double ada = i * xint;
      
      a[i*stride1+j*stride2] = (1.0 - ada)*a[i*stride1]
        + ada*a[i*stride1+nelem*stride2]
        + (1.0 - tse)*a[j*stride2]
        + tse*a[j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[0]
        - (1.0 - tse)*ada*a[nelem*stride1]
        - tse*(1.0 - ada)*a[nelem*stride2]
        - tse*ada*a[nelem*(stride1+stride2)];
      a[nelem+1+i*stride1+j*stride2] = (1.0 - ada)*a[nelem+1+i*stride1]
        + ada*a[nelem+1+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[nelem+1+j*stride2]
        + tse*a[nelem+1+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[nelem+1+0]
        - (1.0 - tse)*ada*a[nelem+1+nelem*stride1]
        - tse*(1.0 - ada)*a[nelem+1+nelem*stride2]
        - tse*ada*a[nelem+1+nelem*(stride1+stride2)];
      a[2*(nelem+1)+i*stride1+j*stride2] = (1.0 - ada)*a[2*(nelem+1)+i*stride1]
        + ada*a[2*(nelem+1)+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[2*(nelem+1)+j*stride2]
        + tse*a[2*(nelem+1)+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[2*(nelem+1)+0]
        - (1.0 - tse)*ada*a[2*(nelem+1)+nelem*stride1]
        - tse*(1.0 - ada)*a[2*(nelem+1)+nelem*stride2]
        - tse*ada*a[2*(nelem+1)+nelem*(stride1+stride2)];
    }
  }
}

void build_coords(const int nelem, double *&coords) 
{
  double ttime0, ttime1, utime1, stime1;
  long imem, rmem;
  print_time(false, ttime0, utime1, stime1, imem, rmem);

    // allocate the memory
  int numv = nelem+1;
  int numv_sq = numv*numv;
  int tot_numv = numv*numv*numv;
  coords = (double*) malloc(3*tot_numv*sizeof(double));

// use FORTRAN-like indexing
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))
  int idx;
  double scale1, scale2, scale3;
    // use these to prevent optimization on 1-scale, etc (real map wouldn't have
    // all these equal)
  scale1 = LENGTH/nelem;
  scale2 = LENGTH/nelem;
  scale3 = LENGTH/nelem;

#ifdef REALTFI
    // use a real TFI xform to compute coordinates
    // compute edges
    // i (stride=1)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,nelem)], nelem, scale1, 1);
    // j (stride=numv)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,nelem)], nelem, scale1, numv);
    // k (stride=numv^2)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,nelem,0)], nelem, scale1, numv_sq);

    // compute faces
    // i=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, numv, numv_sq);
  compute_face(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv, numv_sq);
    // j=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv_sq);
  compute_face(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1, numv_sq);
    // k=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv);
  compute_face(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1, numv);

    // initialize corner indices
  int i000 = VINDEX(0,0,0);
  int ia00 = VINDEX(nelem,0,0);
  int i0t0 = VINDEX(0,nelem,0);
  int iat0 = VINDEX(nelem,nelem,0);
  int i00g = VINDEX(0,0,nelem);
  int ia0g = VINDEX(nelem,0,nelem);
  int i0tg = VINDEX(0,nelem,nelem);
  int iatg = VINDEX(nelem,nelem,nelem);
  double cX, cY, cZ;
  int adaInts = nelem;
  int tseInts = nelem;
  int gammaInts = nelem;
  
  
  for (int i=1; i < nelem; i++) {
    for (int j=1; j < nelem; j++) {
      for (int k=1; k < nelem; k++) {
        idx = VINDEX(i,j,k);
        double tse = i*scale1;
        double ada = j*scale2;
        double gamma = k*scale3;
        double tm1 = 1.0 - tse;
        double am1 = 1.0 - ada;
        double gm1 = 1.0 - gamma;

        cX = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cY = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cZ = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        double *ai0k = &coords[VINDEX(k,0,i)];
        double *aiak = &coords[VINDEX(k,adaInts,i)];
        double *a0jk = &coords[VINDEX(k,j,0)];
        double *atjk = &coords[VINDEX(k,j,tseInts)];
        double *aij0 = &coords[VINDEX(0,j,i)];
        double *aijg = &coords[VINDEX(gammaInts,j,i)];
  
        coords[VINDEX(i,j,k)] = (   am1*ai0k[0] 
                                    + ada*aiak[0] 
                                    + tm1*a0jk[0] 
                                    + tse*atjk[0]
                                    + gm1*aij0[0] 
                                    + gamma*aijg[0] )/2.0 - cX/2.0;

        coords[nelem+1+VINDEX(i,j,k)] = (   am1*ai0k[nelem+1] 
                                            + ada*aiak[nelem+1] 
                                            + tm1*a0jk[nelem+1] 
                                            + tse*atjk[nelem+1]
                                            + gm1*aij0[nelem+1] 
                                            + gamma*aijg[nelem+1] )/2.0 - cY/2.0;

        coords[2*(nelem+1)+VINDEX(i,j,k)] = (   am1*ai0k[2*(nelem+1)] 
                                                + ada*aiak[2*(nelem+1)] 
                                                + tm1*a0jk[2*(nelem+1)] 
                                                + tse*atjk[2*(nelem+1)]
                                                + gm1*aij0[2*(nelem+1)] 
                                                + gamma*aijg[2*(nelem+1)] )/2.0 - cZ/2.0;
      }
    }
  }
  

#else
  for (int i=0; i < numv; i++) {
    for (int j=0; j < numv; j++) {
      for (int k=0; k < numv; k++) {
        idx = VINDEX(i,j,k);
          // blocked coordinate ordering
        coords[idx] = i*scale1;
        coords[tot_numv+idx] = j*scale2;
        coords[2*tot_numv+idx] = k*scale3;
      }
    }
  }
#endif

  print_time(false, ttime1, utime1, stime1, imem, rmem);
  std::cout << "TSTTbinding/MOAB: TFI time = " << ttime1-ttime0 << " sec" 
            << std::endl;
}

void build_connect(const int nelem, const int vstart, int *&connect) 
{
    // allocate the memory
  int nume_tot = nelem*nelem*nelem;
  connect = (int*) malloc(8*nume_tot*sizeof(int));

  int vijk;
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int idx = 0;
  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = vstart+VINDEX(i,j,k);
        connect[idx++] = vijk;
        connect[idx++] = vijk+1;
        connect[idx++] = vijk+1+numv;
        connect[idx++] = vijk+numv;
        connect[idx++] = vijk+numv*numv;
        connect[idx++] = vijk+1+numv*numv;
        connect[idx++] = vijk+1+numv+numv*numv;
        connect[idx++] = vijk+numv+numv*numv;
        assert(i <= numv*numv*numv);
      }
    }
  }
}
