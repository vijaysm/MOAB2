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

// needed to get the proper size for handles
typedef void * Entity_Handle;

using namespace std;

#define ARRAY_PTR(array, type) reinterpret_cast<type*>(array._get_ior()->d_firstElement)
#define HANDLE_ARRAY_PTR(array) reinterpret_cast<Entity_Handle*>(array._get_ior()->d_firstElement)
#define ARRAY_SIZE(array) (array._is_nil() ? 0 : array.upper(0) - array.lower(0) + 1)
#define CHECK_SIZE(array, size)  \
          if (array._is_nil() || ARRAY_SIZE(array) == 0) array = array.create1d(size); \
          else if (ARRAY_SIZE(array) < size) { \
           cerr << "Array passed in is non-zero but too short." << endl; \
           assert(false); }


double LENGTH = 1.0;

// forward declare some functions
void query_mesh(TSTT::Mesh &mesh, const int nelem);
void print_time(const bool print_em, double &tot_time, double &utime, double &stime);
void build_connect(const int nelem, const int vstart, int *&connect);
void testB(TSTT::Mesh &mesh, 
           const int nelem, SIDL::array<double> &coords,
           const int *connect);
void testC(TSTT::Mesh &mesh, const int nelem, SIDL::array<double> &coords);
void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride);
void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2);
void build_coords(const int nelem, SIDL::array<double> &coords);
void build_connect(const int nelem, const int vstart, int *&connect);

int main( int argc, char *argv[] )
{
  // Check command line arg
  std::string filename;

  if (argc < 2) {
    cerr << "Usage: " << argv[0] << " <intervals_per_side>" << std::endl;
    return 1;
  }

  int nelem;
  sscanf(argv[1], "%d", &nelem);

    // initialize the data in native format

    // pre-build the coords array
  double ttime0, utime1, stime1, ttime1;
  double *coords;
  print_time(false, ttime0, utime1, stime1);
  build_coords(nelem, coords);
  print_time(false, ttime1, utime1, stime1);
  assert(NULL != coords);
  std::cout << "MOAB TFI time = " << ttime1-ttime0 << " seconds." << std::endl;

  int *connect = NULL;
  build_connect(nelem, 1, connect);

    // test B: create mesh using bulk interface

    // create an implementation
  MBMesh *mesh = new MBMesh();
  
    // construct the mesh
  testB(*mesh, nelem, coords, connect);

  delete mesh;

    // test C: create mesh using individual interface
  mesh = new MBMesh();
  testC(*mesh, nelem, coords);
  
  delete mesh;
}

void testB(MBMesh &mesh, 
           const int nelem, double *&coords,
           const int *connect) 
{
  double utime1, stime1, ttime1, utime2, stime2, ttime2, ttime3;
  
  print_time(false, ttime1, utime1, stime1);
  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  int num_elems = nelem*nelem*nelem;
  
    // create vertices as a block
  Entity_Handle *sidl_vertices, *dum_handles;
  CHECK_SIZE(dum_handles, num_verts);
  
  createEntities(mesh, TSTT::POINT, dum_handles, sidl_vertices);
  } catch (TSTT::Error err) {
    cerr << "Couldn't create vertices in bulk call" << endl;
    return;
  }

    // set vertex coordinates
  try {
    mesh_mm.setVertexCoordinates(sidl_vertices, TSTT::INTERLEAVED,
                                 coords);
  } catch (TSTT::Error err) {
    cerr << "Couldn't set vertex positions in bulk call" << endl;
    return;
  }

    // need to explicitly fill connectivity array, since we don't know
    // the format of entity handles
  SIDL::array<Entity_Handle> sidl_connect;
  CHECK_SIZE(sidl_connect, 8*num_elems);
  Entity_Handle *sidl_connect_ptr = HANDLE_ARRAY_PTR(sidl_connect);
  int nconnect = 8 * num_elems;
  Entity_Handle *sidl_vertices_ptr = HANDLE_ARRAY_PTR(sidl_vertices);
  
  for (int i = 0; i < nconnect; i++) {
      // use connect[i]-1 because we used starting vertex index (vstart) of 1
    assert(connect[i]-1 < num_verts);
    sidl_connect_ptr[i] = sidl_vertices_ptr[connect[i]-1];
  }
  
    // create the entities
  SIDL::array<Entity_Handle> new_hexes;
  
  try {
    mesh_mm.createEntities(TSTT::HEXAHEDRON, sidl_connect, new_hexes);
  } catch (TSTT::Error err) {
    cerr << "Couldn't create hex elements in bulk call" << endl;
    return;
  }

  print_time(false, ttime2, utime2, stime2);

  query_mesh(mesh, nelem);
  
  print_time(false, ttime3, utime2, stime2);

  std::cout << "TSTT/MOAB ucd blocked: construction = " << ttime2-ttime1 << " sec, query = " 
            << ttime3-ttime2 << " seconds." << std::endl;
  
}

void testC(TSTT::Mesh &mesh, 
           const int nelem, SIDL::array<double> &coords) 
{
  double utime1, stime1, ttime1, utime2, stime2, ttime2, ttime3;
  print_time(false, ttime1, utime1, stime1);

  CAST_INTERFACE(mesh, mesh_mm, ModifiableMesh);

    // need some dimensions
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int num_verts = numv*numv*numv;
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))

    // array to hold vertices created individually
  SIDL::array<Entity_Handle> sidl_vertices;
  CHECK_SIZE(sidl_vertices, num_verts);

    // temporary array to hold (single) vertices
  SIDL::array<Entity_Handle> tmp_vertices, dum_handles;
  CHECK_SIZE(dum_handles, 1);

    // temporary array to hold vertex positions for single vertex
  SIDL::array<double> tmp_coords;
  CHECK_SIZE(tmp_coords, 3);
  double *dum_coords = ARRAY_PTR(tmp_coords, double);

    // get direct pointer to coordinate array
  double *coords_ptr = ARRAY_PTR(coords, double);
  
  for (int i = 0; i < num_verts; i++) {

      // create the vertex
    try {
      mesh_mm.createEntities(TSTT::POINT, dum_handles, tmp_vertices);
    } catch (TSTT::Error err) {
      cerr << "Couldn't create vertex in individual call" << endl;
      return;
    }

    // set vertex coordinates
    dum_coords[0] = coords_ptr[i]; 
    dum_coords[1] = coords_ptr[num_verts+i]; 
    dum_coords[2] = coords_ptr[2*num_verts+i];
    try {
      mesh_mm.setVertexCoordinates(tmp_vertices, TSTT::BLOCKED,
                                   tmp_coords);
    } catch (TSTT::Error err) {
      cerr << "Couldn't set vertex position in individual call" << endl;
      return;
    }

      // assign into permanent vertex array
    sidl_vertices.set(i, tmp_vertices.get(0));
  }
  
  SIDL::array<Entity_Handle> tmp_conn;
  CHECK_SIZE(tmp_conn, 8);
  SIDL::array<Entity_Handle> new_hex;
  CHECK_SIZE(new_hex, 1);

    // get vertex array pointer for reading into tmp_conn
  Entity_Handle *tmp_sidl_vertices = HANDLE_ARRAY_PTR(sidl_vertices);
  
  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        int vijk = VINDEX(i,j,k);
        tmp_conn.set(0, tmp_sidl_vertices[vijk]);
        tmp_conn.set(1, tmp_sidl_vertices[vijk+1]);
        tmp_conn.set(2, tmp_sidl_vertices[vijk+1+numv]);
        tmp_conn.set(3, tmp_sidl_vertices[vijk+numv]);
        tmp_conn.set(4, tmp_sidl_vertices[vijk+numv*numv]);
        tmp_conn.set(5, tmp_sidl_vertices[vijk+1+numv*numv]);
        tmp_conn.set(6, tmp_sidl_vertices[vijk+1+numv+numv*numv]);
        tmp_conn.set(7, tmp_sidl_vertices[vijk+numv+numv*numv]);
        
          // create the entity
  
        try {
          mesh_mm.createEntities(TSTT::HEXAHEDRON, tmp_conn, new_hex);
        } catch (TSTT::Error err) {
          cerr << "Couldn't create hex element in individual call" << endl;
          return;
        }
      }
    }
  }

  print_time(false, ttime2, utime2, stime2);

  query_mesh(mesh, nelem);
  
  print_time(false, ttime3, utime2, stime2);

  std::cout << "TSTT/MOAB ucd indiv: construction = " << ttime2-ttime1 << " sec, query = " 
            << ttime3-ttime2 << " seconds." << std::endl;
  
}

void query_mesh(TSTT::Mesh &mesh, const int nelem)
{
  int nelem_tot = nelem*nelem*nelem;

  CAST_INTERFACE(mesh, mesh_ceq, CoreEntityQuery);
  CAST_INTERFACE(mesh, mesh_aesq, AdvancedEntitySetQuery);
  SIDL::array<Entity_Handle> all_hexes;

    // get all the hex elements
  try {
    mesh_ceq.entitysetGetEntities(0, TSTT::REGION, TSTT::HEXAHEDRON, all_hexes);
  } catch (TSTT::Error err) {
    cerr << "Couldn't get all hex elements in query_mesh" << endl;
    return;
  }

  try {
      // set up some tmp arrays and array ptrs
    Entity_Handle *all_hexes_ptr = HANDLE_ARRAY_PTR(all_hexes);  

    SIDL::array<int> dum_offsets;
    SIDL::array<Entity_Handle> dum_connect, dum_entity;
    CHECK_SIZE(dum_entity, 1);
    TSTT::StorageOrder order = TSTT::BLOCKED;
    SIDL::array<double> dum_coords;
    CHECK_SIZE(dum_coords, 24);
    double *dum_coords_ptr = ARRAY_PTR(dum_coords, double);

      // now loop over elements
    int num_hexes = ARRAY_SIZE(all_hexes);
    for (int i = 0; i < num_hexes; i++) {
      dum_entity.set(0, all_hexes_ptr[i]);

        // get the connectivity of this element; will allocate space on 1st iteration,
        // but will have correct size on subsequent ones
      mesh_aesq.entityGetAdjacencies(dum_entity, TSTT::VERTEX,
                                     dum_connect, dum_offsets);
                                          

        // get vertex coordinates; ; will allocate space on 1st iteration,
        // but will have correct size on subsequent ones
      mesh_ceq.entityGetVertexCoordinates(dum_connect, order,
                                          dum_coords);

      double centroid[3] = {0.0, 0.0, 0.0};
      if (order == TSTT::BLOCKED) {
        for (int j = 0; j < 8; j++) {
          centroid[0] += dum_coords_ptr[3*j];
          centroid[1] += dum_coords_ptr[3*j+1];
          centroid[2] += dum_coords_ptr[3*j+2];
        }
      }
      else {
        for (int j = 0; j < 8; j++) {
          centroid[0] += dum_coords_ptr[j];
          centroid[1] += dum_coords_ptr[8+j];
          centroid[2] += dum_coords_ptr[16+j];
        }
      }
    }
  } catch (TSTT::Error err) {
    cerr << "Problem getting connectivity or vertex coords." << endl;
    return;
  }
}

void print_time(const bool print_em, double &tot_time, double &utime, double &stime) 
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
#ifndef LINUX
 std::cout << "Max resident set size = " << r_usage.ru_maxrss*4096 << " bytes" << std::endl;
 std::cout << "Int resident set size = " << r_usage.ru_idrss << std::endl;
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

void build_coords(const int nelem, SIDL::array<double> &coords) 
{
    // allocate the memory
  int numv = nelem+1;
  int numv_sq = numv*numv;
  int tot_numv = numv*numv*numv;
  CHECK_SIZE(coords, 3*tot_numv);
  double *coords_ptr = ARRAY_PTR(coords, double);

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
  compute_edge(&coords_ptr[VINDEX(0,0,0)], nelem, scale1, 1);
  compute_edge(&coords_ptr[VINDEX(0,nelem,0)], nelem, scale1, 1);
  compute_edge(&coords_ptr[VINDEX(0,0,nelem)], nelem, scale1, 1);
  compute_edge(&coords_ptr[VINDEX(0,nelem,nelem)], nelem, scale1, 1);
    // j (stride=numv)
  compute_edge(&coords_ptr[VINDEX(0,0,0)], nelem, scale1, numv);
  compute_edge(&coords_ptr[VINDEX(nelem,0,0)], nelem, scale1, numv);
  compute_edge(&coords_ptr[VINDEX(0,0,nelem)], nelem, scale1, numv);
  compute_edge(&coords_ptr[VINDEX(nelem,0,nelem)], nelem, scale1, numv);
    // k (stride=numv^2)
  compute_edge(&coords_ptr[VINDEX(0,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords_ptr[VINDEX(nelem,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords_ptr[VINDEX(0,nelem,0)], nelem, scale1, numv_sq);
  compute_edge(&coords_ptr[VINDEX(nelem,nelem,0)], nelem, scale1, numv_sq);

    // compute faces
    // i=0, nelem
  compute_face(&coords_ptr[VINDEX(0,0,0)], nelem, scale1, numv, numv_sq);
  compute_face(&coords_ptr[VINDEX(nelem,0,0)], nelem, scale1, numv, numv_sq);
    // j=0, nelem
  compute_face(&coords_ptr[VINDEX(0,0,0)], nelem, scale1, 1, numv_sq);
  compute_face(&coords_ptr[VINDEX(0,nelem,0)], nelem, scale1, 1, numv_sq);
    // k=0, nelem
  compute_face(&coords_ptr[VINDEX(0,0,0)], nelem, scale1, 1, numv);
  compute_face(&coords_ptr[VINDEX(0,0,nelem)], nelem, scale1, 1, numv);

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

        cX = gm1 *   (am1*(tm1*coords_ptr[i000] + tse*coords_ptr[i0t0])  +
                      ada*(tm1*coords_ptr[ia00] + tse*coords_ptr[iat0])) +
          gamma * (am1*(tm1*coords_ptr[i00g] + tse*coords_ptr[i0tg])  +
                   ada*(tm1*coords_ptr[ia0g] + tse*coords_ptr[iatg]));

        cY = gm1 *   (am1*(tm1*coords_ptr[i000] + tse*coords_ptr[i0t0])  +
                      ada*(tm1*coords_ptr[ia00] + tse*coords_ptr[iat0])) +
          gamma * (am1*(tm1*coords_ptr[i00g] + tse*coords_ptr[i0tg])  +
                   ada*(tm1*coords_ptr[ia0g] + tse*coords_ptr[iatg]));

        cZ = gm1 *   (am1*(tm1*coords_ptr[i000] + tse*coords_ptr[i0t0])  +
                      ada*(tm1*coords_ptr[ia00] + tse*coords_ptr[iat0])) +
          gamma * (am1*(tm1*coords_ptr[i00g] + tse*coords_ptr[i0tg])  +
                   ada*(tm1*coords_ptr[ia0g] + tse*coords_ptr[iatg]));

        double *ai0k = &coords_ptr[VINDEX(k,0,i)];
        double *aiak = &coords_ptr[VINDEX(k,adaInts,i)];
        double *a0jk = &coords_ptr[VINDEX(k,j,0)];
        double *atjk = &coords_ptr[VINDEX(k,j,tseInts)];
        double *aij0 = &coords_ptr[VINDEX(0,j,i)];
        double *aijg = &coords_ptr[VINDEX(gammaInts,j,i)];
  
        coords_ptr[VINDEX(i,j,k)] = (   am1*ai0k[0] 
                                    + ada*aiak[0] 
                                    + tm1*a0jk[0] 
                                    + tse*atjk[0]
                                    + gm1*aij0[0] 
                                    + gamma*aijg[0] )/2.0 - cX/2.0;

        coords_ptr[nelem+1+VINDEX(i,j,k)] = (   am1*ai0k[nelem+1] 
                                            + ada*aiak[nelem+1] 
                                            + tm1*a0jk[nelem+1] 
                                            + tse*atjk[nelem+1]
                                            + gm1*aij0[nelem+1] 
                                            + gamma*aijg[nelem+1] )/2.0 - cY/2.0;

        coords_ptr[2*(nelem+1)+VINDEX(i,j,k)] = (   am1*ai0k[2*(nelem+1)] 
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
        coords_ptr[idx] = i*scale1;
        coords_ptr[tot_numv+idx] = j*scale2;
        coords_ptr[2*tot_numv+idx] = k*scale3;
      }
    }
  }
#endif
}

void build_connect(const int nelem, const int vstart, int *&connect) 
{
    // allocate the memory
  int nume_tot = nelem*nelem*nelem;
  connect = new int[8*nume_tot];

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
