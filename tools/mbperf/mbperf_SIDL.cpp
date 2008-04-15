/** iMesh Mesh Interface brick mesh performance test
 * 
 * This program tests iMesh mesh interface functions used to create a square structured
 * mesh.  Boilerplate taken from iMesh mesh interface test in MOAB and performance test in MOAB
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
#include <time.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <assert.h>

#include "iBase.hh"
#include "iMesh.hh"

#define CAST_INTERFACE(var_in, var_out, iface) \
          iBase::iface var_out; \
          try {var_out = var_in;}\
          catch(iBase::Error) {\
            cerr << "Error: current interface doesn't support iface." << endl; \
            return;}

#define CAST_MINTERFACE(var_in, var_out, iface) \
          iMesh::iface var_out; \
          try {var_out = var_in;}\
          catch(iBase::Error) {\
            cerr << "Error: current interface doesn't support iface." << endl; \
            return;}

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
void query_vert_to_elem(iMesh::Mesh &mesh);
void query_elem_to_vert(iMesh::Mesh &mesh);
void print_time(const bool print_em, double &tot_time, double &utime, double &stime,
                double &mem);
void build_connect(const int nelem, const int vstart, int *&connect);
void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem);
void testB(iMesh::Mesh &mesh, 
           const int nelem, sidl::array<double> &coords,
           const int *connect);
void testC(iMesh::Mesh &mesh, const int nelem, sidl::array<double> &coords);
void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride);
void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2);
void build_coords(const int nelem, sidl::array<double> &coords);
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

    // pre-build the coords array
  sidl::array<double> coords;
  build_coords(nelem, coords);
  assert(NULL != coords);

  int *connect = NULL;
  build_connect(nelem, 1, connect);

    // create an implementation
  iMesh::Mesh mesh;
  try {
    iMesh::Factory::newMesh("", mesh);
  } catch (iBase::Error err) 
  {
    std::cerr << "Trouble creating imesh instance." << std::endl;
    return 1;
  }
  
  
  switch (which_test) {
    case 'B':
        // test B: create mesh using bulk interface
      testB(mesh, nelem, coords, connect);
      break;
      
    case 'C':
    // test C: create mesh using individual interface
      testC(mesh, nelem, coords);
      break;
  }
  
  return 0;
}

void testB(iMesh::Mesh &mesh, 
           const int nelem, sidl::array<double> &coords,
           const int *connect) 
{
  double utime, stime, ttime0, ttime1, ttime2, ttime3, mem;
  
  print_time(false, ttime0, utime, stime, mem);
  std::cout << "Ready to read model into MOAB; memory = " << mem/1.0e6 << " MB." << std::endl;

  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  int num_elems = nelem*nelem*nelem;
  
    // create vertices as a block
  CAST_MINTERFACE(mesh, mesh_arrmod, ArrMod);
  sidl::array<Entity_Handle> sidl_vertices, dum_handles;
  CHECK_SIZE(dum_handles, num_verts);
  int sidl_vertices_size;
  
  try {
    mesh_arrmod.createVtxArr(num_verts, iBase::StorageOrder_BLOCKED,
                             coords, 3*num_verts,
                             sidl_vertices, sidl_vertices_size);
  } catch (iBase::Error err) {
    cerr << "Couldn't create vertices in bulk call" << endl;
    return;
  }

    // need to explicitly fill connectivity array, since we don't know
    // the format of entity handles
  sidl::array<Entity_Handle> sidl_connect;
  int sidl_connect_size = 8 * num_elems;
  CHECK_SIZE(sidl_connect, 8*num_elems);
  Entity_Handle *sidl_connect_ptr = HANDLE_ARRAY_PTR(sidl_connect);

  Entity_Handle *sidl_vertices_ptr = HANDLE_ARRAY_PTR(sidl_vertices);
  for (int i = 0; i < sidl_connect_size; i++) {
      // use connect[i]-1 because we used starting vertex index (vstart) of 1
    assert(connect[i]-1 < num_verts);
    sidl_connect_ptr[i] = sidl_vertices_ptr[connect[i]-1];
  }
  
    // create the entities
  sidl::array<Entity_Handle> new_hexes;
  int new_hexes_size;
  sidl::array<iBase::CreationStatus> status;
  int status_size;
  
  try {
    mesh_arrmod.createEntArr(iMesh::EntityTopology_HEXAHEDRON, sidl_connect, 
                             sidl_connect_size, new_hexes, new_hexes_size,
                             status, status_size);
  } catch (iBase::Error err) {
    cerr << "Couldn't create hex elements in bulk call" << endl;
    return;
  }

  print_time(false, ttime1, utime, stime, mem);
  std::cout << "Read model into MOAB; memory = " << mem/1.0e6 << " MB." << std::endl;

    // query the mesh 2 ways
  query_elem_to_vert(mesh);

  print_time(false, ttime2, utime, stime, mem);
  std::cout << "After E-v query; memory = " << mem/1.0e6 << " MB." << std::endl;

  query_vert_to_elem(mesh);
  
  print_time(false, ttime3, utime, stime, mem);
  std::cout << "After v-E query; memory = " << mem/1.0e6 << " MB." << std::endl;

  std::cout << "iMesh/MOAB ucd blocked: nelem, construct, e_to_v query, v_to_e query = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << " seconds" 
            << std::endl;
}

void testC(iMesh::Mesh &mesh, 
           const int nelem, sidl::array<double> &coords) 
{
  double utime, stime, ttime0, ttime1, ttime2, ttime3, mem;
  print_time(false, ttime0, utime, stime, mem);
  std::cout << "Ready to read data into MOAB; memory = " << mem/1.0e6 << " MB." << std::endl;

  CAST_MINTERFACE(mesh, mesh_arrmod, ArrMod);

    // need some dimensions
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int num_verts = numv*numv*numv;
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))

    // array to hold vertices created individually
  sidl::array<Entity_Handle> sidl_vertices;
  int sidl_vertices_size = num_verts;
  CHECK_SIZE(sidl_vertices, num_verts);

    // temporary array to hold vertex positions for single vertex
  sidl::array<double> tmp_coords;
  int tmp_coords_size = 3;
  CHECK_SIZE(tmp_coords, 3);
  double *dum_coords = ARRAY_PTR(tmp_coords, double);

    // get direct pointer to coordinate array
  double *coords_ptr = ARRAY_PTR(coords, double);
  
  for (int i = 0; i < num_verts; i++) {

      // temporary array to hold (single) vertices
    sidl::array<Entity_Handle> tmp_vertices;
    int tmp_vertices_size = 0;

      // create the vertex
    dum_coords[0] = coords_ptr[i]; 
    dum_coords[1] = coords_ptr[num_verts+i]; 
    dum_coords[2] = coords_ptr[2*num_verts+i];
    try {
      mesh_arrmod.createVtxArr(1, iBase::StorageOrder_BLOCKED,
                               tmp_coords, tmp_coords_size, 
                               tmp_vertices, tmp_vertices_size);
    } catch (iBase::Error err) {
      cerr << "Couldn't create vertex in individual call" << endl;
      return;
    }

      // assign into permanent vertex array
    sidl_vertices.set(i, tmp_vertices.get(0));
  }
  
    // get vertex array pointer for reading into tmp_conn
  Entity_Handle *tmp_sidl_vertices = HANDLE_ARRAY_PTR(sidl_vertices);
  
  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {

        sidl::array<Entity_Handle> tmp_conn;
        int tmp_conn_size = 8;
        CHECK_SIZE(tmp_conn, 8);

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
  
        sidl::array<Entity_Handle> new_hex;
        int new_hex_size = 0;
        sidl::array<iBase::CreationStatus> status;
        int status_size = 0;

        try {
          mesh_arrmod.createEntArr(iMesh::EntityTopology_HEXAHEDRON, 
                                   tmp_conn, tmp_conn_size,
                                   new_hex, new_hex_size,
                                   status, status_size);
        } catch (iBase::Error err) {
          cerr << "Couldn't create hex element in individual call" << endl;
          return;
        }
      }
    }
  }

  print_time(false, ttime1, utime, stime, mem);
  std::cout << "Read data into MOAB; memory  = " << mem/1.0e6 << " MB." << std::endl;

    // query the mesh 2 ways
  query_elem_to_vert(mesh);

  print_time(false, ttime2, utime, stime, mem);
  std::cout << "After E-v query; memory  = " << mem/1.0e6 << " MB." << std::endl;

  query_vert_to_elem(mesh);
  
  print_time(false, ttime3, utime, stime, mem);
  std::cout << "After v-E query; memory  = " << mem/1.0e6 << " MB." << std::endl;

  std::cout << "TSTT/MOAB ucd indiv: nelem, construct, e_to_v query, v_to_e query = " 
            << nelem << ", "
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << " seconds" 
            << std::endl;
}

void query_elem_to_vert(iMesh::Mesh &mesh)
{
  sidl::array<Entity_Handle> all_hexes;
  int all_hexes_size;
  CAST_MINTERFACE(mesh, mesh_ent, Entity);

    // get all the hex elements
  try {
    mesh.getEntities(0, iBase::EntityType_REGION, iMesh::EntityTopology_HEXAHEDRON, 
                     all_hexes, all_hexes_size);
  } catch (iBase::Error err) {
    cerr << "Couldn't get all hex elements in query_mesh" << endl;
    return;
  }

  try {
      // set up some tmp arrays and array ptrs
    Entity_Handle *all_hexes_ptr = HANDLE_ARRAY_PTR(all_hexes);  


      // now loop over elements
    for (int i = 0; i < all_hexes_size; i++) {
      sidl::array<int> dum_offsets;
      sidl::array<Entity_Handle> dum_connect;
      int dum_connect_size = 0;
        // get the connectivity of this element; will allocate space on 1st iteration,
        // but will have correct size on subsequent ones
      mesh_ent.getEntAdj(all_hexes_ptr[i], iBase::EntityType_VERTEX,
                     dum_connect, dum_connect_size);
                                          

        // get vertex coordinates; ; will allocate space on 1st iteration,
        // but will have correct size on subsequent ones
      sidl::array<double> dum_coords;
      int dum_coords_size = 0;
      iBase::StorageOrder order = iBase::StorageOrder_UNDETERMINED;
      mesh.getVtxArrCoords(dum_connect, dum_connect_size, order,
                           dum_coords, dum_coords_size);

      assert(24 == dum_coords_size && ARRAY_SIZE(dum_coords) == 24);
      double *dum_coords_ptr = ARRAY_PTR(dum_coords, double);
      double centroid[3] = {0.0, 0.0, 0.0};
      if (order == iBase::StorageOrder_BLOCKED) {
        for (int j = 0; j < 8; j++) {
          centroid[0] += dum_coords_ptr[j];
          centroid[1] += dum_coords_ptr[8+j];
          centroid[2] += dum_coords_ptr[16+j];
          centroid[0] += dum_coords.get(j);
          centroid[1] += dum_coords.get(8+j);
          centroid[2] += dum_coords.get(16+j);
        }
      }
      else {
        for (int j = 0; j < 8; j++) {
          centroid[0] += dum_coords_ptr[3*j];
          centroid[1] += dum_coords_ptr[3*j+1];
          centroid[2] += dum_coords_ptr[3*j+2];
        }
      }
    }
  } catch (iBase::Error err) {
    cerr << "Problem getting connectivity or vertex coords." << endl;
    return;
  }
}

void query_vert_to_elem(iMesh::Mesh &mesh)
{
  sidl::array<Entity_Handle> all_verts;
  int all_verts_size;
  CAST_MINTERFACE(mesh, mesh_ent, Entity);

    // get all the vertices elements
  try {
    mesh.getEntities(0, iBase::EntityType_VERTEX, 
                     iMesh::EntityTopology_POINT, all_verts, all_verts_size);
  } catch (iBase::Error err) {
    cerr << "Couldn't get all vertices in query_vert_to_elem" << endl;
    return;
  }

  try {
      // set up some tmp arrays and array ptrs
    Entity_Handle *all_verts_ptr = HANDLE_ARRAY_PTR(all_verts);  

      // now loop over vertices
    for (int i = 0; i < all_verts_size; i++) {
      sidl::array<Entity_Handle> dum_hexes;
      int dum_hexes_size;

        // get the connectivity of this element; will have to allocate space on every
        // iteration, since size can vary
      mesh_ent.getEntAdj(all_verts_ptr[i], iBase::EntityType_REGION,
                     dum_hexes, dum_hexes_size);
    }
  } catch (iBase::Error err) {
    cerr << "Problem getting connectivity or vertex coords." << endl;
    return;
  }
}

void print_time(const bool print_em, double &tot_time, double &utime, double &stime, double &mem) 
{
  get_time_mem(tot_time, utime, stime, mem);
  
  if (print_em) {
    std::cout << "User, system, total time = " << utime << ", " << stime 
              << ", " << tot_time << std::endl;
    std::cout << "Total memory = " << mem / 1.0e6 << " MB." << std::endl;
  }
}

void get_time_mem(double &tot_time, double &user_time,
                  double &sys_time, double &tot_mem) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  user_time = (double)r_usage.ru_utime.tv_sec +
    ((double)r_usage.ru_utime.tv_usec/1.e6);
  sys_time = (double)r_usage.ru_stime.tv_sec +
    ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = user_time + sys_time;
  tot_mem = 0;
  if (0 != r_usage.ru_maxrss) {
    tot_mem = r_usage.ru_idrss; 
  }
  else {
      // this machine doesn't return rss - try going to /proc
      // print the file name to open
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    
    close(file_ptr);
    file_str[file_len] = '\0';
      // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    int num_fields = sscanf(file_str, 
                            "%d " // pid
                            "%s " // comm
                            "%c " // state
                            "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                            "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                            "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                            "%u %u " // timeout, itrealvalue
                            "%d " // starttime
                            "%u %u", // vsize, rss
                            &dum_int, 
                            dum_str, 
                            dum_str, 
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, 
                            &dum_int,
                            &vm_size, &rss);
    if (num_fields == 24)
      tot_mem = ((double)vm_size);
  }
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

void build_coords(const int nelem, sidl::array<double> &coords) 
{
  double ttime0, ttime1, utime1, stime1, mem1;
  print_time(false, ttime0, utime1, stime1, mem1);
  std::cout << "Before TFI, total memory = " << mem1 / 1.0e6 << " MB." << std::endl;

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
    // initialize positions of 8 corners
  coords_ptr[VINDEX(0,0,0)] = coords_ptr[VINDEX(0,0,0)+nelem+1] = coords_ptr[VINDEX(0,0,0)+2*(nelem+1)] = 0.0;
  coords_ptr[VINDEX(0,nelem,0)] = coords_ptr[VINDEX(0,nelem,0)+2*(nelem+1)] = 0.0; coords_ptr[VINDEX(0,nelem,0)+nelem+1] = LENGTH;
  coords_ptr[VINDEX(0,0,nelem)] = coords_ptr[VINDEX(0,0,nelem)+nelem+1] = 0.0; coords_ptr[VINDEX(0,0,nelem)+2*(nelem+1)] = LENGTH;
  coords_ptr[VINDEX(0,nelem,nelem)] = 0.0; coords_ptr[VINDEX(0,nelem,nelem)+nelem+1] = coords_ptr[VINDEX(0,nelem,nelem)+2*(nelem+1)] = LENGTH;
  coords_ptr[VINDEX(nelem,0,0)] = LENGTH; coords_ptr[VINDEX(nelem,0,0)+nelem+1] = coords_ptr[VINDEX(nelem,0,0)+2*(nelem+1)] = 0.0;
  coords_ptr[VINDEX(nelem,0,nelem)] = coords_ptr[VINDEX(nelem,0,nelem)+2*(nelem+1)] = LENGTH; coords_ptr[VINDEX(nelem,0,nelem)+nelem+1] = 0.0;
  coords_ptr[VINDEX(nelem,nelem,0)] = coords_ptr[VINDEX(nelem,nelem,0)+nelem+1] = LENGTH; coords_ptr[VINDEX(nelem,nelem,0)+2*(nelem+1)] = 0.0;
  coords_ptr[VINDEX(nelem,nelem,nelem)] = coords_ptr[VINDEX(nelem,nelem,nelem)+nelem+1] = coords_ptr[VINDEX(nelem,nelem,nelem)+2*(nelem+1)] = LENGTH;

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
  print_time(false, ttime1, utime1, stime1, mem1);
  std::cout << "MOAB: TFI time = " << ttime1-ttime0 << " sec, memory = " 
            << mem1/1.0e6 << " MB." << std::endl;
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
