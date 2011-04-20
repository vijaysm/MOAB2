#include <time.h>
#include <assert.h>
#include <iostream>
#include <sstream>
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

#define PRINT_SEQUENCE_COUNT

#ifdef PRINT_SEQUENCE_COUNT
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#  include "EntitySequence.hpp"
#  ifdef MB_ENTITY_SEQUENCE_HPP
#    include "EntitySequenceManager.hpp"
#  else
#    include "SequenceManager.hpp"
#  endif
#endif

using namespace moab;

  // constants
const bool dump_mesh = false;        //!< write mesh to vtk file
const int default_intervals = 25;    //!< defaul interval count for cubic structured hex mesh
const int default_query_count = 100; //!< number of times to do each query set
const int default_order[] = {0, 1, 2};
const int default_create[] = {0,1};
const int default_delete[] = {0,10,30,50,70,90};
#define ARRSIZE(A) (sizeof(A)/sizeof(A[0]))

  // input parameters
long numSideInt, numVert, numElem;   //!< total counts;
int queryCount;                      //!< number of times to do each query set

  // misc globals
Core mb_core;                         //!< moab instance
Interface& mb = mb_core;              //!< moab instance
EntityHandle vertStart, elemStart; //!< first handle
ReadUtilIface *readTool = 0;       
long* queryVertPermutation = 0;      //!< pupulated by init(): "random" order for vertices
long* queryElemPermutation = 0;      //!< pupulated by init(): "random" order for elements

//! Generate random permutation of values in [0,count-1]
long* permutation( long count )
{
  srand( count );
  long* array = new long[count];
  for (long i = 0; i < count; ++i)
    array[i] = i;
  
  for (long i = 0; i < count; ++i) {
    long r = rand();
    if (count > RAND_MAX) {
      r += RAND_MAX * rand();
      if (count/RAND_MAX > RAND_MAX) {
        long t = (long)RAND_MAX * rand();
        r += (long)RAND_MAX * t;
      }
    }
    std::swap( array[i], array[r%count] );
  }
  
  return array;
}

//! Initialize global variables
void init() 
{
  ErrorCode rval = mb.query_interface( readTool );
  if (rval || !readTool) {
    assert(false);
    abort();
  }
 
  queryVertPermutation = permutation( numVert );
  queryElemPermutation = permutation( numElem );
}


void create_vertices_single( ); //!< create vertices one at a time
void create_vertices_block( );  //!< create vertices in block using ReadUtilIface
void create_elements_single( ); //!< create elements one at a time
void create_elements_block( );  //!< create elements in block using ReadUtilIface
 
void forward_order_query_vertices(int percent); //!< calculate mean of all vertex coordinates
void reverse_order_query_vertices(int percent); //!< calculate mean of all vertex coordinates
void  random_order_query_vertices(int percent); //!< calculate mean of all vertex coordinates

void forward_order_query_elements(int percent); //!< check all element connectivity for valid vertex handles
void reverse_order_query_elements(int percent); //!< check all element connectivity for valid vertex handles
void  random_order_query_elements(int percent);  //!< check all element connectivity for valid vertex handles

void forward_order_query_element_verts(int percent); //!< calculate centroid
void reverse_order_query_element_verts(int percent); //!< calculate centroid
void  random_order_query_element_verts(int percent); //!< calculate centroid

void forward_order_delete_vertices( int percent ); //!< delete x% of vertices
void reverse_order_delete_vertices( int percent ); //!< delete x% of vertices
void  random_order_delete_vertices( int percent ); //!< delete x% of vertices

void forward_order_delete_elements( int percent ); //!< delete x% of elements
void reverse_order_delete_elements( int percent ); //!< delete x% of elements
void  random_order_delete_elements( int percent ); //!< delete x% of elements

void create_missing_vertices( int percent ); //!< re-create deleted vertices
void create_missing_elements( int percent ); //!< re-create deleted elements

#ifdef PRINT_SEQUENCE_COUNT
unsigned get_number_sequences( EntityType type );
#endif

/* Build arrays of function pointers, indexed by the order the entities are traversed in */

typedef void (*naf_t)();
typedef void (*iaf_t)(int);

iaf_t query_verts[3] = { &forward_order_query_vertices,
                         &reverse_order_query_vertices,
                         & random_order_query_vertices };

iaf_t query_elems[3] = { &forward_order_query_elements,
                         &reverse_order_query_elements,
                         & random_order_query_elements };

iaf_t query_elem_verts[3] = { &forward_order_query_element_verts,
                              &reverse_order_query_element_verts,
                              & random_order_query_element_verts };

iaf_t delete_verts[3] = { &forward_order_delete_vertices,
                          &reverse_order_delete_vertices,
                          & random_order_delete_vertices };

iaf_t delete_elems[3] = { &forward_order_delete_elements,
                          &reverse_order_delete_elements,
                          & random_order_delete_elements };

const char* order_strs[] = { "Forward", "Reverse", "Random" };

//! Coordinates for ith vertex in structured hex mesh
inline void vertex_coords( long vert_index, double& x, double& y, double& z );
//! Connectivity for ith hex in structured hex mesh
inline void element_conn( long elem_index, EntityHandle conn[8] );
//! True if passed index is one of the x% to be deleted
inline bool deleted_vert( long index, int percent );
//! True if passed index is one of the x% to be deleted
inline bool deleted_elem( long index, int percent );
//! if (deleted_vert(index,percent)) delete vertex
inline void delete_vert( long index, int percent );
//! if (deleted_elem(index,percent)) delete element
inline void delete_elem( long index, int percent );

//! print usage and exit
void usage() {
  std::cerr << "Usage: seqperf [-i <intervals>] [-o <order>] [-d <percent>] [-b|-s] [-q <count>]" << std::endl;
  std::cerr << " -i specify size of cubic structured hex mesh in intervals.  Default: " << default_intervals << std::endl;
  std::cerr << " -o one of \"forward\", \"reverse\", or \"random\".  May be specified multiple times.  Default is all." << std::endl;
  std::cerr << " -d percent of entities to delete.  May be specified multiple times.  Default is {";
  for (unsigned i = 0; i < ARRSIZE(default_delete)-1; ++i)
    std::cerr << default_delete[i] << ",";
  std::cerr << default_delete[ARRSIZE(default_delete)-1] << "}" << std::endl;
  std::cerr << " -b block creation of mesh" << std::endl;
  std::cerr << " -s single entity mesh creation" << std::endl;
  std::cerr << " -q number of times to repeat queries.  Default: " << default_query_count << std::endl;
  exit(1);
}

//! convert CPU time to string
std::string ts( clock_t t )
{
  std::ostringstream s;
  s << ((double)t)/CLOCKS_PER_SEC << 's';
  return s.str();
}

//! run function, printing time spent 
void TIME( const char* str, void (*func)() )
{ 
  std::cout << str << "... " << std::flush; 
  clock_t t = clock();
  (*func)();
  std::cout << ts(clock() - t) << std::endl;
}

//! run function query_repeat times, printing time spent 
void TIME_QRY( const char* str, void (*func)(int percent), int percent )
{ 
  std::cout << str << "... " << std::flush; 
  clock_t t = clock();
  for (int i = 0; i < queryCount; ++i)
    (*func)(percent);
  std::cout << ts(clock() - t) << std::endl;
}

//! run function with integer argument, printing time spent 
void TIME_DEL( const char* str, void (*func)(int), int percent )
{ 
  std::cout << str << "... " << std::flush; 
  clock_t t = clock();
  (*func)(percent);
  std::cout << ts(clock() - t) << std::endl;
}

//! call MB::delete_mesh().  function so can be passed to TIME
void delete_mesh()
{
  mb.delete_mesh();
}

//! Run a single combination of test parameters
void do_test( int create_mode, //!< 0 == single, 1 == block
              int order,       //!< 0 == forward, 1 == reverse, 2 == random
              int percent )    //!< percent of entities to delete
{
  clock_t t = clock();
  if (create_mode) {
    std::cout << "Block Entity Creation (all entities in single block of memory)" << std::endl;
    TIME ( "  Creating initial vertices", create_vertices_block );
    TIME ( "  Creating initial elements", create_elements_block );
    if (dump_mesh && !percent && mb.write_file( "seqperf.vtk" ) == MB_SUCCESS) 
      std::cout << "Wrote mesh to file: seqperf.vtk" << std::endl;
  }
  else {
    std::cout << "Single Entity Creation (entities grouped in memory blocks of constant size)" << std::endl;
    TIME ( "  Creating initial vertices", create_vertices_single );
    TIME ( "  Creating initial elements", create_elements_single );
  }
  
  std::cout << order_strs[order] <<
    " order with deletion of " << percent << "% of vertices and elements" << std::endl;

  TIME_DEL( "  Deleting elements", delete_elems[order], percent );
  TIME_DEL( "  Deleting vertices", delete_verts[order], percent );
  
  int num_vert = 0;
  int num_elem = 0;
  mb.get_number_entities_by_type( 0, MBVERTEX, num_vert );
  mb.get_number_entities_by_type( 0, MBHEX,    num_elem );
  std::cout << "  " << num_vert << " vertices and " << num_elem << " elements remaining" << std::endl;
#ifdef PRINT_SEQUENCE_COUNT
  std::cout << "  " << get_number_sequences(MBVERTEX) << " vertex sequences and "
            << get_number_sequences(MBHEX) << " element sequences." << std::endl;
#endif

  TIME_QRY( "  Querying vertex coordinates", query_verts[order], percent );
  TIME_QRY( "  Querying element connectivity", query_elems[order], percent );
  TIME_QRY( "  Querying element coordinates", query_elem_verts[order], percent );

  TIME_DEL( "  Re-creating vertices", create_missing_vertices, percent );
  TIME_DEL( "  Re-creating elements", create_missing_elements, percent );

  TIME( "  Clearing mesh instance", delete_mesh );
  
  std::cout << "Total time for test: " << ts(clock()-t) << std::endl << std::endl;
}

void parse_order( const char* str, std::vector<int>& list )
{
  if (str[0] == 'f') {
    if (strncmp( str, "forward", strlen(str) ))
      usage();
    list.push_back( 0 );
  }
  else if (str[0] != 'r') 
    usage();
  else if (str[1] == 'e') {
    if (strncmp( str, "reverse", strlen(str) ))
      usage();
    list.push_back( 0 );
  }
  else {
    if (strncmp( str, "random", strlen(str) ))
      usage();
    list.push_back( 0 );
  }
}

void parse_percent( const char* str, std::vector<int>& list )
{
  char* endptr;
  long p = strtol( str, &endptr, 0 );
  if (!endptr || *endptr || p < 0 || p > 100)
    usage();
  
  list.push_back((int)p);
}

int parse_positive_int( const char* str )
{
  char* endptr;
  long p = strtol( str, &endptr, 0 );
  if (!endptr || *endptr || p < 1)
    usage();
  int result = p;
  if (p != (long)result) // overflow
    usage();
  
  return result;
}

void check_default( std::vector<int>& list, const int* array, size_t array_len )
{
  if (list.empty())
    std::copy( array, array+array_len, std::back_inserter(list) );
}
    
  

int main( int argc, char* argv[] )
{
    // Parse arguments
  std::vector<int> createList, orderList, deleteList;
  numSideInt = default_intervals;
  queryCount = default_query_count;
  
  for (int i = 1; i < argc; ++i) {
      // check that arg is a '-' followed by a single character
    if (argv[i][0] != '-' || argv[i][1] == '\0' || argv[i][2] != '\0')
      usage();
  
    const char flag = argv[i][1];
    switch (flag) {
      case 'b': createList.push_back( 1 ); break;
      case 's': createList.push_back( 0 ); break;
      default: 
        if (++i == argc)
          usage();
        switch (flag) {
          case 'i': 
            numSideInt = parse_positive_int( argv[i] );
            break;
          case 'o':
            parse_order( argv[i], orderList );
            break;
          case 'd':
            parse_percent( argv[i], deleteList );
            break;
          case 'q':
            queryCount = parse_positive_int( argv[i] );
            break;
          default:
            usage();
        }
    }
  }
  check_default( createList, default_create, ARRSIZE(default_create) );
  check_default( orderList, default_order, ARRSIZE(default_order) );
  check_default( deleteList, default_delete, ARRSIZE(default_delete) );
  
    // Do some initialization.
  
  int numSideVert = numSideInt + 1;
  numVert = numSideVert * numSideVert * numSideVert;
  numElem = numSideInt * numSideInt * numSideInt;
  if (numVert / numSideVert / numSideVert != numSideVert) // overflow
    usage();
  init();
  
    // Echo input args
  
  std::cout << numSideInt << "x" << numSideInt << "x" << numSideInt << " hex grid: " 
            << numElem << " elements and " << numVert << " vertices" << std::endl;
  
    // Run tests
  
  std::vector<int>::const_iterator i,j,k;
  clock_t t = clock();
  for (std::vector<int>::iterator i = createList.begin(); i != createList.end(); ++i) {
    for (std::vector<int>::iterator j = deleteList.begin(); j != deleteList.end(); ++j) {
      for (std::vector<int>::iterator k = orderList.begin(); k != orderList.end(); ++k) {
        do_test( *i, *k, *j );
      }
    }
  }
  
    // Clean up
  
  std::cout << "TOTAL: " << ts(clock() - t) << std::endl << std::endl;
  delete [] queryVertPermutation;
  delete [] queryElemPermutation;
  return 0;
}


inline void vertex_coords( long vert_index, double& x, double& y, double& z )
{
  const long vs = numSideInt + 1;
  x = vert_index % vs;
  y = (vert_index / vs) % vs;
  z = (vert_index / vs / vs);
}

inline long vert_index( long x, long y, long z )
{
  const long vs = numSideInt + 1;
  return x + vs * (y + vs * z);
}

inline void element_conn( long elem_index, EntityHandle conn[8] )
{
  const long x = elem_index % numSideInt;
  const long y = (elem_index / numSideInt) % numSideInt;
  const long z = (elem_index / numSideInt / numSideInt);
  conn[0] = vertStart + vert_index(x  ,y  ,z  );
  conn[1] = vertStart + vert_index(x+1,y  ,z  );
  conn[2] = vertStart + vert_index(x+1,y+1,z  );
  conn[3] = vertStart + vert_index(x  ,y+1,z  );
  conn[4] = vertStart + vert_index(x  ,y  ,z+1);
  conn[5] = vertStart + vert_index(x+1,y  ,z+1);
  conn[6] = vertStart + vert_index(x+1,y+1,z+1);
  conn[7] = vertStart + vert_index(x  ,y+1,z+1);
}

inline bool deleted_vert( long index, int percent )
{
  return index%(numSideInt+1) >= (numSideInt+1)*(100-percent) / 100;
}

inline bool deleted_elem( long index, int percent )
{
  return index % numSideInt + 1 >= (numSideInt+1)*(100-percent) / 100;
}


void create_vertices_single( )
{
  double coords[3];
  vertex_coords( 0, coords[0], coords[1], coords[2] );
  ErrorCode rval = mb.create_vertex( coords, vertStart );
  assert(!rval);
  
  EntityHandle h;
  for (long i = 1; i < numVert; ++i) {
    vertex_coords( i, coords[0], coords[1], coords[2] );
    rval = mb.create_vertex( coords, h );
    assert(!rval);
    assert(h - vertStart == (EntityHandle)i);
  }
}

void create_vertices_block( )
{
  std::vector<double*> arrays;
  ErrorCode rval = readTool->get_node_coords( 3, numVert, 0, vertStart, arrays );
  if (rval || arrays.size() != 3) {
    assert(false);
    abort();
  }
  double *x = arrays[0], *y = arrays[1], *z = arrays[2];
  assert( x && y && z );
  
  for (long i = 0; i < numVert; ++i) 
    vertex_coords( i, *x++, *y++, *z++ );
}

void create_elements_single( )
{
  EntityHandle conn[8];
  element_conn( 0, conn );
  ErrorCode rval = mb.create_element( MBHEX, conn, 8, elemStart );
  if (rval) {
    assert(false);
    abort();
  }
  
  EntityHandle h;
  for (long i = 1; i < numElem; ++i) {
    element_conn( i, conn );
    rval = mb.create_element( MBHEX, conn, 8, h );
    assert(!rval);
    assert(h - elemStart == (EntityHandle)i);
  }
}


void create_elements_block( )
{
  EntityHandle* conn = 0;
  ErrorCode rval = readTool->get_element_connect( numElem, 8, MBHEX, 0, elemStart, conn );
  if (rval && !conn) {
    assert(false);
    abort();
  }
  
  for (long i = 0; i < numElem; ++i) 
    element_conn( i, conn + 8*i );
}
 
void forward_order_query_vertices(int percent)
{
  ErrorCode r;
  double coords[3];
  long x, y, z;
  const long vert_per_edge = numSideInt + 1;
  const long deleted_x = (numSideInt+1)*(100-percent) / 100;
  EntityHandle h = vertStart;
  for (z = 0; z < vert_per_edge; ++z) {
    for (y = 0; y < vert_per_edge; ++y) {
      for (x = 0; x < deleted_x; ++x, ++h) {
        r = mb.get_coords( &h, 1, coords );
        if (MB_SUCCESS != r) {
          assert(false);
          abort();
        }
      }
      h += (vert_per_edge - deleted_x);
    }
  }
}

void reverse_order_query_vertices(int percent)
{
  ErrorCode r;
  double coords[3];
  long x, y, z;
  const long vert_per_edge = numSideInt + 1;
  const long deleted_x = (numSideInt+1)*(100-percent) / 100;
  EntityHandle h = vertStart + numVert - 1;;
  for (z = vert_per_edge-1; z >= 0; --z) {
    for (y = vert_per_edge-1; y >= 0; --y) {
      h -= (vert_per_edge - deleted_x);
      for (x = deleted_x-1; x >= 0; --x, --h) {
        r = mb.get_coords( &h, 1, coords );
        assert(MB_SUCCESS == r);
      }
    }
  }
}

void random_order_query_vertices(int percent)
{
  ErrorCode r;
  EntityHandle h;
  double coords[3];
  for (long i = 0; i < numVert; ++i) {
    if (!deleted_vert(queryVertPermutation[i],percent)) {
      h = vertStart + queryVertPermutation[i];
      r = mb.get_coords( &h, 1, coords );
      assert(MB_SUCCESS == r);
    }
  }
}

void forward_order_query_elements(int percent)
{
  ErrorCode r;
  const EntityHandle* conn;
  int len;
  long x, y, z;
  const long elem_per_edge = numSideInt;
  const long deleted_x = (numSideInt+1)*(100-percent) / 100 - 1;
  EntityHandle h = elemStart;
  for (z = 0; z < elem_per_edge; ++z) {
    for (y = 0; y < elem_per_edge; ++y) {
      for (x = 0; x < deleted_x; ++x, ++h) {
        r = mb.get_connectivity( h, conn, len );
        assert(MB_SUCCESS == r);
        assert(conn && 8 == len);
      }
      h += (elem_per_edge - deleted_x);
    }
  }
}

void reverse_order_query_elements(int percent)
{
  ErrorCode r;
  const EntityHandle* conn;
  int len;
  long x, y, z;
  const long elem_per_edge = numSideInt;
  const long deleted_x = (numSideInt+1)*(100-percent) / 100 - 1;
  EntityHandle h = elemStart + numElem - 1;;
  for (z = elem_per_edge-1; z >= 0; --z) {
    for (y = elem_per_edge-1; y >= 0; --y) {
      h -= (elem_per_edge - deleted_x);
      for (x = deleted_x-1; x >= 0; --x, --h) {
        r = mb.get_connectivity( h, conn, len );
        assert(MB_SUCCESS == r);
        assert(conn && 8 == len);
      }
    }
  }
}

void  random_order_query_elements(int percent)
{
  ErrorCode r;
  const EntityHandle* conn;
  int len;
  for (long i = 0; i < numElem; ++i) {
    if (!deleted_elem( queryElemPermutation[i], percent )) {
      r = mb.get_connectivity( elemStart + queryElemPermutation[i], conn, len );
      assert(MB_SUCCESS == r);
      assert(conn && 8 == len);
    }
  }
}

/*
static double hex_centroid( double coords[24], double cent[3] )
{
  double a[3], b[3], c[3], vol;
  cent[0] += 0.125*(coords[0] + coords[3] + coords[6] + coords[ 9] + coords[12] + coords[15] + coords[18] + coords[21]);
  cent[1] += 0.125*(coords[1] + coords[4] + coords[7] + coords[10] + coords[13] + coords[16] + coords[19] + coords[22]);
  cent[2] += 0.125*(coords[2] + coords[5] + coords[8] + coords[11] + coords[14] + coords[17] + coords[20] + coords[23]);
  a[0] = coords[0] + coords[3] + coords[15] + coords[12] - coords[ 9] - coords[6] - coords[18] - coords[21];
  a[1] = coords[1] + coords[4] + coords[16] + coords[13] - coords[10] - coords[7] - coords[19] - coords[22];
  a[2] = coords[2] + coords[5] + coords[17] + coords[14] - coords[11] - coords[8] - coords[20] - coords[23];
  b[0] = coords[0] + coords[ 9] + coords[21] + coords[12] - coords[3] - coords[6] - coords[18] - coords[15];
  b[1] = coords[1] + coords[10] + coords[22] + coords[13] - coords[4] - coords[7] - coords[19] - coords[16];
  b[2] = coords[2] + coords[11] + coords[23] + coords[14] - coords[5] - coords[8] - coords[20] - coords[17];
  c[0] = coords[0] + coords[3] + coords[6] + coords[ 9] - coords[12] - coords[15] - coords[18] - coords[21];
  c[1] = coords[1] + coords[4] + coords[7] + coords[10] - coords[13] - coords[16] - coords[19] - coords[22];
  c[2] = coords[2] + coords[5] + coords[8] + coords[11] - coords[14] - coords[17] - coords[20] - coords[23];
  vol = c[0]*(a[1]*b[2] - a[2]*b[1]) + c[1]*(a[2]*b[0] - a[0]*b[2]) + c[2]*(a[0]*b[1] - a[1]*b[0]);
  return (1./64.) * vol;
}
*/

void forward_order_query_element_verts(int percent)
{
  ErrorCode r;
  const EntityHandle* conn;
  int len;
  long x, y, z;
  double coords[24];
  const long elem_per_edge = numSideInt;
  const long deleted_x = (numSideInt+1)*(100-percent) / 100 - 1;
  EntityHandle h = elemStart;
  for (z = 0; z < elem_per_edge; ++z) {
    for (y = 0; y < elem_per_edge; ++y) {
      for (x = 0; x < deleted_x; ++x, ++h) {
        r = mb.get_connectivity( h, conn, len );
        assert(MB_SUCCESS == r);
        assert(conn && 8 == len);
        r = mb.get_coords( conn, len, coords );
        assert(MB_SUCCESS == r );
      }
      h += (elem_per_edge - deleted_x);
    }
  }
}

void reverse_order_query_element_verts(int percent)
{
  ErrorCode r;
  const EntityHandle* conn;
  int len;
  long x, y, z;
  double coords[24];
  const long elem_per_edge = numSideInt;
  const long deleted_x = (numSideInt+1)*(100-percent) / 100 - 1;
  EntityHandle h = elemStart + numElem - 1;;
  for (z = elem_per_edge-1; z >= 0; --z) {
    for (y = elem_per_edge-1; y >= 0; --y) {
      h -= (elem_per_edge - deleted_x);
      for (x = deleted_x-1; x >= 0; --x, --h) {
        r = mb.get_connectivity( h, conn, len );
        assert(MB_SUCCESS == r);
        assert(conn && 8 == len);
        r = mb.get_coords( conn, len, coords );
        assert(MB_SUCCESS == r );
      }
    }
  }
}

void  random_order_query_element_verts(int percent)
{
  ErrorCode r;
  const EntityHandle* conn;
  int len;
  double coords[24];
  for (long i = 0; i < numElem; ++i) {
    if (!deleted_elem( queryElemPermutation[i], percent )) {
      r = mb.get_connectivity( elemStart + queryElemPermutation[i], conn, len );
      assert(MB_SUCCESS == r);
      assert(conn && 8 == len);
      r = mb.get_coords( conn, len, coords );
      assert(MB_SUCCESS == r );
    }
  }
}

void forward_order_delete_vertices( int percent )
{
  for (long i = 0; i < numVert; ++i)
    delete_vert( i, percent );
}

void reverse_order_delete_vertices( int percent )
{
  for (long i = numVert-1; i >= 0; --i)
    delete_vert( i, percent );
}

void  random_order_delete_vertices( int percent )
{
  for (long i = 0; i < numVert; ++i)
    delete_vert( queryVertPermutation[i], percent );
}

void forward_order_delete_elements( int percent )
{
  for (long i = 0; i < numElem; ++i)
    delete_elem( i, percent );
}

void reverse_order_delete_elements( int percent )
{
  for (long i = numElem-1; i >= 0; --i)
    delete_elem( i, percent );
}

void  random_order_delete_elements( int percent )
{
  for (long i = 0; i < numElem; ++i)
    delete_elem( queryElemPermutation[i], percent );
}


void create_missing_vertices( int percent )
{
  EntityHandle h;
  ErrorCode rval;
  double coords[3];
  for (long i = 0; i < numVert; ++i)
    if (deleted_vert( i, percent )) {
      vertex_coords( i, coords[0], coords[1], coords[2] );
      rval = mb.create_vertex( coords, h );
      assert(!rval);
    }
}

void create_missing_elements( int percent )
{
  EntityHandle h;
  ErrorCode rval;
  EntityHandle conn[8];
  for (long i = 0; i < numElem; ++i)
    if (deleted_elem( i, percent )) {
      element_conn( i, conn );
      rval = mb.create_element( MBHEX, conn, 8, h );
      if (rval) {
        assert(false);
        abort();
      }
    }
}

inline void delete_vert( long index, int percent )
{
  if (deleted_vert(index, percent)) {
    EntityHandle h = index + vertStart;
    ErrorCode rval = mb.delete_entities( &h, 1 );
    if (rval) {
      assert(false);
      abort();
    }
  }
}

inline void delete_elem( long index, int percent )
{
  if (deleted_elem(index, percent)) {
    EntityHandle h = index + elemStart;
    ErrorCode rval = mb.delete_entities( &h, 1 );
    if (rval) {
      assert(false);
      abort();
    }
  }
}

#ifdef PRINT_SEQUENCE_COUNT
unsigned get_number_sequences( EntityType type )
{
#ifdef MB_ENTITY_SEQUENCE_HPP
  return mb_core.sequence_manager()->entity_map(type)->size();
#else
  return mb_core.sequence_manager()->entity_map(type).get_sequence_count();
#endif
}
#endif
