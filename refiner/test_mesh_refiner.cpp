#include "MBCore.hpp"
#include "MBEdgeSizeSimpleImplicit.hpp"
#include "MBSimplexTemplateRefiner.hpp"
#include "MBInterface.hpp"
#include "MBParallelConventions.h"

#ifdef USE_MPI
#include "MBParallelComm.hpp"
#include "mpi.h"
#endif // USE_MPI

#include <iostream>
#include <map>

template< int _n >
class MBSplitVertexIndex
{
public:
  MBSplitVertexIndex() { }
  MBSplitVertexIndex( const MBEntityHandle* src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src[i]; std::sort( this->handles, this->handles + _n ); }
  MBSplitVertexIndex( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; }
  MBSplitVertexIndex& operator = ( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; return *this; }

  virtual bool operator < ( const MBSplitVertexIndex<_n>& other ) const
    {
    for ( int i = 0; i < _n; ++ i )
      if ( this->handles[i] < other.handles[i] )
        return true;
      else if ( this->handles[i] > other.handles[i] )
        return false;
    return true;
    }

  MBEntityHandle handles[_n];
};

class MBSplitVerticesBase
{
public:
  MBSplitVerticesBase( MBInterface* m )
    {
    this->mesh = m;
    }
  virtual bool find_or_create( const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle ) = 0;
  MBInterface* mesh;
};

template< int _n >
class MBSplitVertices : public std::map<MBSplitVertexIndex<_n>,MBEntityHandle>, public MBSplitVerticesBase
{
public:
  typedef std::map<MBSplitVertexIndex<_n>,MBEntityHandle> MapType;
  typedef typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator MapIteratorType;

  MBSplitVertices( MBInterface* m )
    : MBSplitVerticesBase( m )
    {
    }
  virtual bool find_or_create( const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle )
    {
    MapIteratorType it = this->find( MBSplitVertexIndex<_n>( split_src ) );
    if ( it == this->end() )
      {
      if ( this->mesh->create_vertex( coords, vert_handle ) != MB_SUCCESS )
        {
        return false;
        }
      return true;
      }
    vert_handle = it->second;
    return false;
    }
};


class MBTestOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  MBInterface* mesh;
  bool input_is_output;
  std::vector<MBSplitVerticesBase*> split_vertices;
  std::vector<MBEntityHandle> elem_vert;
  MBEdgeSizeEvaluator* edge_size_evaluator;

  MBTestOutputFunctor( MBInterface* imesh, MBInterface* omesh, MBEdgeSizeEvaluator* size_eval )
    {
    this->mesh = omesh;
    this->input_is_output = ( imesh == omesh );
    this->edge_size_evaluator = size_eval;

    this->split_vertices.resize( 4 );
    this->split_vertices[0] = 0; // Vertices (0-faces) cannot be split
    this->split_vertices[1] = new MBSplitVertices<1>( this->mesh );
    this->split_vertices[2] = new MBSplitVertices<2>( this->mesh );
    this->split_vertices[3] = new MBSplitVertices<3>( this->mesh );
    }

  ~MBTestOutputFunctor()
    {
    for ( int i = 0; i < 4; ++ i )
      delete this->split_vertices[i];
    }

  void print_vert_crud( MBEntityHandle vout, int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags )
    {
    std::cout << "+ {";
    for ( int i = 0; i < nvhash; ++ i )
      std::cout << " " << vhash[i];
    std::cout << " } -> " << vout << " ";

    std::cout << "[ " << vcoords[0];
    for ( int i = 1; i < 6; ++i )
      std::cout << ", " << vcoords[i];
    std::cout << " ] ";

    double* x = (double*)vtags;
    int* m = (int*)( (char*)vtags + 2 * sizeof( double ) );
    std::cout << "< " << x[0]
              << ", " << x[1];
    for ( int i = 0; i < 4; ++i )
      std::cout << ", " << m[i];
    std::cout << " >\n";
    }

  void assign_tags( MBEntityHandle vhandle, const void* vtags )
    {
    if ( ! vhandle )
      return; // Ignore bad vertices

    int num_tags = this->edge_size_evaluator->get_number_of_vertex_tags();
    MBTag tag_handle;
    int tag_offset;
    for ( int i = 0; i < num_tags; ++i )
      {
      this->edge_size_evaluator->get_output_vertex_tag( i, tag_handle, tag_offset );
      this->mesh->tag_set_data( tag_handle, &vhandle, 1, vtags );
      }
    }

  virtual MBEntityHandle operator () ( MBEntityHandle vhash, const double* vcoords, const void* vtags )
    {
    if ( this->input_is_output )
      { // Don't copy the original vertex!
      this->print_vert_crud( vhash, 1, &vhash, vcoords, vtags );
      return vhash;
      }
    MBEntityHandle vertex_handle;
    if ( this->mesh->create_vertex( vcoords + 3, vertex_handle ) != MB_SUCCESS )
      {
      std::cerr << "Could not insert mid-edge vertex!\n";
      }
    this->assign_tags( vertex_handle, vtags );
    this->print_vert_crud( vertex_handle, 1, &vhash, vcoords, vtags );
    return vertex_handle;
    }

  virtual MBEntityHandle operator () ( int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags )
    {
    MBEntityHandle vertex_handle;
    if ( nvhash == 1 )
      {
      vertex_handle = (*this)( *vhash, vcoords, vtags );
      }
    else if ( nvhash < 4 )
      {
      bool newly_created = this->split_vertices[nvhash]->find_or_create( vhash, vcoords, vertex_handle );
      if ( newly_created )
        {
        this->assign_tags( vertex_handle, vtags );
        }
      if ( ! vertex_handle )
        {
        std::cerr << "Could not insert mid-edge vertex!\n";
        }
      this->print_vert_crud( vertex_handle, nvhash, vhash, vcoords, vtags );
      }
    else
      {
      vertex_handle = 0;
      std::cerr << "Not handling splits on faces with " << nvhash << " corners yet.\n";
      }
    return vertex_handle;
    }

  virtual void operator () ( MBEntityHandle h )
    {
    std::cout << h << " ";
    this->elem_vert.push_back( h );
    }

  virtual void operator () ( MBEntityType etyp )
    {
    MBEntityHandle elem_handle;
    if ( this->mesh->create_element( etyp, &this->elem_vert[0], this->elem_vert.size(), elem_handle ) == MB_FAILURE )
      {
      std::cerr << " *** ";
      }
    this->elem_vert.clear();
    std::cout << "---------> " << elem_handle << " ( " << etyp << " )\n\n";
    }
};

int TestMeshRefiner( int argc, char* argv[] )
{
  int nprocs, rank;
#ifdef USE_MPI
  int err = MPI_Init( &argc, &argv );
  err = MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  std::cout << "Rank: " << ( rank + 1 ) << " of: " << nprocs << "\n";
#else // USE_MPI
  nprocs = 1;
  rank = 0;
#endif // USE_MPI

  bool input_is_output = ( argc > 1 && ! strcmp( argv[1], "-new-mesh" ) ) ? false : true;
  MBInterface* imesh = new MBCore;
  MBInterface* omesh = input_is_output ? imesh : new MBCore;
  MBEdgeSizeSimpleImplicit* eval = new MBEdgeSizeSimpleImplicit( imesh, omesh );

  double coords[][6] = {
    {  0. ,  0.0,  0. ,  0. ,  0.0,  0.  }, // 0
    { -1. ,  0.0,  0.5, -1. ,  0.0,  0.5 }, // 1
    { -0.5, -1.0,  0.5, -0.5, -1.0,  0.5 }, // 2
    {  0. ,  0. ,  1. ,  0. ,  0. ,  1.  }, // 3
    {  0.5,  0.5,  0.5,  0.5,  0.5,  0.5 }, // 4
    {  0.5, -0.5,  0.5,  0.5, -0.5,  0.5 }  // 5
  };

  double default_floatular[2] = {
     38.7,   104. };
  double floatular_values[][2] = {
    {  38.7,   104.,    }, // 0
    {   3.141,   2.718, }, // 1
    { 123.456, 789.012, }, // 2
    {   0.,      1.,    }, // 3
    {   1.,      0.,    }, // 4
    {  -1.,      1.,    }  // 5
  };

  int default_intular[4] = {
     7, 11, 24,  7 };
  int intular_values[][4] = {
    {  7, 11, 24,  7, }, // 0
    { 10,  4, 10, 20, }, // 1
    {  1,  2,  3,  4, }, // 2
    { 13, 17, 19, 23, }, // 3
    {  3,  3,  0,  6, }, // 4
    {  5,  4,  3,  2  }  // 5
  };

  int default_gid[] = { -1 };
  int gid_values[] = {
    0, 1, 2, 3, 4, 5, // vertices
    6, 7, 8, 9        // tetrahedra
  };

  // List of vertices resident on each node.
  int proc_nodes[4][4] = {
    { 0, 2, 3, 1 },
    { 0, 1, 3, 4 },
    { 0, 4, 3, 5 },
    { 0, 5, 3, 2 },
  };

  // List of nodes with a copy of each vertex (first entry is owner)
  int node_procs[6][4] = {
    { 0,  1,  2,  3 }, // 0
    { 0,  1, -1, -1 }, // 1
    { 0,  3, -1, -1 }, // 2
    { 0,  1,  2,  3 }, // 3
    { 1,  2, -1, -1 }, // 4
    { 2,  3, -1, -1 }  // 5
  };

  eval->set_ratio( 2. );

  MBEntityHandle node_handles[4];
  MBEntityHandle tet_handle;
  MBEntityHandle tri_handle;

  MBTag tag_floatular;
  imesh->tag_create( "floatular", 2 * sizeof( double ), MB_TAG_DENSE, MB_TYPE_DOUBLE, tag_floatular, default_floatular );

  MBTag tag_intular;
  imesh->tag_create( "intular", 4 * sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_intular, default_intular );

  MBTag tag_gid;
  imesh->tag_create( PARALLEL_GID_TAG_NAME, sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_gid, default_gid );

  MBTag tag_spr;
  imesh->tag_create( PARALLEL_SHARED_PROCS_TAG_NAME, 4 * sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_spr, default_gid );

  void const* iptrs[4];
  void const* fptrs[4];
  void const* gptrs[4];
  void const* sptrs[4];
  for ( int i = 0; i < 4; ++ i )
    {
    imesh->create_vertex( coords[proc_nodes[rank][i]], node_handles[i] );
    iptrs[i] = (void const*) intular_values[proc_nodes[rank][i]];
    fptrs[i] = (void const*) floatular_values[proc_nodes[rank][i]];
    gptrs[i] = (void const*) &gid_values[proc_nodes[rank][i]];
    sptrs[i] = (void const*) &node_procs[proc_nodes[rank][i]];
    }

  imesh->tag_set_data( tag_floatular, node_handles, 4, fptrs, 0 );
  eval->add_vertex_tag( tag_floatular );

  imesh->tag_set_data( tag_intular, node_handles, 4, iptrs, 0 );
  eval->add_vertex_tag( tag_intular );

  imesh->tag_set_data( tag_gid, node_handles, 4, gptrs, 0 );
  imesh->tag_set_data( tag_spr, node_handles, 4, sptrs, 0 );

  imesh->create_element( MBTET, node_handles, 4, tet_handle );
  imesh->tag_set_data( tag_gid, &tet_handle, 1, gid_values + 6 + rank );
  imesh->list_entities( 0, 1 );

  MBSimplexTemplateRefiner eref( imesh );
  MBTestOutputFunctor* ofunc = new MBTestOutputFunctor( imesh, omesh, eval );
  eref.set_edge_size_evaluator( eval );
  eref.set_output_functor( ofunc );
  eval->create_output_tags();
  eref.refine_entity( tet_handle );

  ofunc->mesh->list_entities( 0, 1 );

  if ( ! ofunc->input_is_output )
    delete ofunc->mesh;
  delete imesh;

#ifdef USE_MPI
  err = MPI_Barrier( MPI_COMM_WORLD );
  err = MPI_Finalize();
#endif // USE_MPI
  return 0;
}

int main( int argc, char* argv[] )
{
  return TestMeshRefiner( argc, argv );
}
