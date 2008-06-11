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

class MBTestOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  typedef std::vector<MBEntityHandle> node_list_t;
  typedef std::pair<MBEntityHandle,MBEntityHandle> node_pair_t;
  typedef std::map<MBEntityHandle,MBEntityHandle> node_hash_t;
  typedef std::map<node_pair_t,MBEntityHandle> edge_hash_t;

  MBInterface* mesh;
  bool input_is_output;
  node_hash_t node_hash;
  edge_hash_t edge_hash;
  node_list_t elem_vert;

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
    else if ( nvhash == 2 )
      {
      node_pair_t pr;
      if ( vhash[0] < vhash[1] )
        {
        pr.first = vhash[0];
        pr.second = vhash[1];
        }
      else
        {
        pr.first = vhash[1];
        pr.second = vhash[0];
        }
      edge_hash_t::iterator it = this->edge_hash.find( pr );
      if ( it == this->edge_hash.end() )
        {
        if ( this->mesh->create_vertex( vcoords + 3, vertex_handle ) != MB_SUCCESS )
          {
          std::cerr << "Could not insert mid-edge vertex!\n";
          }
        this->edge_hash[pr] = vertex_handle;
        }
      else
        {
        vertex_handle = it->second;
        }
      this->print_vert_crud( vertex_handle, 2, vhash, vcoords, vtags );
      }
    else
      {
      vertex_handle = -1;
      std::cerr << "Not handling mid-face vertices yet.\n";
      // FIXME: Handle face midpoint.
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

  MBInterface* iface = new MBCore;
  MBEdgeSizeSimpleImplicit* eval = new MBEdgeSizeSimpleImplicit( iface );

  double p0[6] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
  double p1[6] = { 0.5, 0.0, 0.0,  0.5, 0.0, 0.0 };
  double p2[6] = { 1.0, 0.0, 0.0,  1.0, 0.0, 0.0 };
  double p3[6] = { 0.6, 2.0, 0.0,  0.6, 2.0, 0.0 };

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
  iface->tag_create( "floatular", 2 * sizeof( double ), MB_TAG_DENSE, MB_TYPE_DOUBLE, tag_floatular, default_floatular );

  MBTag tag_intular;
  iface->tag_create( "intular", 4 * sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_intular, default_intular );

  MBTag tag_gid;
  iface->tag_create( PARALLEL_GID_TAG_NAME, sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_gid, default_gid );

  void const* iptrs[4];
  void const* fptrs[4];
  void const* gptrs[4];
  for ( int i = 0; i < 4; ++ i )
    {
    iface->create_vertex( coords[proc_nodes[rank][i]], node_handles[i] );
    iptrs[i] = (void const*) intular_values[proc_nodes[rank][i]];
    fptrs[i] = (void const*) floatular_values[proc_nodes[rank][i]];
    gptrs[i] = (void const*) &gid_values[proc_nodes[rank][i]];
    }

  iface->tag_set_data( tag_floatular, node_handles, 4, fptrs, 0 );
  eval->add_vertex_tag( tag_floatular );

  iface->tag_set_data( tag_intular, node_handles, 4, iptrs, 0 );
  eval->add_vertex_tag( tag_intular );

  iface->tag_set_data( tag_gid, node_handles, 4, gptrs, 0 );

  iface->create_element( MBTET, node_handles, 4, tet_handle );
  iface->tag_set_data( tag_gid, &tet_handle, 1, gid_values + 6 + rank );
  iface->list_entities( 0, 1 );

  MBSimplexTemplateRefiner eref( iface );
  MBTestOutputFunctor* ofunc = new MBTestOutputFunctor;
  ofunc->input_is_output = ( argc > 1 && ! strcmp( argv[1], "-new-mesh" ) ) ? false : true;
  ofunc->mesh = ofunc->input_is_output ? iface : new MBCore;
  eref.set_edge_size_evaluator( eval );
  eref.set_output_functor( ofunc );
  eref.refine_entity( tet_handle );

  ofunc->mesh->list_entities( 0, 1 );

  if ( ! ofunc->input_is_output )
    delete ofunc->mesh;
  delete iface;

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
