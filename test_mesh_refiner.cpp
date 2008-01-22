#include "MBCore.hpp"
#include "MBEdgeSizeSimpleImplicit.hpp"
#include "MBSimplexTemplateRefiner.hpp"
#include "MBInterface.hpp"

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

  virtual void operator () ( const double* vcoords, const void* vtags, MBEntityHandle* vhash )
    {
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
    std::cout << " > { ";
    MBEntityHandle* vin = vhash;
    while ( *vin )
      {
      std::cout << *vin << " ";
      ++ vin;
      }
    std::cout << "}\n";

    MBEntityHandle vertex_handle;
    if ( ! vhash[1] )
      {
      if ( this->input_is_output )
        { // Don't copy the original vertex!
        vertex_handle = vhash[0];
        }
      else
        {
        node_hash_t::iterator it = this->node_hash.find( vhash[0] );
        if ( it == this->node_hash.end() )
          {
          if ( this->mesh->create_vertex( vcoords + 3, vertex_handle ) != MB_SUCCESS )
            {
            std::cerr << "Could not insert corner vertex!\n";
            }
          this->node_hash[vhash[0]] = vertex_handle;
          }
        else
          {
          vertex_handle = it->second;
          }
        }
      }
    else if ( ! vhash[2] )
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
      }
    else
      {
      std::cerr << "Not handling mid-face vertices yet.\n";
      // FIXME: Handle face midpoint.
      }
    std::cout << "        -> " << vertex_handle << "\n";
    this->elem_vert.push_back( vertex_handle );
    }

  virtual void operator () ( MBEntityType etyp )
    {
    std::cout << "---------- " << etyp << "\n\n";
    MBEntityHandle elem_handle;
    this->mesh->create_element( etyp, &this->elem_vert[0], this->elem_vert.size(), elem_handle );
    this->elem_vert.clear();
    }
};

int TestMeshRefiner( int argc, char* argv[] )
{
  MBInterface* iface = new MBCore;
  MBEdgeSizeSimpleImplicit* eval = new MBEdgeSizeSimpleImplicit( iface );

  double p0[6] = { 0.0, 0.0, 0.0,  0.0, 0.0, 0.0 };
  double p1[6] = { 0.5, 0.0, 0.0,  0.5, 0.0, 0.0 };
  double p2[6] = { 1.0, 0.0, 0.0,  1.0, 0.0, 0.0 };
  double p3[6] = { 0.6, 2.0, 0.0,  0.6, 2.0, 0.0 };

  double default_floatular[] = { 38.7, 104. };
  double floatular_values[] = { 38.7, 104., 
                                123.456, 789.012, 
                                0., 1. };

  int default_intular[] = { 7, 11, 24, 7 };
  int intular_values[] = {  7, 11, 24, 7, 
                            1, 2, 3, 4, 
                            13, 17, 19, 23 };

  if ( eval->evaluate_edge( p0, 0, p1, 0, p2, 0 ) )
    {
    return 1;
    }

  eval->set_ratio( 2. );
  if ( ! eval->evaluate_edge( p0, 0, p1, 0, p2, 0 ) )
    {
    return 1;
    }

  MBEntityHandle node_handles[3];
  MBEntityHandle tri_handle;

  MBTag tag_floatular;
  iface->tag_create( "floatular", 2 * sizeof( double ), MB_TAG_DENSE, MB_TYPE_DOUBLE, tag_floatular, default_floatular );

  MBTag tag_intular;
  iface->tag_create( "intular", 4 * sizeof( int ), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_intular, default_intular );

  iface->create_vertex( p0, node_handles[0] );
  iface->create_vertex( p2, node_handles[1] );
  iface->create_vertex( p3, node_handles[2] );

  iface->tag_set_data( tag_floatular, node_handles, 3, (void*)floatular_values );
  eval->add_vertex_tag( tag_floatular );

  iface->tag_set_data( tag_intular, node_handles, 3, (void*)intular_values );
  eval->add_vertex_tag( tag_intular );

  iface->create_element( MBTRI, node_handles, 3, tri_handle );

  MBSimplexTemplateRefiner eref( iface );
  MBTestOutputFunctor* ofunc = new MBTestOutputFunctor;
  ofunc->input_is_output = ( argc > 1 && ! strcmp( argv[1], "-new-mesh" ) ) ? false : true;
  ofunc->mesh = ofunc->input_is_output ? iface : new MBCore;
  eref.set_edge_size_evaluator( eval );
  eref.set_output_functor( ofunc );
  eref.refine_entity( tri_handle );

  if ( ! ofunc->input_is_output )
    delete ofunc->mesh;
  delete iface;
  return 0;
}

int main( int argc, char* argv[] )
{
  return TestMeshRefiner( argc, argv );
}
