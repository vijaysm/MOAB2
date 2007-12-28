#include "MBCore.hpp"
#include "MBEdgeSizeSimpleImplicit.hpp"
#include "MBSimplexTemplateRefiner.hpp"
#include "MBInterface.hpp"

#include <iostream>

class MBTestOutputFunctor : public MBEntityRefinerOutputFunctor
{
  virtual void operator () ( const double* vcoords, const void* vtags )
    {
    std::cout << "[ " << vcoords[0];
    for ( int i = 1; i < 6; ++i )
      std::cout << ", " << vcoords[i];
    std::cout << " ]\n";
    }

  virtual void operator () ( MBEntityType etyp )
    {
    std::cout << "---------- " << etyp << "\n\n";
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

  iface->create_vertex( p0, node_handles[0] );
  iface->create_vertex( p2, node_handles[1] );
  iface->create_vertex( p3, node_handles[2] );

  iface->create_element( MBTRI, node_handles, 3, tri_handle );

  MBSimplexTemplateRefiner eref( iface );
  MBTestOutputFunctor* ofunc = new MBTestOutputFunctor;
  eref.set_edge_size_evaluator( eval );
  eref.set_output_functor( ofunc );
  eref.refine_entity( tri_handle );

  delete iface;
  return 0;
}

int main( int argc, char* argv[] )
{
  return TestMeshRefiner( argc, argv );
}
