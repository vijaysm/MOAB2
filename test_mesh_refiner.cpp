#include "MBCore.hpp"
#include "MBEdgeSizeSimpleImplicit.hpp"
#include "MBInterface.hpp"

int TestMeshRefiner( int argc, char* argv[] )
{
  MBInterface* iface = new MBCore;
  MBEdgeSizeSimpleImplicit eval( iface );
  double p0[6] = { 0.1, 0.0, 0.0,  0.1, 0.0, 0.0 };
  double p1[6] = { 0.6, 0.0, 0.0,  0.6, 0.0, 0.0 };
  double p2[6] = { 1.1, 0.0, 0.0,  1.1, 0.0, 0.0 };
  if ( eval.evaluate_edge( p0, 0, p1, 0, p2, 0 ) )
    {
    return 1;
    }

  eval.set_ratio( 2. );
  if ( ! eval.evaluate_edge( p0, 0, p1, 0, p2, 0 ) )
    {
    return 1;
    }

  delete iface;
  return 0;
}

int main( int argc, char* argv[] )
{
  return TestMeshRefiner( argc, argv );
}
