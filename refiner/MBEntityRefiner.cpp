#include "MBEntityRefiner.hpp"

#include "MBInterface.hpp"

MBEntityRefiner::MBEntityRefiner( MBInterface* parentMesh )
{  
  this->mesh = parentMesh;
  this->edge_size_evaluator = 0;
}

MBEntityRefiner::~MBEntityRefiner()
{
}

bool MBEntityRefiner::set_edge_size_evaluator( MBEdgeSizeEvaluator* ese )
{
  if ( ! ese ) return false;

  this->edge_size_evaluator = ese;

  return true;
}


