#include "MBMeshRefiner.hpp"

#include "MBInterface.hpp"

MBMeshRefiner::MBMeshRefiner( MBInterface* parentMesh )
{  
  this->mesh = parentMesh;
  this->entity_refiner = 0;
}

MBMeshRefiner::~MBMeshRefiner()
{
}

bool MBMeshRefiner::set_entity_refiner( MBEntityRefiner* er )
{
  if ( ! er ) return false;

  this->entity_refiner = er;

  return true;
}

bool MBMeshRefiner::refine_mesh()
{
  return false;
}

