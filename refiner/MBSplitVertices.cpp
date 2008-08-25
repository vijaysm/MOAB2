#include "MBSplitVertices.hpp"
#include "MBRefinerTagManager.hpp"

#include "MBParallelConventions.h"

MBSplitVerticesBase::MBSplitVerticesBase( MBRefinerTagManager* tag_mgr )
{
  this->tag_manager = tag_mgr;
  this->mesh_out = tag_mgr->get_output_mesh();
}

MBSplitVerticesBase::~MBSplitVerticesBase()
{
}

