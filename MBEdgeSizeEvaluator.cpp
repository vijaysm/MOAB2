#include "MBEdgeSizeEvaluator.hpp"

#include "MBInterface.hpp"

#include <assert.h>

MBEdgeSizeEvaluator::MBEdgeSizeEvaluator( MBInterface* parentMesh )
{
  assert( parentMesh );

  this->mesh = parentMesh;
  this->reset_vertex_tags();
}

MBEdgeSizeEvaluator::~MBEdgeSizeEvaluator()
{
}

void MBEdgeSizeEvaluator::reset_vertex_tags()
{
  this->vertexSize = 0;
  this->vertexTags.clear();
}

int MBEdgeSizeEvaluator::add_vertex_tag( MBTag tag_handle )
{
  int offset = this->vertexSize; // old size is offset of tag being added
  int tagSize;
  MBTagType tagType;
  if ( this->mesh->tag_get_size( tag_handle, tagSize ) != MB_SUCCESS )
    return -1;

  if ( this->mesh->tag_get_type( tag_handle, tagType ) != MB_SUCCESS )
    return -1;

  if ( tagType == MB_TAG_BIT )
    {
    // Pad any bit tags to a size in full bytes.
    tagSize = ( tagSize % 8 ? 1 : 0 ) + ( tagSize / 8 );
    }

  // Now pad so that the next tag will be word-aligned:
  while ( tagSize % sizeof(int) )
    ++tagSize;

  this->vertexSize += tagSize;

  this->vertexTags.push_back( std::pair< MBTag, int >( tag_handle, offset ) );
  return offset;
}
