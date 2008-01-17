#include "MBEdgeSizeEvaluator.hpp"

#include "MBInterface.hpp"

#include <assert.h>

/// Construct an evaluator.
MBEdgeSizeEvaluator::MBEdgeSizeEvaluator( MBInterface* parentMesh )
{
  assert( parentMesh );

  this->mesh = parentMesh;
  this->reset_vertex_tags();
}

/// Destruction is virtual so subclasses may clean up after refinement.
MBEdgeSizeEvaluator::~MBEdgeSizeEvaluator()
{
}

/**\fn bool MBEdgeSizeEvaluator::evaluate_edge( \
  *         const double* p0, const void* t0, double* p1, void* t1, const double* p2, const void* t2 )
  *\brief Returns true if the edge \a p0 - \a p2 should be subdivided, false otherwise.
  *
  * The arguments \a p0, \a p1, and \a p2 are all pointers to arrays of 6 doubles each
  * while the arguments \a t0, \a t1, and \a t2 are all pointers to arrays of tag data
  * defined at the corresponding point. While the endpoints \a p0 and \a p2 are
  * immutable, the mid-edge point coordinates \a p1 and tag data \a t1 may be altered by
  * evaluate_edge(). Altered values will be ignored if evaluate_edge() returns false.
  * Be careful to ensure that all calls to evaluate_edge() perform identical modifications
  * given identical input values!
  *
  * A list of tags passed in \a t0, \a t1, and \a t2 is stored in the vertex_tags member.
  * The vertex_size member stores the total length of data associated with each pointer (in bytes).
  * Subclasses may access vertex_tags and vertexSize directly; the refiner uses public methods to
  * populate vertex_tags before evaluate_edge() is called.
  */

/// Clear the list of tag values that will appear past the vertex coordinates in \a p0, \a p1, and \a p2.
void MBEdgeSizeEvaluator::reset_vertex_tags()
{
  this->vertex_size = 0;
  this->vertex_tags.clear();
}

/** Add a tag to the list of tag values that will appear past the vertex coordinates.
  * The return value is the offset into each vertex coordinate pointer (\a p0, \a p1, \a p2) where the
  * tag value(s) will be stored.
  */
int MBEdgeSizeEvaluator::add_vertex_tag( MBTag tag_handle )
{
  int offset = this->vertex_size; // old size is offset of tag being added
  int tag_size;
  MBTagType tagType;
  if ( this->mesh->tag_get_size( tag_handle, tag_size ) != MB_SUCCESS )
    return -1;

  if ( this->mesh->tag_get_type( tag_handle, tagType ) != MB_SUCCESS )
    return -1;

  if ( tagType == MB_TAG_BIT )
    {
    // Pad any bit tags to a size in full bytes.
    tag_size = ( tag_size % 8 ? 1 : 0 ) + ( tag_size / 8 );
    }

  // Now pad so that the next tag will be word-aligned:
  while ( tag_size % sizeof(int) )
    ++tag_size;

  this->vertex_size += tag_size;

  this->vertex_tags.push_back( std::pair< MBTag, int >( tag_handle, offset ) );
  return offset;
}

/**\fn int MBEdgeSizeEvaluator::get_vertex_tag_size()
  *\brief Return the number of bytes to allocate for tag data per point.
  */

/**\fn int MBEdgeSizeEvaluator::get_number_of_vertex_tags() const
  *\brief Return the number of tags that will be output with each new vertex.
  */

/**\brief Return the tag handle and its offset in the array of tag data of each vertex.
  *
  * @param[in] i An index into the list of tags for the vertex.
  * @param[out] tag The tag handle for the $i$-th vertex tag.
  * @param[out] byte_offset The offset (in bytes) of the start of this tag's data in a vertex tag record.
  */
void MBEdgeSizeEvaluator::get_vertex_tag( int i, MBTag& tag, int& byte_offset )
{
  std::vector< std::pair< MBTag, int > >::iterator it = this->vertex_tags.begin() + i;
  tag = it->first;
  byte_offset = it->second;
}
