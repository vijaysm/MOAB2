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
  * A list of tags passed in \a t0, \a t1, and \a t2 is stored in the vertexTags member.
  * The vertexSize member stores the total length of data associated with each pointer (in bytes).
  * Subclasses may access vertexTags and vertexSize directly; the refiner uses public methods to
  * populate vertexTags before evaluate_edge() is called.
  */

/// Clear the list of tag values that will appear past the vertex coordinates in \a p0, \a p1, and \a p2.
void MBEdgeSizeEvaluator::reset_vertex_tags()
{
  this->vertexSize = 0;
  this->vertexTags.clear();
}

/** Add a tag to the list of tag values that will appear past the vertex coordinates.
  * The return value is the offset into each vertex coordinate pointer (\a p0, \a p1, \a p2) where the
  * tag value(s) will be stored.
  */
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

/**\fn int MBEdgeSizeEvaluator::get_vertex_tag_size()
  *\brief Return the number of bytes to allocate for tag data per point.
  */

/**\brief Given endpoint coordinates and tag values plus midpoint coordinates, compute midpoint tag values.
  *
  * Normally, this function will be invoked by the MBEntityRefiner before evaluate_edge is called.
  * However, if evaluate_edge() changes the parametric coordinates of the midpoint,
  * it should call evaluate_tags_at_midpoint() again to update any tag values;
  * that is why this function is a member of MBEdgeSizeEvaluator and not MBEntityRefiner.
  *
  * @param[in] c0 Pointer to endpoint 0 coordinates. The parametric coordinates (3) are followed by world coordinates (3).
  * @param[in] t0 Pointer to endpoint 0 tag values.
  * @param[in] cm Pointer to midpoint coordinates. The parametric coordinates (3) are followed by world coordinates (3).
  * @param[out] tm Pointer to midpoint tag values.
  * @param[in] c1 Pointer to endpoint 1 coordinates. The parametric coordinates (3) are followed by world coordinates (3).
  * @param[in] t1 Pointer to endpoint 1 tag values.
  */
void MBEdgeSizeEvaluator::evaluate_tags_at_midpoint( const double* c0, const void* t0, 
                                                     const double* cm, void* tm, 
                                                     const double* c1, const void* t1 ) const
{
}

