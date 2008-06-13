#include "MBRefinerTagManager.hpp"

#include "MBInterface.hpp"

#include <iostream>
#include <assert.h>

/// Construct an evaluator.
MBRefinerTagManager::MBRefinerTagManager( MBInterface* in_mesh, MBInterface* out_mesh )
{
  assert( in_mesh );
  if ( ! out_mesh )
    out_mesh = in_mesh;

  this->input_mesh = in_mesh;
  this->output_mesh = out_mesh;
  this->reset_vertex_tags();
}

/// Destruction is virtual so subclasses may clean up after refinement.
MBRefinerTagManager::~MBRefinerTagManager()
{
}

/**\fn bool MBRefinerTagManager::evaluate_edge( \
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
  * A list of tags passed in \a t0, \a t1, and \a t2 is stored in the \a input_vertex_tags member.
  * (for tag handles defined on the input mesh) and the \a output_vertex_tags (for the same tag handles
  * defined on the output mesh).
  * The vertex_size member stores the total length of data associated with each pointer (in bytes).
  * Subclasses may access input_vertex_tags, output_vertex_tags, and vertex_size directly;
  * the refiner uses public methods to set them before any entities are evaluated for subdivision.
  * The output_vertex_tags member is populated when the refiner calls create_output_tags().
  * When the input mesh and output mesh pointers are identical, this simply copies input_vertex_tags
  * to output_vertex_tags.
  * When the pointers are distinct, tags are created on the output mesh.
  */

/// Clear the list of tag values that will appear past the vertex coordinates in \a p0, \a p1, and \a p2.
void MBRefinerTagManager::reset_vertex_tags()
{
  this->vertex_size = 0;
  this->input_vertex_tags.clear();
  this->output_vertex_tags.clear();
}

/** Add a tag to the list of tag values that will appear past the vertex coordinates.
  * The return value is the offset into each vertex coordinate pointer (\a p0, \a p1, \a p2) where the
  * tag value(s) will be stored.
  */
int MBRefinerTagManager::add_vertex_tag( MBTag tag_handle )
{
  int offset = this->vertex_size; // old size is offset of tag being added
  int tag_size;
  MBTagType tagType;
  if ( this->input_mesh->tag_get_size( tag_handle, tag_size ) != MB_SUCCESS )
    return -1;

  if ( this->input_mesh->tag_get_type( tag_handle, tagType ) != MB_SUCCESS )
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

  this->input_vertex_tags.push_back( std::pair< MBTag, int >( tag_handle, offset ) );
  return offset;
}

/**\fn int MBRefinerTagManager::get_vertex_tag_size()
  *\brief Return the number of bytes to allocate for tag data per point.
  */

/**\fn int MBRefinerTagManager::get_number_of_vertex_tags() const
  *\brief Return the number of tags that will be output with each new vertex.
  */

/**\brief Populate the list of output tags to match the list of input tags.
  *
  * When the input mesh and output mesh pointers are identical, this simply copies the list of input tags.
  * When the two meshes are distinct, the corresponding tags are created on the output mesh.
  */
void MBRefinerTagManager::create_output_tags()
{
  if ( this->input_mesh == this->output_mesh )
    {
    this->output_vertex_tags = this->input_vertex_tags;
    return;
    }

  std::pair< MBTag, int > tag_rec;
  std::vector< std::pair< MBTag, int > >::iterator it;
  std::vector< char > tag_default;
  std::string tag_name;
  MBTagType tag_type;
  MBDataType tag_data_type;
  int tag_size;
  for ( it = this->input_vertex_tags.begin(); it != this->output_vertex_tags.end(); ++ it )
    {
    MBTag tag_in = it->first;
    tag_rec.second = it->second;
    this->input_mesh->tag_get_name( tag_in, tag_name );
    this->input_mesh->tag_get_size( tag_in, tag_size );
    this->input_mesh->tag_get_type( tag_in, tag_type );
    this->input_mesh->tag_get_data_type( tag_in, tag_data_type );
    this->input_mesh->tag_get_default_value( tag_in, (void*) &tag_default[0] );
    tag_default.resize( tag_size );
    MBErrorCode res = this->output_mesh->tag_create(
      tag_name.c_str(), tag_size, tag_type, tag_data_type, tag_rec.first, (void*) &tag_default[0], true );
    if ( res == MB_FAILURE )
      {
      std::cerr
        << "Could not create output tag name: \"" << tag_name.c_str() << "\" type: "
        << tag_type << " data type: " << tag_data_type << "\n";
      }
    else
      {
      this->output_vertex_tags.push_back( tag_rec );
      }
    }
}

/**\brief Return the tag handle and its offset in the array of tag data of each vertex.
  *
  * @param[in] i An index into the list of tags for the vertex.
  * @param[out] tag The tag handle on the input mesh for the $i$-th vertex tag.
  * @param[out] byte_offset The offset (in bytes) of the start of this tag's data in a vertex tag record.
  */
void MBRefinerTagManager::get_input_vertex_tag( int i, MBTag& tag, int& byte_offset )
{
  std::vector< std::pair< MBTag, int > >::iterator it = this->input_vertex_tags.begin() + i;
  tag = it->first;
  byte_offset = it->second;
}

/**\brief Return the tag handle and its offset in the array of tag data of each vertex.
  *
  * @param[in] i An index into the list of tags for the vertex.
  * @param[out] tag The tag handle on the output mesh for the $i$-th vertex tag.
  * @param[out] byte_offset The offset (in bytes) of the start of this tag's data in a vertex tag record.
  */
void MBRefinerTagManager::get_output_vertex_tag( int i, MBTag& tag, int& byte_offset )
{
  std::vector< std::pair< MBTag, int > >::iterator it = this->output_vertex_tags.begin() + i;
  tag = it->first;
  byte_offset = it->second;
}
