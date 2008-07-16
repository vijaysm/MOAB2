#ifndef MB_SPLITVERTICES_HPP
#define MB_SPLITVERTICES_HPP

#include "MBTypes.h"
#include "MBProcessSet.hpp"
#include "MBTagConventions.hpp"

#include <map>
#include <vector>

class MBRefinerTagManager;

/** An array of existing vertex handles used as a key in a dictionary of new vertices.
  */
template< int _n >
class MBSplitVertexIndex
{
public:
  MBSplitVertexIndex() { }
  MBSplitVertexIndex( const MBEntityHandle* src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src[i]; std::sort( this->handles, this->handles + _n ); }
  MBSplitVertexIndex( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; this->process_set = src.process_set; }
  MBSplitVertexIndex& operator = ( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; this->process_set = src.process_set; return *this; }

  void set_common_processes( const MBProcessSet& procs )
    { this->process_set = procs; }
  MBProcessSet& common_processes()
    { return this->process_set; }
  const MBProcessSet& common_processes() const
    { return this->process_set; }

  bool operator < ( const MBSplitVertexIndex<_n>& other ) const
    {
    // Ignore the process set. Only program errors lead to mismatched process sets with identical handles.
    for ( int i = 0; i < _n; ++ i )
      if ( this->handles[i] < other.handles[i] )
        return true;
      else if ( this->handles[i] > other.handles[i] )
        return false;
    return false;
    }

  MBEntityHandle handles[_n];
  MBProcessSet process_set;
};

/** A non-templated base class that the template subclasses all share.
  *
  * All methods that need to be accessed by other classes should be
  * declared by the base class so that no knowledge of template parameters
  * is required.
  */
class MBSplitVerticesBase
{
public:
  MBSplitVerticesBase( MBRefinerTagManager* tag_mgr );
  virtual ~MBSplitVerticesBase();

  virtual bool find_or_create(
    const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
    std::map<MBProcessSet,int>& proc_partition_counts ) = 0;

  virtual void assign_global_ids( std::map<MBProcessSet,MBEntityHandle>& gids ) = 0;

  /// Determine which processes will contain an output vertex given the split vertices defining it.
  void update_partition_counts( int num, const MBEntityHandle* split_src, std::map<MBProcessSet,int>& proc_partition_counts );

  /// Prepare to compute the processes on which a new split-vertex will live.
  void begin_vertex_procs();

  /// Call this for each existing corner vertex used to define a split-vertex.
  void add_vertex_procs( MBEntityHandle vert_in );

  /// Call this once after all the add_vertex_procs() calls for a split-vertex to prepare queues for the second stage MPI send. 
  void end_vertex_procs();

  MBInterface* mesh_in; // Input mesh. Needed to determine tag values on split_src verts
  MBInterface* mesh_out; // Output mesh. Needed for new vertex set in vert_handle
  MBRefinerTagManager* tag_manager;
  std::vector<int> shared_procs_val; // Used to hold procs sharing an input vert.
  MBProcessSet current_shared_procs; // Holds process list as it is being accumulated
  MBProcessSet common_shared_procs; // Holds intersection of several shared_procs_vals.
  int rank; // This process' rank.
  bool first_vertex; // True just after begin_vertex_procs() is called.
};

/** A map from a set of pre-existing entities to a new mesh entity.
  *
  * This is used as a dictionary to determine whether a new vertex should be
  * created on the given n-simplex (n being the template parameter) or whether
  * it has already been created as part of the refinement of a neighboring entity.
  */
template< int _n >
class MBSplitVertices : public std::map<MBSplitVertexIndex<_n>,MBEntityHandle>, public MBSplitVerticesBase
{
public:
  typedef std::map<MBSplitVertexIndex<_n>,MBEntityHandle> MapType;
  typedef typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator MapIteratorType;

  MBSplitVertices( MBRefinerTagManager* tag_mgr );
  virtual ~MBSplitVertices();
  virtual bool find_or_create(
    const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
    std::map<MBProcessSet,int>& proc_partition_counts );
  virtual void assign_global_ids( std::map<MBProcessSet,MBEntityHandle>& gids );
};

// ------------------------- Template member definitions ----------------------
template< int _n >
MBSplitVertices<_n>::MBSplitVertices( MBRefinerTagManager* tag_mgr )
  : MBSplitVerticesBase( tag_mgr )
{
  this->shared_procs_val.resize( _n * MAX_SHARING_PROCS );
}

template< int _n >
MBSplitVertices<_n>::~MBSplitVertices()
{
}

template< int _n >
bool MBSplitVertices<_n>::find_or_create(
  const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
  std::map<MBProcessSet,int>& proc_partition_counts )
{
  MBSplitVertexIndex<_n> key( split_src );
  MapIteratorType it = this->find( key );
  if ( it == this->end() )
    {
    this->update_partition_counts( _n, split_src, proc_partition_counts );
    key.set_common_processes( this->common_shared_procs );
    if ( this->mesh_out->create_vertex( coords, vert_handle ) != MB_SUCCESS )
      {
      return false;
      }
    (*this)[key] = vert_handle;
    return true;
    }
  vert_handle = it->second;
  return false;
}

template< int _n >
void MBSplitVertices<_n>::assign_global_ids( std::map<MBProcessSet,MBEntityHandle>& gids )
{
  MBTag tag_gid;
  int zero = 0;
  MBErrorCode result = this->mesh_out->tag_create(
    GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag_gid, &zero, true );
  if ( result != MB_SUCCESS && result != MB_ALREADY_ALLOCATED )
    return;

  typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator it;
  for ( it = this->begin(); it != this->end(); ++ it )
    {
    int gid = gids[it->first.process_set] ++;
    this->mesh_out->tag_set_data( tag_gid, &it->second, 1, &gid );
    }
}

#endif /* MB_SPLITVERTICES_HPP */
