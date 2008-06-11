#include "MBMeshRefiner.hpp"

#include "MBEdgeSizeEvaluator.hpp"
#include "MBEntityRefiner.hpp"
#include "MBInterface.hpp"

#ifdef USE_MPI
#include <mpi.h>
#else // USE_MPI
typedef int MPI_Comm;
#endif // USE_MPI

#include <map>

class MBRefinerPartition;

struct MBEdgeIndex
{
  MBEntityHandle corners[2];
};

struct MBFaceIndex
{
  MBEntityHandle corners[3];
};

/// A record of data to be sent to and received from a remote process in order to merge a refined mesh.
class MBRefinerPartition
{
public:
  /// Rank of the local process
  int local_process;
  /// Rank of the remote process to which this data belongs.
  int remote_process;
  /// A list of locally-owned vertices that are shared with the remote process associated with this object
  std::vector<MBEntityHandle> vertex_list_local;
  /// A list of remotely-owned vertices that are shared with the remote process associated with this object
  std::vector<MBEntityHandle> vertex_list_remote;
  /// A list of locally-owned global vertex IDs that are shared with the remote process associated with this object
  std::vector<int> global_ids_local;
  /// A list of remotely-owned global vertex IDs that are shared with the remote process associated with this object
  std::vector<int> global_ids_remote;
  /// Queue asynchronous sends and receives to merge data. You are responsible for calling MPI_Barrier() afterwards.
  void queue_merge( MPI_Comm communicator )
    {
#ifdef USE_MPI
    MPI_Request rq;
    int sz;
    sz = this->vertex_list_local.size();
    if ( sz )
      {
      MPI_Isend( (void*)&global_ids_local[0], sz, MB_MPI_ENTITY_HANDLE_TYPE,
        this->remote, this->local_rank, communicator, 0 /*request*/ );
      }
    sz = this->vertex_list_remote.size();
    if ( sz )
      {
      MPI_Irecv( (void*)&global_ids_remote[0], sz, MB_MPI_ENTITY_HANDLE_TYPE,
        this->remote, this->local_rank, communicator, 0 /*request*/ );
      }
#endif // USE_MPI
    }
  /// Perform the merge using the data. You must have called queue_merge() followed by MPI_Barrier() previous to this function.
  void perform_merge()
    {
    // Now the remote global IDs have been properly transferred from the remote process... assign them to the vertices

    }
};

class MBMeshRefinerOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  /// Mesh to hold output (may/may not be same as input mesh)
  MBInterface* mesh;
  /// Edge size evaluator is required because it has a list of output tags
  MBEdgeSizeEvaluator* edge_size_evaluator;
  /// Hash table of newly-created vertices at edge midpoints
  std::map<MBEdgeIndex, MBEntityHandle> edge_vertices;
  /// Hash table of newly-created vertices interior to triangle faces
  std::map<MBFaceIndex, MBEntityHandle> face_vertices;
  /// Storage allocated for interprocess communication of vertex global IDs
  std::map<int,MBRefinerPartition*> partitions;
  /// Number of tags defined on vertices
  int num_tags;
  /// The accumulated vertices that are used to create a new element when complete.
  std::vector<MBEntityHandle> vertex_accum;
  /**\brief True if the output mesh (this->mesh) is also the input mesh.
    *
    * This indicates whether pre-existing vertices should be inserted into the new mesh or not.
    */
  bool output_is_input;

  MBMeshRefinerOutputFunctor( MBInterface* m, bool oii, MBEdgeSizeEvaluator* es )
    {
    this->mesh = m;
    this->output_is_input = oii;
    this->edge_size_evaluator = es;
    }
  virtual ~MBMeshRefinerOutputFunctor() { }
  /// Insert a new vertex into the mesh and set its tag values. This has no effect when output_is_input is true.
  virtual MBEntityHandle operator () ( MBEntityHandle hash, const double* vcoords, const void* vtags )
    {
    if ( this->output_is_input )
      {
      // Don't insert vertices that already exist in the output
      return hash;
      }

    MBEntityHandle vertex_handle;
    this->mesh->create_vertex( vcoords + 3, vertex_handle );
    MBTag tag_handle;
    int tag_offset;
    for ( int i = 0; i < num_tags; ++i )
      {
      this->edge_size_evaluator->get_vertex_tag( i, tag_handle, tag_offset );
      this->mesh->tag_set_data( tag_handle, &vertex_handle, 1, (char*)vtags + tag_offset );
      }
    return vertex_handle;
    }
  /// Insert a new n-way hashed vertex into the mesh and set its tag values.
  virtual MBEntityHandle operator () ( int nhash, MBEntityHandle* hash, const double* vcoords, const void* vtags )
    {
    // First, see if we've already created this vertex
    MBEntityHandle hv;
    if ( this->find_hashed_vertex( hv, nhash, hash ) )
      {
      return hv;
      }

    this->mesh->create_vertex( vcoords + 3, hv );
    MBTag tag_handle;
    int tag_offset;
    for ( int i = 0; i < num_tags; ++i )
      {
      this->edge_size_evaluator->get_vertex_tag( i, tag_handle, tag_offset );
      this->mesh->tag_set_data( tag_handle, &hv, 1, (char*)vtags + tag_offset );
      }
    return hv;
    }
  /// Accumulate a vertex for use in the connectivity record of an element.
  virtual void operator () ( MBEntityHandle hash )
    {
    this->vertex_accum.push_back( hash );
    }
  /// Create an element from all the accumulated vertices since the last element was created.
  virtual void operator () ( MBEntityType etyp )
    {
    if ( ! this->vertex_accum.size() )
      return; // Ignore creation of vertex-less entities

    MBEntityHandle elem_handle;
    this->mesh->create_element( etyp, &this->vertex_accum[0], this->vertex_accum.size(), elem_handle );
    this->vertex_accum.clear();
    }
  bool find_hashed_vertex( MBEntityHandle& vout, int nhash, MBEntityHandle* hash )
    {
    return false;
    }
  void hash_vertex( MBEntityHandle vout, int nhash, MBEntityHandle* hash )
    {
    }
};

MBMeshRefiner::MBMeshRefiner( MBInterface* parent_mesh )
{  
  this->mesh = parent_mesh;
  this->entity_refiner = 0;
}

MBMeshRefiner::~MBMeshRefiner()
{
  if ( this->entity_refiner )
    delete this->entity_refiner;
}

bool MBMeshRefiner::set_entity_refiner( MBEntityRefiner* er )
{
  if ( ! er || er == this->entity_refiner )
    return false;

  this->entity_refiner = er;

  return true;
}

bool MBMeshRefiner::refine_mesh()
{
  return false;
}

