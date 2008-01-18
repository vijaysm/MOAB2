#include "MBMeshRefiner.hpp"

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
  /// Hash table of newly-created vertices at edge midpoints
  std::map<MBEdgeIndex, MBEntityHandle> edge_vertices;
  /// Hash table of newly-created vertices interior to triangle faces
  std::map<MBFaceIndex, MBEntityHandle> face_vertices;
  /// Storage allocated for interprocess communication of vertex global IDs
  std::map<int,MBRefinerPartition*> partitions;

  virtual ~MBMeshRefinerOutputFunctor() { }
  virtual void operator () ( const double* vcoords, const void* vtags )
    {
    // Assume that vtags contains the proper global ID for the vertex (assigned by the TagAssigner)
    }
  virtual void operator () ( MBEntityType etyp )
    {
    }
};

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

