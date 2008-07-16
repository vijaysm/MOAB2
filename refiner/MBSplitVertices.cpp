#include "MBSplitVertices.hpp"
#include "MBRefinerTagManager.hpp"

MBSplitVerticesBase::MBSplitVerticesBase( MBRefinerTagManager* tag_mgr )
{
  this->tag_manager = tag_mgr;
  this->mesh_in  = tag_mgr->get_input_mesh();
  this->mesh_out = tag_mgr->get_output_mesh();
  this->shared_procs_val.resize( MAX_SHARING_PROCS );
  MBParallelComm* ipcomm = MBParallelComm::get_pcomm( this->mesh_in, 0 );
  this->rank = ipcomm ? ipcomm->proc_config().proc_rank() : 0;
}

MBSplitVerticesBase::~MBSplitVerticesBase()
{
}

/// Determine which processes will contain an output vertex given the split vertices defining it.
void MBSplitVerticesBase::update_partition_counts(
  int num, const MBEntityHandle* split_src, std::map<MBProcessSet,int>& proc_partition_counts )
{
  this->begin_vertex_procs();
  for ( int i = 0; i < num; ++ i )
    {
    this->add_vertex_procs( split_src[i] );
    }
  this->end_vertex_procs();
  proc_partition_counts[this->common_shared_procs]++;
}

/// Prepare to compute the processes on which a new split-vertex will live.
void MBSplitVerticesBase::begin_vertex_procs()
{
  this->first_vertex = true;
  this->common_shared_procs.clear();
}

/// Call this for each existing corner vertex used to define a split-vertex.
void MBSplitVerticesBase::add_vertex_procs( MBEntityHandle vert_in )
{
  int stat;
  bool got = false;
  this->current_shared_procs.clear();
  stat = this->mesh_in->tag_get_data(
    this->tag_manager->shared_proc(), &vert_in, 1, &this->shared_procs_val[0] );
  if ( stat == MB_SUCCESS && this->shared_procs_val[0] != -1 )
    {
    got = true;
    std::cout << " s" << this->rank << " s" << this->shared_procs_val[0] << " | ";
    this->shared_procs_val[1] = -1;
    }
  stat = this->mesh_in->tag_get_data(
    this->tag_manager->shared_procs(), &vert_in, 1, &this->shared_procs_val[0] );
  if ( stat == MB_SUCCESS && this->shared_procs_val[0] != -1 )
    {
    got = true;
    int i;
    for ( i = 0; i < MAX_SHARING_PROCS && this->shared_procs_val[i] != -1; ++ i )
      std::cout << " m" << this->shared_procs_val[i];
    std::cout << " | ";
    }
  if ( got )
    {
    this->current_shared_procs.set_process_members( this->shared_procs_val );
    this->current_shared_procs.set_process_member( this->rank );
    if ( this->first_vertex )
      {
      this->common_shared_procs.unite( this->current_shared_procs );
      this->first_vertex = false;
      }
    else
      {
      this->common_shared_procs.intersect( this->current_shared_procs );
      }
    }
  else
    {
    std::cout << " not shared | ";
    }
}

/// Call this once after all the add_vertex_procs() calls for a split-vertex to prepare queues for the second stage MPI send. 
void MBSplitVerticesBase::end_vertex_procs()
{
  std::cout << "    Common procs " << this->common_shared_procs;
  std::cout << "\n";
  // FIXME: Here is where we add the vertex to the appropriate queues.
}

