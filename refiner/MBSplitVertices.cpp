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

MBEntitySource::MBEntitySource( int nc, MBRefinerTagManager* tag_mgr )
{
  this->tag_manager = tag_mgr;
  this->mesh_out = tag_mgr->get_output_mesh();
  this->num_corners = nc;
}

MBEntitySource::~MBEntitySource()
{
}

bool MBEntitySource::create_element(
  MBEntityType etyp, int nconn, const MBEntityHandle* elem_verts, MBEntityHandle& elem_handle,
  std::map<MBProcessSet,int>& proc_partition_counts )
{
  // Get the global IDs of the input vertices
  int stat;
  proc_partition_counts[this->tag_manager->get_element_procs()]++;
  if ( this->mesh_out->create_element( etyp, elem_verts, nconn, elem_handle ) != MB_SUCCESS )
    {
    return false;
    }
  this->push_back( MBEntitySourceRecord( this->num_corners, elem_handle, this->tag_manager->get_element_procs() ) );
  this->tag_manager->set_sharing( elem_handle, this->tag_manager->get_element_procs() );
  return true;
}

void MBEntitySource::assign_global_ids( std::map<MBProcessSet,int>& gids )
{
  std::vector<MBEntityHandle> adjacencies;
  adjacencies.resize( this->num_corners );
  std::vector<MBEntitySourceRecord>::iterator it;
  int stat;
  for ( it = this->begin(); it != this->end(); ++ it )
    {
    int num_nodes;
    const MBEntityHandle* conn;
    this->mesh_out->get_connectivity( it->handle, conn, num_nodes );
    stat = this->tag_manager->get_output_gids( this->num_corners, conn, it->ids );
    std::sort( it->ids.begin(), it->ids.end() );
    }
  std::sort( this->begin(), this->end() );
  for ( it = this->begin(); it != this->end(); ++ it )
    {
    int gid = gids[it->process_set] ++;
    this->tag_manager->set_gid( it->handle, gid );
    std::cout << "Assigning entity: " << it->handle << " GID: " << gid << "\n";
    }
}

