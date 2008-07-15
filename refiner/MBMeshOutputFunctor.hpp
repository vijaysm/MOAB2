#ifndef MB_MESHOUTPUTFUNCTOR_HPP
#define MB_MESHOUTPUTFUNCTOR_HPP

#include "MBTypes.h"
#include "MBEntityRefiner.hpp"
#include "MBParallelComm.hpp"

#include <iostream>
#include <set>
#include <map>
#include <iterator>
#include <algorithm>

#include <string.h>

/**\brief Represent a set of processes using a bit vector.
  *
  * This is used by the mesh refiner when determining where to record
  * split vertices so that labeling can be inferred across process
  * boundaries without communicating anything other than the number of
  * entities in a given partition.
  */
class MBProcessSet
{
public:
  enum
    {
    SHARED_PROC_BYTES = (MAX_SHARING_PROCS / 8 + (MAX_SHARING_PROCS % 8 ? 1 : 0))
    };

  MBProcessSet()
    {
    this->clear();
    }

  MBProcessSet( const unsigned char* psetbits )
    {
    for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
      this->processes[i] = psetbits[i];
    }

  ~MBProcessSet()
    {
    }

  void unite( const MBProcessSet& other )
    {
    for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
      {
      this->processes[i] |= other.processes[i];
      }
    }

  void intersect( const MBProcessSet& other )
    {
    for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
      {
      this->processes[i] &= other.processes[i];
      }
    }

  void clear()
    {
    memset( this->processes, 0, SHARED_PROC_BYTES );
    }

  void set_process_member( int i )
    {
    int byte = i / 8;
    int bitmask = 1 << ( i % 8 );
    this->processes[byte] |= bitmask;
    }

  void set_process_members( const std::vector<int>& procs )
    {
    for ( std::vector<int>::const_iterator it = procs.begin(); it != procs.end() && *it != -1; ++ it )
      {
      this->set_process_member( *it );
      }
    }

  bool is_process_member( int i ) const
    {
    int byte = i / 8;
    int bitmask = 1 << ( i % 8 );
    return ( this->processes[byte] & bitmask ) ? true : false;
    }

  const unsigned char* data() const
    {
    return this->processes;
    }

  bool operator < ( const MBProcessSet& other ) const
    {
    for ( int i = 0; i < SHARED_PROC_BYTES; ++ i )
      {
      if ( this->processes[i] < other.processes[i] )
        return true;
      else if ( this->processes[i] > other.processes[i] )
        return false;
      }
    return false; // equality
    }

  friend std::ostream& operator << ( std::ostream& os, const MBProcessSet& pset )
    {
    for ( int i = 0; i < MAX_SHARING_PROCS; ++ i )
      {
      os << ( pset.is_process_member( i ) ? "1" : "0" );
      }
    return os;
    }

protected:
  unsigned char processes[SHARED_PROC_BYTES];
};

template< int _n >
class MBSplitVertexIndex
{
public:
  MBSplitVertexIndex() { }
  MBSplitVertexIndex( const MBEntityHandle* src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src[i]; std::sort( this->handles, this->handles + _n ); }
  MBSplitVertexIndex( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; }
  MBSplitVertexIndex& operator = ( const MBSplitVertexIndex<_n>& src )
    { for ( int i = 0; i < _n; ++ i ) this->handles[i] = src.handles[i]; return *this; }

  bool operator < ( const MBSplitVertexIndex<_n>& other ) const
    {
    for ( int i = 0; i < _n; ++ i )
      if ( this->handles[i] < other.handles[i] )
        return true;
      else if ( this->handles[i] > other.handles[i] )
        return false;
    return false;
    }

  MBEntityHandle handles[_n];
};

class MBSplitVerticesBase
{
public:
  MBSplitVerticesBase( MBRefinerTagManager* tag_mgr )
    {
    this->tag_manager = tag_mgr;
    this->mesh_in  = tag_mgr->get_input_mesh();
    this->mesh_out = tag_mgr->get_output_mesh();
    this->shared_procs_val.resize( MAX_SHARING_PROCS );
    MBParallelComm* ipcomm = MBParallelComm::get_pcomm( this->mesh_in, 0 );
    this->rank = ipcomm ? ipcomm->proc_config().proc_rank() : 0;
    }
  virtual ~MBSplitVerticesBase() { }
  virtual bool find_or_create(
    const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
    std::map<MBProcessSet,int>& proc_partition_counts ) = 0;

  /// Determine which processes will contain an output vertex given the split vertices defining it.
  void update_partition_counts( int num, const MBEntityHandle* split_src, std::map<MBProcessSet,int>& proc_partition_counts )
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
  void begin_vertex_procs()
    {
    this->first_vertex = true;
    this->common_shared_procs.clear();
    }

  /// Call this for each existing corner vertex used to define a split-vertex.
  void add_vertex_procs( MBEntityHandle vert_in )
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
  void end_vertex_procs()
    {
    std::cout << "    Common procs " << this->common_shared_procs;
    std::cout << "\n";
    // FIXME: Here is where we add the vertex to the appropriate queues.
    }

  MBInterface* mesh_in; // Input mesh. Needed to determine tag values on split_src verts
  MBInterface* mesh_out; // Output mesh. Needed for new vertex set in vert_handle
  MBRefinerTagManager* tag_manager;
  std::vector<int> shared_procs_val; // Used to hold procs sharing an input vert.
  MBProcessSet current_shared_procs; // Holds process list as it is being accumulated
  MBProcessSet common_shared_procs; // Holds intersection of several shared_procs_vals.
  int rank; // This process' rank.
  bool first_vertex; // True just after begin_vertex_procs() is called.
};

template< int _n >
class MBSplitVertices : public std::map<MBSplitVertexIndex<_n>,MBEntityHandle>, public MBSplitVerticesBase
{
public:
  typedef std::map<MBSplitVertexIndex<_n>,MBEntityHandle> MapType;
  typedef typename std::map<MBSplitVertexIndex<_n>,MBEntityHandle>::iterator MapIteratorType;

  MBSplitVertices( MBRefinerTagManager* tag_mgr )
    : MBSplitVerticesBase( tag_mgr )
    {
    this->shared_procs_val.resize( _n * MAX_SHARING_PROCS );
    }
  virtual ~MBSplitVertices() { }
  virtual bool find_or_create(
    const MBEntityHandle* split_src, const double* coords, MBEntityHandle& vert_handle,
    std::map<MBProcessSet,int>& proc_partition_counts )
    {
    MBSplitVertexIndex<_n> key( split_src );
    MapIteratorType it = this->find( key );
    if ( it == this->end() )
      {
      if ( this->mesh_out->create_vertex( coords, vert_handle ) != MB_SUCCESS )
        {
        return false;
        }
      (*this)[key] = vert_handle;
      this->update_partition_counts( _n, split_src, proc_partition_counts );
      return true;
      }
    vert_handle = it->second;
    return false;
    }
};

class MBMeshOutputFunctor : public MBEntityRefinerOutputFunctor
{
public:
  MBInterface* mesh_in;
  MBInterface* mesh_out;
  bool input_is_output;
  std::vector<MBSplitVerticesBase*> split_vertices;
  std::vector<MBEntityHandle> elem_vert;
  MBRefinerTagManager* tag_manager;
  MBEntityHandle destination_set;
  std::map<MBProcessSet,int> proc_partition_counts;

  MBMeshOutputFunctor( MBRefinerTagManager* tag_mgr )
    {
    this->mesh_in  = tag_mgr->get_input_mesh();
    this->mesh_out = tag_mgr->get_output_mesh();
    this->input_is_output = ( this->mesh_in == this->mesh_out );
    this->tag_manager = tag_mgr;
    this->destination_set = 0; // don't place output entities in a set by default.

    this->split_vertices.resize( 4 );
    this->split_vertices[0] = 0; // Vertices (0-faces) cannot be split
    this->split_vertices[1] = new MBSplitVertices<1>( this->tag_manager );
    this->split_vertices[2] = new MBSplitVertices<2>( this->tag_manager );
    this->split_vertices[3] = new MBSplitVertices<3>( this->tag_manager );
    }

  ~MBMeshOutputFunctor()
    {
    for ( int i = 0; i < 4; ++ i )
      delete this->split_vertices[i];
    }

  void print_vert_crud( MBEntityHandle vout, int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags )
    {
    std::cout << "+ {";
    for ( int i = 0; i < nvhash; ++ i )
      std::cout << " " << vhash[i];
    std::cout << " } -> " << vout << " ";

    std::cout << "[ " << vcoords[0];
    for ( int i = 1; i < 6; ++i )
      std::cout << ", " << vcoords[i];
    std::cout << " ] ";

#if 0
    double* x = (double*)vtags;
    int* m = (int*)( (char*)vtags + 2 * sizeof( double ) );
    std::cout << "< " << x[0]
              << ", " << x[1];
    for ( int i = 0; i < 4; ++i )
      std::cout << ", " << m[i];
#endif // 0

    std::cout << " >\n";
    }

  void print_partition_counts( MBParallelComm* comm )
    {
    int lnparts = this->proc_partition_counts.size();
    std::vector<unsigned char> lpdefns;
    std::vector<int> lpsizes;
    lpdefns.resize( MBProcessSet::SHARED_PROC_BYTES * lnparts );
    lpsizes.resize( lnparts );
    std::cout << "**** Partition Counts ****\n";
    int i = 0;
    std::map<MBProcessSet,int>::iterator it;
    for ( it = this->proc_partition_counts.begin(); it != this->proc_partition_counts.end(); ++ it, ++ i )
      {
      for ( int j = 0; j < MBProcessSet::SHARED_PROC_BYTES; ++ j )
        lpdefns[MBProcessSet::SHARED_PROC_BYTES * i + j] = it->first.data()[j];
      lpsizes[i] = it->second;
      std::cout << "Partition " << it->first << ": " << it->second << "\n";
      }

    if ( ! comm )
      return;

    std::vector<int> nparts;
    std::vector<int> dparts;
    unsigned long prank = comm->proc_config().proc_rank();
    unsigned long psize = comm->proc_config().proc_size();
    nparts.resize( psize );
    dparts.resize( psize + 1 );
    MPI_Allgather( &lnparts, 1, MPI_INT, &nparts[0], 1, MPI_INT, comm->proc_config().proc_comm() );
    unsigned long ndefs = 0;
    for ( int rank = 1; rank <= psize; ++ rank )
      {
      dparts[rank] = nparts[rank - 1] + dparts[rank - 1];
      std::cout << "Proc " << rank << ": " << nparts[rank-1] << " partitions, offset: " << dparts[rank] << "\n";
      }
    std::vector<unsigned char> part_defns;
    std::vector<int> part_sizes;
    part_defns.resize( MBProcessSet::SHARED_PROC_BYTES * dparts[psize] );
    part_sizes.resize( dparts[psize] );
    MPI_Allgatherv(
      &lpsizes[0], lnparts, MPI_INT,
      &part_sizes[0], &nparts[0], &dparts[0], MPI_INT, comm->proc_config().proc_comm() );
    for ( int rank = 0; rank < psize; ++ rank )
      {
      nparts[rank] *= MBProcessSet::SHARED_PROC_BYTES;
      dparts[rank] *= MBProcessSet::SHARED_PROC_BYTES;
      }
    MPI_Allgatherv(
      &lpdefns[0], MBProcessSet::SHARED_PROC_BYTES * lnparts, MPI_UNSIGNED_CHAR,
      &part_defns[0], &nparts[0], &dparts[0], MPI_UNSIGNED_CHAR, comm->proc_config().proc_comm() );
    for ( int i = 0; i < dparts[psize]; ++ i )
      {
      MBProcessSet pset( &part_defns[MBProcessSet::SHARED_PROC_BYTES * i] );
      std::map<MBProcessSet,int>::iterator it = this->proc_partition_counts.find( pset );
      if ( it != this->proc_partition_counts.end() )
        {
        std::cout << "Partition " << pset << ( it->second == part_sizes[i] ? " matches" : " broken" ) << ".\n";
        }
      else
        {
        this->proc_partition_counts[pset] = part_sizes[i];
        }
      }
    std::map<MBProcessSet,int>::iterator pcit;
    for ( pcit = this->proc_partition_counts.begin(); pcit != this->proc_partition_counts.end(); ++ pcit )
      {
      std::cout << "Partition " << pcit->first << ": " << pcit->second << " #\n";
      }
    }

  void assign_tags( MBEntityHandle vhandle, const void* vtags )
    {
    if ( ! vhandle )
      return; // Ignore bad vertices

    int num_tags = this->tag_manager->get_number_of_vertex_tags();
    MBTag tag_handle;
    int tag_offset;
    for ( int i = 0; i < num_tags; ++i )
      {
      this->tag_manager->get_output_vertex_tag( i, tag_handle, tag_offset );
      this->mesh_out->tag_set_data( tag_handle, &vhandle, 1, vtags );
      }
    }

  virtual MBEntityHandle operator () ( MBEntityHandle vhash, const double* vcoords, const void* vtags )
    {
    if ( this->input_is_output )
      { // Don't copy the original vertex!
      this->print_vert_crud( vhash, 1, &vhash, vcoords, vtags );
      return vhash;
      }
    MBEntityHandle vertex_handle;
    if ( this->mesh_out->create_vertex( vcoords + 3, vertex_handle ) != MB_SUCCESS )
      {
      std::cerr << "Could not insert mid-edge vertex!\n";
      }
    this->assign_tags( vertex_handle, vtags );
    this->print_vert_crud( vertex_handle, 1, &vhash, vcoords, vtags );
    return vertex_handle;
    }

  virtual MBEntityHandle operator () ( int nvhash, MBEntityHandle* vhash, const double* vcoords, const void* vtags )
    {
    MBEntityHandle vertex_handle;
    if ( nvhash == 1 )
      {
      vertex_handle = (*this)( *vhash, vcoords, vtags );
      }
    else if ( nvhash < 4 )
      {
      bool newly_created = this->split_vertices[nvhash]->find_or_create(
        vhash, vcoords, vertex_handle, this->proc_partition_counts );
      if ( newly_created )
        {
        this->assign_tags( vertex_handle, vtags );
        }
      if ( ! vertex_handle )
        {
        std::cerr << "Could not insert mid-edge vertex!\n";
        }
      this->print_vert_crud( vertex_handle, nvhash, vhash, vcoords, vtags );
      }
    else
      {
      vertex_handle = 0;
      std::cerr << "Not handling splits on faces with " << nvhash << " corners yet.\n";
      }
    return vertex_handle;
    }

  virtual void operator () ( MBEntityHandle h )
    {
    std::cout << h << " ";
    this->elem_vert.push_back( h );
    }

  virtual void operator () ( MBEntityType etyp )
    {
    MBEntityHandle elem_handle;
    if ( this->mesh_out->create_element( etyp, &this->elem_vert[0], this->elem_vert.size(), elem_handle ) == MB_FAILURE )
      {
      std::cerr << " *** ";
      }
    this->elem_vert.clear();
    std::cout << "---------> " << elem_handle << " ( " << etyp << " )\n\n";
    }
};

#endif // MB_MESHOUTPUTFUNCTOR_HPP
