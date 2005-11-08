
#undef DEBUG

#ifdef DEBUG
#  include <stdio.h>
#  include <stdarg.h>
#endif

#ifndef HDF5_FILE
#  error Attempt to compile WriteHDF5Parallel with HDF5 support disabled
#endif

#include <stdlib.h>
#include <string.h>

#include <vector>
#include <set>
#include <map>
#include <utility>

#include <mpi.h>

#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include <H5FDmpi.h>
#include <H5FDmpio.h>

#include "mhdf.h"

#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"

#include "WriteHDF5Parallel.hpp"


#ifdef DEBUG
#  define START_SERIAL                     \
     for (int _x = 0; _x < numProc; ++_x) {\
       MPI_Barrier( MPI_COMM_WORLD );      \
       if (_x != myRank) continue     
#  define END_SERIAL                       \
     }                                     \
     MPI_Barrier( MPI_COMM_WORLD )
#else
#  define START_SERIAL
#  define END_SERIAL
#endif


#define DEBUG_OUT_STREAM stdout

#ifndef DEBUG
static void printdebug( const char*, ... ) {}
#else
static void printdebug( const char* fmt, ... )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  fprintf( DEBUG_OUT_STREAM, "[%d] ", rank );
  va_list args;
  va_start( args, fmt );
  vfprintf( DEBUG_OUT_STREAM, fmt, args );
  va_end( args );
  fflush( DEBUG_OUT_STREAM );
}
#endif

#ifndef DEBUG
void WriteHDF5Parallel::printrange( MBRange& ) {}
#else
void WriteHDF5Parallel::printrange( MBRange& r )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MBEntityType type = MBMAXTYPE;
  for (MBRange::const_pair_iterator i = r.pair_begin(); i != r.pair_end(); ++i)
  {
    MBEntityHandle a, b;
    a = (*i).first;
    b = (*i).second;
    MBEntityType mytype = iFace->type_from_handle(a);
    if (mytype != type)
    {
      type = mytype;
      fprintf(DEBUG_OUT_STREAM, "%s[%d]  %s", type == MBMAXTYPE ? "" : "\n", rank, MBCN::EntityTypeName( type ) );
    }
    unsigned long id1 = iFace->id_from_handle( a );
    unsigned long id2 = iFace->id_from_handle( b );
    if (id1 == id2)
      fprintf(DEBUG_OUT_STREAM, " %lu", id1 );
    else
      fprintf(DEBUG_OUT_STREAM, " %lu-%lu", id1, id2 );
  }
  fprintf(DEBUG_OUT_STREAM, "\n");
  fflush( DEBUG_OUT_STREAM );
}
#endif


#ifndef DEBUG
static void print_type_sets( MBInterface* , int , int , MBRange& ) {}
#else
static void print_type_sets( MBInterface* iFace, int myRank, int numProc, MBRange& sets )
{
  MBTag gid, did, bid, sid, nid, iid;
  iFace->tag_get_handle( GLOBAL_ID_TAG_NAME, gid ); 
  iFace->tag_get_handle( GEOM_DIMENSION_TAG_NAME, did );
  iFace->tag_get_handle( MATERIAL_SET_TAG_NAME, bid );
  iFace->tag_get_handle( DIRICHLET_SET_TAG_NAME, nid );
  iFace->tag_get_handle( NEUMANN_SET_TAG_NAME, sid );
  iFace->tag_get_handle( PARALLEL_INTERFACE_TAG_NAME, iid );
  MBRange typesets[10];
  const char* typenames[] = {"Block", "Sideset", "NodeSet", "Vertex", "Curve", "Surface", "Volume", "Body", "Interfaces", "Other"};
  for (MBRange::iterator riter = sets.begin(); riter != sets.end(); ++riter)
  {
    unsigned dim, id, proc[2], oldsize;
    if (MB_SUCCESS == iFace->tag_get_data(bid, &*riter, 1, &id)) 
      dim = 0;
    else if (MB_SUCCESS == iFace->tag_get_data(sid, &*riter, 1, &id))
      dim = 1;
    else if (MB_SUCCESS == iFace->tag_get_data(nid, &*riter, 1, &id))
      dim = 2;
    else if (MB_SUCCESS == iFace->tag_get_data(did, &*riter, 1, &dim)) {
      id = 0;
      iFace->tag_get_data(gid, &*riter, 1, &id);
      dim += 3;
    }
    else if (MB_SUCCESS == iFace->tag_get_data(iid, &*riter, 1, proc)) {
      assert(proc[0] == (unsigned)myRank || proc[1] == (unsigned)myRank);
      id = proc[proc[0] == (unsigned)myRank];
      dim = 8;
    }
    else {
      id = *riter;
      dim = 9;
    }

    oldsize = typesets[dim].size();
    typesets[dim].insert( id );
    assert( typesets[dim].size() - oldsize == 1 );  
  }
  for (int ii = 0; ii < 10; ++ii)
  {
    char num[16];
    std::string line(typenames[ii]);
    sprintf(num, "(%u):", typesets[ii].size());
    line += num;
    for (MBRange::const_pair_iterator piter = typesets[ii].pair_begin();
         piter != typesets[ii].pair_end(); ++piter)
    {
      sprintf(num," %d", (*piter).first);
      line += num;
      if ((*piter).first != (*piter).second) {
        sprintf(num,"-%d", (*piter).second);
        line += num;
      }
    }

    printdebug ("%s\n", line.c_str());
  }
  printdebug("Total: %u\n", sets.size());
}
#endif

#ifdef NDEBUG
#  define assert(A)
#else
#  define assert(A) if (!(A)) do_assert(__FILE__, __LINE__, #A)
   static void do_assert( const char* file, int line, const char* condstr )
   {
     int rank;
     MPI_Comm_rank( MPI_COMM_WORLD, &rank );
     fprintf( DEBUG_OUT_STREAM, "[%d] Assert(%s) failed at %s:%d\n", rank, condstr, file, line );
     fflush( DEBUG_OUT_STREAM );
     abort();
   }
#endif


void range_remove( MBRange& from, const MBRange& removed )
{
  
/* The following should be more efficient, but isn't due
   to the inefficient implementation of MBRange::erase(iter,iter)
  MBRange::const_iterator s, e, n = from.begin();
  for (MBRange::const_pair_iterator p = removed.pair_begin();
       p != removed.pair_end(); ++p)
  {
    e = s = MBRange::lower_bound(n, from.end(), (*p).first);
    e = MBRange::lower_bound(s, from.end(), (*p).second);
    if (e != from.end() && *e == (*p).second)
      ++e;
    n = from.erase( s, e );
  }
*/

  if (removed.size())
    from = from.subtract(removed);
}


WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface,
                                      const char** multiproc_tags )
  : WriteHDF5(iface)
{
  if (multiproc_tags)
  {
    while (*multiproc_tags)
    {
      multiProcSetTags.push_back( *multiproc_tags );
      ++multiproc_tags;
    }
  }
  else
  {
    multiProcSetTags.push_back(  MATERIAL_SET_TAG_NAME );
    multiProcSetTags.push_back( DIRICHLET_SET_TAG_NAME );
    multiProcSetTags.push_back(   NEUMANN_SET_TAG_NAME );
  }  
  std::sort(multiProcSetTags.begin(), multiProcSetTags.end());
}


MBErrorCode WriteHDF5Parallel::gather_interface_meshes()
{
  MBRange range;
  MBErrorCode result;
  MBTag iface_tag, geom_tag;
  int i, proc_pair[2];
  
  START_SERIAL;
  printdebug( "Pre-interface mesh:\n");
  printrange(nodeSet.range);
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    printrange(eiter->range);
  printrange(setSet.range);
  
    // Allocate space for remote mesh data
  remoteMesh.resize( numProc );
  
    // Get tag handles
  result = iFace->tag_get_handle( PARALLEL_INTERFACE_TAG_NAME, iface_tag );
  if (MB_SUCCESS != result) return result;
  result = iFace->tag_get_handle( PARALLEL_GEOM_TOPO_TAG_NAME, geom_tag );
  if (MB_SUCCESS != result) return result;
  
  
    // Get interface mesh sets
  result = iFace->get_entities_by_type_and_tag( 0,
                                                MBENTITYSET,
                                                &iface_tag,
                                                0,
                                                1,
                                                range );
  if (MB_SUCCESS != result) return result;
  
  
    // Populate lists of interface mesh entities
  for (MBRange::iterator iiter = range.begin(); iiter != range.end(); ++iiter)
  {
    result = iFace->tag_get_data( iface_tag, &*iiter, 1, proc_pair );
    if (MB_SUCCESS != result) return result;
    const int remote_proc = proc_pair[0];
    
      // Get list of all entities in interface and 
      // the subset of that list that are meshsets.
    MBRange entities, sets;
    result = iFace->get_entities_by_handle( *iiter, entities );
    if (MB_SUCCESS != result) return result;
    result = iFace->get_entities_by_type( *iiter, MBENTITYSET, sets );
    if (MB_SUCCESS != result) return result;

      // Put any non-meshset entities in the list directly.
    range_remove( entities, sets );
    remoteMesh[remote_proc].merge( entities );
    remoteMesh[remote_proc].insert( *iiter );
    
    for (MBRange::iterator siter = sets.begin(); siter != sets.end(); ++siter)
    {
        // For current parallel meshing code, root processor owns
        // all curve and geometric vertex meshes.  
      int dimension;
      result = iFace->tag_get_data( geom_tag, &*siter, 1, &dimension );
      if (result == MB_SUCCESS && dimension < 2)
        continue;
        
        // Put entities in list for appropriate processor.
      remoteMesh[remote_proc].insert( *siter );
      entities.clear();
      result = iFace->get_entities_by_handle( *siter, entities );
      if (MB_SUCCESS != result) return result;
      remoteMesh[remote_proc].merge( entities );
    }
  }
  
    // For current parallel meshing code, root processor owns
    // all curve and geometric vertex meshes.  Find them and
    // allocate them appropriately.
  MBRange curves_and_verts;
  MBTag tags[] = { geom_tag, geom_tag };
  int value_ints[] = { 0, 1 };
  const void* values[] = {value_ints, value_ints + 1};
  result = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET,
                                                tags, values, 2,
                                                curves_and_verts, 
                                                MBInterface::UNION );
                                                assert(MB_SUCCESS == result);
  MBRange edges, nodes;
  for (MBRange::iterator riter = curves_and_verts.begin();
       riter != curves_and_verts.end(); ++riter)
  {
    result = iFace->get_entities_by_type( *riter, MBVERTEX, nodes ); assert(MB_SUCCESS == result);
    result = iFace->get_entities_by_type( *riter, MBEDGE, edges ); assert(MB_SUCCESS == result);
  }
  std::list<ExportSet>::iterator eiter = exportList.begin();
  for ( ; eiter != exportList.end() && eiter->type != MBEDGE; ++eiter );
  
  if (myRank == 0)
  {
    nodeSet.range.merge( nodes );
    setSet.range.merge(curves_and_verts);
    #ifdef DEBUG
    setSet.range.sanity_check();
    #endif
    eiter->range.merge( edges );
  } 
  else
  {
    MBRange diff = nodeSet.range.intersect( nodes );
    range_remove( nodeSet.range, diff );
    remoteMesh[0].merge( diff );
    
    diff = eiter->range.intersect( edges );
    range_remove( eiter->range, diff );
    remoteMesh[0].merge( diff );
    
    diff = setSet.range.intersect( curves_and_verts );
    range_remove( setSet.range, diff );
    #ifdef DEBUG
    setSet.range.sanity_check();
    #endif
    remoteMesh[0].merge( diff );
  }
  edges.merge(nodes);
  edges.merge(curves_and_verts);
  for (i = 1; i < numProc; i++)
  {
    MBRange diff = edges.intersect( remoteMesh[i] );
    range_remove(remoteMesh[i], diff);
    remoteMesh[0].merge(diff);
  }
  
  
  
    // For all remote mesh entities, remove them from the
    // lists of local mesh to be exported and give them a 
    // junk file Id of 1.  Need to specify a file ID greater
    // than zero so the code that gathers adjacencies and 
    // such doesn't think that the entities aren't being
    // exported.
  for (i = 0; i < numProc; i++)
  {
    if (i == myRank) continue;
    
    MBRange& range = remoteMesh[i];
    
    range_remove( nodeSet.range, range );
    range_remove( setSet.range, range );
    #ifdef DEBUG
    setSet.range.sanity_check();
    #endif
    for (std::list<ExportSet>::iterator eiter = exportList.begin();
         eiter != exportList.end(); ++eiter )
      range_remove( eiter->range, range );
    
    int id = 1;
    for (MBRange::iterator riter = remoteMesh[i].begin(); 
         riter != remoteMesh[i].end() && iFace->type_from_handle(*riter) != MBENTITYSET; 
         ++riter)
    {
      result = iFace->tag_set_data( idTag, &*riter, 1, &id );
      if (MB_SUCCESS != result) return result;
    }
    
    //allRemoteMesh.merge( range );
  }
  
  printdebug("Remote mesh:\n");
  for (int ii = 0; ii < numProc; ++ii)
  {
    printdebug("  proc %d : %d\n", ii, remoteMesh[ii].size());
    printrange( remoteMesh[ii] );
  }

  printdebug( "Post-interface mesh:\n");
  printrange(nodeSet.range);
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    printrange(eiter->range);
  printrange(setSet.range);

  END_SERIAL;
  
  return MB_SUCCESS;
}



MBErrorCode WriteHDF5Parallel::create_file( const char* filename,
                                            bool overwrite,
                                            std::vector<std::string>& qa_records,
                                            int dimension )
{
  MBErrorCode rval;
  int result;
  mhdf_Status status;
    
  result = MPI_Comm_rank( MPI_COMM_WORLD, &myRank );
  assert(MPI_SUCCESS == result);
  result = MPI_Comm_size( MPI_COMM_WORLD, &numProc );
  assert(MPI_SUCCESS == result);
  
  rval = gather_interface_meshes();
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create actual file and write meta info ***************/

  if (myRank == 0)
  {
      // create the file
    const char* type_names[MBMAXTYPE];
    memset( type_names, 0, MBMAXTYPE * sizeof(char*) );
    for (MBEntityType i = MBEDGE; i < MBENTITYSET; ++i)
      type_names[i] = MBCN::EntityTypeName( i );
   
    filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, &status );
    if (!filePtr)
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
    rval = write_qa( qa_records );
    if (MB_SUCCESS != rval) return rval;
  }
  
  
     /**************** Create node coordinate table ***************/
 
  rval = create_node_table( dimension );
  if (MB_SUCCESS != rval) return rval;
  
  
    /**************** Create element tables ***************/

  rval = negotiate_type_list();
  if (MB_SUCCESS != rval) return rval;
  rval = create_element_tables();
  if (MB_SUCCESS != rval) return rval;
  

    /**************** Comminucate all remote IDs ***********************/
  
  rval = communicate_remote_ids( MBVERTEX );
  for (std::list<ExportSet>::iterator ex_itor = exportList.begin(); 
       ex_itor != exportList.end(); ++ex_itor)
  {
    rval = communicate_remote_ids( ex_itor->type );
    assert(MB_SUCCESS == rval);
  }
  
  
    /**************** Create adjacency tables *********************/
  
  rval = create_adjacency_tables();
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create meshset tables *********************/
  
  rval = create_meshset_tables();
  if (MB_SUCCESS != rval) return rval;
  
  
    /* Need to write tags for shared sets this proc is responsible for */
  
  MBRange parallel_sets;
  for (std::list<ParallelSet>::const_iterator psiter = parallelSets.begin();
       psiter != parallelSets.end(); ++psiter)
    if (psiter->description)
      parallel_sets.insert( psiter->handle );
  
  setSet.range.merge( parallel_sets );
  rval = gather_tags();
  if (MB_SUCCESS != rval)
    return rval;
  range_remove( setSet.range, parallel_sets );   
  #ifdef DEBUG
  setSet.range.sanity_check();
  #endif
  

    /**************** Create tag data *********************/
  
  std::list<SparseTag>::iterator tag_iter;
  sort_tags_by_name();
  const int num_tags = tagList.size();
  std::vector<int> tag_offsets(num_tags), tag_counts(num_tags);
  std::vector<int>::iterator tag_off_iter = tag_counts.begin();
  for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter, ++tag_off_iter)
    *tag_off_iter = tag_iter->range.size();
  
  printdebug("Exchanging tag data for %d tags.\n", num_tags);
  std::vector<int> proc_tag_offsets(num_tags*numProc);
  result = MPI_Gather( &tag_counts[0], num_tags, MPI_INT,
                 &proc_tag_offsets[0], num_tags, MPI_INT,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    tag_counts[i] = 0;
    int next_offset = 0;
    for (int j = 0; j < numProc; j++)
    {
      int count = proc_tag_offsets[i + j*num_tags];
      proc_tag_offsets[i + j*num_tags] = next_offset;
      next_offset += count;
      tag_counts[i] += count;
    }

    if (0 == myRank)
    {
      rval = create_tag( tag_iter->tag_id, next_offset );
      assert(MB_SUCCESS == rval);
      printdebug( "Creating table of size %d for tag 0x%lx\n", (int)next_offset, (unsigned long)tag_iter->tag_id);
    }
  }
  
  result = MPI_Bcast( &tag_counts[0], num_tags, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
  result = MPI_Scatter( &proc_tag_offsets[0], num_tags, MPI_INT,
                             &tag_offsets[0], num_tags, MPI_INT,
                             0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);


  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    tag_iter->offset = tag_offsets[i];
    tag_iter->write = tag_counts[i] > 0;
  }

  #ifdef DEBUG
  START_SERIAL;  
  printdebug("Tags: %16s %8s %8s %8s\n", "Name", "Count", "Offset", "Handle");

  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    std::string name;
    iFace->tag_get_name( tag_iter->tag_id, name );
    printdebug("      %16s %8d %8d %8lx\n", name.c_str(), tag_counts[i], tag_offsets[i], (unsigned long)tag_iter->tag_id );
  }
  END_SERIAL;  
  #endif
  
  /************** Close serial file and reopen parallel *****************/
  
  if (0 == myRank)
  {
    mhdf_closeFile( filePtr, &status );
  }
  
  unsigned long junk;
  hid_t hdf_opt = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_fapl_mpio( hdf_opt, MPI_COMM_WORLD, MPI_INFO_NULL );
  filePtr = mhdf_openFileWithOpt( filename, 1, &junk, hdf_opt, &status );
  if (!filePtr)
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return MB_FAILURE;
  }
  
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5Parallel::create_node_table( int dimension )
{
  int result;
  mhdf_Status status;
 
    // gather node counts for each processor
  std::vector<int> node_counts(numProc);
  int num_nodes = nodeSet.range.size();
  result = MPI_Gather( &num_nodes, 1, MPI_INT, &node_counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // create node data in file
  long first_id;
  if (myRank == 0)
  {
    int total = 0;
    for (int i = 0; i < numProc; i++)
      total += node_counts[i];
      
    hid_t handle = mhdf_createNodeCoords( filePtr, dimension, total, &first_id, &status );
    if (mhdf_isError( &status ))
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
    mhdf_closeData( filePtr, handle, &status );
 }
    
    // send id offset to every proc
  result = MPI_Bcast( &first_id, 1, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  nodeSet.first_id = (id_t)first_id;
   
      // calculate per-processor offsets
  if (myRank == 0)
  {
    int prev_size = node_counts[0];
    node_counts[0] = 0;
    for (int i = 1; i < numProc; ++i)
    {
      int mysize = node_counts[i];
      node_counts[i] = node_counts[i-1] + prev_size;
      prev_size = mysize;
    }
  }
  
    // send each proc it's offset in the node table
  int offset;
  result = MPI_Scatter( &node_counts[0], 1, MPI_INT, 
                        &offset, 1, MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  nodeSet.offset = offset;
  
  writeUtil->assign_ids( nodeSet.range, idTag, (id_t)(nodeSet.first_id + nodeSet.offset) );

  return MB_SUCCESS;
}



struct elemtype {
  int mbtype;
  int numnode;
  
  elemtype( int vals[2] ) : mbtype(vals[0]), numnode(vals[1]) {}
  elemtype( int t, int n ) : mbtype(t), numnode(n) {}
  
  bool operator==( const elemtype& other ) const
  {
    return mbtype == other.mbtype &&
            (mbtype == MBPOLYGON ||
             mbtype == MBPOLYHEDRON ||
             mbtype == MBENTITYSET ||
             numnode == other.numnode);
  }
  bool operator<( const elemtype& other ) const
  {
    if (mbtype > other.mbtype)
      return false;
   
    return mbtype < other.mbtype ||
           (mbtype != MBPOLYGON &&
            mbtype != MBPOLYHEDRON &&
            mbtype != MBENTITYSET &&
            numnode < other.numnode);
  }
  bool operator!=( const elemtype& other ) const
    { return !this->operator==(other); }
};


MBErrorCode WriteHDF5Parallel::negotiate_type_list()
{
  int result;
  
  exportList.sort();
  
    // Get number of types each processor has
  int num_types = 2*exportList.size();
  std::vector<int> counts(numProc);
  result = MPI_Gather( &num_types, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Get list of types on this processor
  std::vector<int> my_types(num_types);
  std::vector<int>::iterator viter = my_types.begin();
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
       eiter != exportList.end(); ++eiter)
  {
    *viter = eiter->type;      ++viter;
    *viter = eiter->num_nodes; ++viter;
  }

  #ifdef DEBUG
  START_SERIAL;
  printdebug( "Local Element Types:\n");
  viter = my_types.begin();
  while (viter != my_types.end())
  {
    int type = *viter; ++viter;
    int count = *viter; ++viter;
    printdebug("  %s : %d\n", MBCN::EntityTypeName((MBEntityType)type), count);
  }
  END_SERIAL;
  #endif

    // Get list of types from each processor
  std::vector<int> displs(numProc + 1);
  displs[0] = 0;
  for (int i = 1; i <= numProc; ++i)
    displs[i] = displs[i-1] + counts[i-1];
  int total = displs[numProc];
  std::vector<int> alltypes(total);
  result = MPI_Gatherv( &my_types[0], my_types.size(), MPI_INT,
                        &alltypes[0], &counts[0], &displs[0], MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Merge type lists
  std::list<elemtype> type_list;
  std::list<elemtype>::iterator liter;
  for (int i = 0; i < numProc; ++i)
  {
    int* proc_type_list = &alltypes[displs[i]];
    liter = type_list.begin();
    for (int j = 0; j < counts[i]; j += 2)
    {
      elemtype type( &proc_type_list[j] );
        // skip until insertion spot
      for (; liter != type_list.end() && *liter < type; ++liter);
      
      if (liter == type_list.end() || *liter != type)
        liter = type_list.insert( liter, type );
    }
  }
  
    // Send total number of types to each processor
  total = type_list.size();
  result = MPI_Bcast( &total, 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Send list of types to each processor
  std::vector<int> intlist(total * 2);
  viter = intlist.begin();
  for (liter = type_list.begin(); liter != type_list.end(); ++liter)
  {
    *viter = liter->mbtype;  ++viter;
    *viter = liter->numnode; ++viter;
  }
  result = MPI_Bcast( &intlist[0], 2*total, MPI_INT, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

  #ifdef DEBUG
  START_SERIAL;
  printdebug( "Global Element Types:\n");
  viter = intlist.begin();
  while (viter != intlist.end())
  {
    int type = *viter; ++viter;
    int count = *viter; ++viter;
    printdebug("  %s : %d\n", MBCN::EntityTypeName((MBEntityType)type), count);
  }
  END_SERIAL;
  #endif
  
    // Insert missing types into exportList, with an empty
    // range of entities to export.
  std::list<ExportSet>::iterator ex_iter = exportList.begin();
  viter = intlist.begin();
  for (int i = 0; i < total; ++i)
  {
    int mbtype = *viter; ++viter;
    int numnode = *viter; ++viter;
    while (ex_iter != exportList.end() && ex_iter->type < mbtype)
      ++ex_iter;
    
    bool equal = ex_iter != exportList.end() && ex_iter->type == mbtype;
    if (equal && mbtype != MBPOLYGON && mbtype != MBPOLYHEDRON)
    {
      while (ex_iter != exportList.end() && ex_iter->num_nodes < numnode)
        ++ex_iter;
        
      equal = ex_iter != exportList.end() && ex_iter->num_nodes == numnode;
    }
    
    if (!equal)
    {
      ExportSet insert;
      insert.type = (MBEntityType)mbtype;
      insert.num_nodes = numnode;
      insert.first_id = 0;
      insert.offset = 0;
      insert.poly_offset = 0;
      insert.adj_offset = 0;
      ex_iter = exportList.insert( ex_iter, insert );
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::create_element_tables()
{
  int result;
  MBErrorCode rval;
  std::list<ExportSet>::iterator ex_iter;
  std::vector<long>::iterator viter;
  
    // Get number of each element type from each processor
  const int numtypes = exportList.size();
  std::vector<long> my_counts(numtypes);
  std::vector<long> counts(numtypes * numProc + numtypes);
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { *viter = ex_iter->range.size(); ++viter; }
  
  result = MPI_Gather( &my_counts[0], numtypes, MPI_LONG,
                       &counts[0],    numtypes, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Convert counts to offsets
  for (int i = 0; i < numtypes; i++) 
  {
    long prev = 0;
    for (int j = 0; j <= numProc; j++)
    {
      long tmp = counts[j*numtypes + i];
      counts[j*numtypes+i] = prev;
      prev += tmp;
    }
  }
  
    // Send offsets to each processor
  result = MPI_Scatter( &counts[0],    numtypes, MPI_LONG,
                        &my_counts[0], numtypes, MPI_LONG,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Update store offsets in ExportSets
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    ex_iter->offset = (id_t)*(viter++);
  
    // If polygons or polyhedra, send calculate offsets for each
  std::vector<int> perproc(numProc+1);
  ExportSet *poly[] = {0,0};
  int polycount[2];
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    if (ex_iter->type == MBPOLYGON)
    {
      assert(!poly[0]);
      poly[0] = &*ex_iter;
    }
    else if(ex_iter->type == MBPOLYHEDRON)
    {
      assert(!poly[1]);
      poly[1] = &*ex_iter;
    }
  }
  for (int i = 0; i < 2; i++)
  {
    ExportSet* ppoly = poly[i];
    if (!ppoly)
      continue;
  
    int count;
    rval = writeUtil->get_poly_array_size( ppoly->range.begin(),
                                           ppoly->range.end(),
                                           count );
    assert(MB_SUCCESS == rval);
    result = MPI_Gather( &count, 1, MPI_INT, &perproc[0], 1, MPI_INT, 0, MPI_COMM_WORLD );
    assert(MPI_SUCCESS == result);
    
    int prev = 0;
    for (int j = 1; j <= numProc; j++)
    {
      int tmp = perproc[j];
      perproc[j] = prev;
      prev += tmp;
    }
                                           
    polycount[i] = perproc[numProc];
    result = MPI_Scatter( &perproc[0], 1, MPI_INT, &count, 1, MPI_INT, 0, MPI_COMM_WORLD );
    assert(MPI_SUCCESS == result);
    ppoly->poly_offset = count;
  }
  
    // Create element tables
  std::vector<long> start_ids(numtypes);
  if (myRank == 0)
  {
    viter = start_ids.begin();
    long* citer = &counts[numtypes * numProc];
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    {
      switch(ex_iter->type) {
      case MBPOLYGON:
        rval = create_poly_tables( MBPOLYGON,
                                   *citer,
                                   polycount[0],
                                   *viter );
        break;
      case MBPOLYHEDRON:
        rval = create_poly_tables( MBPOLYHEDRON,
                                   *citer,
                                   polycount[1],
                                   *viter );
        break;
      default:
        rval = create_elem_tables( ex_iter->type,
                                   ex_iter->num_nodes,
                                   *citer,
                                   *viter );
      }
      assert(MB_SUCCESS == rval);
      ++citer;
      ++viter;
    }
  }
  
    // send start IDs to each processor
  result = MPI_Bcast( &start_ids[0], numtypes, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Assign IDs to local elements
  viter = start_ids.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    ex_iter->first_id = *(viter++);
    id_t myfirst = (id_t)(ex_iter->first_id + ex_iter->offset);
    rval = writeUtil->assign_ids( ex_iter->range, idTag, myfirst );
    assert(MB_SUCCESS == rval);
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode WriteHDF5Parallel::create_adjacency_tables()
{
  MBErrorCode rval;
  mhdf_Status status;
  int i, j, result;
  const int numtypes = exportList.size();
  std::vector<long>::iterator viter;
  std::list<ExportSet>::iterator ex_iter;
  std::vector<long> local(numtypes), all(numProc * numtypes + numtypes);
  
    // Get adjacency counts for local processor
  viter = local.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    id_t num_adj;
    rval = count_adjacencies( ex_iter->range, num_adj );
    assert (MB_SUCCESS == rval);
    *viter = num_adj; ++viter;
  }
  
    // Send local adjacency counts to root processor
  result = MPI_Gather( &local[0], numtypes, MPI_LONG,
                       &all[0],   numtypes, MPI_LONG, 
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Convert counts to offsets
  for (i = 0; i < numtypes; i++) 
  {
    long prev = 0;
    for (j = 0; j <= numProc; j++)
    {
      long tmp = all[j*numtypes + i];
      all[j*numtypes+i] = prev;
      prev += tmp;
    }
  }
  
    // For each element type for which there is no adjacency data,
    // send -1 to all processors as the offset
  for (i = 0; i < numtypes; ++i)
    if (all[numtypes*numProc+i] == 0)
      for (j = 0; j < numProc; ++j)
        all[j*numtypes+i] = -1;
  
    // Send offsets back to each processor
  result = MPI_Scatter( &all[0],   numtypes, MPI_LONG,
                        &local[0], numtypes, MPI_LONG,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Record the adjacency offset in each ExportSet
  viter = local.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { ex_iter->adj_offset = *viter; ++viter; }
  
    // Create data tables in file
  if (myRank == 0)
  {
    viter = all.begin() + (numtypes * numProc);
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter, ++viter)
    {
      if (!*viter) 
        continue;
      
      hid_t handle = mhdf_createAdjacency( filePtr,
                                           ex_iter->name(),
                                           *viter,
                                           &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        return MB_FAILURE;
      }
      mhdf_closeData( filePtr, handle, &status );
    }
  }

  return MB_SUCCESS;
}


struct RemoteSetData {
  MBTag handle;
  MBRange range;
  std::vector<int> counts, displs, all_values, local_values;
};

MBErrorCode WriteHDF5Parallel::get_remote_set_data( const char* tagname,
                                                    RemoteSetData& data,
                                                    long& offset )
{
  MBErrorCode rval;
  int i, result;
  MBRange::iterator riter;
  
  rval = iFace->tag_get_handle( tagname, data.handle );
  if (rval != MB_SUCCESS && rval != MB_TAG_NOT_FOUND)
    return rval;

  printdebug("Negotiating multi-proc meshsets for tag: \"%s\"\n", tagname);

    // Get sets with tag, or leave range empty if the tag
    // isn't defined on this processor.
  if (rval != MB_TAG_NOT_FOUND)
  {
    rval = iFace->get_entities_by_type_and_tag( 0, 
                                                MBENTITYSET, 
                                                &data.handle,
                                                0,
                                                1,
                                                data.range );
    if (rval != MB_SUCCESS) return rval;
    data.range = data.range.intersect( setSet.range );
    range_remove( setSet.range, data.range );
    #ifdef DEBUG
    setSet.range.sanity_check();
    #endif
  }
  
  printdebug("Found %d meshsets with \"%s\" tag.\n", data.range.size(), tagname );

    // Exchange number of sets with tag between all processors
  data.counts.resize(numProc);
  int count = data.range.size();
  result = MPI_Allgather( &count,          1, MPI_INT, 
                          &data.counts[0], 1, MPI_INT,
                          MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

    // Exchange tag values for sets between all processors
  data.displs.resize(numProc+1);
  data.displs[0] = 0;
  for (i = 1; i <= numProc; i++)
    data.displs[i] = data.displs[i-1] + data.counts[i-1];
  int total = data.displs[numProc];
  data.all_values.resize(total);
  data.local_values.resize(count);
  std::vector<int>::iterator viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
    rval = iFace->tag_get_data( data.handle, &*riter, 1, &*viter ); ++viter;
    assert(MB_SUCCESS == rval);
  }
  result = MPI_Allgatherv( &data.local_values[0], count, MPI_INT,
                           &data.all_values[0], &data.counts[0], &data.displs[0], MPI_INT,
                           MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

    // Find sets that span multple processors and update appropriately.
    // The first processor (sorted by MPI rank) that contains a given set
    // will be responsible for writing the set description.  All multi-
    // processor sets will be written at the beginning of the set tables.
    // Processors will write set contents/children for a given set in
    // the order of their MPI rank.
    //
    // Identify which meshsets will be managed by this processor and
    // the corresponding offset in the set description table. 
  std::set<int> tag_values;
  for (i = 0; i < total; ++i)
  {
    int id = 0;
    if (tag_values.insert(data.all_values[i]).second)
    {
      id = (int)++offset;
    }

    const unsigned int values_offset = (unsigned)(i - data.displs[myRank]);
    if (values_offset < (unsigned)count)
    {
      riter = data.range.begin();
      riter += values_offset;
      rval = iFace->tag_set_data( idTag, &*riter, 1, &id );
      assert(MB_SUCCESS == rval);
      //if (!id) allRemoteMesh.insert(*riter);
    }
  }
  
  return MB_SUCCESS;
}
  


MBErrorCode WriteHDF5Parallel::create_meshset_tables()
{
  MBErrorCode rval;
  int result, i;
  long total_offset = 0;
  MBRange::const_iterator riter;

  START_SERIAL;
  print_type_sets( iFace, myRank, numProc, setSet.range );
  END_SERIAL;

    // Gather data about multi-processor meshsets - removes sets from setSet.range
  std::vector<RemoteSetData> remote_set_data( multiProcSetTags.size() );
  for (i = 0; i< (int)multiProcSetTags.size(); i++)
  {
    rval = get_remote_set_data( multiProcSetTags[i].c_str(),
                                remote_set_data[i],
                                total_offset ); assert(MB_SUCCESS == rval);
  }

  START_SERIAL;
  print_type_sets( iFace, myRank, numProc, setSet.range );
  END_SERIAL;

    // Gather counts from each proc
  std::vector<long> set_offsets(numProc + 1);
  long local_count = setSet.range.size();
  result = MPI_Gather( &local_count,    1, MPI_LONG,
                       &set_offsets[0], 1, MPI_LONG,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  for (i = 0; i <= numProc; i++)
  {
    long tmp = set_offsets[i];
    set_offsets[i] = total_offset;
    total_offset += tmp;
  }
  
    // Send each proc its offsets
  long sets_offset;
  result = MPI_Scatter( &set_offsets[0], 1, MPI_LONG,
                        &sets_offset,    1, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setSet.offset = (id_t)(sets_offset);

    // Create the table
  long total_count_and_start_id[2] = { set_offsets[numProc], 0 };
  if (myRank == 0 && total_count_and_start_id[0] > 0)
  {
    rval = create_set_meta( (id_t)total_count_and_start_id[0], total_count_and_start_id[1] );
    assert (MB_SUCCESS == rval);
  }
  
    // Send totals to all procs.
  result = MPI_Bcast( total_count_and_start_id, 2, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setSet.first_id = total_count_and_start_id[1];
  writeSets = total_count_and_start_id[0] > 0;

  START_SERIAL;  
  printdebug("Non-shared sets: %ld local, %ld global, offset = %ld, first_id = %ld\n",
    local_count, total_count_and_start_id[0], sets_offset, total_count_and_start_id[1] );
  END_SERIAL;
  
    // Not writing any sets??
  if (!writeSets)
    return MB_SUCCESS;
  
    // Assign set IDs
  writeUtil->assign_ids( setSet.range, idTag, (id_t)(setSet.first_id + setSet.offset) );
  for (i = 0; i < (int)remote_set_data.size(); ++i)
    fix_remote_set_ids( remote_set_data[i], setSet.first_id );
  
  
    // Communicate sizes for remote sets
  long data_offsets[2] = { 0, 0 };
  for (i = 0; i < (int)remote_set_data.size(); ++i)
  {
    rval = negotiate_remote_set_contents( remote_set_data[i], data_offsets ); 
    assert(MB_SUCCESS == rval);
  }
  remote_set_data.clear();
  
    // Exchange IDs for remote/adjacent sets not shared between procs
  //rval = communicate_remote_ids( MBENTITYSET ); assert(MB_SUCCESS == rval);
  
    // Communicate counts for local sets
  long data_counts[2];
  rval = count_set_size( setSet.range, rangeSets, data_counts[0], data_counts[1] );
  if (MB_SUCCESS != rval) return rval;
  std::vector<long> set_counts(2*numProc);
  result = MPI_Gather( data_counts,    2, MPI_LONG,
                       &set_counts[0], 2, MPI_LONG,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  for (i = 0; i < 2*numProc; ++i)
  {
    long tmp = set_counts[i];
    set_counts[i] = data_offsets[i%2];
    data_offsets[i%2] += tmp;
  }
  long all_counts[] = {data_offsets[0], data_offsets[1]};
  result = MPI_Scatter( &set_counts[0], 2, MPI_LONG,
                        data_offsets,   2, MPI_LONG,
                        0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setContentsOffset = data_offsets[0];
  setChildrenOffset = data_offsets[1];
  
    // Create set contents and set children tables
  if (myRank == 0)
  {
    rval = create_set_tables( all_counts[0], all_counts[1] );
    if (MB_SUCCESS != rval) return rval;
  }
  
    // Send totals to all processors
  result = MPI_Bcast( all_counts, 2, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  writeSetContents = all_counts[0] > 0;
  writeSetChildren = all_counts[1] > 0;

  START_SERIAL;  
  printdebug("Non-shared set contents: %ld local, %ld global, offset = %ld\n",
    data_counts[0], all_counts[0], data_offsets[0] );
  printdebug("Non-shared set children: %ld local, %ld global, offset = %ld\n",
    data_counts[1], all_counts[1], data_offsets[1] );
  END_SERIAL;
  
  return MB_SUCCESS;
}

void WriteHDF5Parallel::remove_remote_entities( MBRange& range )
{
  MBRange result;
  result.merge( range.intersect( nodeSet.range ) );
  result.merge( range.intersect( setSet.range ) );  
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
  {
    result.merge( range.intersect( eiter->range ) );
  }
  range = result;
}

void WriteHDF5Parallel::remove_remote_entities( std::vector<MBEntityHandle>& vect )
{
  MBRange intrsct;
  for (std::vector<MBEntityHandle>::const_iterator iter = vect.begin();
       iter != vect.end(); ++iter)
    intrsct.insert(*iter);
  remove_remote_entities( intrsct );
  
  unsigned int read, write;
  for (read = write = 0; read < vect.size(); ++read)
  {
    if (intrsct.find(vect[read]) != intrsct.end())
    {
      if (read != write)
        vect[write] = vect[read];
      ++write;
    }
  }
  if (write != vect.size())
    vect.resize(write);
}



MBErrorCode WriteHDF5Parallel::negotiate_remote_set_contents( RemoteSetData& data,
                                                              long* offsets /* long[2] */ )
{
  unsigned i;
  MBErrorCode rval;
  MBRange::const_iterator riter;
  int result;
  const unsigned count = data.range.size();
  const unsigned total = data.all_values.size();
  std::vector<int>::iterator viter;

    // Calculate counts for each meshset
  std::vector<long> local_sizes(2*count);
  std::vector<long>::iterator sizes_iter = local_sizes.begin();
  MBRange tmp_range;
  std::vector<MBEntityHandle> child_list;
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
      // Count contents
    *sizes_iter = 0;
    tmp_range.clear();
    rval = iFace->get_entities_by_handle( *riter, tmp_range );
    remove_remote_entities( tmp_range );
    assert (MB_SUCCESS == rval);
    for (MBRange::iterator iter = tmp_range.begin(); iter != tmp_range.end(); ++iter)
    {
      int id = 0;
      rval = iFace->tag_get_data( idTag, &*iter, 1, &id );
      if (rval != MB_TAG_NOT_FOUND && rval != MB_SUCCESS)
        { assert(0); return MB_FAILURE; }
      if (id > 0)
        ++*sizes_iter;
    }
    ++sizes_iter;
    
      // Count children
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_child_meshsets( *riter, child_list );
    remove_remote_entities( child_list );
    assert (MB_SUCCESS == rval);
    for (std::vector<MBEntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
    {
      int id = 0;
      rval = iFace->tag_get_data( idTag, &*iter, 1, &id );
      if (rval != MB_TAG_NOT_FOUND && rval != MB_SUCCESS)
        { assert(0); return MB_FAILURE; }
      if (id > 0)
        ++*sizes_iter;
    }
    ++sizes_iter;
  }
  
    // Exchange sizes for sets between all processors.
  std::vector<long> all_sizes(2*total);
  std::vector<int> counts(numProc), displs(numProc);
  for (i = 0; i < (unsigned)numProc; i++)
    counts[i] = 2 * data.counts[i];
  displs[0] = 0;
  for (i = 1; i < (unsigned)numProc; i++)
    displs[i] = displs[i-1] + counts[i-1];
  result = MPI_Allgatherv( &local_sizes[0], 2*count, MPI_LONG,
                           &all_sizes[0], &counts[0], &displs[0], MPI_LONG,
                           MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);

  
    // Update information in-place in the array from the Allgatherv.
    
    // Change the corresponding sizes for the first instance of a tag
    // value such that it ends up being the total size of the set.
    // Change the size to -1 for the later instances of a tag value.
    //
    // For the sets that this processor has, update the offsets at
    // which the set data is to be written.  Store the offset of the data
    // on this processor for the set *relative* to the start of the
    // data of *the set*.
  std::vector<long> local_offsets(2*count);
  std::map<int,int> tagsort;  // Map of {tag value, index of first set w/ value}
  for (i = 0; i < total; ++i)
  {
    const std::map<int,int>::iterator p = tagsort.find( data.all_values[i] );
    const unsigned r = (unsigned)(i - data.displs[myRank]);  // offset in "local" array
    
      // If this is the first instance of this tag value, 
      // then the processor with this instance is responsible
      // for writing the tag description
    if ( p == tagsort.end() )  
    {
      tagsort.insert( std::make_pair(data.all_values[i], i) );
        // If within the range for this processor, save offsets
      if (r < (unsigned)count) 
      {
        local_offsets[2*r] = local_offsets[2*r+1] = 0;
      }
    }
      // Otherwise update the total size in the table
      // for the processor that is responsible for writing
      // the data and mark the data for the current processor
      // with a -1.
    else 
    {
        // If within the range for this processor, save offsets
      int j = p->second;
      if (r < (unsigned)count) 
      {
        local_offsets[2*r  ] = all_sizes[2*j  ];
        local_offsets[2*r+1] = all_sizes[2*j+1];
      }
      
      all_sizes[2*j  ] += all_sizes[2*i  ];
      all_sizes[2*j+1] += all_sizes[2*i+1];
      all_sizes[2*i  ]  = all_sizes[2*i+1] = -1;
    }
  }  

  
    // Store the total size of each set (rather than the
    // number of entities local to this processor) in the
    // local_sizes array for each meshset.  Only need this
    // for the sets this processor is writing the description
    // for, but it's easier to get it for all of them.
  sizes_iter = local_sizes.begin();
  viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter ); ++viter;
    assert( p != tagsort.end() );
    int j = 2 * p->second;
    *sizes_iter = all_sizes[j  ]; ++sizes_iter;
    *sizes_iter = all_sizes[j+1]; ++sizes_iter;
  }
  
    // Now calculate the offset of the data for each (entire, parallel) set in
    // the set contents and set children tables.
  for (i = 0; i < all_sizes.size(); ++i)
  {
    if (all_sizes[i] >= 0)
    {
      int j = i % 2;              // contents or children list ?
      long tmp = offsets[j];      // save current, running offset
      offsets[j] += all_sizes[i]; // next set's offset is current plus the size of this set
      all_sizes[i] = tmp;         // size of this set is running offset.
    }
  }
  
    // Local offsets for this processor are stored as values relative to the
    // start of each set's data.  Convert them to offsets relative to the
    // start of all the set data.
  sizes_iter = local_offsets.begin();
  viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter ); ++viter;
    assert( p != tagsort.end() );
    int j = 2 * p->second;
    *sizes_iter += all_sizes[j  ]; ++sizes_iter;
    *sizes_iter += all_sizes[j+1]; ++sizes_iter;
  }

#ifdef DEBUG  
START_SERIAL; if (counts[myRank]) {
std::string name;
iFace->tag_get_name( data.handle, name );
printdebug("Remote set data for tag \"%s\"\n", name.c_str() );
printdebug("    tag value         owner local offsets   total counts\n");
for (unsigned d = 0; d < (unsigned)counts[myRank]; ++d) {
if (0==(d%2)) {
printdebug("%13d %13s %13d %13d\n", data.all_values[(d+displs[myRank])/2], all_sizes[d+displs[myRank]] < 0 ? "no" : "yes", local_offsets[d], local_sizes[d] );
} else {
printdebug("                            %13d %13d\n", local_offsets[d], local_sizes[d] );
} 
}
} END_SERIAL;
#endif
  
    // Store each parallel meshset in the list
  sizes_iter = local_sizes.begin();
  std::vector<long>::iterator offset_iter = local_offsets.begin();
  std::vector<long>::iterator all_iter = all_sizes.begin() + displs[myRank];
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
    ParallelSet info;
    info.handle = *riter;
    info.contentsOffset = *offset_iter; ++offset_iter;
    info.childrenOffset = *offset_iter; ++offset_iter;
    info.contentsCount = *sizes_iter; ++sizes_iter;
    info.childrenCount = *sizes_iter; ++sizes_iter;
    info.description = *all_iter >= 0; all_iter += 2;
    parallelSets.push_back( info );
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::fix_remote_set_ids( RemoteSetData& data, long first_id )
{
  const id_t id_diff = (id_t)(first_id - 1);
  id_t file_id;
  MBErrorCode rval;

  for (MBRange::iterator iter = data.range.begin(); iter != data.range.end(); ++iter)
  {
    rval = iFace->tag_get_data( idTag, &*iter, 1, &file_id );
    assert( MB_SUCCESS == rval );
    file_id += id_diff;
    rval = iFace->tag_set_data( idTag, &*iter, 1, &file_id );
    assert( MB_SUCCESS == rval );
  }
  
  return MB_SUCCESS;
}   


MBErrorCode WriteHDF5Parallel::write_shared_set_descriptions( hid_t table )
{
  const id_t start_id = setSet.first_id;
  MBErrorCode rval;
  mhdf_Status status;
  
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    if (!iter->description)
      continue;  // handled by a different processor
    
      // Get offset in table at which to write data
    int file_id;
    rval = iFace->tag_get_data( idTag, &(iter->handle), 1, &file_id );
    file_id -= start_id;
    
      // Get flag data
    unsigned int flags;
    rval = iFace->get_meshset_options( iter->handle, flags );
    assert( MB_SUCCESS == rval );
      
      // Write the data
    long data[3] = { iter->contentsOffset + iter->contentsCount - 1, 
                     iter->childrenOffset + iter->childrenCount - 1, 
                     flags };
    mhdf_writeSetMeta( table, file_id, 1, H5T_NATIVE_LONG, data, &status );
    if (mhdf_isError(&status))
      printdebug("Meshset %d : %s\n", ID_FROM_HANDLE(iter->handle), mhdf_message(&status));
    assert( !mhdf_isError( &status ) );
  }

  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_contents( hid_t table )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_entities_by_handle( iter->handle, handle_list );
    assert( MB_SUCCESS == rval );
    remove_remote_entities( handle_list );
    
    id_list.clear();
    for (unsigned int i = 0; i < handle_list.size(); ++i)
    {
      int id;
      rval = iFace->tag_get_data( idTag, &handle_list[i], 1, &id );
      assert( MB_SUCCESS == rval );
      if (id > 0)
        id_list.push_back(id);
    }
    
    if (id_list.empty())
      continue;
    
    mhdf_writeSetData( table, 
                       iter->contentsOffset, 
                       id_list.size(),
                       id_type,
                       &id_list[0],
                       &status );
    assert(!mhdf_isError(&status));
  }
  
  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_children( hid_t table )
{
  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  printdebug("Writing %d parallel sets.\n", parallelSets.size());
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_child_meshsets( iter->handle, handle_list );
    assert( MB_SUCCESS == rval );
    remove_remote_entities( handle_list );
    
    id_list.clear();
    for (unsigned int i = 0; i < handle_list.size(); ++i)
    {
      int id;
      rval = iFace->tag_get_data( idTag, &handle_list[i], 1, &id );
      assert( MB_SUCCESS == rval );
      if (id > 0)
        id_list.push_back(id);
    }
    
    if (!id_list.empty())
    {
      mhdf_writeSetChildren( table, 
                             iter->childrenOffset, 
                             id_list.size(),
                             id_type,
                             &id_list[0],
                             &status );
      assert(!mhdf_isError(&status));
    }
  }

  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::write_finished()
{
  parallelSets.clear();
  return WriteHDF5::write_finished();
}


class TagNameCompare {
  MBInterface* iFace;
  std::string name1, name2;
public:
  TagNameCompare( MBInterface* iface ) : iFace(iface) {}
  bool operator() (const WriteHDF5::SparseTag& t1, 
                   const WriteHDF5::SparseTag& t2);
};
bool TagNameCompare::operator() (const WriteHDF5::SparseTag& t1, 
                                 const WriteHDF5::SparseTag& t2)
{
  MBErrorCode rval;
  rval = iFace->tag_get_name( t1.tag_id, name1 );
  rval = iFace->tag_get_name( t2.tag_id, name2 );
  return name1 < name2;
}  

void WriteHDF5Parallel::sort_tags_by_name( )
{
  tagList.sort( TagNameCompare( iFace ) );
}


MBErrorCode WriteHDF5Parallel::communicate_remote_ids( MBEntityType type )
{
  int result;
  MBErrorCode rval;

    // Get the export set for the specified type
  ExportSet* export_set = 0;
  if (type == MBVERTEX)
    export_set = &nodeSet;
  else if(type == MBENTITYSET)
    export_set = &setSet;
  else
  {
    for (std::list<ExportSet>::iterator esiter = exportList.begin();
         esiter != exportList.end(); ++esiter)
      if (esiter->type == type)
      {
        export_set = &*esiter;
        break;
      }
  }
  assert(export_set != NULL);
  
    // Get the ranges in the set
  std::vector<unsigned long> myranges;
  MBRange::const_pair_iterator p_iter( export_set->range.pair_begin() );
  const MBRange::const_pair_iterator p_end( export_set->range.pair_end() );
  for ( ; p_iter != p_end; ++p_iter)
  {
    myranges.push_back( (*p_iter).first );
    myranges.push_back( (*p_iter).second );
  }

  START_SERIAL;
  printdebug("%s ranges to communicate:\n", MBCN::EntityTypeName(type));
  for (unsigned int xx = 0; xx != myranges.size(); xx+=2)
    printdebug("  %lu - %lu\n", myranges[xx], myranges[xx+1] );
  END_SERIAL;
  
    // Communicate the number of ranges and the start_id for
    // each processor.
  std::vector<int> counts(numProc), offsets(numProc), displs(numProc);
  int mycount = myranges.size();
  int mystart = export_set->first_id + export_set->offset;
  result = MPI_Allgather( &mycount, 1, MPI_INT, &counts[0], 1, MPI_INT, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  result = MPI_Allgather( &mystart, 1, MPI_INT, &offsets[0], 1, MPI_INT, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Communicate the ranges 
  displs[0] = 0;
  for (int i = 1; i < numProc; ++i)
    displs[i] = displs[i-1] + counts[i-1];
  std::vector<unsigned long> allranges( displs[numProc-1] + counts[numProc-1] );
  result = MPI_Allgatherv( &myranges[0], myranges.size(), MPI_UNSIGNED_LONG,
                           &allranges[0], &counts[0], &displs[0],
                           MPI_UNSIGNED_LONG, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
  MBTag global_id_tag;
  iFace->tag_get_handle( PARALLEL_GLOBAL_ID_TAG_NAME, global_id_tag );
  
    // Set file IDs for each communicated entity
    
    // For each processor
  for (int proc = 0; proc < numProc; ++proc)
  {
    if (proc == myRank)
      continue;
    
      // Get data for corresponding processor
    const int offset = offsets[proc];
    const int count = counts[proc];
    const unsigned long* const ranges = &allranges[displs[proc]];
    
      // For each geometry meshset in the interface
    MBRange::iterator r_iter = MBRange::lower_bound( remoteMesh[proc].begin(),
                                                     remoteMesh[proc].end(),
                                                     CREATE_HANDLE(type,0,result) );
    MBRange::iterator r_stop = MBRange::lower_bound( r_iter,
                                                     remoteMesh[proc].end(),
                                                     CREATE_HANDLE(type+1,0,result) );
    for ( ; r_iter != r_stop; ++r_iter)
    {
      MBEntityHandle entity = *r_iter;

        // Get handle on other processor
      MBEntityHandle global;
      rval = iFace->tag_get_data( global_id_tag, &entity, 1, &global );
      assert(MB_SUCCESS == rval);

        // Find corresponding fileid on other processor.
        // This could potentially be n**2, but we will assume that
        // the range list from each processor is short (typically 1).
      int j, steps = 0;
      unsigned long low, high;
      for (j = 0; j < count; j += 2)
      {
        low = ranges[j];
        high = ranges[j+1];
        if (low <= global && high >= global)
          break;
        steps += (high - low) + 1;
      }
      if (j >= count) {
      printdebug("*** handle = %u, type = %d, id = %d, proc = %d\n",
      (unsigned)global, (int)(iFace->type_from_handle(global)), (int)(iFace->id_from_handle(global)), proc);
      for (int ii = 0; ii < count; ii+=2) 
      printdebug("***  %u to %u\n", (unsigned)ranges[ii], (unsigned)ranges[ii+1] );
      }
      assert(j < count);
      int fileid = offset + steps + (global - low);
      rval = iFace->tag_set_data( idTag, &entity, 1, &fileid );
      assert(MB_SUCCESS == rval);
    } // for(r_iter->range)
  } // for(each processor)
  
  return MB_SUCCESS;
}
