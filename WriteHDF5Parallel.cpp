
#define DEBUG

#include <assert.h>
#ifdef DEBUG
#  include <stdio.h>
#  include <stdarg.h>
#endif

#include <vector>
#include <set>
#include <map>
#include <utility>

#include <mpi.h>

#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include <H5FDmpio.h>

#include "mhdf.h"

#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"

#include "WriteHDF5Parallel.hpp"

void WriteHDF5Parallel::printdebug( const char* fmt, ... )
{
#ifdef DEBUG
  fprintf( stderr, "[%d] ", myRank );
  va_list args;
  va_start( args, fmt );
  vfprintf( stderr, fmt, args );
  va_end( args );
  fflush( stderr );
#endif
}

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


  // These tags are used to describe meshsets that span
  // multiple processors during a parallel write.  
  //
  // The presence of the HDF5_PARALLE_SET_TAG on a meshset
  // indicates that the set spans multiple processors.  The
  // data of the tag is a pair of unsigned long values
  // indicating the offset into the contents and children
  // tables, respectively, at which this processor is to 
  // write it's members in the set.
  //
  // The presence of the HDF5_PARALLEL_SET_META_TAG on a
  // meshset indicates that this processor is responsible
  // for writing the set description.  The data for the tag
  // is the tripple of values (stored as unsigned long) 
  // that are to be written to the set description table
  // ({last contents index, last child index, flags}).  The
  // implementation may choose to create this tag for all
  // parallel meshsets on the root processor.  The data is
  // stored in a tag because a) it needs to be stored somewhere
  // before the set is written and b) the normal code as 
  // written in this file need to be aware of parallel details
  // such as whether or not it's running on the root processor.
#define HDF5_PARALLEL_SET_TAG "h5set_offset"
#define HDF5_PARALLEL_SET_META_TAG "h5set_size"

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
  int proc_pair[2];
  
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
    entities.subtract( sets );    
    remoteMesh[remote_proc].merge( entities );
    remoteMesh[remote_proc].insert( *iiter );
    
    for (MBRange::iterator siter = sets.begin(); siter != sets.end(); ++siter)
    {
        // For current parallel meshing code, root processor owns
        // all curve and geometric vertex meshes.  
      int dimension;
      result = iFace->tag_get_data( geom_tag, &*siter, 1, &dimension );
      const int tmp_remote_proc = 
        (result == MB_SUCCESS && dimension < 2) ? 0 : remote_proc;
      
        // Put entities in list for appropriate processor.
      remoteMesh[tmp_remote_proc].insert( *siter );
      entities.clear();
      result = iFace->get_entities_by_handle( *siter, entities );
      if (MB_SUCCESS != result) return result;
      remoteMesh[tmp_remote_proc].merge( entities );
    }
  }
  
    // For all remote mesh entities, remove them from the
    // lists of local mesh to be exported and give them a 
    // junk file Id of 1.  Need to specify a file ID greater
    // than zero so the code that gathers adjacencies and 
    // such doesn't think that the entities aren't being
    // exported.
  for (int i = 0; i < numProc; i++)
  {
    if (i == myRank) continue;
    
    MBRange& range = remoteMesh[i];
    
    nodeSet.range.subtract( range );
    setSet.range.subtract( range );
    for (std::list<ExportSet>::iterator eiter = exportList.begin();
         eiter != exportList.end(); ++eiter )
      eiter->range.subtract( range );
    
    int id = 1;
    for (MBRange::iterator riter = remoteMesh[i].begin(); 
         riter != remoteMesh[i].end(); ++riter)
    {
      result = iFace->tag_set_data( idTag, &*riter, 1, &id );
      if (MB_SUCCESS != result) return result;
    }
  }
      
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
  
  rval = gather_tags();
  if (MB_SUCCESS != rval)
    return rval;
  
    /**************** Create actual file and write meta info ***************/

  if (myRank == 0)
  {
      // create the file
    const char* type_names[MBMAXTYPE];
    bzero( type_names, MBMAXTYPE * sizeof(char*) );
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
   
  
    /**************** Create tag data *********************/
  
  std::list<SparseTag>::iterator tag_iter;
  sort_tags_by_name();
  const int num_tags = tagList.size();
  int* tag_offsets = new int[num_tags];
  int* tag_off_iter = tag_offsets;
  for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter, ++tag_off_iter)
    *tag_off_iter = tag_iter->range.size();
  
  int* proc_tag_offsets = 0;
  if (myRank == 0)
    proc_tag_offsets = new int[num_tags*numProc];
  result = MPI_Gather( tag_offsets, num_tags, MPI_INT,
                  proc_tag_offsets, num_tags, MPI_INT,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
  if (myRank == 0)
  {
    tag_iter = tagList.begin();
    for (int i = 0; i < num_tags; ++i, ++tag_iter)
    {
      int next_offset = 0;
      for (int j = 0; j < numProc; j++)
      {
        int count = proc_tag_offsets[i*numProc + j];
        proc_tag_offsets[i*numProc+j] = next_offset;
        next_offset += count;
      }
      
      rval = create_tag( tag_iter->tag_id, next_offset );
      assert(MB_SUCCESS == rval);
    }
  }
  
  result = MPI_Scatter( proc_tag_offsets, num_tags, MPI_INT,
                             tag_offsets, num_tags, MPI_INT,
                             0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  delete [] proc_tag_offsets;
  
  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
    tag_iter->offset = tag_offsets[i];
  
  
  /************** Close serial file and reopen parallel *****************/
  
  if (0 == myRank)
  {
    delete [] tag_offsets;
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
  int result;
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
  
    // End local adjacency counts to root processor
  result = MPI_Gather( &local[0], numtypes, MPI_LONG,
                       &all[0],   numtypes, MPI_LONG, 
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Convert counts to offsets
  for (int i = 0; i < numtypes; i++) 
  {
    long prev = 0;
    for (int j = 0; j <= numProc; j++)
    {
      long tmp = all[j*numtypes + i];
      all[j*numtypes+i] = prev;
      prev += tmp;
    }
  }
  
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

MBErrorCode WriteHDF5Parallel::create_meshset_tables()
{
  MBErrorCode rval;
  int result, id = 1;
  long prev[3];
    
    // Give all sets to be exported an ID of 1 as a mark
    // so we know which connected sets are to be exported,
    // etc.  Will set it to the actual file ID later.
  for (MBRange::iterator riter = setSet.range.begin(); 
       riter != setSet.range.end(); ++riter)
  {
    rval = iFace->tag_set_data( idTag, &*riter, 1, &id );
    if (MB_SUCCESS != rval) return rval;
  }
  
  rval = negotiate_shared_meshsets(prev);
  if (MB_SUCCESS != rval) return rval;

    // Gather counts from each proc
  std::vector<long> set_offsets(3*numProc+3);
  long counts[3];
  counts[0] = setSet.range.size();
  rval = count_set_size( setSet.range, rangeSets, counts[1], counts[2] );
  if (MB_SUCCESS != rval) return rval;

printdebug( "Local set counts: { %ld %ld %ld }\n", counts[0], counts[1], counts[2] );
  
  result = MPI_Gather( counts,          3, MPI_LONG, 
                       &set_offsets[0], 3, MPI_LONG,
                       0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Convert counts to offsets
  for (int i = 0; i <= numProc; i++)
  {
    long* curr = &set_offsets[3*i];
    for (int j = 0; j < 3; j++)
    {
      int tmp = curr[j];
      curr[j] = prev[j];
      prev[j] += tmp;
    }
  }
  
    // Send each proc its offsets
  result = MPI_Scatter( &set_offsets[0], 3, MPI_LONG,
                        counts,          3, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setSet.offset = (id_t)(counts[0]);
  setContentsOffset = counts[1];
  setChildrenOffset = counts[2];

printdebug( "Local set offsets: { %ld %ld %ld }\n", counts[0], counts[1], counts[2] );
  
    // Create the tables
    
    // The last three entries in the vector contain
    // the total table sizes.  Reuse the previous entry
    // to hold the start Id.
  long* totals = &set_offsets[3*numProc-1];
  if (myRank == 0 && totals[1])
  {
printdebug( "Total set counts: { %ld %ld %ld }\n", totals[1], totals[2], totals[3] );
    rval = create_set_tables( totals[1], totals[2], totals[3], totals[0] );
    if (MB_SUCCESS != rval) return rval;
  }
  
    // Send totals to all processors
  result = MPI_Bcast( totals, 4, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  setSet.first_id  = totals[0];
  writeSets        = totals[1] > 0;
  writeSetContents = totals[2] > 0;
  writeSetChildren = totals[3] > 0;
  id_t myfirst = (id_t)(setSet.first_id + setSet.offset);
  rval = writeUtil->assign_ids( setSet.range, idTag, myfirst );
  assign_parallel_set_ids();
  
  return rval;
}

MBErrorCode WriteHDF5Parallel::negotiate_shared_meshsets( long* offsets )
{
  MBErrorCode rval = MB_SUCCESS;
  offsets[0] =  0; // row at which to write first set description
  offsets[1] = -1; // last used row in set contents table
  offsets[2] = -1; // last used row in set children table
  for (std::vector<std::string>::iterator niter = multiProcSetTags.begin();
       niter != multiProcSetTags.end(); ++niter)
  {
    rval = negotiate_shared_meshsets( niter->c_str(), offsets );
    assert(MB_SUCCESS == rval);
  }
  offsets[1]++;
  offsets[2]++;
  
printdebug( "parallel set counts : { %ld %ld %ld }\n", offsets[0], offsets[1], offsets[2] );

  return rval;
}

MBErrorCode WriteHDF5Parallel::negotiate_shared_meshsets( const char* tagname,
                                                          long* offsets )
{
  int i;
  MBErrorCode rval;
  MBRange::const_iterator riter;
  int result;

  MBTag handle;
  rval = iFace->tag_get_handle( tagname, handle );
  if (rval != MB_SUCCESS && rval != MB_TAG_NOT_FOUND)
    return rval;

    // Get sets with tag, or leave range empty if the tag
    // isn't defined on this processor.
  MBRange range;
  if (rval != MB_TAG_NOT_FOUND)
  {
    rval = iFace->get_entities_by_type_and_tag( 0, 
                                                MBENTITYSET, 
                                                &handle,
                                                0,
                                                1,
                                                range );
    if (rval != MB_SUCCESS) return rval;
  }

    // Exchange number of sets with tag between all processors
  std::vector<int> counts(numProc);
  int count = range.size();
  result = MPI_Allgather( &count,     1, MPI_INT, 
                          &counts[0], 1, MPI_INT,
                          MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  
    // Exchange tag values for sets between all processors
  std::vector<int> displs(numProc+1);
  displs[0] = 0;
  for (i = 1; i <= numProc; i++)
    displs[i] = displs[i-1] + counts[i-1];
  int total = displs[numProc];
  std::vector<int> all_values(total);
  std::vector<int> local_values(count);
  std::vector<int>::iterator viter = local_values.begin();
  for (riter = range.begin(); riter != range.end(); ++riter)
  {
    rval = iFace->tag_get_data( handle, &*riter, 1, &*viter ); ++viter;
    assert(MB_SUCCESS == rval);
  }
  result = MPI_Allgatherv( &local_values[0], count, MPI_INT,
                           &all_values[0], &counts[0], &displs[0], MPI_INT,
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
  int meta_offset = (int)offsets[0];
  std::set<int> tag_values;
  riter = range.begin();
  for (i = 0; i < total; ++i)
  {
    const unsigned int values_offset = (unsigned)(i - displs[myRank]);
    if (tag_values.insert(all_values[i]).second)
    {
      meta_offset++;
      if (values_offset < (unsigned)count)
      {
        rval = iFace->tag_set_data( idTag, &*riter, 1, &meta_offset );
        assert(MB_SUCCESS == rval);
      }
    }
    else
    {
      const int zero = 0;
      if (values_offset < (unsigned)count)
      {
        rval = iFace->tag_set_data( idTag, &*riter, 1, &zero );
        assert(MB_SUCCESS == rval);
      }
    }
  }
  offsets[0] = meta_offset;
  
    // Calculate counts for each meshset
  std::vector<long> local_sizes(2*count);
  std::vector<long>::iterator sizes_iter = local_sizes.begin();
  MBRange tmp_range;
  std::vector<MBEntityHandle> child_list;
  for (riter = range.begin(); riter != range.end(); ++riter)
  {
      // Count contents
    *sizes_iter = 0;
    tmp_range.clear();
    rval = iFace->get_entities_by_handle( *riter, tmp_range );
    assert (MB_SUCCESS == rval);
    for (MBRange::iterator iter = tmp_range.begin(); iter != tmp_range.end(); ++iter)
    {
      int id = 0;
      rval = iFace->tag_get_data( idTag, &*iter, 1, &id );
      if (rval != MB_TAG_NOT_FOUND && rval != MB_SUCCESS)
        { assert(0); return MB_FAILURE; }
      if (id)
        ++*sizes_iter;
    }
    ++sizes_iter;
    
      // Count children
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_child_meshsets( *riter, child_list );
    assert (MB_SUCCESS == rval);
    for (std::vector<MBEntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
    {
      int id = 0;
      rval = iFace->tag_get_data( idTag, &*iter, 1, &id );
      if (rval != MB_TAG_NOT_FOUND && rval != MB_SUCCESS)
        { assert(0); return MB_FAILURE; }
      if (id)
        ++*sizes_iter;
    }
    ++sizes_iter;
  }
  
    // Exchange sizes for sets between all processors.
  std::vector<long> all_sizes(2*total);
  std::vector<int> counts2(numProc), displs2(numProc);
  for (i = 0; i < numProc; i++)
    counts2[i] = 2 * counts[i];
  displs2[0] = 0;
  for (i = 1; i < numProc; i++)
    displs2[i] = displs2[i-1] + counts[i-1];
  result = MPI_Allgatherv( &local_sizes[0], 2*count, MPI_LONG,
                           &all_sizes[0], &counts2[0], &displs2[0], MPI_LONG,
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
  std::vector<long> local_offsets;
  std::map<int,int> tagsort;  // Map of {tag value, index of first set w/ value}
  for (i = 0; i < total; ++i)
  {
    const std::map<int,int>::iterator p = tagsort.find( all_values[i] );
    const unsigned r = (unsigned)(i - displs[myRank]);  // offset in "local" array
    
      // If this is the first instance of this tag value, 
      // then the processor with this instance is responsible
      // for writing the tag description
    if ( p == tagsort.end() )  
    {
      tagsort.insert( std::make_pair(all_values[i], i) );
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
        local_offsets[r+1] = all_sizes[2*j  ];
        local_offsets[r+2] = all_sizes[2*j+1];
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
  viter = local_values.begin();
  for (riter = range.begin(); riter != range.end(); ++riter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter ); ++viter;
    assert( p != tagsort.end() );
    int j = 2 * p->second;
    *sizes_iter = all_sizes[j  ]; ++sizes_iter;
    *sizes_iter = all_sizes[j+1]; ++sizes_iter;
  }
  
    // Now calculate the offset of the data for each set in
    // the set contents and set children tables.
  for (std::map<int,int>::iterator map_iter = tagsort.begin();
       map_iter != tagsort.end(); ++map_iter)
  {
    int j = 2 * map_iter->second;
    all_sizes[j  ] = (offsets[1] += all_sizes[j  ]);
    all_sizes[j+1] = (offsets[2] += all_sizes[j+1]);
  }
  
    // Convert the offsets relative to the beginning of the
    // set data to absolute offsets in the table for the
    // data to be written by this processor
  sizes_iter = local_offsets.begin();
  viter = local_values.begin();
  for (riter = range.begin(); riter != range.end(); ++riter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter ); ++viter;
    assert( p != tagsort.end() );
    int j = 2 * p->second;
    *sizes_iter = all_sizes[j  ]; ++sizes_iter;
    *sizes_iter = all_sizes[j+1]; ++sizes_iter;
  }
  
    // Store each parallel meshset in the list
  sizes_iter = local_sizes.begin();
  std::vector<long>::iterator offset_iter = local_offsets.begin();
  std::vector<long>::iterator all_iter = all_sizes.begin() + displs2[myRank];
  for (riter = range.begin(); riter != range.end(); ++riter)
  {
    ParallelSet info;
    info.handle = *riter;
    info.contentsOffset = *offset_iter; ++offset_iter;
    info.childrenOffset = *offset_iter; ++offset_iter;
    info.contentsCount = *sizes_iter; ++sizes_iter;
    info.childrenCount = *sizes_iter; ++sizes_iter;
    info.description = *all_iter > 0; all_iter += 2;
    parallelSets.push_back( info );

      // Remove parallel meshsets from global list to
      // be written by non-parallel export code.
    setSet.range.erase( *riter );
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::assign_parallel_set_ids()
{
  const id_t start_id = setSet.first_id;
  int file_id;
  MBErrorCode rval;
  
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    rval = iFace->tag_get_data( idTag, &(iter->handle), 1, &file_id );
    assert( MB_SUCCESS == rval );
    file_id += start_id - 1;
    rval = iFace->tag_set_data( idTag, &(iter->handle), 1, &file_id );
    assert( MB_SUCCESS == rval );
  }
  return MB_SUCCESS;
}    
  

MBErrorCode WriteHDF5Parallel::write_shared_set_descriptions( hid_t table )
{
  const id_t start_id = setSet.first_id;
  MBErrorCode rval;
  mhdf_Status status;
  MBRange my_sets;
  
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
    unsigned long data[3] = { iter->contentsCount, iter->childrenCount, flags };
    mhdf_writeSetMeta( table, file_id, 1, H5T_NATIVE_LONG, data, &status );
    assert( !mhdf_isError( &status ) );
    
      // Store handle because will need to write tag data later.
    my_sets.insert( iter->handle );
  }
  
    // Put the parallel meshsets this proc is responsible for back
    // in the global list.  The serial code has already finished
    // writing the set data.  Need to make sure tag data gets written.
  setSet.range.merge( my_sets );
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
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list );
    assert( MB_SUCCESS == rval );
    
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
  
  for( std::list<ParallelSet>::iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_child_meshsets( iter->handle, handle_list );
    assert( MB_SUCCESS == rval );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list );
    assert( MB_SUCCESS == rval );
    
    if (id_list.empty())
      continue;
    
    mhdf_writeSetChildren( table, 
                           iter->childrenOffset, 
                           id_list.size(),
                           id_type,
                           &id_list[0],
                           &status );
    assert(!mhdf_isError(&status));
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
  assert( MB_SUCCESS == rval );
  rval = iFace->tag_get_name( t2.tag_id, name2 );
  assert( MB_SUCCESS == rval );
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
  
    // Communicate the number of ranges and the start_id for
    // each processor.
  std::vector<int> counts(numProc), offsets(numProc), displs(numProc);
  int mycount = myranges.size();
  int mystart = export_set->first_id + export_set->offset;
  result = MPI_Allgather( &mycount, 1, MPI_INT, &counts[0], numProc, MPI_INT, MPI_COMM_WORLD );
  assert(MPI_SUCCESS == result);
  result = MPI_Allgather( &mystart, 1, MPI_INT, &offsets[0], numProc, MPI_INT, MPI_COMM_WORLD );
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
    for (MBRange::iterator r_iter = remoteMesh[proc].begin();
         r_iter != remoteMesh[proc].end(); ++r_iter)
    {
        // Get the entities of the passed type in the geomtery meshset
      MBRange range;
      const MBEntityType tmp_type = iFace->type_from_handle( *r_iter );
      if (tmp_type == MBENTITYSET)
      {
        rval = iFace->get_entities_by_type( *r_iter, type, range );
        assert(MB_SUCCESS == rval);
      }
      else if (tmp_type == type)
      {
        range.insert( *r_iter );
      }
      
        // Use communicated processor information to find
        // corresponding fileid.
      for (MBRange::iterator sr_iter = range.begin();
           sr_iter != range.end(); ++sr_iter)
      {
        MBEntityHandle entity = *sr_iter;
        
          // Get handle on other processor
        MBEntityHandle global;
        rval = iFace->tag_get_data( global_id_tag, &entity, 1, &global );
        assert(MB_SUCCESS == rval);
        
          // Find corresponding fileid on other processor.
          // This could potentially be n**2, but we will assume that
          // the range list from each processor is short (typically 1).
        int j, steps = 0;
        unsigned long low, high;
        for (j = 0; j < count; ++j)
        {
          low = ranges[2*count];
          high = ranges[2*count+1];
          if (low >= global && high <= global)
            break;
          steps += (high - low) + 1;
        }
        assert(j < count);
        int fileid = offset + steps + (global - low);
        rval = iFace->tag_set_data( idTag, &entity, 1, &fileid );
        assert(MB_SUCCESS == rval);
      } // for(range)
    } // for(i_iter->range)
  } // for(interfaceList)
  
  return MB_SUCCESS;
}
