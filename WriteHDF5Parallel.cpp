#include <assert.h>
#include "mhdf.h"
#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include <H5FDmpio.h>
#include "MBInterface.hpp"
#include "MBInternals.hpp"
#include "MBTagConventions.hpp"
#include "WriteHDF5Parallel.hpp"
#include "MBParallelConventions.h"



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
             mbtype == MBMESHSET ||
             numnode == other.numnode);
  }
  bool operator<( const elemtype& other ) const
  {
    if (mbtype > other.mbtype)
      return false;
   
    return mbtype < other.mbtype ||
           (mbtype != MBPOLYGON &&
            mbtype != MBPOLYHEDRON &&
            mbtype != MBMESHSET &&
            numnode < other.numnode);
  }
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
                                      const char*[] multiproc_tags )
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
  multiProcSetTags.sort();
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
    result = iFace->tag_get_data( iface_tag, *&iiter, 1, proc_pair );
    if (MB_SUCCESS != result) return result;
    const int remote_proc = proc_pair[0];
    
      // Get list of all entities in interface and 
      // the subset of that list that are meshsets.
    MBRange entities, sets;
    result = iFace->get_entities_by_handle( *iiter, entities );
    if (MB_SUCCESS != result) return result;
    result = iFace->get_entities_by_type( *iiter, MBMESHSET, sets );
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
    
    nodeSet.range.subract( range );
    setSet.range.subract( range );
    for (std::list<ExportSet>::iterator eiter = exportList.begin();
         eiter != exportList.end(); ++eiter )
      eiter->range.subtract( range );
    
    int id = 1;
    for (MBRange::iterator riter = remoteList[type].begin(); 
         riter != remoteList[type].end(); ++riter)
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
  assert(result);
  result = MPI_Comm_size( MPI_COMM_WORLD, &numProc );
  assert(result);
  
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
  for (ex_itor = exportList.begin(); ex_itor != exportList.end(); ++ex_itor)
  {
    rval = communicate_remote_ids( ex_itor->type );
    assert(MB_SUCCESS == rval);
  }
  
  
    /**************** Create adjacency tables *********************/
  
  rval = create_adjcency_tables();
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
  if (rank == 0)
    proc_tag_offsets = new int[num_tags*numproc];
  result = MPI_Gather( tag_offsets, num_tags, MPI_INT,
                  proc_tag_offsets, num_tags, MPI_INT,
                       0, MPI_COMM_WORLD );
  assert(result);
  
  if (rank == 0)
  {
    tag_iter = tagList.begin();
    for (int i = 0; i < num_tags; ++i, ++tag_iter)
    {
      int next_offset = 0;
      for (int j = 0; j < numproc; j++)
      {
        int count = proc_tag_offsets[i*numproc + j];
        proc_tag_offsets[i*numproc+j] = next_offset;
        next_offset += count;
      }
      
      rval = create_tag( tag_iter->tag_id, next_offset );
      assert(MB_SUCCESS == rval);
    }
  }
  
  result = MPI_Scatter( proc_tag_offsets, num_tags, MPI_INT,
                             tag_offsets, num_tags, MPI_INT,
                             0, MPI_COMM_WORLD );
  assert(result);
  delete [] proc_tag_offsets;
  
  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
    tag_iter->offset = tag_offsets[i];
  
  
  /************** Close serial file and reopen parallel *****************/
  
  if (0 == rank)
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
 
    // gather node counts for each processor
  std::vector<int> node_counts(numProc);
  int num_nodes = nodeSet.range.size();
  result = MPI_Gather( &num_nodes, 1, MPI_INT, &node_counts[0], 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(result);
  
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
  assert(result);
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
  result = MPI_Scatter( node_counts, 1, MPI_INT, 
                        &offset, 1, MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(result);
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
  assert(result);
  
    // Get list of types on this processor
  std::vector<int> my_types(num_types);
  std::vector<int>::iterator viter = my_types.begin();
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
       eiter != exportList.end(); ++eiter)
  {
    *viter = eiter->type;
    ++viter;
    *viter = eiter->num_nodes;
    ++viter;
  }
  
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
  assert(result);
  
    // Merge type lists
  std::list<elemtype> type_list;
  std::list<elemtype>::iterator liter;
  for (int i = 0; i < numProc; ++i)
  {
    int num_types = counts[i] / 2;
    int* typelist = &alltypes[displs[i]];
    liter = type_list.begin();
    for (int j = 0; j < counts[i]; j += 2)
    {
      int* type = &typelist[j];
        // skip until insertion spot
      for (; liter != typelist.end() && *liter < type; ++liter);
      
      if (liter == typelist.end() || *liter != type)
        liter = typeList.insert( liter, type );
    }
  }
  
    // Send total number of types to each processor
  total = type_list.size();
  result = MPI_Bcast( &total, 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(result);
  
    // Send list of types to each processor
  std::vector<int> intlist(total * 2);
  viter = intlist.begin();
  for (liter = type_list.begin(); liter != type_list.end(); ++liter)
  {
    *viter = liter->mbtype;
    ++viter;
    *viter = liter->numnode;
    ++viter;
  }
  result = MPI_Bcast( &intlist[0], 2*total, MPI_INT, 0, MPI_COMM_WORLD );
  assert(result);
  
    // Insert missing types into exportList, with an empty
    // range of entities to export.
  std::list<ExportSet>::iterator ex_iter = exportList.begin();
  viter = intlist;
  for (int i = 0; i < total; ++i)
  {
    int mbtype = *viter; ++viter;
    int numnode = *viter; ++viter;
    while (ex_iter != exportList.end() && ex_iter->type < mbtype)
      ++ex_iter;
    
    bool equal = ex_iter != exportList.end() && ex_iter->type == mbtype;
    if (equal && mbtype != MBPOLYGON && mbtype != MBPOLYHEDRON)
    {
      while (ex_iter != exportList.end() && ex_iter->num_nodes < numnodes)
        ++ex_iter;
        
      equal = ex_iter != exportList.end() && ex_iter->num_nodes == numnodes;
    }
    
    if (!equal)
    {
      ExportSet insert;
      insert.type = mbtype;
      insert.num_nodes = numnode;
      insert.first_id = 0;
      insert.offset = 0;
      insert.poly_offset = 0;
      insert.adj_offset = 0;
      insert.type2 = 0;
      ex_iter = exportList.insert( ex_iter, insert );
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::create_element_tables()
{
  int result;
  MBErrorCode rval;
  std::list<ExportSet>::interator ex_iter;
  std::vector<long>::iterator viter;
  
    // Get number of each element type from each processor
  const int numtypes = exportList.size();
  std::vector<long> my_counts(numtypes);
  std::vector<long> counts(numtypes * numProc + numtypes);
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { *viter = ex_iter->range()->size(); ++viter; }
  
  result = MPI_Gather( &my_counts[0], numtypes, MPI_LONG,
                       &counts[0],    numtypes, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(result);
  
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
  assert(result);
  
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
    else if(ex_iter->type == POLYHEDRON)
    {
      asesrt(!poly[1]);
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
    assert(result);
    
    int prev = 0;
    for (int j = 1; j <= numProc; j++)
    {
      int tmp = perproc[j];
      perproc[j] = prev;
      prev += tmp;
    }
                                           
    polycount[i] = perproc[numProc];
    result = MPI_Scatter( &perproc[0], 1, MPI_INT, &count, 1, MPI_INT, 0, MPI_COMM_WORLD );
    assert(result);
    ppoly->poly_offset = count;
  }
  
    // Write element tables
  std::vector<long> start_ids(numtypes);
  std::vector<mhdf_ElemHandle> mhdf_types(numtypes);
  if (myRank == 0)
  {
    viter = start_ids.begin();
    long* citer = &counts[numtypes * numProc];
    std::vector<mhdf_ElemHandle>::iterator hiter = mhdf_types.begin();
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    {
      switch(ex_iter->type) {
      case MBPOLYGON:
        rval = create_poly_tables( MBPOLYGON,
                                   *citer,
                                   polycount[0],
                                   *hiter,
                                   *viter );
        break;
      case MBPOLYHEDRON:
        rval = create_poly_tables( MBPOLYHEDRON,
                                   *citer,
                                   polycount[1],
                                   *hiter,
                                   *viter );
        break;
      default:
        rval = create_elem_tables( ex_iter->type,
                                   ex_iter->num_nodes,
                                   *citer,
                                   *hiter,
                                   *viter );
      }
      assert(MB_SUCCESS == rval);
      ex_iter->type2 = *hiter;
      ++citer;
      ++hiter;
      ++viter;
    }
  }
  
    // send start IDs to each processor
  result = MPI_Bcast( &start_ids[0], numtypes, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(result);
  
    // Assign IDs to local elements
  viter = start_ids.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    ex_iter->first_id = *(viter++);
    id_t myfirst = (id_t)(ex_iter->first_id + ex_iter->offset)
    rval = writeUtil->assign_ids( ex_iter->range, idTag, myfirst );
    assert(MB_SUCCESS == rval);
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode WriteHDF5Parallel::create_adjacency_tables()
{
  MBErrorCode rval;
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
    rval = count_adjcencies( ex_iter->range, num_adj );
    assert (MB_SUCCESS == rval);
    *viter = num_adj; ++viter;
  }
  
    // End local adjacency counts to root processor
  result = MPI_Gather( &local[0], numtypes, MPI_LONG,
                       &all[0],   numtypes, MPI_LONG, 
                       0, MPI_COMM_WORLD );
  
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
  
    // Record the adjacency offset in each ExportSet
  viter = local.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { ex_iter->adj_offset = *viter; ++viter; }
  
    // Create data tables in file
  if (myRank == 0)
  {
    viter = &all[numtypes * numProc];
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter, ++viter)
    {
      if (!*viter) 
        continue;
      
      hid_t handle = mhdf_createAdjacency( filePtr,
                                           ex_iter->type2,
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
  int result;
  long prev[3];
  
  rval = negotiate_shared_meshsets(prev);
  if (MB_SUCCESS != rval) return rval;

    // Gather counts from each proc
  std::vector<long> set_offsets(3*numProc+3);
  long counts[3];
  counts[1] = setSet.range.size();
  rval = count_set_size( setSet.range, rangeSets, counts[1], counts[2] );
  if (MB_SUCCESS != rval) return rval;
  
  result = MPI_Gather( counts,          3, MPI_LONG, 
                       &set_offests[0], 3, MPI_LONG,
                       0, MPI_COMM_WORLD );
  assert( result );
  
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
  assert(result);
  setSet.offset = (id_t)(counts[0]);
  setContentsOffset = counts[1];
  setChildrenOFfset = counts[2];
  
    // Create the tables
    
    // The last three entries in the vector contain
    // the total table sizes.  Reuse the previous entry
    // to hold the start Id.
  long* totals = &set_offsets[3*numProc-1];
  if (myRank == 0 && totals[0])
  {
    rval = create_set_tables( totals[1], totals[2], totals[3], totals[0] );
    if (MB_SUCCESS != rval) return rval;
  }
  
    // Send totals to all processors
  result = MPI_Bcast( totals, 4, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(result);
  setSet.first_id  = totals[0];
  writeSets        = totals[1] > 0;
  writeSetContents = totals[2] > 0;
  writeSetChildren = totals[3] > 0;

  return MB_SUCCESS;
}

static int int_sort( const void* ptr1, const void* ptr2 )
  { return *(const int*)ptr1 - *(const int*)ptr2; }

MBErrorCode WriteHDF5Parallel::negotiate_shared_meshsets( long* offsets )
{
  MBErrorCode rval;
  bzero( offsets, 3*sizeof(long) );
  for (std::vector<string>::iterator niter = multiProcSetTags.begin();
       niter != multiProcSetTags.end(); ++niter)
  {
    rval = negotiate_shared_meshests( niter->c_str(), offsets );
    assert(MB_SUCCESS == rval);
  }
  return rval;
}

MBErrorCode WriteHDF5Parallel::negotiate_shared_meshests( const char* tagname,
                                                          long* offsets )
{
  int i;
  MBErrorCode rval;
  MBRange::const_iterator riter;

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
                                                MBMESHSET, 
                                                &handle,
                                                0,
                                                1,
                                                range );
    if (rval != MB_SUCCESS) return rval;
  }

    // Exchange number of sets with tag between all processors
  std::vector<int> counts(numProc);
  const int count = 3*range.size();
  result = MPI_Allgather( &count,     1, MPI_INT, 
                          &counts[0], 1, MPI_INT,
                          0, MPI_COMM_WORLD );
  assert(result);

    // Exchange tag values and sizes for sets between all processors
    // Send three integers for each set with the tag: {tag value, 
    // num contents, num children}.  Values are interleaved in array.
  std::vector<int> displs(numProc+1);
  displs[0] = 0;
  for (i = 0; i <= numProc; i++)
    displs[i] = displs[i-1] + counts[i-1];
  int total = displs[numProc];
  std::vector<int> all(total);
  std::vector<int> local(count);
  std::vector<int>::iterator viter = local.begin();
  for (riter = range.begin(); riter != range.end(); ++riter)
  {
    rval = iFace->tag_get_data( handle, &*riter, 1, &*viter ); ++viter;
    assert(MB_SUCCESS == rval);
    rval = iFace->get_number_entities_by_handle( *riter, *viter ); ++viter;
    assert(MB_SUCCESS == rval);
    rval = iFace->num_child_meshsets( *riter, &*viter ); ++viter;
    assert(MB_SUCCESS == rval);
  }      
  result = MPI_Allgatherv( &local[0], count, MPI_INT,
                           &all[0], &counts[0], &displs[0], MPI_INT, 
                           0, MPI_COMM_WORLD );
  assert(result);

    // Find sets that span multple processors and update appropriately.
    // The first processor (sorted by MPI rank) that contains a given set
    // will be responsible for writing the set description.  All multi-
    // processor sets will be written at the beginning of the set tables.
    // Processors will write set contents/children for a given set in
    // the order of their MPI rank.
    //
    // Update information in-place in the array from the Allgatherv.
    // Change the corresponding sizes for the first instance of a tag
    // value such that it ends up being the total size of the set.
    // Change the size to -1 for the later instances of a tag value.
    //
    // For the sets that this processor has, update the offsets at
    // which the set data is to be written.  Re-use the array of 
    // local counts ("local") to store this data.
  std::map<int,int> tagsort;  // Map of {tag value, index of first set w/ value}
  for (i = 0; i < total; i+=3)
  {
      // Update count of data table offset
    offsets[1] += all[i+1]; // update count of set contents
    offsets[2] += all[i+2]; // update count of set children
    
    const std::map<int,int>::iterator p = tagsort.find( all[i] );
    const int r = i - displs[myRank];  // offset in "local" array
    
      // If this is the first instance of this tag value, 
      // then the processor with this instance is responsible
      // for writing the tag description
    long first_meta_offset = -1;
    if ( p == tagsort.end() )  
    {
      tagsort.insert( make_pair(all[i], i) );
        // If within the range for this processor, save offsets
      if (r >= 0 && r < count) 
      {
        local[r+1] = local[r+2] = 0;
        if (first_meta_offset < 0)
          first_meta_offset = offsets[0];
      }
      offsets[0]++; // update number of unique sets
    }
      // Otherwise update the total size in the table
      // for the processor that is responsible for writing
      // the data and mark the data for the current processor
      // with a -1.
    else 
    {
        // If within the range for this processor, save offsets
      int j = p->second;
      if (r >= 0 && r < count) 
      {
        local[r+1] = all[j+1];
        local[r+2] = all[j+2];
      }
      
      all[j+1] += all[i+1];
      all[j+2] += all[i+2];
      all[i+1] = -1;
    }
  }
  
    // And information to list of parallel sets.
  riter = range.begin();
  viter = local.begin();
  std::vector<int> aiter = all.begin() + displs[myRank];
  while (viter != local.end())
  {
    assert( riter != range.end() );
    
    ParallelSet info;
    info.handle = *riter; ++riter;
    ++viter;
    info.dataOffset = *viter; ++viter;
    info.childOffset = *viter; ++viter;
    ++aiter;
    if (*aiter >=0)
    {
      info.metaOffset = first_meta_offset++;
      info.dataCount = *aiter; ++aiter;
      info.childCount = *aiter; ++aiter;
    }
    else
    {
      info.metaOffset = -1;
      aiter += 2;
    }
    
    parallelSets.push_back( info );
  }
  
  return MB_SUCCESS;
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


MBErrorCode WriteHDF5Parallel::communicate_remote_ids( MBEnittyType type )
{
  int rank, numproc, result;
  MBErrorCode rval;
  result = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  assert(result);
  result = MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  assert(result);


    // Get the export set for the specified type
  ExportSet* export_set = 0;
  if (type == MBVERTEX)
    export_set = &nodeSet;
  else if(type == MBMESHSET)
    export_set = &setSet;
  else
  {
    for (std::list<ExportSet>::iterator esiter = exportList.begin();
         esiter != exportList.end(); ++esiter)
      if (esiter->type == type)
      {
        export_set = &*iter;
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
    myranges.push_back( p_iter->first );
    myranges.push_back( p_iter->second );
  }
  
    // Communicate the number of ranges and the start_id for
    // each processor.
  std::vector<int> counts(numproc), offsets(numproc), displs(numproc);
  const int mycount = myranges.size();
  const int mystart = export_set->first_id + export_set->offset;
  result = MPI_Allgather( &mycount, 1, MPI_INT, &counts[0], numproc, MPI_INT, MPI_COMM_WORLD );
  assert(result);
  result = MPI_Allgather( &mystart, 1, MPI_INT, &offsets[0], numproc, MPI_INT, MPI_COMM_WORLD );
  assert(result);
  
    // Communicate the ranges 
  displs[0] = 0;
  for (int i = 1; i < numproc; ++i)
    displs[i] = displs[i-1] + counts[i-1];
  std::vector<unsigned long> allranges( displs[numproc-1] + counts[numproc-1] );
  result = MPI_Allgatherv( &myranges[0], myranges.size(), MPI_UNSIGNED_LONG,
                           &allranges[0], &counts[0], &displs[0],
                           MPI_UNSIGNED_LONG, MPI_COMM_WORLD );
  assert(result);
  
  MBTag global_id_tag;
  iFace->tag_get_handle( PARALLEL_GLOBAL_ID_TAG_NAME, global_id_tag );
  
    // Set file IDs for each communicated entity
    
    // For each interface meshset
  for (std::list<Interface>::iterator i_iter = interfaceList.begin();
       i_iter != iterfaceList.end(); ++i_iter)
  {
      // Get data for corresponding processor
    const int proc = i_iter->rank;
    const int offset = offsets[proc];
    const int count = counts[proc];
    const unsigned long* const ranges = &allranges[displs[proc]];
    
      // For each geometry meshset in the interface
    for (MBRange::iterator r_iter = i_iter->range.begin();
         r_iter != i_iter->range.end(); ++r_iter)
    {
        // Get the entities of the passed type in the geomtery meshset
      MBRange range;
      if (iFace->type_from_id(*r_iter) == MBENTITYSET)
      {
        rval = iFace->get_entities_by_type( *r_iter, type, range );
        assert(MB_SUCCESS == rval);
      }
      else
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
        int steps = 0;
        for (int j = 0; j < count; ++j)
        {
          const unsigned long low = ranges[2*count];
          const unsigned long high = ranges[2*count+1];
          
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


MBErrorCode WriteHDF5Parallel::write_parallel_sets()
{
  mhdf_Status status;
  MBErrorCode rval;
  long first_id, meta_size, data_size, child_size;
  hid_t set_table = 0, content_table = 0, child_table = 0;
  setSet.type2 = mhdf_set_type_handle();
  
  /* Get this from the file because we if working in parallel, 
     we may need to do the collective open/closes even if there
     is no set data on this processor */
  int do_sets, do_contents, do_children;
  do_sets = mhdf_haveSets( filePtr, &do_contents, &do_children, &status );
  CHK_MHDF_ERR_0(status);
  
  /* If no sets, just return success */
  if (!do_sets)
    return MB_SUCCESS;
  
  /* Write set description table and set contents table */
  
  /* Create the table */
  set_table = mhdf_openSetMeta( filePtr, &meta_size, &first_id, &status );
  CHK_MHDF_ERR_0(status);
  
  if (do_contents)
  {
    content_table = mhdf_openSetData( filePtr, &data_size, &status );
    CHK_MHDF_ERR_1(status, set_table);
  }
  
    
  MBRange set_contents;
  MBRange::const_iterator iter = sets.begin();
  MBRange::const_iterator comp = rangeSets.begin();
  const MBRange::const_iterator end = sets.end();
  long set_data[3];
  long set_offset = setSet.offset;
  long content_offset = setContentsOffset;
  long child_offset = setChildrenOffset;
  unsigned long flags;
  std::vector<id_t> id_list;
  std::vector<MBEntityHandle> handle_list;
  for ( ; iter != end; ++iter )
  {
    rval = get_set_info( *iter, data_size, child_size, flags );
    CHK_MB_ERR_2C(rval, set_table, do_contents, content_table, status);
    
    if (*iter == *comp)
    {
      set_contents.clear();
      id_list.clear();
      
      rval = iFace->get_entities_by_handle( *iter, set_contents, false );
      CHK_MB_ERR_2C(rval, set_table, do_contents, content_table, status);

      int length;
      rval = range_to_id_list( set_contents, &id_list, length );
      CHK_MB_ERR_2C(rval, set_table, do_contents, content_table, status);

      assert (length < data_size);
      flags |= mhdf_SET_RANGE_BIT;
      data_size = length;
      ++comp;
    }
    else
    {
      id_list.clear();
      handle_list.clear();
      
      rval = iFace->get_entities_by_handle( *iter, handle_list, false );
      CHK_MB_ERR_2C(rval, set_table, do_contents, content_table, status);
      
      rval = vector_to_id_list( handle_list, id_list );
      CHK_MB_ERR_2C(rval, set_table, do_contents, content_table, status);
    }
    
    child_offset += child_size;
    set_data[0] = content_offset + data_size - 1;
    set_data[1] = child_offset - 1;
    set_data[2] = flags;

    mhdf_writeSetMeta( set_table, set_offset++, 1L, H5T_NATIVE_LONG, set_data, &status );
    CHK_MHDF_ERR_2C(status, set_table, do_contents, content_table );
    
    if (id_list.size())
    {
      mhdf_writeSetData( content_table, 
                         content_offset,
                         id_list.size(),
                         id_type,
                         &id_list[0],
                         &status );
      CHK_MHDF_ERR_2C(status, set_table, do_contents, content_table );
      content_offset += data_size;
    }
  }
  
  mhdf_closeData( filePtr, set_table, &status );
  mhdf_closeData( filePtr, content_table, &status );
  if (!do_children)
    return MB_SUCCESS;
      
  
  /* Write set children */
  
  child_offset = setChildrenOffset;
  child_table = mhdf_openSetChildren( filePtr, &child_size, &status );
  CHK_MHDF_ERR_0(status);
   
  for (iter = sets.begin(); iter != end; ++iter)
  {
    handle_list.clear();
    id_list.clear();
    rval = iFace->get_child_meshsets( *iter, handle_list, 1 );
    CHK_MB_ERR_1(rval, child_table, status);

    if (handle_list.size() == 0)
      continue;

    rval = vector_to_id_list( handle_list, id_list );
    CHK_MB_ERR_1(rval, child_table, status);


    mhdf_writeSetChildren( child_table, 
                           child_offset, 
                           id_list.size(), 
                           id_type, 
                           &id_list[0], 
                           &status );
    CHK_MHDF_ERR_1(status, child_table);
    child_offset += id_list.size();
  }

  mhdf_closeData( filePtr, child_table, &status );
  return MB_SUCCESS;
}
  
