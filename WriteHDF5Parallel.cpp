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

struct elemdata {
  int mbtype;
  int numnode;
  int count;
};

MBErrorCode WriteHDF5Parallel::gather_interface_meshes()
{
  MBRange range;
  MBErrorCode result;
  MBTag iface_tag;
  int proc_pair[2];
  int rank, numproc;

  int init;
  MPI_Initialized( &init );
  assert( init );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  
    // Get interface mesh sets
  result = iFace->tag_get_handle( PARALLEL_INTERFACE_TAG_NAME, iface_tag );
  if (MB_SUCCESS != result) return result;
  result = iFace->get_entities_by_type_and_tag( 0,
                                                MBENTITYSET,
                                                &iface_tag,
                                                0,
                                                1,
                                                range );
  if (MB_SUCCESS != result) return result;
  
    // Populate lists with interface mesh sets
  for (MBRange::iterator iter = range.begin(); iter != range.end(); ++iter)
  {
    Interface iface;
    iface.set = *iter;
    result = iFace->tag_get_data( iface_tag, &iface.set, 1, proc_pair );
    if (proc_pair[0] == rank)
    {
      iface.rank = proc_pair[1];
      iface.remote = false;
      interfaceList.push_back( iface );
    }
    else
    {
      assert( proc_pair[1] == rank );
      iface.rank = proc_pair[0];
      iface.remote = true;
      interfaceList.push_back( iface );
    }
  }
  //std::sort( interfaceList );
  interfaceList.sort();
  
    // Make the contents of each interface meshset unique 
    // (not intersecting with others)
  for (std::list<Interface>::iterator iter = interfaceList.begin();
       iter != interfaceList.end(); ++iter)
  {
    result = iFace->get_entities_by_type( iter->set, MBENTITYSET, iter->range );
    if (MB_SUCCESS != result) return result;
    
    for (std::list<Interface>::iterator iter2 = interfaceList.begin();
         iter2 != iter; ++iter2)
      iter->range = iter->range.subtract( iter2->range );
  }
  
    // Get lists of all non-local mesh
  MBTag geom_tag;
  result = iFace->tag_get_handle( PARALLEL_GEOM_TOPO_TAG_NAME, geom_tag );
  if (MB_SUCCESS != result) return result;
  
  MBRange remoteList[MBPOLYGON], temp, remoteSets;
  for (std::list<Interface>::iterator iter = interfaceList.begin();
       iter != interfaceList.end(); ++iter)
  {
    if (!iter->remote)
      continue;
    
    remoteSets.insert( iter->set );
    for (MBRange::iterator riter = iter->range.begin(); 
         riter != iter->range.end(); ++riter)
    {
      MBEntityHandle set = *riter;
      assert( iFace->type_from_handle( set ) == MBENTITYSET );
      remoteSets.insert( set );
 
      int dimension;
      result = iFace->tag_get_data( geom_tag, &set, 1, &dimension );
      if (MB_SUCCESS != result) return result;
    
      for (MBEntityType type = MBVERTEX; 
           type <= MBPOLYGON && MBCN::Dimension(type) <= dimension; 
           ++type)
      {
        temp.clear();
        result = iFace->get_entities_by_type( set, type, temp );
        if (MB_SUCCESS != result) return result;
      
        remoteList[type].merge( temp );
      }
    }
  }
  
    // Remove remote mesh entities from export sets 
    // (not writing them on this processor.)
  for (std::list<ExportSet>::iterator iter = exportList.begin();
       iter != exportList.end(); ++iter)
    if (iter->type <= MBPOLYGON)
      iter->range.subtract( remoteList[iter->type] );
  nodeSet.range.subtract( remoteList[MBVERTEX] );
  setSet.range.subtract( remoteSets );
  
    // Give all remote entities an junk ID of 1.
    // We need to do this so that the code that gathers
    // adjacencies doesn't think these aren't being exported.
    // Will assign them proper IDs later.
  int id = 1;
  for (MBEntityType type = MBVERTEX; type <= MBPOLYGON; ++type)
  {
    for (MBRange::iterator riter = remoteList[type].begin(); 
         riter != remoteList[type].end(); ++riter)
    {
      MBEntityHandle ent = *riter;
      result = iFace->tag_set_data( idTag, &ent, 1, &id );
      if (MB_SUCCESS != result) return result;
    }
  }
  for (MBRange::iterator riter = remoteSets.begin(); 
       riter != remoteSets.end(); ++riter)
  {
    MBEntityHandle ent = *riter;
    result = iFace->tag_set_data( idTag, &ent, 1, &id );
    if (MB_SUCCESS != result) return result;
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
  
  rval = gather_interface_meshes();
  if (MB_SUCCESS != rval) return rval;
  
  rval = gather_tags();
  if (MB_SUCCESS != rval)
    return rval;
    
  int rank, numproc;
  result = MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  assert(result);
  result = MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  
  
    /**************** Create actual file and write meta info ***************/

  if (rank == 0)
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
 
    // gather node counts for each processor
  int* node_counts = new int[numproc];
  int num_nodes = nodeSet.range.size();
  result = MPI_Gather( &num_nodes, 1, MPI_INT, node_counts, 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(result);
  
    // create node data in file
  long first_id;
  if (rank == 0)
  {
    int total = 0;
    for (int i = 0; i < numproc; i++)
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
  if (rank == 0)
  {
    int prev_size = node_counts[0];
    node_counts[0] = 0;
    for (int i = 1; i < numproc; ++i)
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
  delete [] node_counts;
  
  writeUtil->assign_ids( nodeSet.range, idTag, (id_t)(nodeSet.first_id + nodeSet.offset) );

  
  
    /**************** Create element tables ***************/

    // Get element types and counts from each processor
  
  //std::sort( exportList );
  exportList.sort();
  int num_types = exportList.size();
  int max_num_types;
  result = MPI_Allreduce( &num_types, &max_num_types, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD );
  assert(result);
  
  elemdata* type_array = new elemdata[max_num_types];
  bzero( type_array, sizeof(elemdata)*max_num_types );
  elemdata* typeiter = type_array;
  for (std::list<ExportSet>::iterator iter = exportList.begin();
       iter != exportList.end(); ++iter)
  {
    typeiter->mbtype = iter->type;
    typeiter->numnode = iter->num_nodes;
    typeiter->count = iter->range.size();
    if (typeiter->mbtype == MBPOLYGON ||
        typeiter->mbtype == MBPOLYHEDRON)
    {
      rval = writeUtil->get_poly_array_size( iter->range.begin(), 
                                             iter->range.end(),
                                             typeiter->numnode );
      assert(MB_SUCCESS == rval);
    }
    ++typeiter;
  }
  
  elemdata* all_type_array = new elemdata[max_num_types*numproc];
  result = MPI_Gather( type_array,     3*max_num_types, MPI_INT,
                       all_type_array, 3*max_num_types, MPI_INT,
                       0, MPI_COMM_WORLD );
  assert(result);
  delete [] type_array;
  
  mhdf_ElemHandle* type2_array = 0;
  int* type_offset_list = 0;
  const int offset_list_size = 2 * num_types;
  if (rank == 0)
  {
      // break type array into per-proc arrays
    elemdata** iterlist = new elemdata*[numproc];
    elemdata** endlist = new elemdata*[numproc];
    elemdata* iteriter = all_type_array;
    for (int i = 0; i < numproc; i++)
    {
      iterlist[i] = iteriter;
      iteriter += numproc;
      endlist[i] = iteriter;
      
      while (endlist[i] > iterlist[i] && endlist[i-1]->mbtype == 0)
        --(endlist[i]);
    }
  
      // merge into one array
    elemdata* result_type_array = new elemdata[max_num_types*numproc];
    int result_index = 0, first = 0;
    while (first < numproc)  // while still a nonempty list
    {
        // find which type goes next in merged array
      int i, mintype, minnode;
      mintype = iterlist[first]->mbtype;
      minnode = iterlist[first]->numnode;
      for (i = first + 1; i < numproc; ++i)
      {
        if (iterlist[i] == endlist[i])
          ;
        else if (iterlist[i]->mbtype < mintype)
        {
          mintype = iterlist[i]->mbtype;
          minnode = iterlist[i]->numnode;
        }
        else if (iterlist[i]->mbtype == mintype &&
                 iterlist[i]->numnode < minnode)
        {
          minnode = iterlist[i]->numnode;
        }
      }
      
      result_type_array[result_index].mbtype = mintype;
      result_type_array[result_index].numnode = minnode;
      result_type_array[result_index].count = 0;
      if (mintype == MBPOLYGON || mintype == MBPOLYHEDRON)
        result_type_array[result_index].numnode = 0;
      
        // add values for all procs with the type
      for (i = first; i < numproc; i++)
      {
        if (iterlist[i] == endlist[i] ||
            iterlist[i]->mbtype != mintype)
          continue;
        
        if (mintype == MBPOLYGON || mintype == MBPOLYHEDRON)
        {
          result_type_array[result_index].count += iterlist[i]->count;
          result_type_array[result_index].numnode += iterlist[i]->numnode;
          ++(iterlist[i]);
        }
        else if (iterlist[i]->numnode == minnode)
        {
          result_type_array[result_index].count += iterlist[i]->count;
          ++(iterlist[i]);
        }
      }
      
        // iterate
      ++result_index;
      while (first < numproc && iterlist[first] == endlist[first])
        ++first;
    } // while (first < numproc)
    
    num_types = result_index;
    
      // create element tables in file
    type2_array = new mhdf_ElemHandle[num_types];
    for (int i = 0; i < num_types; i++)
    {
      if (result_type_array[i].mbtype == MBPOLYGON ||
          result_type_array[i].mbtype == MBPOLYHEDRON)
      {
        rval = create_poly_tables( (MBEntityType)result_type_array[i].mbtype,
                                   result_type_array[i].count,
                                   result_type_array[i].numnode,
                                   type2_array[i], 
                                   first_id );
      }
      else
      {
        rval = create_elem_tables( (MBEntityType)result_type_array[i].mbtype,
                                   result_type_array[i].numnode,
                                   result_type_array[i].count,
                                   type2_array[i],
                                   first_id );
      }
      
      if (MB_SUCCESS != rval)
        return rval;
        
        // re-use count field to communicate the first ID back to others.
      result_type_array[i].count = (int)first_id;
    }
    
      // create list of offsets for each processor
      
    type_offset_list = new int[offset_list_size * numproc];
    iteriter = all_type_array;
    for (int i = 0; i < numproc; i++)
    {
      iterlist[i] = iteriter;
      iteriter += numproc;
    }
    for (int i = 0; i < num_types; i++)
    {
      int prev_off = 0;
      int prev_poff = 0;
      for (int j = 0; j < numproc; j++)
      {
        int count = 0;
        int pcount = 0;
        if (iterlist[j] < endlist[j] && 
            iterlist[j]->mbtype == result_type_array[i].mbtype)
        {
          if (iterlist[j]->mbtype == MBPOLYGON || 
              iterlist[j]->mbtype == MBPOLYHEDRON)
          {
            count = iterlist[j]->count;
            pcount = iterlist[j]->numnode;
            ++(iterlist[j]);
          }
          else if (iterlist[j]->numnode == result_type_array[i].numnode)
          {
            count = iterlist[j]->count;
            ++(iterlist[j]);
          }
        }
        
        type_offset_list[i*offset_list_size + j*2] = prev_off;
        type_offset_list[i*offset_list_size + j*2+1] = prev_poff;
        prev_off += count;
        prev_poff += pcount;
      }
    }
    
    delete [] endlist;
    delete [] iterlist;
    delete [] all_type_array;
  }
  
    // send number of types to each processor
  result = MPI_Bcast( &num_types, 1, MPI_INT, 0, MPI_COMM_WORLD );
  assert(result);
  
    // send type array to each processor
  elemdata* result_type_array = 0;
  if (rank != 0)
  {
    result_type_array = new elemdata[num_types];
  }
  result = MPI_Bcast( result_type_array, 3*num_types, MPI_INT, 0, MPI_COMM_WORLD );
  assert(result);
  
    // send offsets to each processor
  int* offsets = new int[offset_list_size];
  result = MPI_Scatter( type_offset_list, offset_list_size, MPI_INT,
                        offsets,          offset_list_size, MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(result);
  delete [] type_offset_list;
  
    // fill in missing types and assign IDs for each type
  std::list<ExportSet>::iterator ex_itor = exportList.begin();
  for (int i = 0; i < num_types; ++i)
  {
    if (result_type_array[i].mbtype == MBPOLYGON ||
        result_type_array[i].mbtype == MBPOLYHEDRON)
    {
      if (ex_itor->type != result_type_array[i].mbtype)
      {
        ExportSet tmp;
        tmp.type = (MBEntityType)result_type_array[i].mbtype;
        tmp.num_nodes = 0;
        ex_itor = exportList.insert( ex_itor, tmp );
      }
      ex_itor->first_id = result_type_array[i].count;
      ex_itor->offset = offsets[2*i];
      ex_itor->poly_offset = offsets[2*i+1];
    }
    else
    {
      if (ex_itor->type != result_type_array[i].mbtype ||
          ex_itor->num_nodes != result_type_array[i].numnode)
      {
        ExportSet tmp;
        tmp.type = (MBEntityType)result_type_array[i].mbtype;
        tmp.num_nodes = result_type_array[i].numnode;
        ex_itor = exportList.insert( ex_itor, tmp );
      }
      ex_itor->first_id = result_type_array[i].count;
      ex_itor->offset = offsets[2*i];
      ex_itor->poly_offset = 0;
    }
  
    writeUtil->assign_ids( ex_itor->range, idTag,
                           (id_t)(ex_itor->first_id + ex_itor->offset));
  }
  
  delete [] result_type_array;
  

    /**************** Comminucate all remote IDs ***********************/
  
  rval = communicate_remote_ids();
  if (MB_SUCCESS != rval)
    return rval;
  
  
    /**************** Create adjacency tables *********************/
  
    // Count adjacencies for each element type
  int *type_adj_list = new int[num_types];
  int *adj_iter = type_adj_list;
  for (std::list<ExportSet>::iterator iter = exportList.begin();
       iter != exportList.end(); ++iter)
  {
    id_t num_adj;
    rval = count_adjacencies( iter->range, num_adj );
    if (MB_SUCCESS != rval) return rval;
  
    *adj_iter = num_adj;
    ++adj_iter;
  }
  
    // Gather counts on root processor
  int* proc_adj_counts = 0;
  if (rank == 0)
  {
    proc_adj_counts = new int[numproc * num_types];
  }
  result = MPI_Gather( type_adj_list,   num_types, MPI_INT,
                       proc_adj_counts, num_types, MPI_INT,
                       0, MPI_COMM_WORLD );
  assert(result);
  
    // create adjacency data tables in file
  if (rank == 0)
  {
    for (int i = 0; i < num_types; ++i)
    {
      int next_offset = 0;
      for (int j = 0; j < numproc; ++j)
      {
        int count = proc_adj_counts[i*num_types + j];
        proc_adj_counts[i*num_types + j] = next_offset;
        next_offset += count;
      }
      
      if (next_offset == 0)
        continue;
      
      hid_t handle = mhdf_createAdjacency( filePtr,
                                           type2_array[i],
                                           next_offset,
                                           &status );
      if (mhdf_isError( &status ))
      {
        writeUtil->report_error( "%s\n", mhdf_message( &status ) );
        return MB_FAILURE;
      }
      mhdf_closeData( filePtr, handle, &status );
    }
    delete [] type2_array;
  }
  
    // scatter adjacency table offsets to all processors
  result = MPI_Scatter( proc_adj_counts, num_types, MPI_INT,
                        type_adj_list,   num_types, MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(result);
  delete [] proc_adj_counts;
  
    // update offset lists on each proc
  adj_iter = type_adj_list;
  for (std::list<ExportSet>::iterator iter = exportList.begin();
       iter != exportList.end(); ++iter, ++adj_iter)
    iter->adj_offset = *adj_iter;
  
  
    /**************** Create meshset tables *********************/
  
  int* proc_set_offsets = 0;
  if (rank == 0)
    proc_set_offsets = new int[3*numproc];
  
  long tmp1, tmp2;
  rval = count_set_size( setSet.range, rangeSets, tmp1, tmp2 );
  if (MB_SUCCESS != rval) return rval;
  int set_offsets[] = { setSet.range.size(), (int)tmp1, (int)tmp2 };
  result = MPI_Gather( proc_set_offsets, 3, MPI_INT,
                            set_offsets, 3, MPI_INT,
                        0, MPI_COMM_WORLD );
  assert(result);
  
  if (rank == 0)
  {
    int next[] = { 0, 0, 0};
    for (int i = 0; i < numproc; i++)
    {
      int* array = proc_set_offsets + 3*i;
      for (int j = 0; j < 3; j++)
      {
        int count = array[j];
        array[j] = next[j];
        next[j] += count;
      }
    }
    
    rval = create_set_tables( next[0], next[1], next[2], first_id );
    if (MB_SUCCESS != rval) return rval;
  }
  result = MPI_Bcast( &first_id, 1, MPI_LONG, 0, MPI_COMM_WORLD );
  assert(result);
  
  result = MPI_Scatter( proc_set_offsets, 3, MPI_INT,
                             set_offsets, 3, MPI_INT,
                            0, MPI_COMM_WORLD );
  assert(result);
  delete [] proc_set_offsets;
  
  setSet.offset = set_offsets[0];
  setContentsOffset = set_offsets[1];
  setChildrenOffset = set_offsets[2];
  
   
  
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

// statics initialized in communicate_remote_ids().
static MBInterface* sorterIface = 0;
static MBTag sorterTag = 0;
static int qsortCompare( const void* p1, const void* p2 )
{
  MBEntityHandle handles[] = { *(const MBEntityHandle*)p1, 
                               *(const MBEntityHandle*)p2 };
  MBEntityHandle ids[2];
  MBErrorCode rval = sorterIface->tag_get_data( sorterTag, handles, 2, ids );
  assert (MB_SUCCESS == rval);
  return ids[0] < ids[1] ? -1 : ids[0] > ids[1] ? 1 : 0;
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

MBErrorCode WriteHDF5Parallel::start_send( Interface& iface )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  
    // Get handles for mesh entities in interface
  MBRange range;
  mesh_from_geom_sets( iface.range, range );
    // Copy range into array buffer
  iface.buffer.resize( range.size() );
  std::copy( iface.range.begin(), iface.range.end(), iface.buffer.begin() );
    // Sort array by "global" identifiers
  qsort( &(iface.buffer[0]), iface.buffer.size(), sizeof(MBEntityHandle), &qsortCompare );
    // **FIX ME** need to handle 64-bits, but is it possible to
    //            have 32-bit IDs and 64-bit handles?
  assert( sizeof(MBEntityHandle) == sizeof(int) );
    // Change buffer to contain IDs
  for (std::vector<MBEntityHandle>::iterator iter = iface.buffer.begin();
       iter != iface.buffer.end(); ++iter)
  {
    MBEntityHandle handle = *iter;
    MBErrorCode rval = iFace->tag_get_data( idTag, &handle, 1, &*iter );
    assert( MB_SUCCESS != rval && 0 != *iter );
  }
    // Send buffer
  int result = MPI_Isend( &(iface.buffer[0]), iface.buffer.size(), MPI_INT,
                          iface.rank, rank, MPI_COMM_WORLD, &iface.req );
  assert(result);
  return MB_SUCCESS; 
}

MBErrorCode WriteHDF5Parallel::start_recv( Interface& iface )
{
    // Get handles for mesh entities in interface
  MBRange range;
  mesh_from_geom_sets( iface.range, range );

  iface.buffer.resize( range.size() );
  int result = MPI_Irecv( &(iface.buffer[0]), iface.buffer.size(), MPI_INT,
                          iface.rank, iface.rank, MPI_COMM_WORLD, &iface.req );
  assert(result);
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::finish_recv( Interface& iface )
{
    // Get handles for mesh entities in interface
  MBRange range;
  mesh_from_geom_sets( iface.range, range );
    // Copy range into array 
  std::vector<MBEntityHandle> list( range.size() );
  std::copy( iface.range.begin(), iface.range.end(), list.begin() );
    // Sort array by "global" identifiers
  qsort( &(list[0]), list.size(), sizeof(MBEntityHandle), &qsortCompare );
    // **FIX ME** need to handle 64-bits, but is it possible to
    //            have 32-bit IDs and 64-bit handles?
  assert( sizeof(MBEntityHandle) == sizeof(int) );
    // Set IDs for entities
  MBErrorCode rval = iFace->tag_set_data( idTag,
                                          &(list[0]),
                                          list.size(),
                                          &(iface.buffer[0]) );
  assert(MB_SUCCESS == rval);
  return rval;
}

MBErrorCode WriteHDF5Parallel::communicate_remote_ids()
{
    // set up statics for use in qsort compare function
  sorterIface = iFace;
  iFace->tag_get_handle( PARALLEL_GLOBAL_ID_TAG_NAME, sorterTag );
  
  std::list<Interface>::iterator iter;
  const std::list<Interface>::iterator end = interfaceList.end();
  
    // Start all communications
  for (iter = interfaceList.begin(); iter != end; ++iter)
  {
    MBErrorCode rval = iter->remote ? start_recv( *iter ) : start_send( *iter );
    assert( MB_SUCCESS == rval );
    if (MB_SUCCESS != rval) return rval;
  }
  
    // Handle communications as they finish
  std::vector<MPI_Request> req_list( interfaceList.size() );
  while (interfaceList.size())
  {
      // Get array of pending requests
    req_list.resize( interfaceList.size() );
    std::vector<MPI_Request>::iterator req_iter = req_list.begin();
    for (iter = interfaceList.begin(); iter != end; ++iter, ++req_iter)
      *req_iter = iter->req;
    
      // Wait for one to finish
    int index;
    int result = MPI_Waitany( req_list.size(),
                              &(req_list[0]),
                              &index,
                              MPI_STATUS_IGNORE );
    assert(result);
    MPI_Request req = req_list[index];
    
      // Get the interface that completed
    for (iter = interfaceList.begin(); iter != interfaceList.end(); ++iter)
      if (iter->req == req)
        break;
    assert(iter != interfaceList.end());
    
      // Handle finished request
    if (iter->remote)
    {
      MBErrorCode rval = finish_recv( *iter );
      assert (MB_SUCCESS == rval);
      if (MB_SUCCESS != rval) return rval;
    }
    interfaceList.erase( iter );
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::mesh_from_geom_sets( const MBRange& sets, MBRange& mesh )
{
  MBRange temp;
  for (MBRange::const_iterator iter = sets.begin(); iter != sets.end(); ++iter)
  {
    temp.clear();
    MBErrorCode rval = iFace->get_entities_by_handle( *iter, temp );
    if (MB_SUCCESS != rval) { assert(0); return rval; }
    mesh.merge( temp );
  }
  return MB_SUCCESS;
}

