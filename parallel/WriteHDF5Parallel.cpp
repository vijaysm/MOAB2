
#undef DEBUG

#ifdef DEBUG
#  include <stdio.h>
#  include <stdarg.h>
#endif

#include <stdio.h>
#include <time.h>

#ifndef HDF5_FILE
#  error Attempt to compile WriteHDF5Parallel with HDF5 support disabled
#endif

#include <stdlib.h>
#include <string.h>

#include <vector>
#include <set>
#include <map>
#include <utility>
#include <iostream>
#include <sstream>

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
#include "MBParallelComm.hpp"
#include "MBCN.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBRange.hpp"

#include "WriteHDF5Parallel.hpp"


#define CHECK_MPI(A) do { if (MPI_SUCCESS != (A)) { \
  writeUtil->report_error( "MPI Failure at " __FILE__ ":%d\n", __LINE__ ); \
  return MB_FAILURE; } } while(false)

#define CHECK_HDF(A) do { if (mhdf_isError(&(A))) { \
  writeUtil->report_error( "MHDF Failure at " __FILE__ ":%d : %s\n", \
    __LINE__, mhdf_message(&(A)) ); \
  return MB_FAILURE; } } while(false)


#ifdef VALGRIND
#  include <valgrind/memcheck.h>
#else
#  ifndef VALGRIND_CHECK_MEM_IS_DEFINED
#    define VALGRIND_CHECK_MEM_IS_DEFINED
#  endif
#  ifndef VALGRIND_CHECK_MEM_IS_ADDRESSABLE
#    define VALGRIND_CHECK_MEM_IS_ADDRESSABLE
#  endif
#  ifndef VALGRIND_MAKE_MEM_UNDEFINED
#    define VALGRIND_MAKE_MEM_UNDEFINED
#  endif
#endif

template <typename T> inline 
void VALGRIND_MAKE_VEC_UNDEFINED( std::vector<T>& v ) {
    VALGRIND_MAKE_MEM_UNDEFINED( &v[0], v.size() * sizeof(T) );
}

#define TPRINT(A)
//#define TPRINT(A) tprint( (A) )
static void tprint(const char* A) 
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  char buffer[128]; 
  sprintf(buffer,"%02d: %6.2f: %s\n", rank, (double)clock()/CLOCKS_PER_SEC, A);
  fputs( buffer, stderr ); 
}

#ifdef DEBUG
#  define START_SERIAL                     \
     for (unsigned _x = 0; _x < myPcomm->proc_config().proc_size(); ++_x) {\
       MPI_Barrier( myPcomm->proc_config().proc_comm() );      \
       if (_x != myPcomm->proc_config().proc_rank()) continue     
#  define END_SERIAL                       \
     }                                     \
     MPI_Barrier( myPcomm->proc_config().proc_comm() )
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


#ifdef NDEBUG
#  undef assert
#  define assert
#else
#  undef assert
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


#ifndef DEBUG
void WriteHDF5Parallel::printrange( MBRange& ) {}
#else
void WriteHDF5Parallel::printrange( MBRange& r )
{
  int rank;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MBEntityType type = MBMAXTYPE;
  for (MBRange::const_pair_iterator i = r.const_pair_begin(); i != r.const_pair_end(); ++i)
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
static void print_type_sets( MBInterface* iFace, int rank, int size, MBRange& sets )
{
  MBTag gid, did, bid, sid, nid;
  iFace->tag_get_handle( GLOBAL_ID_TAG_NAME, gid ); 
  iFace->tag_get_handle( GEOM_DIMENSION_TAG_NAME, did );
  iFace->tag_get_handle( MATERIAL_SET_TAG_NAME, bid );
  iFace->tag_get_handle( DIRICHLET_SET_TAG_NAME, nid );
  iFace->tag_get_handle( NEUMANN_SET_TAG_NAME, sid );
  MBRange typesets[10];
  const char* typenames[] = {"Block", "Sideset", "NodeSet", "Vertex", "Curve", "Surface", "Volume", "Body", "Other"};
  for (MBRange::iterator riter = sets.begin(); riter != sets.end(); ++riter)
  {
    unsigned dim, id, oldsize;
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
    else {
      id = *riter;
      dim = 9;
    }

    oldsize = typesets[dim].size();
    typesets[dim].insert( id );
    assert( typesets[dim].size() - oldsize == 1 );  
  }
  for (int ii = 0; ii < 9; ++ii)
  {
    char num[16];
    std::string line(typenames[ii]);
    if (typesets[ii].empty())
      continue;
    sprintf(num, "(%lu):",(unsigned long)typesets[ii].size());
    line += num;
    for (MBRange::const_pair_iterator piter = typesets[ii].const_pair_begin();
         piter != typesets[ii].const_pair_end(); ++piter)
    {
      sprintf(num," %lx", (unsigned long)(*piter).first);
      line += num;
      if ((*piter).first != (*piter).second) {
        sprintf(num,"-%lx", (unsigned long)(*piter).second);
        line += num;
      }
    }

    printdebug ("%s\n", line.c_str());
  }
  printdebug("Total: %lu\n", (unsigned long)sets.size());
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
  {
    MBRange tmp = subtract( from, removed);
    from.swap( tmp );
  }
}

MBWriterIface* WriteHDF5Parallel::factory( MBInterface* iface )
  { return new WriteHDF5Parallel( iface ); }

void WriteHDF5Parallel::MultiProcSetTags::add( const std::string& name )
  { list.push_back( Data(name) ); }

void WriteHDF5Parallel::MultiProcSetTags::add( const std::string& filter, 
                                               const std::string& data )
  { list.push_back( Data(filter,data) ); }

void WriteHDF5Parallel::MultiProcSetTags::add( const std::string& filter, 
                                               int filterval,
                                               const std::string& data )
  { list.push_back( Data(filter,data,filterval) ); }


WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface)
    : WriteHDF5(iface), myPcomm(NULL), pcommAllocated(false)
{
  multiProcSetTags.add(  MATERIAL_SET_TAG_NAME );
  multiProcSetTags.add( DIRICHLET_SET_TAG_NAME );
  multiProcSetTags.add(   NEUMANN_SET_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 0, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 1, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 2, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 3, GLOBAL_ID_TAG_NAME );
}

WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface,
                                      const std::vector<std::string>& tag_names )
  : WriteHDF5(iface), myPcomm(NULL), pcommAllocated(false)
{
  for(std::vector<std::string>::const_iterator i = tag_names.begin();
      i != tag_names.end(); ++i)
    multiProcSetTags.add( *i );
}

WriteHDF5Parallel::WriteHDF5Parallel( MBInterface* iface,
                                      const MultiProcSetTags& set_tags )
    : WriteHDF5(iface), multiProcSetTags(set_tags), myPcomm(NULL), 
      pcommAllocated(false)
{}

WriteHDF5Parallel::~WriteHDF5Parallel()
{
  if (pcommAllocated && myPcomm) 
    delete myPcomm;
}

// The parent WriteHDF5 class has ExportSet structs that are
// populated with the entities to be written, grouped by type
// (and for elements, connectivity length).  This function:
//  o determines which entities are to be written by a remote processor
//  o removes those entities from the ExportSet structs in WriteMesh
//  o puts them in the 'interfaceMesh' array of MBRanges in this class
//  o sets their file Id to '1'
MBErrorCode WriteHDF5Parallel::gather_interface_meshes()
{
  MBErrorCode result;
  
  //START_SERIAL;
  printdebug( "Pre-interface mesh:\n");
  printrange(nodeSet.range);
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    printrange(eiter->range);
  printrange(setSet.range);
  
  MBRange iface_sets = myPcomm->interface_sets();

    // Populate lists of interface mesh entities
  MBRange tmpset;
  for (MBRange::iterator ifit = iface_sets.begin(); ifit != iface_sets.end(); ifit++) {
    int owner;
    result = myPcomm->get_owner(*ifit, owner);
    if (MB_SUCCESS != result || -1 == owner) return result;

    tmpset.clear();
    result = iFace->get_entities_by_handle(*ifit, tmpset, true);
    if (MB_SUCCESS != result) return result;
    interfaceMesh[owner].merge( tmpset );
  }
  
    // interfaceMesh currently contains interface entities for the 
    // entire mesh.  We now need to a) remove from the sets of entities
    // and handles that this proc doesn't own and b) remove from interfaceMesh
    // any handles for entities that we aren't going to write (on any proc.)
  
    // First move handles of non-owned entities from lists of entities
    // that this processor will write to the 'nonowned' list.
    
  MBRange nonowned;
  result = myPcomm->filter_pstatus( nodeSet.range, PSTATUS_NOT_OWNED, PSTATUS_AND, -1, &nonowned);
  if (MB_SUCCESS != result)
    return result;
  nodeSet.range = subtract( nodeSet.range, nonowned);
  
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
       eiter != exportList.end(); ++eiter ) {
    tmpset.clear();
    result = myPcomm->filter_pstatus( eiter->range, PSTATUS_NOT_OWNED, PSTATUS_AND, -1, &tmpset);
    if (MB_SUCCESS != result)
      return result;
    eiter->range = subtract( eiter->range,  tmpset );
    nonowned.merge(tmpset);
  }
  
    // Now remove from interfaceMesh any entities that are not
    // in 'nonowned' because we aren't writing those entities
    // (on any processor.)
  for (proc_iter i = interfaceMesh.begin(); i != interfaceMesh.end(); ++i)
    if (i->first != myPcomm->proc_config().proc_rank())
      i->second = intersect( nonowned,  i->second );
  
    // For the 'interfaceMesh' list for this processor, just remove
    // entities we aren't writing.
  tmpset.clear();
  tmpset.merge( nodeSet.range );
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
       eiter != exportList.end(); ++eiter ) 
    tmpset.merge( eiter->range );
  MBRange& my_remote_mesh = interfaceMesh[myPcomm->proc_config().proc_rank()];
  my_remote_mesh = intersect( my_remote_mesh,  tmpset );
  
    // print some debug output summarizing what we've accomplished
  printdebug("Remote mesh:\n");
  for (proc_iter i = interfaceMesh.begin(); i != interfaceMesh.end(); ++i)
  {
    printdebug("  proc %u : %d\n", i->first, i->second.size());
    printrange( i->second );
  }

  printdebug( "Post-interface mesh:\n");
  printrange(nodeSet.range);
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    printrange(eiter->range);
  printrange(setSet.range);

  //END_SERIAL;
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::parallel_create_file( const char* filename,
                                            bool overwrite,
                                            const std::vector<std::string>& qa_records,
                                            const MBTag* user_tag_list,
                                            int user_tag_count,
                                            int dimension,
                                            int pcomm_no)
{
  MBErrorCode rval;
  int result;
  mhdf_Status status;

  myPcomm = MBParallelComm::get_pcomm(iFace, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new MBParallelComm(iFace);
    pcommAllocated = true;
  }

TPRINT("Gathering interface meshes");
  rval = gather_interface_meshes();
  if (MB_SUCCESS != rval) return rval;

    /**************** get tag names for sets likely to be shared ***********/
TPRINT("Getting shared entity sets");
  rval = get_sharedset_tags();
  if (MB_SUCCESS != rval) return rval;
  

    /**************** Create actual file and write meta info ***************/

  if (myPcomm->proc_config().proc_rank() == 0)
  {
      // create the file
    const char* type_names[MBMAXTYPE];
    memset( type_names, 0, MBMAXTYPE * sizeof(char*) );
    for (MBEntityType i = MBEDGE; i < MBENTITYSET; ++i)
      type_names[i] = MBCN::EntityTypeName( i );
   
TPRINT("call mhdf_createFile");
    filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, &status );
    if (!filePtr)
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return MB_FAILURE;
    }
    
TPRINT("call write_qa");
    rval = write_qa( qa_records );
    if (MB_SUCCESS != rval) return rval;
  }
  
  
     /**************** Create node coordinate table ***************/
 
TPRINT("creating node table");
  rval = create_node_table( dimension );
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create element tables ***************/

TPRINT("negotiating element types");
  rval = negotiate_type_list();
  if (MB_SUCCESS != rval) return rval;
TPRINT("creating element table");
  rval = create_element_tables();
  if (MB_SUCCESS != rval) return rval;
  
  
    /*************** Exchange file IDs *****************/

TPRINT("communicating file ids");
  rval = exchange_file_ids();
  if (MB_SUCCESS != rval) return rval;
 

    /**************** Create adjacency tables *********************/
  
TPRINT("creating adjacency table");
  rval = create_adjacency_tables();
  if (MB_SUCCESS != rval) return rval;
  
    /**************** Create meshset tables *********************/
  
TPRINT("creating meshset table");
  rval = create_meshset_tables();
  if (MB_SUCCESS != rval) return rval;
  
  
    /* Need to write tags for shared sets this proc is responsible for */
  
  MBRange parallel_sets;
  for (std::list<ParallelSet>::const_iterator psiter = parallelSets.begin();
       psiter != parallelSets.end(); ++psiter)
    if (psiter->description)
      parallel_sets.insert( psiter->handle );
  
  setSet.range.merge( parallel_sets );
  rval = gather_tags( user_tag_list, user_tag_count );
  if (MB_SUCCESS != rval)
    return rval;
  range_remove( setSet.range, parallel_sets );   
  

    /**************** Create tag data *********************/
  
  std::list<SparseTag>::iterator tag_iter;
  sort_tags_by_name();
  const int num_tags = tagList.size();
  
    // Construct vector (tag_counts) containing a pair of values for each
    // tag, where the first value in the pair is the number of entities on
    // this processor for which the tag has been set.  The second value is
    // zero for normal tags.  For variable-length tags it is the total number
    // of tag values set for all entities on this processor.
  std::vector<unsigned long> tag_offsets(2*num_tags), tag_counts(2*num_tags);
  std::vector<unsigned long>::iterator tag_off_iter = tag_counts.begin();
  for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
    int s;
    assert(tag_off_iter < tag_counts.begin() + tag_counts.size());
    *tag_off_iter = tag_iter->range.size();
    ++tag_off_iter;
    if (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s )) {
      unsigned long total_len;
      rval = get_tag_data_length( *tag_iter, total_len );
      if (MB_SUCCESS != rval)
        return rval;
      
      assert(tag_off_iter < tag_counts.begin() + tag_counts.size());
      *tag_off_iter = total_len;
      assert(total_len == *tag_off_iter);
    }
    else {
      *tag_off_iter = 0;
    }
    ++tag_off_iter;
  }
  
    // Populate proc_tag_offsets on root processor with the values from
    // tag_counts on each processor.
TPRINT("communicating tag metadata");
  printdebug("Exchanging tag data for %d tags.\n", num_tags);
  std::vector<unsigned long> proc_tag_offsets(2*num_tags*myPcomm->proc_config().proc_size());
  VALGRIND_CHECK_MEM_IS_DEFINED( &tag_counts[0], 2*num_tags*sizeof(long) );
  VALGRIND_MAKE_VEC_UNDEFINED( proc_tag_offsets );
  result = MPI_Gather( &tag_counts[0], 2*num_tags, MPI_UNSIGNED_LONG,
                 &proc_tag_offsets[0], 2*num_tags, MPI_UNSIGNED_LONG,
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Calculate the total counts over all processors (tag_counts)
    // and the offset at which each processor should begin writing
    // its data (proc_tag_offsets).  Both lists contain a pair of
    // values, where the first is the number of entities and the second
    // is either unused for fixed-length tags or the total data table
    // size for variable-length tags.
  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    assert((unsigned)(2*i+1) < tag_counts.size());
    tag_counts[2*i] = tag_counts[2*i+1] = 0;
    unsigned long next_offset = 0;
    unsigned long next_var_len_offset = 0;
    for (unsigned int j = 0; j < myPcomm->proc_config().proc_size(); j++)
    {
      unsigned idx = 2*i + j*2*num_tags;
      assert(idx < proc_tag_offsets.size());
      unsigned long count = proc_tag_offsets[idx];
      proc_tag_offsets[idx] = next_offset;
      next_offset += count;
      assert((unsigned)(2*i) < tag_counts.size());
      tag_counts[2*i] += count;
      
      ++idx;
      assert(idx < proc_tag_offsets.size());
      count = proc_tag_offsets[idx];
      proc_tag_offsets[idx] = next_var_len_offset;
      next_var_len_offset += count;
      assert((unsigned)(2*i + 1) < tag_counts.size());
      tag_counts[2*i + 1] += count;
    }

    if (0 == myPcomm->proc_config().proc_rank())
    {
      rval = create_tag(tag_iter->tag_id, next_offset, next_var_len_offset);
      if (MB_SUCCESS != rval)
        return rval;
      printdebug( "Creating table of size %lu for tag 0x%lx\n", 
                  next_var_len_offset ? next_var_len_offset : next_offset, 
                  (unsigned long)tag_iter->tag_id );
    }
  }
  
    // Send total counts to all processors.  This is necessary because all 
    // processors need to know if we are not writing anything for the tag (count == 0).  
  if (myPcomm->rank() == 0) {
    VALGRIND_CHECK_MEM_IS_DEFINED( &tag_counts[0], 2*num_tags*sizeof(long) );
  }
  result = MPI_Bcast( &tag_counts[0], 2*num_tags, MPI_UNSIGNED_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Send to each processor its per-tag offset values.
  if (myPcomm->rank() == 0) {
    VALGRIND_CHECK_MEM_IS_DEFINED(  &proc_tag_offsets[0], proc_tag_offsets.size()*sizeof(long) );
  }
  result = MPI_Scatter( &proc_tag_offsets[0], 2*num_tags, MPI_UNSIGNED_LONG,
                             &tag_offsets[0], 2*num_tags, MPI_UNSIGNED_LONG,
                             0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);


  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    assert((unsigned)(2*i) < tag_counts.size());
    tag_iter->offset = tag_offsets[2*i];
    tag_iter->write = tag_counts[2*i] > 0;
    tag_iter->varDataOffset = tag_offsets[2*i + 1];
  }

  #ifdef DEBUG
  START_SERIAL;  
  printdebug("Tags: %12s %8s %8s %8s %8s %8s\n", "Name", "Count", "Offset", "Var Off", "Var Len", "Handle");

  tag_iter = tagList.begin();
  for (int i = 0; i < num_tags; ++i, ++tag_iter)
  {
    std::string name;
    iFace->tag_get_name( tag_iter->tag_id, name );
    printdebug("%18s %8lu %8lu %8lu %8lu 0x%7lx\n", name.c_str(), tag_counts[2*i], tag_offsets[2*i], tag_offsets[2*i+1], tag_counts[2*i+1], (unsigned long)tag_iter->tag_id );
  }
  END_SERIAL;  
  #endif
  
  /************** Close serial file and reopen parallel *****************/
  
  if (0 == myPcomm->proc_config().proc_rank())
  {
    mhdf_closeFile( filePtr, &status );
  }
  
TPRINT("(re)opening file in parallel mode");
  unsigned long junk;
  hid_t hdf_opt = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_fapl_mpio( hdf_opt, myPcomm->proc_config().proc_comm(), MPI_INFO_NULL );
  filePtr = mhdf_openFileWithOpt( filename, 1, &junk, hdf_opt, &status );
  H5Pclose( hdf_opt );
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
  std::vector<long> node_counts(myPcomm->proc_config().proc_size());
  long num_nodes = nodeSet.range.size();
  VALGRIND_CHECK_MEM_IS_DEFINED( &num_nodes, sizeof(long) );
  result = MPI_Gather( &num_nodes, 1, MPI_LONG, &node_counts[0], 1, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // create node data in file
  long first_id;
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    int total = 0;
    for (unsigned int i = 0; i < myPcomm->proc_config().proc_size(); i++)
      total += node_counts[i];
      
    hid_t handle = mhdf_createNodeCoords( filePtr, dimension, total, &first_id, &status );
    CHECK_HDF(status);
    mhdf_closeData( filePtr, handle, &status );
  }
    
    // send id offset to every proc
  result = MPI_Bcast( &first_id, 1, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  nodeSet.first_id = (id_t)first_id;
   
      // calculate per-processor offsets
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    int prev_size = node_counts[0];
    node_counts[0] = 0;
    for (unsigned int i = 1; i < myPcomm->proc_config().proc_size(); ++i)
    {
      int mysize = node_counts[i];
      node_counts[i] = node_counts[i-1] + prev_size;
      prev_size = mysize;
    }
  }
  
    // send each proc it's offset in the node table
  long offset;
  if (myPcomm->rank() == 0) {
    VALGRIND_CHECK_MEM_IS_DEFINED( &node_counts[0], sizeof(long) );
  }
  result = MPI_Scatter( &node_counts[0], 1, MPI_LONG, 
                        &offset, 1, MPI_LONG,
                        0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  nodeSet.offset = offset;

  return assign_ids( nodeSet.range, nodeSet.first_id + nodeSet.offset );
}



struct elemtype {
  int mbtype;
  int numnode;
  
  elemtype( int vals[2] ) : mbtype(vals[0]), numnode(vals[1]) {}
  elemtype( int t, int n ) : mbtype(t), numnode(n) {}
  
  bool operator==( const elemtype& other ) const
  {
    return mbtype == other.mbtype &&
            (mbtype == MBENTITYSET ||
             numnode == other.numnode);
  }
  bool operator<( const elemtype& other ) const
  {
    if (mbtype > other.mbtype)
      return false;
   
    return mbtype < other.mbtype ||
           (mbtype != MBENTITYSET &&
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
  std::vector<int> counts(myPcomm->proc_config().proc_size());
  result = MPI_Gather( &num_types, 1, MPI_INT, &counts[0], 1, MPI_INT, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Get list of types on this processor
  std::vector<int> my_types(num_types);
  VALGRIND_MAKE_VEC_UNDEFINED( my_types );
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
  std::vector<int> displs(myPcomm->proc_config().proc_size() + 1);
  VALGRIND_MAKE_VEC_UNDEFINED( displs );
  displs[0] = 0;
  for (unsigned long i = 1; i <= myPcomm->proc_config().proc_size(); ++i)
    displs[i] = displs[i-1] + counts[i-1];
  int total = displs[myPcomm->proc_config().proc_size()];
  std::vector<int> alltypes(total);
  VALGRIND_MAKE_VEC_UNDEFINED( alltypes );
  VALGRIND_CHECK_MEM_IS_DEFINED( &my_types[0], my_types.size()*sizeof(int) );
  result = MPI_Gatherv( &my_types[0], my_types.size(), MPI_INT,
                        &alltypes[0], &counts[0], &displs[0], MPI_INT,
                        0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Merge type lists
  std::list<elemtype> type_list;
  std::list<elemtype>::iterator liter;
  for (unsigned int i = 0; i < myPcomm->proc_config().proc_size(); ++i)
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
  result = MPI_Bcast( &total, 1, MPI_INT, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Send list of types to each processor
  std::vector<int> intlist(total * 2);
  VALGRIND_MAKE_VEC_UNDEFINED( intlist );
  viter = intlist.begin();
  for (liter = type_list.begin(); liter != type_list.end(); ++liter)
  {
    *viter = liter->mbtype;  ++viter;
    *viter = liter->numnode; ++viter;
  }
  result = MPI_Bcast( &intlist[0], 2*total, MPI_INT, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);

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
  for (long i = 0; i < total; ++i)
  {
    long mbtype = *viter; ++viter;
    long numnode = *viter; ++viter;
    while (ex_iter != exportList.end() && ex_iter->type < mbtype)
      ++ex_iter;
    
    bool equal = ex_iter != exportList.end() && ex_iter->type == mbtype;
    if (equal)
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
  std::vector<long> counts(numtypes * myPcomm->proc_config().proc_size() + numtypes);
  VALGRIND_MAKE_VEC_UNDEFINED( my_counts );
  VALGRIND_MAKE_VEC_UNDEFINED( counts );
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { *viter = ex_iter->range.size(); ++viter; }
  
  VALGRIND_CHECK_MEM_IS_DEFINED( &my_counts[0], numtypes*sizeof(long) );
  result = MPI_Gather( &my_counts[0], numtypes, MPI_LONG,
                       &counts[0],    numtypes, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Convert counts to offsets
  for (int i = 0; i < numtypes; i++) 
  {
    long prev = 0;
    for (unsigned int j = 0; j <= myPcomm->proc_config().proc_size(); j++)
    {
      long tmp = counts[j*numtypes + i];
      counts[j*numtypes+i] = prev;
      prev += tmp;
    }
  }
  
    // Send offsets to each processor
  if (myPcomm->rank() == 0) {
    VALGRIND_CHECK_MEM_IS_DEFINED( &counts[0], counts.size()*sizeof(long) );
  }
  result = MPI_Scatter( &counts[0],    numtypes, MPI_LONG,
                        &my_counts[0], numtypes, MPI_LONG,
                        0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Update store offsets in ExportSets
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    ex_iter->offset = (id_t)*(viter++);
  
    // Create element tables
  std::vector<long> start_ids(numtypes);
  VALGRIND_MAKE_VEC_UNDEFINED( start_ids );
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    viter = start_ids.begin();
    long* citer = &counts[numtypes * myPcomm->proc_config().proc_size()];
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    {
      rval = create_elem_tables( ex_iter->type,
                                 ex_iter->num_nodes,
                                 *citer,
                                 *viter );
      if (MB_SUCCESS != rval)
        return rval;
      ++citer;
      ++viter;
    }
  }
  
    // send start IDs to each processor
  result = MPI_Bcast( &start_ids[0], numtypes, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Assign IDs to local elements
  viter = start_ids.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    ex_iter->first_id = *(viter++);
    id_t myfirst = (id_t)(ex_iter->first_id + ex_iter->offset);
    rval = assign_ids( ex_iter->range, myfirst );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
  return MB_SUCCESS;
}
  
MBErrorCode WriteHDF5Parallel::create_adjacency_tables()
{
  MBErrorCode rval;
  mhdf_Status status;
  unsigned int j;
  int i, result;

#ifdef WRITE_NODE_ADJACENCIES  
  const int numtypes = exportList.size()+1;
#else
  const int numtypes = exportList.size();
#endif
  std::vector<long>::iterator viter;
  std::list<ExportSet>::iterator ex_iter;
  std::vector<long> local(numtypes), all(myPcomm->proc_config().proc_size() * numtypes + numtypes);
  VALGRIND_MAKE_VEC_UNDEFINED( local );
  VALGRIND_MAKE_VEC_UNDEFINED( all );
  
    // Get adjacency counts for local processor
  viter = local.begin();
  id_t num_adj;
#ifdef WRITE_NODE_ADJACENCIES  
  rval = count_adjacencies( nodeSet.range, num_adj );
  assert (MB_SUCCESS == rval);
  *viter = num_adj; ++viter;
#endif

  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
  {
    rval = count_adjacencies( ex_iter->range, num_adj );
    if (MB_SUCCESS != rval)
      return rval;
    *viter = num_adj; ++viter;
  }
  
    // Send local adjacency counts to root processor
  VALGRIND_CHECK_MEM_IS_DEFINED( &local[0], numtypes*sizeof(long) );
  result = MPI_Gather( &local[0], numtypes, MPI_LONG,
                       &all[0],   numtypes, MPI_LONG, 
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
  
  if (myPcomm->rank() == 0) {
      // Convert counts to offsets
    for (i = 0; i < numtypes; i++) 
    {
      long prev = 0;
      for (unsigned j = 0; j <= myPcomm->proc_config().proc_size(); j++)
      {
        long tmp = all[j*numtypes + i];
        all[j*numtypes+i] = prev;
        prev += tmp;
      }
    }

      // For each element type for which there is no adjacency data,
      // send -1 to all processors as the offset
    for (i = 0; i < numtypes; ++i)
      if (all[numtypes*myPcomm->proc_config().proc_size()+i] == 0)
        for (j = 0; j < myPcomm->proc_config().proc_size(); ++j)
          all[j*numtypes+i] = -1;
  }
  
    // Send offsets back to each processor
  result = MPI_Scatter( &all[0],   numtypes, MPI_LONG,
                        &local[0], numtypes, MPI_LONG,
                        0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
    // Record the adjacency offset in each ExportSet
  viter = local.begin();
#ifdef WRITE_NODE_ADJACENCIES  
  nodeSet.adj_offset = *viter; ++viter;
#endif
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    { ex_iter->adj_offset = *viter; ++viter; }
  
    // Create data tables in file
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    viter = all.begin() + (numtypes * myPcomm->proc_config().proc_size());
#ifdef WRITE_NODE_ADJACENCIES  
    if (*viter) {
      hid_t handle = mhdf_createAdjacency( filePtr, 
                                           mhdf_node_type_handle(),
                                           *viter,
                                           &status );
      CHECK_HDF(status);
      mhdf_closeData( filePtr, handle, &status );
    }
    ++viter;
#endif
    for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter, ++viter)
    {
      if (!*viter) 
        continue;
      
      hid_t handle = mhdf_createAdjacency( filePtr,
                                           ex_iter->name(),
                                           *viter,
                                           &status );
      CHECK_HDF(status);
      mhdf_closeData( filePtr, handle, &status );
    }
  }

  return MB_SUCCESS;
}

/*
MBErrorCode WriteHDF5Parallel::get_interface_set_data( RemoteSetData& data,
                                                       long& offset )
{
  const char* PROC_ID_TAG = "HDF5Writer_Rank";
  MBTag iface_tag, proc_tag;
  MBErrorCode rval;
  
  rval = iFace->tag_get_handle( PARALLEL_INTERFACE_TAG_NAME, iface_tag );
  if (MB_SUCCESS != rval) return rval;
  
  rval = iFace->tag_get_handle( PROC_ID_TAG, proc_tag );
  if (MB_SUCCESS == rval) 
    iFace->tag_delete( proc_tag );
  rval = iFace->tag_create( PROC_ID_TAG, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, proc_tag, 0 );
  if (MB_SUCCESS != rval) return rval;
    
  MBRange interface_sets, sets;
  rval = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET, &iface_tag, 0, 1, interface_sets );
  if (MB_SUCCESS != rval) return rval;
  
  std::vector<int> list;
  for (MBRange::iterator i = interface_sets.begin(); i != interface_sets.end(); ++i)
  {
    int proc_ids[2];
    rval = iFace->tag_get_data( iface_tag, &*i, 1, proc_ids );
    if (MB_SUCCESS != rval) return rval;
    
    sets.clear();
    rval = iFace->get_entities_by_type( *i, MBENTITYSET, sets );
    if (MB_SUCCESS != rval) return rval;
  
    list.clear();
    list.resize( sets.size(), proc_ids[0] );
    rval = iFace->tag_set_data( proc_tag, sets, &list[0] );
    if (MB_SUCCESS != rval) return rval;
  }
  
  return get_remote_set_data( PROC_ID_TAG, PARALLEL_GLOBAL_ID_TAG_NAME, data, offset );
}
*/
  
/** Working data for group of global sets identified by an ID tag */
struct RemoteSetData {
  MBTag data_tag;    //!< The ID tag for matching sets across processors
  MBTag filter_tag;  //!< Optional tag to filter on (e.g. geometric dimension)
  int filter_value;  //!< Value of filter_tag for this group
  MBRange range;     //!< Set handles with data_tag set (and optionally filter_tag == filter_value)
  std::vector<int> counts;       //!< Number of sets with tag on each proc, indexed by MPI rank
  std::vector<int> displs;       //!< Offset in all_values at which the data_tag values for
                                 //!< each processor begin. displs[n] = sum(i from 0 to n-1)(counts[i])
  std::vector<int> all_values;   //!< data_tag values for sets on all processors, 
                                 //!< counts[0] values for proc 0, then counts[1] values for proc 
                                 //!< 1, etc.
  std::vector<int> local_values; //!< data_tag values for sets that exist on this processor
};

static void print_remote_set_data( std::ostream& s, MBInterface* mb,
                            const struct RemoteSetData& data,
                            const char* pfx = "" )
{
  std::string n1,n2;
  if (data.data_tag)
    mb->tag_get_name( data.data_tag, n1 );
  if (data.filter_tag)
    mb->tag_get_name( data.filter_tag, n2 );

  std::cerr << pfx << "data_tag:     " << n1 << std::endl
            << pfx << "filter_tag:   " << n2 << std::endl
            << pfx << "filter_value: " << data.filter_value << std::endl
            << pfx << "range:        " << data.range << std::endl
            << pfx << "counts:       ";
  std::copy( data.counts.begin(), data.counts.end(),
             std::ostream_iterator<int>(std::cerr,",") );
  std::cerr << std::endl << pfx << "displs:       ";
  std::copy( data.displs.begin(), data.displs.end(),
             std::ostream_iterator<int>(std::cerr,",") );
  std::cerr << std::endl << pfx << "all_values:  ";
  std::copy( data.all_values.begin(), data.all_values.end(),
             std::ostream_iterator<int>(std::cerr,",") );
  std::cerr << std::endl << pfx << "local_values:";
  std::copy( data.local_values.begin(), data.local_values.end(),
             std::ostream_iterator<int>(std::cerr,",") );
  std::cerr << std::endl;
}

MBErrorCode WriteHDF5Parallel::get_remote_set_data( 
                        const WriteHDF5Parallel::MultiProcSetTags::Data& tags,
                        RemoteSetData& data, long& offset )
{
  MBErrorCode rval;
  int i;
  int result;
  MBRange::iterator riter;

  rval = iFace->tag_get_handle( tags.filterTag.c_str(), data.filter_tag );
  if (rval != MB_SUCCESS) return rval;
  if (tags.useFilterValue) 
  {
    i = 0;
    iFace->tag_get_size( data.filter_tag, i );
    if (i != sizeof(int)) {
      fprintf(stderr, "Cannot use non-int tag data for filtering remote sets.\n" );
      assert(0);
      return MB_FAILURE;
    }  
    data.filter_value = tags.filterValue;
  }
  else
  {
    data.filter_value = 0;
  }
  
  rval = iFace->tag_get_handle( tags.dataTag.c_str(), data.data_tag );
  if (rval != MB_SUCCESS) return rval;
  i = 0;
  iFace->tag_get_size( data.data_tag, i );
  if (i != sizeof(int)) {
    fprintf(stderr, "Cannot use non-int tag data for matching remote sets.\n" );
    assert(0);
    return MB_FAILURE;
  }  
    

  printdebug("Negotiating multi-proc meshsets for tag: \"%s\"\n", tags.filterTag.c_str());

    // Get sets with tag, or leave range empty if the tag
    // isn't defined on this processor.
  if (rval != MB_TAG_NOT_FOUND)
  {
    MBTag handles[] = { data.filter_tag, data.data_tag };
    const void* values[] = { tags.useFilterValue ? &tags.filterValue : 0, 0 };
    rval = iFace->get_entities_by_type_and_tag( 0, 
                                                MBENTITYSET, 
                                                handles,
                                                values,
                                                2,
                                                data.range );
    if (rval != MB_SUCCESS) return rval;
    MBRange tmp = intersect( data.range,  setSet.range );
    data.range.swap( tmp );
    range_remove( setSet.range, data.range );
  }
  
  printdebug("Found %d meshsets with \"%s\" tag.\n", data.range.size(), tags.filterTag.c_str() );

    // Exchange number of sets with tag between all processors
  data.counts.resize(myPcomm->proc_config().proc_size());
  VALGRIND_MAKE_VEC_UNDEFINED( data.counts );
  int count = data.range.size();
  result = MPI_Allgather( &count,          1, MPI_INT, 
                          &data.counts[0], 1, MPI_INT,
                          myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);

    // Exchange tag values for sets between all processors
  data.displs.resize(myPcomm->proc_config().proc_size()+1);
  VALGRIND_MAKE_VEC_UNDEFINED( data.displs );
  data.displs[0] = 0;
  for (unsigned int j = 1; j <= myPcomm->proc_config().proc_size(); j++)
    data.displs[j] = data.displs[j-1] + data.counts[j-1];
  int total = data.displs[myPcomm->proc_config().proc_size()];
  data.all_values.resize(total);
  VALGRIND_MAKE_VEC_UNDEFINED( data.all_values );
  data.local_values.resize(count);
  VALGRIND_MAKE_VEC_UNDEFINED( data.local_values );
  rval = iFace->tag_get_data( data.data_tag, data.range, &data.local_values[0] );
  if (MB_SUCCESS != rval)
    return rval;
  VALGRIND_CHECK_MEM_IS_DEFINED( &data.local_values[0], count*sizeof(int) );
  result = MPI_Allgatherv( &data.local_values[0], count, MPI_INT,
                           &data.all_values[0], &data.counts[0], &data.displs[0], MPI_INT,
                           myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);


    // Remove from the list any sets that are unique to one processor
  std::vector<int> sorted( data.all_values );
  std::sort( sorted.begin(), sorted.end() );
  int r = 0, w = 0;
  for (unsigned j = 0; j < myPcomm->proc_config().proc_size(); ++j)
  {
    const int start = w;
    for (int i = 0; i < data.counts[j]; ++i)
    {
      std::vector<int>::iterator p 
        = std::lower_bound( sorted.begin(), sorted.end(), data.all_values[r] );
      ++p;
      if (p != sorted.end() && *p == data.all_values[r])
      {
        data.all_values[w] = data.all_values[r];
        ++w;
      }
      ++r;
    }
    data.counts[j] = w - start;
  }
  total = w;
  data.all_values.resize( total );
  r = w = 0;
  for (i = 0; i < count; ++i)
  {
    std::vector<int>::iterator p 
      = std::lower_bound( sorted.begin(), sorted.end(), data.local_values[r] );
    ++p;
    if (p != sorted.end() && *p == data.local_values[r])
    {
      data.local_values[w] = data.local_values[r];
      ++w;
    }
    else
    {
      riter = data.range.begin();
      riter += w;
      setSet.range.insert( *riter );
      data.range.erase( riter );
    }
    ++r;
  }
  count = data.range.size();
  assert( count == data.counts[myPcomm->proc_config().proc_rank()] );
  assert( count == w );
  data.local_values.resize( count );
  sorted.clear(); // release storage
  
    // recalculate displacements for updated counts
  data.displs[0] = 0;
  for (unsigned int j = 1; j <= myPcomm->proc_config().proc_size(); j++)
    data.displs[j] = data.displs[j-1] + data.counts[j-1];
  
    // Find sets that span multple processors and update appropriately.
    // The first processor (sorted by MPI rank) that contains a given set
    // will be responsible for writing the set description.  All multi-
    // processor sets will be written at the beginning of the set tables.
    // Processors will write set contents/children for a given set in
    // the order of their MPI rank.
    //
    // Identify which meshsets will be managed by this processor and
    // the corresponding offset in the set description table. 
  int cpu = 0;
  std::set<int> seen;
  for (i = 0; i < total; ++i)
  {
    if (data.displs[cpu+1] == i)
      ++cpu;

    if (seen.insert(data.all_values[i]).second)
      ++offset;

    std::vector<int>::iterator loc 
      = std::find( data.local_values.begin(), data.local_values.end(), data.all_values[i] );
    if (loc != data.local_values.end()) 
    {
      riter = data.range.begin();
      riter += loc - data.local_values.begin();
      cpuParallelSets[cpu].insert( *riter );
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::set_shared_set_ids( RemoteSetData& data, long& offset )
{
    // Determine File ID for each shared set.
  std::map<int,long> val_id_map; // Map from tag value to file ID for set
  for (size_t i = 0; i < data.all_values.size(); ++i)
  {
    std::map<int,long>::iterator p = val_id_map.find( data.all_values[i] );
    if (p == val_id_map.end())
      val_id_map[data.all_values[i]] = offset++;
  }
  MBRange::const_iterator riter = data.range.begin();
  for (size_t i = 0; i < data.local_values.size(); ++i, ++riter)
  {
    std::map<int,long>::iterator p = val_id_map.find( data.local_values[i] );
    assert( p != val_id_map.end() );
    long id = p->second;
    
    if (idMap.end() == idMap.insert( *riter, id, 1 )) {
      for (unsigned x = 0; x < myPcomm->size(); ++x) {
        MPI_Barrier( myPcomm->proc_config().proc_comm() );      
        if (x != myPcomm->rank()) continue;   

        std::ostringstream s;
        s << "[" << myPcomm->rank() << "] ";
        std::string pfx1 = s.str();
        s << "  ";
        std::string pfx2 = s.str();

        std::cerr << pfx1 << "Duplicate shared set handle or internal accounting error" << std::endl;
        std::cerr << pfx1 << "RemoteSetData:  " << std::endl;
        print_remote_set_data( std::cerr, iFace, data, pfx2.c_str() );

        std::cerr << pfx1 << "val_id_map: " << std::endl;
        for (p = val_id_map.begin(); p != val_id_map.end(); ++p)
          std::cerr << pfx2 << p->first << "->" << p->second << std::endl;

        std::cerr << pfx1 << "idMap: " << std::endl;
        print_id_map( std::cerr, pfx2.c_str() );

        std::cerr << pfx1 << "Failed at: (" << i << ") " << data.local_values[i] 
                  << "->" << id << " for " 
                  << MBCN::EntityTypeName(TYPE_FROM_HANDLE(*riter)) 
                  << " " << ID_FROM_HANDLE(*riter) << std::endl;
        std::cerr.flush();
      }
      
      assert(false);
      return MB_FAILURE;
    }
  }
  
  return MB_SUCCESS;
}
  

MBErrorCode WriteHDF5Parallel::create_meshset_tables()
{
  MBErrorCode rval = MB_SUCCESS;
  int result, i;
  long total_offset = 0;
  MBRange::const_iterator riter;

  START_SERIAL;
  print_type_sets( iFace, myPcomm->proc_config().proc_rank(), myPcomm->proc_config().proc_size(), setSet.range );
  END_SERIAL;

    // Gather data about multi-processor meshsets - removes sets from setSet.range
  std::vector<RemoteSetData> remote_set_data( multiProcSetTags.list.size() );
  for (i = 0; i< (int)multiProcSetTags.list.size(); i++)
  {
    rval = get_remote_set_data( multiProcSetTags.list[i],
                                remote_set_data[i],
                                total_offset ); 
    if (MB_SUCCESS != rval)
      return rval;
  }

  START_SERIAL;
  printdebug("myLocalSets\n");
  print_type_sets( iFace, myPcomm->proc_config().proc_rank(), myPcomm->proc_config().proc_size(), setSet.range );
  END_SERIAL;

    // Gather counts of non-shared sets from each proc
    // to determine total table size.
  std::vector<long> set_offsets(myPcomm->proc_config().proc_size() + 1);
  VALGRIND_MAKE_VEC_UNDEFINED( set_offsets );
  long local_count = setSet.range.size();
  result = MPI_Gather( &local_count,    1, MPI_LONG,
                       &set_offsets[0], 1, MPI_LONG,
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  for (unsigned int j = 0; j <= myPcomm->proc_config().proc_size(); j++)
  {
    long tmp = set_offsets[j];
    set_offsets[j] = total_offset;
    total_offset += tmp;
  }
  
    // Send each proc its offsets in the set description table.
  long sets_offset;
  result = MPI_Scatter( &set_offsets[0], 1, MPI_LONG,
                        &sets_offset,    1, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setSet.offset = (id_t)(sets_offset);

    // Create the set description table
  long total_count_and_start_id[2] = { set_offsets[myPcomm->proc_config().proc_size()], 0 };
  if (myPcomm->proc_config().proc_rank() == 0 && total_count_and_start_id[0] > 0)
  {
    rval = create_set_meta( (id_t)total_count_and_start_id[0], total_count_and_start_id[1] );
    if (MB_SUCCESS != rval)
      return rval;
  }
  
    // Send totals to all procs.
  result = MPI_Bcast( total_count_and_start_id, 2, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setSet.first_id = total_count_and_start_id[1];
  writeSets = total_count_and_start_id[0] > 0;

  START_SERIAL;  
  printdebug("Non-shared sets: %ld local, %ld global, offset = %ld, first_id = %ld\n",
    local_count, total_count_and_start_id[0], sets_offset, total_count_and_start_id[1] );
  printdebug("my Parallel Sets:\n");
  print_type_sets(iFace, myPcomm->proc_config().proc_rank(), myPcomm->proc_config().proc_size(), cpuParallelSets[myPcomm->proc_config().proc_rank()] );
  END_SERIAL;
  
    // Not writing any sets??
  if (!writeSets)
    return MB_SUCCESS;
  
  long start_id = setSet.first_id;
  for (i = 0; i < (int)remote_set_data.size(); ++i)
    set_shared_set_ids( remote_set_data[i], start_id );
  
    // Assign set IDs
  assign_ids( setSet.range, setSet.first_id + setSet.offset );
  
    // Communicate sizes for remote sets
  long data_offsets[3] = { 0, 0, 0 };
  for (i = 0; i < (int)remote_set_data.size(); ++i)
  {
    rval = negotiate_remote_set_contents( remote_set_data[i], data_offsets ); 
    if (MB_SUCCESS != rval)
      return rval;
  }
  remote_set_data.clear();
  
    // Exchange IDs for remote/adjacent sets not shared between procs
  //rval = communicate_remote_ids( MBENTITYSET ); assert(MB_SUCCESS == rval);
  
    // Communicate counts for local sets
  long data_counts[3];
  rval = count_set_size( setSet.range, data_counts[0], data_counts[1], data_counts[2] );
  if (MB_SUCCESS != rval) 
    return rval;
  std::vector<long> set_counts(3*myPcomm->proc_config().proc_size());
  VALGRIND_MAKE_VEC_UNDEFINED( set_counts );
  result = MPI_Gather( data_counts,    3, MPI_LONG,
                       &set_counts[0], 3, MPI_LONG,
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
    
  for (unsigned int j = 0; j < 3*myPcomm->proc_config().proc_size(); ++j)
  {
    long tmp = set_counts[j];
    set_counts[j] = data_offsets[j%3];
    data_offsets[j%3] += tmp;
  }
  long all_counts[] = {data_offsets[0], data_offsets[1], data_offsets[2]};
  result = MPI_Scatter( &set_counts[0], 3, MPI_LONG,
                        data_offsets,   3, MPI_LONG,
                        0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setContentsOffset = data_offsets[0];
  setChildrenOffset = data_offsets[1];
  setParentsOffset = data_offsets[2];
  
    // Create set contents and set children tables
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    rval = create_set_tables( all_counts[0], all_counts[1], all_counts[2] );
    if (MB_SUCCESS != rval) 
      return rval;
  }
  
    // Send totals to all processors
  result = MPI_Bcast( all_counts, 3, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  writeSetContents = all_counts[0] > 0;
  writeSetChildren = all_counts[1] > 0;
  writeSetParents  = all_counts[2] > 0;

  START_SERIAL;  
  printdebug("Non-shared set contents: %ld local, %ld global, offset = %ld\n",
    data_counts[0], all_counts[0], data_offsets[0] );
  printdebug("Non-shared set children: %ld local, %ld global, offset = %ld\n",
    data_counts[1], all_counts[1], data_offsets[1] );
  printdebug("Non-shared set parents: %ld local, %ld global, offset = %ld\n",
    data_counts[2], all_counts[2], data_offsets[2] );
  END_SERIAL;
  
  return MB_SUCCESS;
}

void WriteHDF5Parallel::remove_remote_entities( MBEntityHandle relative,
                                                MBRange& range )
{
  MBRange result;
  result.merge( intersect( range,  nodeSet.range ) );
  result.merge( intersect( range,  setSet.range ) );  
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
  {
    result.merge( intersect( range, eiter->range ) );
  }
  //result.merge( intersect( range, myParallelSets ) );
  MBRange sets;
  int junk;
  sets.merge( MBRange::lower_bound( range.begin(), range.end(), CREATE_HANDLE(MBENTITYSET, 0, junk )), range.end() );
  remove_remote_sets( relative, sets );
  result.merge( sets );
  range.swap(result);
}

void WriteHDF5Parallel::remove_remote_sets( MBEntityHandle relative, 
                                            MBRange& range )
{
  MBRange result( intersect( range,  setSet.range ) );
  //result.merge( intersect( range, yParallelSets ) );
  MBRange remaining( subtract( range, result ) );
  
  for(MBRange::iterator i = remaining.begin(); i != remaining.end(); ++i)
  {
      // Look for the first CPU which knows about both sets.
    proc_iter cpu;
    for (cpu = cpuParallelSets.begin(); cpu != cpuParallelSets.end(); ++cpu)
      if (cpu->second.find(relative) != cpu->second.end() &&
          cpu->second.find(*i) != cpu->second.end())
        break;
      // If we didn't find one, it may indicate a bug.  However,
      // it could also indicate that it is a link to some set that
      // exists on this processor but is not being written, because
      // the caller requested that some subset of the mesh be written.
    //assert(cpu < myPcomm->proc_config().proc_size());
      // If I'm the first set that knows about both, I'll handle it.
    if (cpu->first == myPcomm->proc_config().proc_rank())
      result.insert( *i );
  }
  
  range.swap( result );
}
  
  

void WriteHDF5Parallel::remove_remote_entities( MBEntityHandle relative,
                                                std::vector<MBEntityHandle>& vect )
{
  MBRange intrsct;
  for (std::vector<MBEntityHandle>::const_iterator iter = vect.begin();
       iter != vect.end(); ++iter)
    intrsct.insert(*iter);
  remove_remote_entities( relative, intrsct );
  
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

  

void WriteHDF5Parallel::remove_remote_sets( MBEntityHandle relative,
                                            std::vector<MBEntityHandle>& vect )
{
  MBRange intrsct;
  for (std::vector<MBEntityHandle>::const_iterator iter = vect.begin();
       iter != vect.end(); ++iter)
    intrsct.insert(*iter);
  remove_remote_sets( relative, intrsct );
  
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

// Given a RemoteSetData object describing the set information for a 
// single tag (or tag pair), populate the list of parallel sets
// (this->parallelSets) with the per-entityset data.
MBErrorCode WriteHDF5Parallel::negotiate_remote_set_contents( RemoteSetData& data,
                                                              long* offsets /* long[3] */ )
{
  unsigned i;
  MBErrorCode rval;
  MBRange::const_iterator riter;
  int result;
  const unsigned count = data.range.size();
  const unsigned total = data.all_values.size();
  std::vector<int>::iterator viter, viter2;

    // Calculate counts for each meshset
  std::vector<long> local_sizes(3*count);
  VALGRIND_MAKE_VEC_UNDEFINED( local_sizes );
  std::vector<long>::iterator sizes_iter = local_sizes.begin();
  MBRange tmp_range;
  std::vector<MBEntityHandle> child_list;
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
      // Count contents
    *sizes_iter = 0;
    tmp_range.clear();
    rval = iFace->get_entities_by_handle( *riter, tmp_range );
    if (MB_SUCCESS != rval)
      return rval;
    remove_remote_entities( *riter, tmp_range );
    for (MBRange::iterator iter = tmp_range.begin(); iter != tmp_range.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    ++sizes_iter;
    
      // Count children
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_child_meshsets( *riter, child_list );
    if (MB_SUCCESS != rval)
      return rval;
    remove_remote_sets( *riter, child_list );
    for (std::vector<MBEntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    ++sizes_iter;
    
      // Count parents
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_parent_meshsets( *riter, child_list );
    if (MB_SUCCESS != rval)
      return rval;
    remove_remote_sets( *riter, child_list );
    for (std::vector<MBEntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    ++sizes_iter;
  }
  
    // Exchange sizes for sets between all processors.
  std::vector<long> all_sizes(3*total);
  VALGRIND_MAKE_VEC_UNDEFINED( all_sizes );
  std::vector<int> counts(myPcomm->proc_config().proc_size()), displs(myPcomm->proc_config().proc_size());
  VALGRIND_MAKE_VEC_UNDEFINED( counts );
  VALGRIND_MAKE_VEC_UNDEFINED( displs );
  for (i = 0; i < (unsigned)myPcomm->proc_config().proc_size(); i++)
    counts[i] = 3 * data.counts[i];
  displs[0] = 0;
  for (i = 1; i < (unsigned)myPcomm->proc_config().proc_size(); i++)
    displs[i] = displs[i-1] + counts[i-1];
  VALGRIND_CHECK_MEM_IS_DEFINED( &local_sizes[0], 3*count*sizeof(long) );
  result = MPI_Allgatherv( &local_sizes[0], 3*count, MPI_LONG,
                           &all_sizes[0], &counts[0], &displs[0], MPI_LONG,
                           myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);

  
    // Update information in-place in the array from the Allgatherv.
    
    // Change the corresponding sizes for the first instance of a tag
    // value such that it ends up being the total size of the set.
    // Change the size to -1 for the later instances of a tag value.
    //
    // For the sets that this processor has, update the offsets at
    // which the set data is to be written.  Store the offset of the data
    // on this processor for the set *relative* to the start of the
    // data of *the set*.
  std::vector<long> local_offsets(3*count);
  VALGRIND_MAKE_VEC_UNDEFINED( local_offsets );
  std::map<int,int> tagsort;  // Map of {tag value, index of first set w/ value}
  for (i = 0; i < total; ++i)
  {
    const std::map<int,int>::iterator p = tagsort.find( data.all_values[i] );
    const unsigned r = (unsigned)(i - data.displs[myPcomm->proc_config().proc_rank()]);  // offset in "local" array
    
      // If this is the first instance of this tag value, 
      // then the processor with this instance is responsible
      // for writing the tag description
    if ( p == tagsort.end() )  
    {
      tagsort[data.all_values[i]] = i;
        // If within the range for this processor, save offsets
      if (r < (unsigned)count) 
      {
        local_offsets[3*r] = local_offsets[3*r+1] = local_offsets[3*r+2] = 0;
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
          // the offset for this processor, from the start of the data
          // for this group of sets, is the current total count for the
          // group of sets.
        local_offsets[3*r  ] = all_sizes[3*j  ];  // contents
        local_offsets[3*r+1] = all_sizes[3*j+1];  // children
        local_offsets[3*r+2] = all_sizes[3*j+2];  // parents
      }
      
        // update the total count for the set in the first position in
        // all_sizes at which the set occurs (the one corresponding to
        // the processor that owns the set.)
      all_sizes[3*j  ] += all_sizes[3*i  ]; // contents
      all_sizes[3*j+1] += all_sizes[3*i+1]; // children
      all_sizes[3*j+2] += all_sizes[3*i+2]; // parents
        // set the size to -1 in the positions corresponding to non-owning processor
      all_sizes[3*i  ] = all_sizes[3*i+1] = all_sizes[3*i+2] = -1;
    }
  }  
    
  
    // Store the total size of each set (rather than the
    // number of entities local to this processor) in the
    // local_sizes array for each meshset.  Only need this
    // for the sets this processor is writing the description
    // for, but it's easier to get it for all of them.
  sizes_iter = local_sizes.begin();
  viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter, ++viter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter ); 
    assert( p != tagsort.end() );
    int j = 3 * p->second;
    *sizes_iter = all_sizes[j  ]; ++sizes_iter;  // contents
    *sizes_iter = all_sizes[j+1]; ++sizes_iter;  // children
    *sizes_iter = all_sizes[j+2]; ++sizes_iter;  // parents
  }
  
    // Now calculate the offset of the data for each (entire, parallel) set in
    // the set contents, children and parents tables.  offsets is long[3], and
    // is both input and output of this function.  We increment offsets by the
    // total count (over all processors) for each set such that it contains
    // the next open row in the table.  This will be passed back into this
    // function for the next tag (or tag pair) such that ultimately it will
    // contain the beginning of the non-shared set data in each of the three tables.
    // all_sizes is re-used to store the global offset in each table for each 
    // set with the tag.
  for (i = 0; i < all_sizes.size(); ++i)
  {
    if (all_sizes[i] >= 0) // value is -1 (from above) if not this processor
    {
      int j = i % 3;              // contents, children or parents list ?
      long tmp = offsets[j];      // save current, running offset
      offsets[j] += all_sizes[i]; // next set's offset is current plus the size of this set
      all_sizes[i] = tmp;         // size of this set is running offset.
    }
  }
  
    // Local offsets for this processor are stored as values relative to the
    // start of each set's data.  Convert them to offsets relative to the
    // start of all the set data.  Add the offset *from* the start of the set
    // data (local_offsets) to the offset *of* the start of the set data 
    // (stored in all_sizes in the previous loop) 
  std::vector<long>::iterator offset_iter = local_offsets.begin();
  viter = data.local_values.begin();
  for (riter = data.range.begin(); riter != data.range.end(); ++riter, ++viter)
  {
    const std::map<int,int>::iterator p = tagsort.find( *viter );
    assert( p != tagsort.end() );
    int j = 3 * p->second;
    *offset_iter += all_sizes[j  ]; ++offset_iter; // contents
    *offset_iter += all_sizes[j+1]; ++offset_iter; // children
    *offset_iter += all_sizes[j+2]; ++offset_iter; // parents
  }

#ifdef DEBUG  
START_SERIAL; if (counts[myPcomm->proc_config().proc_rank()]) {
std::string name1, name2;
iFace->tag_get_name( data.data_tag, name1 );
iFace->tag_get_name( data.filter_tag, name2 );
printdebug("Remote set data\n" );
printdebug("    %13s %13s owner local_offsets total_counts\n", name1.c_str(), name2.c_str());
for (unsigned d = 0; d < (unsigned)counts[myPcomm->proc_config().proc_rank()]; ++d) {
switch(d%3) {
  case 0: // data/contents
printdebug("   %13d %13d %5s %13d %12d\n", data.all_values[(d+displs[myPcomm->proc_config().proc_rank()])/3], 
 data.filter_value, 
 all_sizes[d+displs[myPcomm->proc_config().proc_rank()]] < 0 ? "no" : "yes", 
 local_offsets[d], local_sizes[d] );
  break;
  case 1: // children
printdebug("                          (children) %13d %12d\n", local_offsets[d], local_sizes[d] );
  break;
  case 2: // parents
printdebug("                           (parents) %13d %12d\n", local_offsets[d], local_sizes[d] );
  break;
} 
}
} 
END_SERIAL;
#endif
  
    // Store each parallel meshset in the list
  sizes_iter = local_sizes.begin();
  offset_iter = local_offsets.begin();
  std::vector<long>::iterator all_iter = all_sizes.begin() + displs[myPcomm->proc_config().proc_rank()];
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
    ParallelSet info;
    info.handle = *riter;
    info.contentsOffset = *offset_iter; ++offset_iter;
    info.childrenOffset = *offset_iter; ++offset_iter;
    info.parentsOffset = *offset_iter; ++offset_iter;
    info.contentsCount = *sizes_iter; ++sizes_iter;
    info.childrenCount = *sizes_iter; ++sizes_iter;
    info.parentsCount = *sizes_iter; ++sizes_iter;
    info.description = *all_iter >= 0; all_iter += 3;
    parallelSets.push_back( info );
  }
  
  return MB_SUCCESS;
}


MBErrorCode WriteHDF5Parallel::write_shared_set_descriptions( hid_t table )
{
//char buffer[256];
//sprintf(buffer, "write_shared_set_descriptions( %u )", (unsigned)parallelSets.size() );
//TPRINT( buffer );
TPRINT( "write_shared_set_descriptions" );

  const id_t start_id = setSet.first_id;
  MBErrorCode rval;
  mhdf_Status status;
  
  for( std::list<ParallelSet>::const_iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    if (!iter->description)
      continue;  // handled by a different processor
    
      // Get offset in table at which to write data
    id_t file_id = idMap.find( iter->handle );
    assert( file_id && file_id >= start_id );
    file_id -= start_id;
    
      // Get flag data
    unsigned int flags;
    rval = iFace->get_meshset_options( iter->handle, flags );
    if (MB_SUCCESS != rval)
      return rval;
      
      // Write the data
    long data[4] = { iter->contentsOffset + iter->contentsCount - 1, 
                     iter->childrenOffset + iter->childrenCount - 1, 
                     iter->parentsOffset  + iter->parentsCount  - 1,
                     flags };
    mhdf_writeSetMeta( table, file_id, 1, H5T_NATIVE_LONG, data, &status );
    CHECK_HDF(status);
  }

TPRINT( "finished write_shared_set_descriptions" );

  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_contents( hid_t table )
{
TPRINT( "write_shared_set_contents" );

  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  for( std::list<ParallelSet>::const_iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_entities_by_handle( iter->handle, handle_list );
    if (MB_SUCCESS != rval)
      return rval;
    remove_remote_entities( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (id_list.empty())
      continue;
    
    mhdf_writeSetData( table, 
                       iter->contentsOffset, 
                       id_list.size(),
                       id_type,
                       &id_list[0],
                       &status );
    CHECK_HDF(status);
  }
  

TPRINT( "finished write_shared_set_contents" );
  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_children( hid_t table )
{
TPRINT( "write_shared_set_children" );

  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  printdebug("Writing %d parallel sets.\n", parallelSets.size());
  for( std::list<ParallelSet>::const_iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_child_meshsets( iter->handle, handle_list );
    if (MB_SUCCESS != rval)
      return rval;
    remove_remote_sets( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (!id_list.empty())
    {
      mhdf_writeSetParentsChildren( table, 
                                    iter->childrenOffset, 
                                    id_list.size(),
                                    id_type,
                                    &id_list[0],
                                    &status );
      CHECK_HDF(status);
    }
  }

TPRINT( "finished write_shared_set_children" );
  return MB_SUCCESS;
}
    

MBErrorCode WriteHDF5Parallel::write_shared_set_parents( hid_t table )
{
TPRINT( "write_shared_set_parents" );

  MBErrorCode rval;
  mhdf_Status status;
  std::vector<MBEntityHandle> handle_list;
  std::vector<id_t> id_list;
  
  printdebug("Writing %d parallel sets.\n", parallelSets.size());
  for( std::list<ParallelSet>::const_iterator iter = parallelSets.begin();
        iter != parallelSets.end(); ++iter)
  {
    handle_list.clear();
    rval = iFace->get_parent_meshsets( iter->handle, handle_list );
    if (MB_SUCCESS != rval)
      return rval;
    remove_remote_sets( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (!id_list.empty())
    {
      mhdf_writeSetParentsChildren( table, 
                                    iter->parentsOffset, 
                                    id_list.size(),
                                    id_type,
                                    &id_list[0],
                                    &status );
      CHECK_HDF(status);
    }
  }

TPRINT( "finished write_shared_set_parents" );
  return MB_SUCCESS;
}

MBErrorCode WriteHDF5Parallel::write_finished()
{
  parallelSets.clear();
  cpuParallelSets.clear();
  //myParallelSets.clear();
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


MBErrorCode WriteHDF5Parallel::exchange_file_ids()
{
    // create tag to store file IDs
  MBEntityHandle default_val = 0;
  MBTag file_id_tag = 0;
  MBErrorCode rval = iFace->tag_create( "__hdf5_ll_fileid", 
                                        sizeof(MBEntityHandle),
                                        MB_TAG_DENSE,
                                        MB_TYPE_HANDLE,
                                        file_id_tag,
                                        &default_val );
  if (MB_SUCCESS != rval)
    return rval;

    // get file ids for my interface entities
  MBRange::const_iterator i;
  const MBRange& imesh = interfaceMesh[myPcomm->proc_config().proc_rank()];
  std::vector<MBEntityHandle> file_id_vect( imesh.size() );
  VALGRIND_MAKE_VEC_UNDEFINED( file_id_vect );
  std::vector<MBEntityHandle>::iterator j = file_id_vect.begin();
  for (i = imesh.begin(); i != imesh.end(); ++i, ++j) {
    *j = idMap.find( *i );
    if (!*j) {
      iFace->tag_delete( file_id_tag );
      return MB_FAILURE;
    }
  }

printdebug( "Interface entities:\n" );
printrange( interfaceMesh[myPcomm->proc_config().proc_rank()] );
  
    // store file IDs in tag
  rval = iFace->tag_set_data( file_id_tag, imesh, &file_id_vect[0] );
  if (MB_SUCCESS != rval) {
    iFace->tag_delete( file_id_tag );
    return rval;
  }
  
    // do communication
  MBRange dum_range;
  rval = myPcomm->exchange_tags( file_id_tag, dum_range );
  if (MB_SUCCESS != rval) {
    iFace->tag_delete( file_id_tag );
    return rval;
  }
  
    // store file IDs for remote entities
  for (proc_iter p = interfaceMesh.begin(); p != interfaceMesh.end(); ++p) {
    if (p->first == myPcomm->proc_config().proc_rank())
      continue;

    file_id_vect.resize( p->second.size() );
    VALGRIND_MAKE_VEC_UNDEFINED( file_id_vect );
    rval = iFace->tag_get_data( file_id_tag, p->second, &file_id_vect[0] );
    if (MB_SUCCESS != rval) {
      iFace->tag_delete( file_id_tag );
      return rval;
    }
    
    j = file_id_vect.begin();
    for (i = p->second.begin(); i != p->second.end(); ++i, ++j) {
      if (*j == 0 || idMap.insert( *i, *j, 1 ) == idMap.end()) {
         iFace->tag_delete( file_id_tag );
         return MB_FAILURE;
      }
    }
  }
  
  return iFace->tag_delete( file_id_tag );
}


MBErrorCode WriteHDF5Parallel::get_sharedset_tags() 
{
    // get all the sets
  MBRange all_sets;
  MBErrorCode result = iFace->get_entities_by_type_and_tag(0, MBENTITYSET, NULL, NULL, 0,
                                                           all_sets);
  if (MB_SUCCESS != result || all_sets.empty()) return result;
  
    // get all the tags on those sets & test against known exceptions
  std::set<MBTag> all_tags;
  std::vector<MBTag> all_tagsv;
  std::string tag_name;
  
  for (MBRange::iterator rit = all_sets.begin(); rit != all_sets.end(); rit++) {
    all_tagsv.clear();
    result = iFace->tag_get_tags_on_entity(*rit, all_tagsv);
    if (MB_SUCCESS != result) return result;

    for (std::vector<MBTag>::iterator vit = all_tagsv.begin(); vit != all_tagsv.end(); vit++) {
        // don't look at tags already selected
      if (std::find(all_tags.begin(), all_tags.end(), *vit) != all_tags.end()) continue;

        // get name
      result = iFace->tag_get_name(*vit, tag_name);
      if (MB_SUCCESS != result) return result;

        // look for known exclusions
      const char *tag_cstr = tag_name.c_str();
      if (
        !((tag_cstr[0] == '_' && tag_cstr[1] == '_') ||
          tag_name == PARALLEL_SHARED_PROC_TAG_NAME ||
          tag_name == PARALLEL_SHARED_PROCS_TAG_NAME ||
          tag_name == PARALLEL_SHARED_HANDLE_TAG_NAME ||
          tag_name == PARALLEL_SHARED_HANDLES_TAG_NAME ||
          tag_name == PARALLEL_STATUS_TAG_NAME 
          ))
        all_tags.insert(*vit);
    }
  }
  
    // now add the tag names to the list
  for (std::set<MBTag>::iterator sit = all_tags.begin(); sit != all_tags.end(); sit++) {
    result = iFace->tag_get_name(*sit, tag_name);
    if (MB_SUCCESS != result) return result;
    multiProcSetTags.add( tag_name);
  }
    
  return MB_SUCCESS;
}
