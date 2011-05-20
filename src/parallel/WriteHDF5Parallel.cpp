
#undef DEBUG
#undef TIME_DEBUG

#include <stdio.h>
#include <stdarg.h>

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

#include "WriteHDF5Parallel.hpp"

#include <H5Tpublic.h>
#include <H5Ppublic.h>
#include <H5FDmpi.h>
#include <H5FDmpio.h>
#include <H5Spublic.h>

#include "mhdf.h"

#include "moab/Interface.hpp"
#include "Internals.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "moab/ParallelComm.hpp"
#include "moab/CN.hpp"
#include "moab/Range.hpp"

#include "IODebugTrack.hpp"

namespace moab {

#define CHECK_MPI(A) do { if (MPI_SUCCESS != (A)) { \
  writeUtil->report_error( "MPI Failure at " __FILE__ ":%d\n", __LINE__ ); \
  return error(MB_FAILURE); } } while(false)

#define CHECK_HDF(A) do { if (mhdf_isError(&(A))) { \
  writeUtil->report_error( "MHDF Failure at " __FILE__ ":%d : %s\n", \
    __LINE__, mhdf_message(&(A)) ); \
  return error(MB_FAILURE); } } while(false)


#ifdef VALGRIND
#  include <valgrind/memcheck.h>
#else
#  ifndef VALGRIND_CHECK_MEM_IS_DEFINED
#    define VALGRIND_CHECK_MEM_IS_DEFINED(a,b)
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
     fprintf( stderr, "[%d] Assert(%s) failed at %s:%d\n", rank, condstr, file, line );
     fflush( stderr );
     abort();
   }
#endif

class CpuTimer {
private:
  double atBirth, atLast;
public:
  CpuTimer() : atBirth(MPI_Wtime()), atLast(MPI_Wtime()) {}
  double since_birth() { return (atLast = MPI_Wtime()) - atBirth; };
  double elapsed() { double tmp = atLast; return (atLast = MPI_Wtime()) - tmp; }
};

static int my_Gatherv( void* sendbuf, 
                       int sendcount, 
                       MPI_Datatype sendtype,
                       std::vector<unsigned char>& recvbuf,
                       std::vector<int>& recvcounts,
                       int root, 
                       MPI_Comm comm )
{
  int nproc, rank, bytes, err;
  MPI_Comm_size( comm, &nproc );
  MPI_Comm_rank( comm, &rank );
  MPI_Type_size( sendtype, &bytes );
  
  recvcounts.resize( rank == root ? nproc : 0 );
  err = MPI_Gather( &sendcount, 1, MPI_INT, &recvcounts[0], 1, MPI_INT, root, comm );
  if (MPI_SUCCESS != err) return err;
  
  
  std::vector<int> disp(recvcounts.size());
  if (root == rank) {
    disp[0] = 0;
    for (int i = 1; i < nproc; ++i)
      disp[i] = disp[i-1] + recvcounts[i-1];
    recvbuf.resize( bytes * (disp.back() + recvcounts.back()) );
  }
  
  return MPI_Gatherv( sendbuf, sendcount, sendtype,
                      &recvbuf[0], &recvcounts[0], &disp[0],
                      sendtype, root, comm );
}


// This function doesn't do anything useful.  It's just a nice
// place to set a break point to determine why the reader fails.
static inline ErrorCode error( ErrorCode rval )
  { return rval; }

static void print_type_sets( Interface* iFace, DebugOutput* str, Range& sets )
{
  const unsigned VB = 2;
  if (str->get_verbosity() < VB)
    return;

  Tag gid, did, bid, sid, nid;
  iFace->tag_get_handle( GLOBAL_ID_TAG_NAME, gid ); 
  iFace->tag_get_handle( GEOM_DIMENSION_TAG_NAME, did );
  iFace->tag_get_handle( MATERIAL_SET_TAG_NAME, bid );
  iFace->tag_get_handle( DIRICHLET_SET_TAG_NAME, nid );
  iFace->tag_get_handle( NEUMANN_SET_TAG_NAME, sid );
  Range typesets[10];
  const char* typenames[] = {"Block ", "Sideset ", "NodeSet", "Vertex", "Curve", "Surface", "Volume", "Body", "Other"};
  for (Range::iterator riter = sets.begin(); riter != sets.end(); ++riter)
  {
    unsigned dim, id ; //, oldsize;
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

    //oldsize = typesets[dim].size();
    typesets[dim].insert( id );
//    assert( typesets[dim].size() - oldsize == 1 );  
  }
  for (int ii = 0; ii < 9; ++ii)
  {
    char tmp[64];
    sprintf(tmp,"%s (%lu) ",typenames[ii],(unsigned long)typesets[ii].size());
    str->print( VB, tmp, typesets[ii] );
  }
  str->printf(VB,"Total: %lu\n", (unsigned long)sets.size());
}


static void range_remove( Range& from, const Range& removed )
{
  
/* The following should be more efficient, but isn't due
   to the inefficient implementation of Range::erase(iter,iter)
  Range::const_iterator s, e, n = from.begin();
  for (Range::const_pair_iterator p = removed.pair_begin();
       p != removed.pair_end(); ++p)
  {
    e = s = Range::lower_bound(n, from.end(), (*p).first);
    e = Range::lower_bound(s, from.end(), (*p).second);
    if (e != from.end() && *e == (*p).second)
      ++e;
    n = from.erase( s, e );
  }
*/

  if (removed.size())
  {
    Range tmp = subtract( from, removed);
    from.swap( tmp );
  }
}

#define debug_barrier() debug_barrier_line(__LINE__)

void WriteHDF5Parallel::debug_barrier_line(int lineno)
{
  const unsigned threshold = 2;
  static unsigned long count = 0;
  if (dbgOut.get_verbosity() >= threshold) {
    dbgOut.printf( threshold, "*********** Debug Barrier %lu (@%d)***********\n", ++count, lineno);
    MPI_Barrier( myPcomm->proc_config().proc_comm() );
  }
}

WriterIface* WriteHDF5Parallel::factory( Interface* iface )
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

ErrorCode 
WriteHDF5Parallel::MultiProcSetTags::Data::get_sets( Interface* moab,
                                                     Range& sets,
                                                     Tag& id_tag ) const
{
  ErrorCode rval;
  int size;
  const void* values[] = { 0, 0 };
  Tag handles[2];

    // Get first tag handle: this tag is used to filter sets
    // (we only return sets for which this tag is set)

  rval = moab->tag_get_handle( filterTag.c_str(), handles[0] );
  if (MB_TAG_NOT_FOUND == rval)
    return MB_SUCCESS; // return nothing if tag isn't defined on this proc
  else if (MB_SUCCESS != rval)
    return error(rval);

    // Get data tag.  This tag contains IDs used to globally identify
    // sets in the group with the filter_tag.  This might be the same
    // tag as the filter tag.
  
  rval = moab->tag_get_handle( dataTag.c_str(), handles[1] );
  if (MB_TAG_NOT_FOUND == rval)
    return MB_SUCCESS; // return nothing if tag isn't defined on this proc
  else if (MB_SUCCESS != rval)
    return error(rval);
  
  moab->tag_get_size( handles[1], size );
  if (size != (int)sizeof(int)) {
    fprintf(stderr, "Cannot use non-int tag data for matching remote sets.\n" );
    assert(0);
    return error(MB_FAILURE);
  }  
  
    // We can optionally also filter on the value of the filter tag.
    // (e.g. GEOM_DIMENSION)
  
  if (useFilterValue) {
    moab->tag_get_size( handles[0], size );
    if (size != (int)sizeof(int)) {
      fprintf(stderr, "Cannot use non-int tag data for filtering remote sets.\n" );
      assert(0);
      return error(MB_FAILURE);
    }
    values[0] = &filterValue;
  }
  
    // If data tag and filter tag are the same, don't pass the
    // same tag handle in twice (it would work, but would be inefficient)
  int num_tags = (handles[0] == handles[1]) ? 1 : 2;

    // get the tagged entities
  id_tag = handles[1];
  return moab->get_entities_by_type_and_tag( 0,
                                             MBENTITYSET,
                                             handles,
                                             values,
                                             num_tags, 
                                             sets );
}



WriteHDF5Parallel::WriteHDF5Parallel( Interface* iface)
    : WriteHDF5(iface), myPcomm(NULL), pcommAllocated(false), hslabOp(H5S_SELECT_OR)
{
  multiProcSetTags.add(  MATERIAL_SET_TAG_NAME );
  multiProcSetTags.add( DIRICHLET_SET_TAG_NAME );
  multiProcSetTags.add(   NEUMANN_SET_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 0, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 1, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 2, GLOBAL_ID_TAG_NAME );
  multiProcSetTags.add( GEOM_DIMENSION_TAG_NAME, 3, GLOBAL_ID_TAG_NAME );
}

WriteHDF5Parallel::WriteHDF5Parallel( Interface* iface,
                                      const std::vector<std::string>& tag_names )
  : WriteHDF5(iface), myPcomm(NULL), pcommAllocated(false)
{
  for(std::vector<std::string>::const_iterator i = tag_names.begin();
      i != tag_names.end(); ++i)
    multiProcSetTags.add( *i );
}

WriteHDF5Parallel::WriteHDF5Parallel( Interface* iface,
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
//  o passes them back in a Range
ErrorCode WriteHDF5Parallel::gather_interface_meshes(Range& nonowned)
{
  ErrorCode result;
  
  //START_SERIAL;
  dbgOut.print( 3, "Pre-interface mesh:\n");
  dbgOut.print( 3, nodeSet.range );
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    dbgOut.print( 3, eiter->range );
  dbgOut.print( 3, setSet.range );
  
    // Move handles of non-owned entities from lists of entities
    // that this processor will write to the 'nonowned' list.
    
  nonowned.clear();
  result = myPcomm->filter_pstatus( nodeSet.range, PSTATUS_NOT_OWNED, PSTATUS_AND, -1, &nonowned);
  if (MB_SUCCESS != result)
    return error(result);
  nodeSet.range = subtract( nodeSet.range, nonowned );
  
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
       eiter != exportList.end(); ++eiter ) {
    Range tmpset;
    result = myPcomm->filter_pstatus( eiter->range, PSTATUS_NOT_OWNED, PSTATUS_AND, -1, &tmpset);
    if (MB_SUCCESS != result)
      return error(result);
    eiter->range = subtract( eiter->range,  tmpset );
    nonowned.merge(tmpset);
  }

  dbgOut.print( 3, "Post-interface mesh:\n");
  dbgOut.print( 3, nodeSet.range );
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
    dbgOut.print( 3, eiter->range );
  dbgOut.print( 3, setSet.range );

  //END_SERIAL;
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5Parallel::parallel_create_file( const char* filename,
                                            bool overwrite,
                                            const std::vector<std::string>& qa_records,
                                            const Tag* user_tag_list,
                                            int user_tag_count,
                                            int dimension,
                                            int pcomm_no,
                                            double* times)
{
  ErrorCode rval;
  mhdf_Status status;

  myPcomm = ParallelComm::get_pcomm(iFace, pcomm_no);
  if (0 == myPcomm) {
    myPcomm = new ParallelComm(iFace);
    pcommAllocated = true;
  }
  
  dbgOut.set_rank( myPcomm->proc_config().proc_rank() );

  Range nonlocal;
  debug_barrier();
  dbgOut.tprint(1,"Gathering interface meshes\n");
  rval = gather_interface_meshes( nonlocal );
  if (MB_SUCCESS != rval) return error(rval);

    /**************** get tag names for sets likely to be shared ***********/
  //debug_barrier();
  //dbgOut.tprint(1,"Getting shared entity sets\n");
  //rval = get_sharedset_tags();
  //if (MB_SUCCESS != rval) return error(rval);
  

    /**************** Create actual file and write meta info ***************/

  debug_barrier();
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    dbgOut.tprintf(1,"Creating file: %s\n",filename);
    
      // create the file
    const char* type_names[MBMAXTYPE];
    memset( type_names, 0, MBMAXTYPE * sizeof(char*) );
    for (EntityType i = MBEDGE; i < MBENTITYSET; ++i)
      type_names[i] = CN::EntityTypeName( i );
   
    dbgOut.tprint(1,"call mhdf_createFile\n");
    filePtr = mhdf_createFile( filename, overwrite, type_names, MBMAXTYPE, id_type, &status );
    if (!filePtr)
    {
      writeUtil->report_error( "%s\n", mhdf_message( &status ) );
      return error(MB_FAILURE);
    }
    
    dbgOut.tprint(1,"call write_qa\n");
    rval = write_qa( qa_records );
    if (MB_SUCCESS != rval) return error(rval);
  }
  
  
     /**************** Create node coordinate table ***************/
  CpuTimer timer;
  debug_barrier();
  dbgOut.tprint(1,"creating node table\n");
  topState.start("creating node table");
  rval = create_node_table( dimension );
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[CREATE_NODE_TIME] = timer.elapsed();
  
    /**************** Create element tables ***************/

  debug_barrier();
  dbgOut.tprint(1,"negotiating element types\n");
  topState.start("negotiating element types");
  rval = negotiate_type_list();
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[NEGOTIATE_TYPES_TIME] = timer.elapsed();
  dbgOut.tprint(1,"creating element table\n");
  topState.start("creating element tables");
  rval = create_element_tables();
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[CREATE_ELEM_TIME] = timer.elapsed();
  
  
    /*************** Exchange file IDs *****************/

  debug_barrier();
  dbgOut.tprint(1,"communicating file ids\n");
  topState.start("communicating file ids");
  rval = exchange_file_ids( nonlocal );
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[FILEID_EXCHANGE_TIME] = timer.elapsed();
 

    /**************** Create adjacency tables *********************/
  
  debug_barrier();
  dbgOut.tprint(1,"creating adjacency table\n");
  topState.start("creating adjacency tables");
  rval = create_adjacency_tables();
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[CREATE_ADJ_TIME] = timer.elapsed();
  
    /**************** Create meshset tables *********************/
  
  debug_barrier();
  dbgOut.tprint(1,"creating meshset table\n");
  topState.start("creating meshset tables");
  rval = create_meshset_tables(times);
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[CREATE_SET_TIME] = timer.elapsed();
  
  
    /* Need to write tags for shared sets this proc is responsible for */
  
  Range parallel_sets;
  for (std::list<ParallelSet>::const_iterator psiter = parallelSets.begin();
       psiter != parallelSets.end(); ++psiter)
    if (psiter->description)
      parallel_sets.insert( psiter->handle );
  
  setSet.range.merge( parallel_sets );
  rval = gather_tags( user_tag_list, user_tag_count );
  if (MB_SUCCESS != rval)
    return error(rval);
  range_remove( setSet.range, parallel_sets );   
  
    /**************** Create tag data *********************/

  debug_barrier();
  dbgOut.tprint(1,"creating tag tables\n");
  topState.start("creating tag tables");
  rval = create_tag_tables();
  topState.end(rval);
  if (MB_SUCCESS != rval) return error(rval);
  if (times) times[CREATE_TAG_TIME] = timer.elapsed();
    
  /************** Close serial file and reopen parallel *****************/
  
  if (0 == myPcomm->proc_config().proc_rank())
  {
    mhdf_closeFile( filePtr, &status );
  }
  
  MPI_Barrier( myPcomm->proc_config().proc_comm() );
  
  dbgOut.tprint(1,"(re)opening file in parallel mode\n");
  unsigned long junk;
  hid_t hdf_opt = H5Pcreate( H5P_FILE_ACCESS );
  H5Pset_fapl_mpio( hdf_opt, myPcomm->proc_config().proc_comm(), MPI_INFO_NULL );
  filePtr = mhdf_openFileWithOpt( filename, 1, &junk, id_type, hdf_opt, &status );
  H5Pclose( hdf_opt );
  if (!filePtr)
  {
    writeUtil->report_error( "%s\n", mhdf_message( &status ) );
    return error(MB_FAILURE);
  }
  
  if (collectiveIO) {
    dbgOut.print(1,"USING COLLECTIVE IO\n");
    writeProp = H5Pcreate( H5P_DATASET_XFER );
    H5Pset_dxpl_mpio( writeProp, H5FD_MPIO_COLLECTIVE );
  }
  
    /* Test if we can use H5S_APPEND when selecting hyperslabs */
  if (HDF5_can_append_hyperslabs()) {
    dbgOut.print(1,"HDF5 library supports H5Sselect_hyperlsab with H5S_SELECT_APPEND\n");
    hslabOp = H5S_SELECT_APPEND;
  }
  
  dbgOut.tprint(1,"Exiting parallel_create_file\n");
  return MB_SUCCESS;
}


class TagNameCompare {
  Interface* iFace;
  std::string name1, name2;
public:
  TagNameCompare( Interface* iface ) : iFace(iface) {}
  bool operator() (const WriteHDF5::TagDesc& t1, 
                   const WriteHDF5::TagDesc& t2);
};
bool TagNameCompare::operator() (const WriteHDF5::TagDesc& t1, 
                                 const WriteHDF5::TagDesc& t2)
{
  iFace->tag_get_name( t1.tag_id, name1 );
  iFace->tag_get_name( t2.tag_id, name2 );
  return name1 < name2;
}  

struct serial_tag_data {
  TagType storage;
  DataType type;
  int size;
  int name_len;
  char name[sizeof(unsigned long)];
  
  static size_t len( int name_len ) {
    return sizeof(serial_tag_data) + name_len - sizeof(unsigned long);
  }
  size_t len() const { return len( name_len ); }
};

ErrorCode WriteHDF5Parallel::append_serial_tag_data( 
                                         std::vector<unsigned char>& buffer,
                                         const WriteHDF5::TagDesc& tag )
{
  ErrorCode rval;
  
  std::string name;
  rval = iFace->tag_get_name( tag.tag_id, name );
  if (MB_SUCCESS != rval)
    return error(rval);
  
    // get name length padded to preserve alignement
  size_t name_len = name.size();
  if (!name_len) return MB_SUCCESS; // skip tags with no name
    // we pad such that we always terminate with a null character
  name_len += sizeof(unsigned long) - (name_len % sizeof(unsigned long));

    // allocate struct within buffer
  size_t init_size = buffer.size();
  buffer.resize( init_size + serial_tag_data::len(name_len) );
  serial_tag_data* ptr = reinterpret_cast<serial_tag_data*>(&buffer[init_size]);
  
    // populate struct
  rval = iFace->tag_get_type( tag.tag_id, ptr->storage );
  if (MB_SUCCESS != rval) return error(rval);
  rval = iFace->tag_get_data_type( tag.tag_id, ptr->type );
  if (MB_SUCCESS != rval) return error(rval);
  rval = iFace->tag_get_size( tag.tag_id, ptr->size );
  if (MB_VARIABLE_DATA_LENGTH == rval)
    ptr->size = MB_VARIABLE_LENGTH;
  else if (MB_SUCCESS != rval)
    return error(rval);
  ptr->name_len = name_len;
  Range range;
  memset( ptr->name, 0, ptr->name_len );
  memcpy( ptr->name, name.data(), name.size() );
  
  return MB_SUCCESS;
}
   

ErrorCode WriteHDF5Parallel::check_serial_tag_data( 
                               const std::vector<unsigned char>& buffer,
                               std::vector<TagDesc*>* missing,
                               std::vector<TagDesc*>* newlist )
{
  ErrorCode rval;

    // Use 'write_sparse' field as a 'visited' mark
  std::list<TagDesc>::iterator tag_iter;
  if (missing)
    for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter)
      tag_iter->write_sparse = true;

    // Use a set as a temporary for what will ultimately go in
    // newlist because we need to pass back newlist in the order
    // of the tagList member.
  std::set<TagDesc*> newset;

    // Iterate over data from, updating the local list of tags.
    // Be carefull to keep tagList sorted such that in the end all 
    // procs have the same list in the same order.
  std::vector<unsigned char>::const_iterator diter = buffer.begin();
  tag_iter = tagList.begin();
  while (diter < buffer.end()) {
      // Get struct from buffer
    const serial_tag_data* ptr = reinterpret_cast<const serial_tag_data*>(&*diter);
    
      // Find local struct for tag
    std::string name(ptr->name);
    std::string n;
    iFace->tag_get_name( tag_iter->tag_id, n );  // second time we've called, so shouldnt fail
    if (n > name) {
      tag_iter = tagList.begin(); // new proc, start search from beginning
    }
    iFace->tag_get_name( tag_iter->tag_id, n );
    while (n < name) {
      ++tag_iter;
      if (tag_iter == tagList.end())
        break;
      iFace->tag_get_name( tag_iter->tag_id, n );
    }
    if (tag_iter == tagList.end() || n != name) { // new tag
      TagDesc newtag;
      
      if (ptr->size == MB_VARIABLE_DATA_LENGTH) 
        rval = iFace->tag_create_variable_length( name.c_str(), ptr->storage, ptr->type, newtag.tag_id, 0 );
      else
        rval = iFace->tag_create( name.c_str(), ptr->size, ptr->storage, ptr->type, newtag.tag_id, 0 );
      if (MB_SUCCESS != rval)
        return error(rval);
      
      newtag.sparse_offset = 0;
      newtag.var_data_offset = 0;
      newtag.write_sparse = false;
      newtag.max_num_ents = 0;
      newtag.max_num_vals = 0;

      tag_iter = tagList.insert( tag_iter, newtag );
      if (newlist)
        newset.insert( &*tag_iter );
    }
    else { // check that tag is as expected
      DataType type;
      iFace->tag_get_data_type( tag_iter->tag_id, type );
      if (type != ptr->type) {
        writeUtil->report_error("Processes have inconsistent data type for tag \"%s\"", name.c_str() );
        return error(MB_FAILURE);
      }
      int size;
      iFace->tag_get_size( tag_iter->tag_id, size );
      if (type != ptr->type) {
        writeUtil->report_error("Processes have inconsistent size for tag \"%s\"", name.c_str() );
        return error(MB_FAILURE);
      }
      tag_iter->write_sparse = false;
    }
  
      // Step to next variable-length struct.
    diter += ptr->len();
  }

    // now pass back any local tags that weren't in the buffer
  if (missing) {
    for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
      if (tag_iter->write_sparse) {
        tag_iter->write_sparse = false;
        missing->push_back( &*tag_iter );
      }
    }
  }
  
    // be careful to populate newlist in the same, sorted, order as tagList
  if (newlist) {
    for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) 
      if (newset.find(&*tag_iter) != newset.end())
        newlist->push_back(&*tag_iter);
  }
  
  return MB_SUCCESS;
}

static void set_bit( int position, unsigned char* bytes )
{
  int byte = position/8;
  int bit = position%8;
  bytes[byte] |= (((unsigned char)1)<<bit);
}

static bool get_bit( int position, const unsigned char* bytes )
{
  int byte = position/8;
  int bit = position%8;
  return 0 != (bytes[byte] & (((unsigned char)1)<<bit));
}


ErrorCode WriteHDF5Parallel::create_tag_tables()
{  
  std::list<TagDesc>::iterator tag_iter;
  ErrorCode rval;
  int err;
  const int num_proc = myPcomm->proc_config().proc_size();
  const int rank = myPcomm->proc_config().proc_rank();
  const MPI_Comm comm = myPcomm->proc_config().proc_comm();

  subState.start( "negotiating tag list" );

  dbgOut.tprint(1,"communicating tag metadata\n");

  dbgOut.printf(2,"Exchanging tag data for %d tags.\n", (int)tagList.size() ); 
  
    // Sort tagList contents in alphabetical order by tag name
  tagList.sort( TagNameCompare( iFace ) );
  
    // Negotiate total list of tags to write
  
    // Build concatenated list of all tag data
  std::vector<unsigned char> tag_buffer;
  for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
    rval = append_serial_tag_data( tag_buffer, *tag_iter );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // broadcast list from root to all other procs
  unsigned long size = tag_buffer.size();
  err = MPI_Bcast( &size, 1, MPI_UNSIGNED_LONG, 0, comm );
  CHECK_MPI(err);
  tag_buffer.resize(size);
  err = MPI_Bcast( &tag_buffer[0], size, MPI_UNSIGNED_CHAR, 0, comm );
  CHECK_MPI(err);
  
    // update local tag list
  std::vector<TagDesc*> missing;
  rval = check_serial_tag_data( tag_buffer, &missing, 0 );
  if (MB_SUCCESS != rval)
    return error(rval);
  
    // check if we're done (0->done, 1->more, 2+->error)
  int code, lcode = (MB_SUCCESS != rval) ? rval + 2 : missing.empty() ? 0 : 1;
  err = MPI_Allreduce( &lcode, &code, 1, MPI_INT, MPI_MAX, comm );
  CHECK_MPI(err);
  if (code > 1) {
    writeUtil->report_error("Inconsistent tag definitions between procs.");
    return error((ErrorCode)(code-2));
  }  
  
  dbgOut.print(1,"Not all procs had same tag definitions, negotiating...\n");
  
    // if not done...
  if (code) {
      // get tags defined on this proc but not on root proc
    tag_buffer.clear();
    for (size_t i = 0; i < missing.size(); ++i) {
      rval = append_serial_tag_data( tag_buffer, *missing[i] );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
    
      // gather extra tag definitions on root processor
    std::vector<int> junk; // don't care how many from each proc
    assert(rank || tag_buffer.empty()); // must be empty on root
    err = my_Gatherv( &tag_buffer[0], tag_buffer.size(),
                       MPI_UNSIGNED_CHAR, tag_buffer, junk, 0, comm );
    CHECK_MPI(err);
    
      // process serialized tag descriptions on root, and 
    rval = MB_SUCCESS;
    if (0 == rank) {
        // process serialized tag descriptions on root, and 
      std::vector<TagDesc*> newlist;
      rval = check_serial_tag_data( tag_buffer, 0, &newlist );
      tag_buffer.clear();
        // re-serialize a unique list of new tag defintitions
      for (size_t i = 0; MB_SUCCESS == rval && i != newlist.size(); ++i) {
        rval = append_serial_tag_data( tag_buffer, *newlist[i] );
        if (MB_SUCCESS != rval)
          return error(rval);
      }
    }
    
      // broadcast any new tag definitions from root to other procs
    long size = tag_buffer.size();
    if (MB_SUCCESS != rval)
      size = -rval;
    err = MPI_Bcast( &size, 1, MPI_LONG, 0, comm );
    CHECK_MPI(err);
    if (size < 0) {
      writeUtil->report_error("Inconsistent tag definitions between procs.");
      return error((ErrorCode)-size);
    }
    tag_buffer.resize(size);
    err = MPI_Bcast( &tag_buffer[0], size, MPI_UNSIGNED_CHAR, 0, comm );
    CHECK_MPI(err);
    
      // process new tag definitions
    rval = check_serial_tag_data( tag_buffer, 0, 0 );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
  subState.end();
  subState.start("negotiate which element/tag combinations are dense");

    // Figure out for which tag/element combinations we can
    // write dense tag data.  

    // Construct a table of bits,
    // where each row of the table corresponds to a tag
    // and each column to an element group.

    // Two extra, because first is nodes and last is sets.
    // (n+7)/8 is ceil(n/8)
  const int bytes_per_tag = (exportList.size() + 9)/8;
  std::vector<unsigned char> data(bytes_per_tag * tagList.size(), 0);
  std::vector<unsigned char> recv( data.size(), 0 );
  unsigned char* iter = &data[0];
  if (writeTagDense) {
    for (tag_iter = tagList.begin(); tag_iter != tagList.end();
         ++tag_iter, iter += bytes_per_tag) {

      Range tagged;
      rval = get_sparse_tagged_entities( *tag_iter, tagged );
      if (MB_SUCCESS != rval)
        return error(rval);

      int s;
      if (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s )) 
        continue;  

      std::string n;
      iFace->tag_get_name( tag_iter->tag_id, n );  // second time we've called, so shouldnt fail

      // Check if we want to write this tag in dense format even if not
      // all of the entities have a tag value.  The criterion of this
      // is that the tag be dense, have a default value, and have at
      // least 2/3 of the entities tagged.
      bool prefer_dense = false;
      TagType type;
      rval = iFace->tag_get_type( tag_iter->tag_id, type );
      if (MB_SUCCESS != rval)
        return error(rval);
      if (MB_TAG_DENSE == type) {
        const void* defval = 0;
        rval = iFace->tag_get_default_value( tag_iter->tag_id, defval, s );
        if (MB_SUCCESS == rval)
          prefer_dense = true;
      }

      int i = 0;
      if (check_dense_format_tag( nodeSet, tagged, prefer_dense )) {
        set_bit( i, iter );
        dbgOut.printf( 2, "Can write dense data for \"%s\"/Nodes\n", n.c_str());
      }
      std::list<ExportSet>::const_iterator ex_iter = exportList.begin();
      for (++i; ex_iter != exportList.end(); ++i, ++ex_iter) {
        if (check_dense_format_tag( *ex_iter, tagged, prefer_dense )) {
          set_bit( i, iter );
          dbgOut.printf( 2, "Can write dense data for \"%s\"/%s\n", n.c_str(),
            ex_iter->name());
        }
      }
      if (check_dense_format_tag( setSet, tagged, prefer_dense )) {
        set_bit( i, iter );
        dbgOut.printf( 2, "Can write dense data for \"%s\"/Sets\n", n.c_str());
      }
    }

      // Do bit-wise AND of list over all processors (only write dense format
      // if all proccesses want dense format for this group of entities).
    err = MPI_Allreduce( &data[0], &recv[0], data.size(), MPI_UNSIGNED_CHAR,
                         MPI_BAND, myPcomm->proc_config().proc_comm() );
    CHECK_MPI(err);
  } // if (writeTagDense)
  
    // Store initial counts for sparse-formatted tag data.
    // The total number of values to send and receive will be the number of 
    // tags plus the number of var-len tags because we need to negitiate 
    // offsets into two different tables for the var-len tags.
  std::vector<unsigned long> counts;
  
    // Record dense tag/element combinations
  iter = &recv[0];
  const unsigned char* iter2 = &data[0];
  for (tag_iter = tagList.begin(); tag_iter != tagList.end();
       ++tag_iter, iter += bytes_per_tag, iter2 += bytes_per_tag) {

    Range tagged;
    rval = get_sparse_tagged_entities( *tag_iter, tagged );
    if (MB_SUCCESS != rval)
      return error(rval);

    std::string n;
    iFace->tag_get_name( tag_iter->tag_id, n );  // second time we've called, so shouldnt fail

    int i = 0;
    if (get_bit(i, iter)) {
      assert(get_bit(i, iter2));
      tag_iter->dense_list.push_back(nodeSet);
      tagged -= nodeSet.range;
      dbgOut.printf( 2, "Will write dense data for \"%s\"/Nodes\n", n.c_str());
    }
    std::list<ExportSet>::const_iterator ex_iter = exportList.begin();
    for (++i; ex_iter != exportList.end(); ++i, ++ex_iter) {
      if (get_bit(i, iter)) {
        assert(get_bit(i, iter2));
        tag_iter->dense_list.push_back(*ex_iter);
        dbgOut.printf( 2, "WIll write dense data for \"%s\"/%s\n", n.c_str(),
          ex_iter->name());
        tagged -= ex_iter->range;
      }
    }
    if (get_bit(i, iter)) {
      assert(get_bit(i, iter2));
      tag_iter->dense_list.push_back(setSet);
      dbgOut.printf( 2, "Will write dense data for \"%s\"/Sets\n", n.c_str());
      tagged -= setSet.range;
    }

    counts.push_back( tagged.size() );

    int s;
    if (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s )) {
      unsigned long data_len;
      rval = get_tag_data_length( *tag_iter, tagged, data_len );
      assert(MB_SUCCESS == rval);
      if (MB_SUCCESS != rval) return error(rval);
      counts.push_back( data_len );
    }
  }
  

  subState.end();
  subState.start("Negotiate offsets for sparse tag info");

  
    // For each tag, gather counts on root and offsets into tag data
    // tables back to each process.  The total number of values to send
    // and receive will be the number of tags plus the number of var-len
    // tags because we need to negitiate offsets into two different 
    // tables for the var-len tags.

  std::vector<unsigned long> allcnt(counts.size()*num_proc);
  err = MPI_Gather( &counts[0], counts.size(), MPI_UNSIGNED_LONG, 
                    &allcnt[0], counts.size(), MPI_UNSIGNED_LONG,
                    0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(err);
  
  std::vector<unsigned long> maxima(counts.size(),0), sums(counts.size(),0);
  if (0 == myPcomm->proc_config().proc_rank()) {
    for (unsigned i = 0; i < myPcomm->proc_config().proc_size(); ++i) {
      for (unsigned j = 0; j < counts.size(); ++j) {
         if (allcnt[i*counts.size()+j] > maxima[j])
          maxima[j] = allcnt[i*counts.size()+j];
        size_t offset = sums[j];
        sums[j] += allcnt[i*counts.size()+j];
        allcnt[i*counts.size()+j] = offset;
      }
    }
  }
  err = MPI_Scatter( &allcnt[0], counts.size(), MPI_UNSIGNED_LONG,
                     &counts[0], counts.size(), MPI_UNSIGNED_LONG,
                     0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(err);
  err = MPI_Bcast( &maxima[0], maxima.size(), MPI_UNSIGNED_LONG, 0, 
                   myPcomm->proc_config().proc_comm() );
  CHECK_MPI(err);
  
    // Copy values into local structs and if root then create tables
  std::vector<unsigned long>::iterator citer = counts.begin(), miter = maxima.begin();
  for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
    assert(citer != counts.end() && miter != maxima.end());
    tag_iter->sparse_offset = *citer; ++citer;
    tag_iter->max_num_ents = *miter; ++miter;
    tag_iter->write_sparse = (0 != tag_iter->max_num_ents);
    int s;
    if (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s )) {
      assert(citer != counts.end() && miter != maxima.end());
      tag_iter->var_data_offset = *citer; ++citer;
      tag_iter->max_num_vals = *miter; ++miter;
    }
    else {
      tag_iter->var_data_offset = 0;
      tag_iter->max_num_vals = 0;
    }
  }
  
  subState.end();

    // Create tag tables on root process
  if (0 == myPcomm->proc_config().proc_rank()) {
    citer = sums.begin();
    for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
      assert(citer != sums.end());
      unsigned long num_ents = *citer; ++citer;
      unsigned long num_val = 0;
      int s;
      if (MB_VARIABLE_DATA_LENGTH == iFace->tag_get_size( tag_iter->tag_id, s )) {
        assert(citer != sums.end());
        num_val = *citer; ++citer;
      }
      dbgOut.printf( 2, "Writing tag description for tag 0x%lx with %lu values\n", 
                     (unsigned long)tag_iter->tag_id, num_val ? num_val : num_ents );

      rval = create_tag( *tag_iter, num_ents, num_val );
      if (MB_SUCCESS != rval)
        return error(rval);
    }
  }
  
  if (dbgOut.get_verbosity() > 1) {
    dbgOut.printf(2,"Tags: %12s %8s %8s %8s %8s %8s\n", "Name", "Count", "Offset", "Var Off", "Max Ent", "Handle");

    for (tag_iter = tagList.begin(); tag_iter != tagList.end(); ++tag_iter) {
      std::string name;
      iFace->tag_get_name( tag_iter->tag_id, name );
      size_t size;
      get_num_sparse_tagged_entities( *tag_iter, size );
      dbgOut.printf(2,"%18s %8lu %8lu %8lu %8lu 0x%7lx\n", name.c_str(), 
        (unsigned long)size, 
        (unsigned long)tag_iter->sparse_offset,
        (unsigned long)tag_iter->var_data_offset,
        (unsigned long)tag_iter->max_num_ents,
        (unsigned long)tag_iter->tag_id );
    }
  }

  return MB_SUCCESS;
}

ErrorCode WriteHDF5Parallel::create_node_table( int dimension )
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
  long first_id_and_max_count[2];
  if (myPcomm->proc_config().proc_rank() == 0)
  {
    int total = 0;
    for (unsigned int i = 0; i < myPcomm->proc_config().proc_size(); i++)
      total += node_counts[i];
    
    nodeSet.total_num_ents = total;
    hid_t handle = mhdf_createNodeCoords( filePtr, dimension, total, &first_id_and_max_count[0], &status );
    CHECK_HDF(status);
    mhdf_closeData( filePtr, handle, &status );
    
    first_id_and_max_count[1] = *std::max_element( node_counts.begin(), node_counts.end() );
  }
    
    // send id offset to every proc
  result = MPI_Bcast( first_id_and_max_count, 2, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  nodeSet.first_id = (id_t)first_id_and_max_count[0];
  nodeSet.max_num_ents = (id_t)first_id_and_max_count[1];
   
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


ErrorCode WriteHDF5Parallel::negotiate_type_list()
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

  dbgOut.print( 2,"Local Element Types:\n");
  viter = my_types.begin();
  while (viter != my_types.end())
  {
    int type = *viter; ++viter;
    int count = *viter; ++viter;
    dbgOut.printf(2, "  %s : %d\n", CN::EntityTypeName((EntityType)type), count);
  }

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

  dbgOut.print(2, "Global Element Types:\n");
  viter = intlist.begin();
  while (viter != intlist.end())
  {
    int type = *viter; ++viter;
    int count = *viter; ++viter;
    dbgOut.printf(2,"  %s : %d\n", CN::EntityTypeName((EntityType)type), count);
  }
  
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
      insert.type = (EntityType)mbtype;
      insert.num_nodes = numnode;
      insert.first_id = 0;
      insert.offset = 0;
      insert.adj_offset = 0;
      ex_iter = exportList.insert( ex_iter, insert );
    }
  }
  
  return MB_SUCCESS;
}

ErrorCode WriteHDF5Parallel::create_element_tables()
{
  int result;
  ErrorCode rval;
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
  if (myPcomm->proc_config().proc_rank() == 0) {
    for (int i = 0; i < numtypes; i++) 
    {
      long prev = 0;
      for (unsigned int j = 0; j < myPcomm->proc_config().proc_size(); j++)
      {
        long tmp = counts[j*numtypes + i];
        if (tmp > my_counts[i]) // put max count in my_counts
          my_counts[i] = tmp;
        counts[j*numtypes+i] = prev;
        prev += tmp;
      }
      counts[myPcomm->proc_config().proc_size()*numtypes + i] = prev;
    }
  }
  
    // Broadcast max count for each type
  result = MPI_Bcast( &my_counts[0], numtypes, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  viter = my_counts.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    ex_iter->max_num_ents = *(viter++);
  
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
      ex_iter->total_num_ents = *citer;
      rval = create_elem_tables( *ex_iter, *viter );
      if (MB_SUCCESS != rval)
        return error(rval);
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
      return error(rval);
  }
  
  return MB_SUCCESS;
}
  
ErrorCode WriteHDF5Parallel::create_adjacency_tables()
{
  ErrorCode rval;
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
      return error(rval);
    *viter = num_adj; ++viter;
  }
  
    // Send local adjacency counts to root processor
  VALGRIND_CHECK_MEM_IS_DEFINED( &local[0], numtypes*sizeof(long) );
  result = MPI_Gather( &local[0], numtypes, MPI_LONG,
                       &all[0],   numtypes, MPI_LONG, 
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
  
      // Record in 'local' the maximum adjacency count for each type
  if (myPcomm->rank() == 0) {
      // Convert counts to offsets
    for (i = 0; i < numtypes; i++) 
    {
      long prev = 0;
      for (j = 0; j < myPcomm->proc_config().proc_size(); j++)
      {
        long tmp = all[j*numtypes + i];
        if (tmp > local[i])
          local[i] = tmp;
        all[j*numtypes+i] = prev;
        prev += tmp;
      }
      all[myPcomm->proc_config().proc_size()*numtypes+i] = prev;
    }

      // For each element type for which there is no adjacency data,
      // send -1 to all processors as the offset.
    for (i = 0; i < numtypes; ++i)
      if (all[numtypes*myPcomm->proc_config().proc_size()+i] == 0)
        for (j = 0; j < myPcomm->proc_config().proc_size(); ++j)
          all[j*numtypes+i] = -1;
  }
  
    // Broadcast maximum count for each type
  result = MPI_Bcast( &local[0], numtypes, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  viter = local.begin();
  for (ex_iter = exportList.begin(); ex_iter != exportList.end(); ++ex_iter)
    ex_iter->max_num_adjs = *(viter++);
  
  
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
  
/** Working data for group of global sets identified by an ID tag */
struct RemoteSetData {
  Range range;     //!< Set handles with data_tag set (and optionally filter_tag == filter_value)
  std::vector<int> counts;       //!< Number of sets with tag on each proc, indexed by MPI rank
  std::vector<int> displs;       //!< Offset in all_values at which the data_tag values for
                                 //!< each processor begin. displs[n] = sum(i from 0 to n-1)(counts[i])
  std::vector<int> all_values;   //!< data_tag values for sets on all processors, 
                                 //!< counts[0] values for proc 0, then counts[1] values for proc 
                                 //!< 1, etc.
  std::vector<int> local_values; //!< data_tag values for sets that exist on this processor
};

/*
static void print_remote_set_data( DebugOutput* ds, Interface* mb,
                            const struct RemoteSetData& data,
                            const char* pfx = "" )
{
  std::string n1,n2;
  if (data.data_tag)
    mb->tag_get_name( data.data_tag, n1 );
  if (data.filter_tag)
    mb->tag_get_name( data.filter_tag, n2 );

  ds->printf( 2, "data_tag:     \"%s\"\n", n1.c_str() );
  ds->printf( 2, "filter_tag:   \"%s\"\n", n2.c_str() );
  ds->printf( 2, "filter_value:  %d\n", data.filter_value );
  ds->print( 2, "range:        ", data.range );
  ds->print( 2, "counts:" );
  for (size_t i = 0; i < data.counts.size(); ++i)
    ds->printf( 2, " %d,", data.counts[i] );
  ds->print( 2, "\ndispls:" );
  for (size_t i = 0; i < data.displs.size(); ++i)
    ds->printf( 2, " %d,", data.displs[i] );
  ds->print( 2, "\nall_values:" );
  for (size_t i = 0; i < data.all_values.size(); ++i)
    ds->printf( 2, " %d,", data.all_values[i] );
  ds->print( 2, "\nlocal_values:" );
  for (size_t i = 0; i < data.local_values.size(); ++i)
    ds->printf( 2, " %d,", data.local_values[i] );
  ds->print( 2, "\n\n" );
}
*/

ErrorCode WriteHDF5Parallel::get_remote_set_data( 
                        const WriteHDF5Parallel::MultiProcSetTags::Data& tags,
                        RemoteSetData& data, long& offset )
{
  ErrorCode rval;
  int i;
  int result;
  Range::iterator riter;
  Tag id_tag;

  rval = tags.get_sets( iFace, data.range, id_tag );
  if (MB_SUCCESS != rval)
    return error(rval);

  if (!data.range.empty()) {
    Range tmp = intersect( data.range,  setSet.range );
    data.range.swap( tmp );
    range_remove( setSet.range, data.range );
  }
  
  dbgOut.printf(1,"Found %lu meshsets with \"%s\" tag.\n", (unsigned long)data.range.size(), tags.filterTag.c_str() );

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
  rval = iFace->tag_get_data( id_tag, data.range, &data.local_values[0] );
  if (MB_SUCCESS != rval)
    return error(rval);
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
    for (i = 0; i < data.counts[j]; ++i)
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

ErrorCode WriteHDF5Parallel::set_shared_set_ids( RemoteSetData& data, long& offset )
{
    // Determine File ID for each shared set.
  std::map<int,long> val_id_map; // Map from tag value to file ID for set
  for (size_t i = 0; i < data.all_values.size(); ++i)
  {
    std::map<int,long>::iterator p = val_id_map.find( data.all_values[i] );
    if (p == val_id_map.end())
      val_id_map[data.all_values[i]] = offset++;
  }
  Range::const_iterator riter = data.range.begin();
  for (size_t i = 0; i < data.local_values.size(); ++i, ++riter)
  {
    std::map<int,long>::iterator p = val_id_map.find( data.local_values[i] );
    assert( p != val_id_map.end() );
    long id = p->second;
    
    if (!idMap.insert( *riter, id, 1 ).second) {
      for (unsigned x = 0; x < myPcomm->size(); ++x) {
        MPI_Barrier( myPcomm->proc_config().proc_comm() );      
        if (x != myPcomm->rank()) continue;   

        std::ostringstream s;
        s << "[" << myPcomm->rank() << "] ";
        std::string pfx1 = s.str();
        s << "  ";
        std::string pfx2 = s.str();

        std::cerr << pfx1 << "Duplicate shared set handle or internal accounting error" << std::endl;
        //std::cerr << pfx1 << "RemoteSetData:  " << std::endl;
        //print_remote_set_data( &dbgOut, iFace, data, pfx2.c_str() );

        std::cerr << pfx1 << "val_id_map: " << std::endl;
        for (p = val_id_map.begin(); p != val_id_map.end(); ++p)
          std::cerr << pfx2 << p->first << "->" << p->second << std::endl;

        std::cerr << pfx1 << "idMap: " << std::endl;
        print_id_map( std::cerr, pfx2.c_str() );

        std::cerr << pfx1 << "Failed at: (" << i << ") " << data.local_values[i] 
                  << "->" << id << " for " 
                  << CN::EntityTypeName(TYPE_FROM_HANDLE(*riter)) 
                  << " " << ID_FROM_HANDLE(*riter) << std::endl;
        std::cerr.flush();
      }
      
      assert(false);
      return error(MB_FAILURE);
    }
  }
  
  return MB_SUCCESS;
}
  

ErrorCode WriteHDF5Parallel::create_meshset_tables(double* times)
{
  ErrorCode rval = MB_SUCCESS;
  int result, i;
  long total_offset = 0;
  Range::const_iterator riter;

  START_SERIAL;
  print_type_sets( iFace, &dbgOut, setSet.range );
  END_SERIAL;
  CpuTimer timer;

    // Gather data about multi-processor meshsets - removes sets from setSet.range
  std::vector<RemoteSetData> remote_set_data( multiProcSetTags.list.size() );
  for (i = 0; i< (int)multiProcSetTags.list.size(); i++)
  {
    subState.start( "Resolving shared sets for tag: ", multiProcSetTags.list[i].filterTag.c_str() );
    rval = get_remote_set_data( multiProcSetTags.list[i],
                                remote_set_data[i],
                                total_offset ); 
    subState.end(rval);
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  if (times) times[RESOLVE_SHARED_SET_TIME] = timer.elapsed();

  subState.start( "Negotiating offsets for local sets" );
  dbgOut.print(2, "myLocalSets\n");
  print_type_sets( iFace, &dbgOut, setSet.range );

    // Gather counts of non-shared sets from each proc
    // to determine total table size.
  std::vector<long> set_offsets(myPcomm->proc_config().proc_size() + 1);
  //VALGRIND_MAKE_VEC_UNDEFINED( set_offsets );
  long local_count = setSet.range.size(), max_count = 0;
  result = MPI_Gather( &local_count,    1, MPI_LONG,
                       &set_offsets[0], 1, MPI_LONG,
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  if (myPcomm->proc_config().proc_rank() == 0) {
    for (unsigned int j = 0; j <= myPcomm->proc_config().proc_size(); j++)
    {
      long tmp = set_offsets[j];
      set_offsets[j] = total_offset;
      total_offset += tmp;
      if (max_count < tmp)
        max_count = tmp;
    }
  }
  
    // Send each proc its offsets in the set description table.
  long sets_offset;
  result = MPI_Scatter( &set_offsets[0], 1, MPI_LONG,
                        &sets_offset,    1, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setSet.offset = (id_t)(sets_offset);
    // Send each proc the max count so that we can do collective IO
  result = MPI_Bcast( &max_count, 1, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setSet.max_num_ents = max_count;

    // Create the set description table
  long total_count_and_start_id[2] = { set_offsets[myPcomm->proc_config().proc_size()], 0 };
  if (myPcomm->proc_config().proc_rank() == 0 && total_count_and_start_id[0] > 0)
  {
    setSet.total_num_ents = total_count_and_start_id[0];
    rval = create_set_meta( total_count_and_start_id[1] );
    if (MB_SUCCESS != rval)
      return error(rval);
  }
  
    // Send totals to all procs.
  result = MPI_Bcast( total_count_and_start_id, 2, MPI_LONG, 0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setSet.first_id = total_count_and_start_id[1];
  writeSets = total_count_and_start_id[0] > 0;

  dbgOut.printf(2,"Non-shared sets: %ld local, %ld global, offset = %ld, first_id = %ld\n",
    local_count, total_count_and_start_id[0], sets_offset, total_count_and_start_id[1] );
  dbgOut.printf(2,"my Parallel Sets:\n");
  print_type_sets(iFace, &dbgOut, cpuParallelSets[myPcomm->proc_config().proc_rank()] );
  
  subState.end();
  if (times) times[LOCAL_SET_OFFSET_TIME] = timer.elapsed();
  
    // Not writing any sets??
  if (!writeSets)
    return MB_SUCCESS;
  
  subState.start( "Negotiating offsets for shared set data" );
  
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
      return error(rval);
  }
  remote_set_data.clear();
  subState.end();
  if (times) times[SHARED_SET_OFFSET_TIME] = timer.elapsed();
  
    // Exchange IDs for remote/adjacent sets not shared between procs
  //rval = communicate_remote_ids( MBENTITYSET ); assert(MB_SUCCESS == rval);
  
    // Communicate counts for local sets
  long data_counts[3];
  rval = count_set_size( setSet.range, data_counts[0], data_counts[1], data_counts[2] );
  if (MB_SUCCESS != rval) 
    return error(rval);
  std::vector<long> set_counts(3*myPcomm->proc_config().proc_size());
  VALGRIND_MAKE_VEC_UNDEFINED( set_counts );
  result = MPI_Gather( data_counts,    3, MPI_LONG,
                       &set_counts[0], 3, MPI_LONG,
                       0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  
  
  long max_counts[3] = { 0, 0, 0 };
  if (myPcomm->proc_config().proc_rank() == 0)
  {
      // calculate offsets for each set, and maximum count for each proccessor
    for (unsigned int j = 0; j < 3*myPcomm->proc_config().proc_size(); ++j)
    {
      long tmp = set_counts[j];
      set_counts[j] = data_offsets[j%3];
      data_offsets[j%3] += tmp;
      if (tmp > max_counts[j%3])
        max_counts[j%3] = tmp;
    }
    
      // Create set contents and set children tables
    rval = create_set_tables( data_offsets[0], data_offsets[1], data_offsets[2] );
    if (MB_SUCCESS != rval) 
      return error(rval);
  }
  
    // Broadcast both max_counts and total table sizes.
    // max_counts is necessary for collective IO because each
    // proc needs to know the total number of writes that must
    // be done collectively while writing the contents of non-shared
    // sets. 
    // Also send the total table sizes, which also accounts for shared
    // sets, so that processors know whether or not any data is to be
    // written for the corresponding table.
  long summary_counts[6] = { max_counts[0], max_counts[1], max_counts[2],
                             data_offsets[0], data_offsets[1], data_offsets[2] };
  result = MPI_Bcast( summary_counts, 6, MPI_LONG, 0, 
                      myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);

  result = MPI_Scatter( &set_counts[0], 3, MPI_LONG,
                        data_offsets,   3, MPI_LONG,
                        0, myPcomm->proc_config().proc_comm() );
  CHECK_MPI(result);
  setContentsOffset = data_offsets[0];
  setChildrenOffset = data_offsets[1];
  setParentsOffset = data_offsets[2];
  maxNumSetContents = summary_counts[0];
  maxNumSetChildren = summary_counts[1];
  maxNumSetParents  = summary_counts[2];
  writeSetContents = summary_counts[3] > 0;
  writeSetChildren = summary_counts[4] > 0;
  writeSetParents  = summary_counts[5] > 0;

  dbgOut.printf(2,"Non-shared set contents: %ld local, %ld global, offset = %ld\n",
    data_counts[0], summary_counts[3], data_offsets[0] );
  dbgOut.printf(2,"Non-shared set children: %ld local, %ld global, offset = %ld\n",
    data_counts[1], summary_counts[4], data_offsets[1] );
  dbgOut.printf(2,"Non-shared set parents: %ld local, %ld global, offset = %ld\n",
    data_counts[2], summary_counts[5], data_offsets[2] );
  if (times) times[LOCAL_SET_OFFSET_TIME] += timer.elapsed();
  
  return MB_SUCCESS;
}

void WriteHDF5Parallel::remove_remote_entities( EntityHandle relative,
                                                Range& range )
{
  Range result;
  result.merge( intersect( range,  nodeSet.range ) );
  result.merge( intersect( range,  setSet.range ) );  
  for (std::list<ExportSet>::iterator eiter = exportList.begin();
           eiter != exportList.end(); ++eiter )
  {
    result.merge( intersect( range, eiter->range ) );
  }
  //result.merge( intersect( range, myParallelSets ) );
  Range sets;
  int junk;
  sets.merge( Range::lower_bound( range.begin(), range.end(), CREATE_HANDLE(MBENTITYSET, 0, junk )), range.end() );
  remove_remote_sets( relative, sets );
  result.merge( sets );
  range.swap(result);
}

void WriteHDF5Parallel::remove_remote_sets( EntityHandle relative, 
                                            Range& range )
{
  Range result( intersect( range,  setSet.range ) );
  //result.merge( intersect( range, yParallelSets ) );
  Range remaining( subtract( range, result ) );
  
  for(Range::iterator i = remaining.begin(); i != remaining.end(); ++i)
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
  
  

void WriteHDF5Parallel::remove_remote_entities( EntityHandle relative,
                                                std::vector<EntityHandle>& vect )
{
  Range intrsct;
  for (std::vector<EntityHandle>::const_iterator iter = vect.begin();
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

  

void WriteHDF5Parallel::remove_remote_sets( EntityHandle relative,
                                            std::vector<EntityHandle>& vect )
{
  Range intrsct;
  for (std::vector<EntityHandle>::const_iterator iter = vect.begin();
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
ErrorCode WriteHDF5Parallel::negotiate_remote_set_contents( RemoteSetData& data,
                                                              long* offsets /* long[3] */ )
{
  unsigned i;
  ErrorCode rval;
  Range::const_iterator riter;
  int result;
  const unsigned count = data.range.size();
  const unsigned total = data.all_values.size();
  std::vector<int>::iterator viter, viter2;

    // Calculate counts for each meshset
  std::vector<long> local_sizes(3*count);
  VALGRIND_MAKE_VEC_UNDEFINED( local_sizes );
  std::vector<long>::iterator sizes_iter = local_sizes.begin();
  Range tmp_range;
  std::vector<EntityHandle> child_list;
  for (riter = data.range.begin(); riter != data.range.end(); ++riter)
  {
      // Count contents
    *sizes_iter = 0;
    tmp_range.clear();
    rval = iFace->get_entities_by_handle( *riter, tmp_range );
    if (MB_SUCCESS != rval)
      return error(rval);
    remove_remote_entities( *riter, tmp_range );
    for (Range::iterator iter = tmp_range.begin(); iter != tmp_range.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    ++sizes_iter;
    
      // Count children
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_child_meshsets( *riter, child_list );
    if (MB_SUCCESS != rval)
      return error(rval);
    remove_remote_sets( *riter, child_list );
    for (std::vector<EntityHandle>::iterator iter = child_list.begin();
         iter != child_list.end(); ++iter)
      if (0 != idMap.find( *iter ))
        ++*sizes_iter;
    ++sizes_iter;
    
      // Count parents
    *sizes_iter = 0;
    child_list.clear();
    rval = iFace->get_parent_meshsets( *riter, child_list );
    if (MB_SUCCESS != rval)
      return error(rval);
    remove_remote_sets( *riter, child_list );
    for (std::vector<EntityHandle>::iterator iter = child_list.begin();
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


ErrorCode WriteHDF5Parallel::write_shared_set_descriptions( hid_t table,
                                              IODebugTrack* dbg_track )
{
  dbgOut.tprint( 1, "write_shared_set_descriptions\n" );

  const id_t start_id = setSet.first_id;
  ErrorCode rval;
  
    // Count number of parallel sets for which we need to write
    // the description
  std::list<ParallelSet>::const_iterator iter;
  size_t count = 0;
  for (iter = parallelSets.begin(); iter != parallelSets.end(); ++iter)
    if (iter->description)
      ++count;
  
    // get a large enough buffer for all of the shared set metadata
  std::vector<mhdf_index_t> tmp_buffer;
  size_t buffer_size = bufferSize / (4*sizeof(mhdf_index_t));
  mhdf_index_t* buffer;
  if (buffer_size >= count) {
    buffer = reinterpret_cast<mhdf_index_t*>(dataBuffer);
  }
  else {
    tmp_buffer.resize(4 * count);
    buffer = &tmp_buffer[0];
  }
  
    // get a buffer for offset locations
  Range indices;
  Range::iterator hint;
  
  mhdf_index_t* buff_iter = buffer;
  for( iter = parallelSets.begin();
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
      return error(rval);
      
      // Write the data
    *buff_iter = iter->contentsOffset + iter->contentsCount - 1; ++buff_iter;
    *buff_iter = iter->childrenOffset + iter->childrenCount - 1; ++buff_iter;
    *buff_iter = iter->parentsOffset  + iter->parentsCount  - 1, ++buff_iter;
    *buff_iter = flags;                                          ++buff_iter;
    hint = indices.insert( hint, file_id+1 );
    
    if (dbg_track) dbg_track->record_io( file_id, 1 );
  }

  herr_t herr = 0;
  hid_t space = H5Dget_space( table );
  hsize_t dims[2] = { count, 4 }, offsets[2] = { 0, 0 };
  hid_t mem;
  H5S_seloper_t op = H5S_SELECT_SET;
  if (count) {
    mem = H5Screate_simple( 2, dims, NULL );
    Range::const_pair_iterator pi;
    dims[1] = 4;
    for (pi = indices.const_pair_begin(); pi != indices.const_pair_end(); ++pi) {
      offsets[0] = pi->first-1;
      dims[0] = pi->second - pi->first + 1;
      herr = H5Sselect_hyperslab( space, op, offsets, 0, dims, 0 );
      if (herr < 0)
        break;
      op = hslabOp;
    }
  }
  else {
    hsize_t one = 1;
    mem = H5Screate_simple( 1, &one, NULL );
    herr = H5Sselect_none( mem );
    if (herr < 0) {
      H5Sclose( space );
      H5Sclose( mem );
      writeUtil->report_error( "H5Sselect_none failed for set contents" );
      return error(MB_FAILURE);
    }
    herr = H5Sselect_none( space );
  }

  if (herr < 0) {
    H5Sclose( mem );
    H5Sclose( space );
    dbgOut.tprint(1,"H5Sselect_elements failed\n");
    writeUtil->report_error("H5Sselect_elements failed");
    return error(MB_FAILURE);
  }
  herr = H5Dwrite( table, MHDF_INDEX_TYPE, mem, space, writeProp, buffer );
  H5Sclose( space );
  H5Sclose( mem );
  if (herr < 0) {
    dbgOut.tprint(1,"Failed to write shared set descriptions\n");
    writeUtil->report_error("H5Dwrite of shared set descriptions failed.");
    return error(MB_FAILURE);
  }

  return MB_SUCCESS;
}

ErrorCode WriteHDF5Parallel::write_shared_set_data( hid_t table,
                                   WriteUtilIface::EntityListType which_data,
                                   IODebugTrack* dbg_track )
{
  herr_t err;
  ErrorCode rval;
  std::vector<EntityHandle> handle_list;
  std::vector<id_t> id_list;
  std::vector<id_t> big_list;
  hid_t data_space = H5Dget_space( table );
  H5Sselect_none( data_space );
  
  std::list<ParallelSet>::const_iterator iter;
  for (iter = parallelSets.begin(); iter != parallelSets.end(); ++iter) {
    handle_list.clear();
    long offset;
    switch (which_data) {
      case WriteUtilIface::CONTENTS:
        rval = iFace->get_entities_by_handle( iter->handle, handle_list );
        offset = iter->contentsOffset;
        break;
      case WriteUtilIface::CHILDREN:
        rval = iFace->get_child_meshsets( iter->handle, handle_list );
        offset = iter->childrenOffset;
        break;
      case WriteUtilIface::PARENTS:
        rval = iFace->get_parent_meshsets( iter->handle, handle_list );
        offset = iter->parentsOffset;
        break;
    }
    if (MB_SUCCESS != rval) {
      H5Sclose( data_space );
      return error(rval);
    }
    remove_remote_entities( iter->handle, handle_list );
    
    id_list.clear();
    vector_to_id_list( handle_list, id_list, true );
    if (id_list.empty())
      continue;
    
    if (dbg_track) 
      dbg_track->record_io( offset, id_list.size() );
      
    hsize_t start = offset;
    hsize_t count = id_list.size();
    err = H5Sselect_hyperslab( data_space, hslabOp, &start, 0, &count, 0 );
    if (err < 0) {
      H5Sclose( data_space );
      writeUtil->report_error( "H5Sselect_hyperslab failed for set contents" );
      return error(MB_FAILURE);
    }
    std::copy( id_list.begin(), id_list.end(), std::back_inserter(big_list) );
  }
  
  hid_t mem;
  if (big_list.empty()) {
//#if H5_VERS_MAJOR > 1 || H5_VERS_MINOR >= 8
//    mem = H5Screate(H5S_NULL); 
//#else
    hsize_t one = 1;
    mem = H5Screate_simple( 1, &one, NULL );
    err = H5Sselect_none( mem );
    if (err < 0) {
      H5Sclose( data_space );
      H5Sclose( mem );
      writeUtil->report_error( "H5Sselect_none failed for set contents" );
      return error(MB_FAILURE);
    }
//#endif
    H5Sselect_none(data_space);
  }
  else {
    hsize_t dim = big_list.size();
    mem = H5Screate_simple( 1, &dim, NULL );
  }
  
  err = H5Dwrite( table, id_type, mem, data_space, writeProp, &big_list[0] );
  H5Sclose( data_space );
  H5Sclose( mem );
  if (err < 0) {
    dbgOut.tprint(1,"Failed to write shared set descriptions\n");
    writeUtil->report_error("H5Dwrite of shared set descriptions failed.");
    return error(MB_FAILURE);
  }
    
  return MB_SUCCESS;    
}

ErrorCode WriteHDF5Parallel::write_shared_set_contents( hid_t table,
                                            IODebugTrack* dbg_track )
{
  dbgOut.tprint(1, "write_shared_set_contents\n" );
  return write_shared_set_data( table, WriteUtilIface::CONTENTS, dbg_track );
}
    

ErrorCode WriteHDF5Parallel::write_shared_set_children( hid_t table,
                                            IODebugTrack* dbg_track )
{
  dbgOut.tprint(1, "write_shared_set_children\n" );
  return write_shared_set_data( table, WriteUtilIface::CHILDREN, dbg_track );
}
    

ErrorCode WriteHDF5Parallel::write_shared_set_parents( hid_t table,
                                            IODebugTrack* dbg_track )
{
  dbgOut.tprint(1, "write_shared_set_parents\n" );
  return write_shared_set_data( table, WriteUtilIface::PARENTS, dbg_track );
}

ErrorCode WriteHDF5Parallel::write_finished()
{
  parallelSets.clear();
  cpuParallelSets.clear();
  //myParallelSets.clear();
  return WriteHDF5::write_finished();
}


ErrorCode WriteHDF5Parallel::exchange_file_ids( const Range& nonlocal )
{
  ErrorCode rval;
  
    // For each entity owned on the interface, write its file id to
    // a tag.  The sets of entities to be written should already contain
    // only owned entities so by intersecting with them we not only
    // filter by entities to be written, but also restrict to entities
    // owned by the proc
    
    // Get list of interface entities
  Range imesh, tmp;
  for (std::list<ExportSet>::reverse_iterator i = exportList.rbegin();
       i != exportList.rend(); ++i) {
    tmp.clear();
    rval = myPcomm->filter_pstatus( i->range, PSTATUS_SHARED, PSTATUS_AND, -1, &tmp );
    if (MB_SUCCESS != rval) return error(rval);
    imesh.merge(tmp);
  }
  tmp.clear();
  rval = myPcomm->filter_pstatus( nodeSet.range, PSTATUS_SHARED, PSTATUS_AND, -1, &tmp );
  if (MB_SUCCESS != rval) return error(rval);
  imesh.merge(tmp);

  
    // create tag to store file IDs
  EntityHandle default_val = 0;
  Tag file_id_tag = 0;
  rval = iFace->tag_create( "__hdf5_ll_fileid", 
                              sizeof(EntityHandle),
                              MB_TAG_DENSE,
                              MB_TYPE_HANDLE,
                              file_id_tag,
                              &default_val );
  if (MB_SUCCESS != rval)
    return error(rval);

  
    // copy file IDs into tag
  std::vector<EntityHandle> file_id_vect( imesh.size() );
  Range::const_iterator i;
  std::vector<EntityHandle>::iterator j = file_id_vect.begin();
  for (i = imesh.begin(); i != imesh.end(); ++i, ++j) {
    *j = idMap.find( *i );
    if (!*j) {
      iFace->tag_delete( file_id_tag );
      return error(MB_FAILURE);
    }
  }
  rval = iFace->tag_set_data( file_id_tag, imesh, &file_id_vect[0] );
  if (MB_SUCCESS != rval) {
    iFace->tag_delete( file_id_tag );
    return error(rval);
  }

  
    // do communication
  rval = myPcomm->exchange_tags( file_id_tag, imesh );
  if (MB_SUCCESS != rval) {
    iFace->tag_delete( file_id_tag );
    return error(rval);
  }
  
    // copy file IDs from tag into idMap for remote entities
  file_id_vect.resize( nonlocal.size() );
  rval = iFace->tag_get_data( file_id_tag, nonlocal, &file_id_vect[0] );
  if (MB_SUCCESS != rval) {
    iFace->tag_delete( file_id_tag );
    return error(rval);
  }
    
  j = file_id_vect.begin();
  for (i = nonlocal.begin(); i != nonlocal.end(); ++i, ++j) {
    if (*j == 0) {
       int owner = -1;
       myPcomm->get_owner( *i, owner );
       const char* name = CN::EntityTypeName(TYPE_FROM_HANDLE(*i));
       int id = ID_FROM_HANDLE(*i);
       writeUtil->report_error( "Process %u did not receive valid id handle "
                                "for shared %s %d owned by process %d",
                                myPcomm->proc_config().proc_rank(),
                                name, id, owner );
       dbgOut.printf(1,"Did not receive valid remote id for "
                                "shared %s %d owned by process %d",
                                name, id, owner );
       iFace->tag_delete( file_id_tag );
       return error(MB_FAILURE);
    }
    else {
      if (!idMap.insert( *i, *j, 1 ).second) {
        iFace->tag_delete( file_id_tag );
        return error(MB_FAILURE);
      }
    }
  }
  
#ifndef NDEBUG
    // check that writer is correct with regards to which entities
    // that it owns by verifying that the file ids that we thought
    // we were sending where not received instead
  file_id_vect.resize( imesh.size() );
  rval = iFace->tag_get_data( file_id_tag, imesh, &file_id_vect[0] );
  if (MB_SUCCESS != rval) {
    iFace->tag_delete( file_id_tag );
    return error(rval);
  }
  int invalid_count = 0;
  j = file_id_vect.begin();
  for (i = imesh.begin(); i != imesh.end(); ++i, ++j) {
    EntityHandle h = idMap.find(*i);
    if (*j != h) {
      ++invalid_count;
      dbgOut.printf(1,"Conflicting owneship for %s %ld\n",
        CN::EntityTypeName(TYPE_FROM_HANDLE(*i)),
        (long)ID_FROM_HANDLE(*i));
    }
  }
  if (invalid_count) {
    writeUtil->report_error("%d entities with conflicting ownership found "
                            "by process %u.  This will result in duplicate "
                            "entities written to file.\n",
                            invalid_count, myPcomm->proc_config().proc_rank() );
    iFace->tag_delete( file_id_tag );
    return error(MB_FAILURE);
  }
#endif   
  
  return iFace->tag_delete( file_id_tag );
}


ErrorCode WriteHDF5Parallel::get_sharedset_tags() 
{
    // get all the sets
  Range all_sets;
  ErrorCode result = iFace->get_entities_by_type_and_tag(0, MBENTITYSET, NULL, NULL, 0,
                                                           all_sets);
  if (MB_SUCCESS != result) return error(result);
  if (all_sets.empty()) return MB_SUCCESS;
  
    // get all the tags on those sets & test against known exceptions
  std::set<Tag> all_tags;
  std::vector<Tag> all_tagsv;
  std::string tag_name;
  
  for (Range::iterator rit = all_sets.begin(); rit != all_sets.end(); rit++) {
    all_tagsv.clear();
    result = iFace->tag_get_tags_on_entity(*rit, all_tagsv);
    if (MB_SUCCESS != result) return error(result);

    for (std::vector<Tag>::iterator vit = all_tagsv.begin(); vit != all_tagsv.end(); vit++) {
        // don't look at tags already selected
      if (std::find(all_tags.begin(), all_tags.end(), *vit) != all_tags.end()) continue;

        // get name
      result = iFace->tag_get_name(*vit, tag_name);
      if (MB_SUCCESS != result) return error(result);

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
  for (std::set<Tag>::iterator sit = all_tags.begin(); sit != all_tags.end(); sit++) {
    result = iFace->tag_get_name(*sit, tag_name);
    if (MB_SUCCESS != result) return error(result);
    multiProcSetTags.add( tag_name);
  }
    
  return MB_SUCCESS;
}


void WriteHDF5Parallel::print_times( const double* times ) const
{
  if (!myPcomm) {
    WriteHDF5::print_times( times );
  }
  else {
    double recv[NUM_TIMES];
    MPI_Reduce( (void*)times, recv, NUM_TIMES, MPI_DOUBLE, MPI_MAX, 0, myPcomm->proc_config().proc_comm() );
    if (0 == myPcomm->proc_config().proc_rank())
      WriteHDF5::print_times( recv );
  }
}


} // namespace moab
