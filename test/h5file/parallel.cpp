#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <unistd.h>
#include <signal.h>
/*#include <assert.h>*/
#include <mpi.h>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "WriteHDF5Parallel.hpp"
#include "FileOptions.hpp"

#include "testdir.h"

int numproc, rank;

#ifndef NDEBUG
#  define assert(A) do_assert(__FILE__, __LINE__, (A), #A)
#else
#  define assert(A)
#endif
void do_assert( const char* file, int line, bool condition, const char* condstr )
{
  if (condition)
    return;
  
  fprintf(stderr, "[%d] Assert(%s) failed at %s:%d\n", rank, condstr, file, line );
  abort();
}

#  define START_SERIAL                     \
     for (int _x = 0; _x < numproc; ++_x) {\
       MPI_Barrier( MPI_COMM_WORLD );      \
       if (_x != rank) continue     
#  define END_SERIAL                       \
     }                                     \
     MPI_Barrier( MPI_COMM_WORLD )

MBCore *iFace = NULL;
MBTag blockTag, geomTag, ifaceTag, idTag, gidTag;

static void printerror( const char* format, ... )
{
  FILE* const ostream = stdout;
  fprintf(ostream, "[%d] ", rank );
  va_list args;
  va_start(args, format);
  vfprintf(ostream, format, args );
  va_end(args);
  fprintf(ostream, "\n" );
  fflush(ostream);
}

      
void create_interface( MBEntityHandle geom, int proc );

volatile int go = 0;

extern "C" void handle_signal( int sig )
{
  go = 1;
}

int main( int argc, char* argv[] )
{
  MBErrorCode rval;
  int i;
  MBRange::iterator riter;
  std::vector<MBEntityHandle>::iterator viter;

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );

  iFace = new MBCore();
  
  char* defaults[] = { "out.blocks.h5m", "../../" TEST_DIR "/h5file/blocks.h5m" };
  bool wait_for_sig = false;
  char** fnames = argv + 1;
  int fname_count = argc - 1;
  if (argc > 1 && !strcmp( argv[i], "-w" )) {
    ++fnames;
    --fname_count;
  }
  if (0 == fname_count) {
    fnames = defaults;
    fname_count = 2;
  }
  else if (1 == fname_count) {
    fprintf(stderr, "Usage: %s [-w] output_file input_file [input_file2 ...]\n", argv[0] );
    return 1;
  }

  
  if (wait_for_sig)
  {
    signal( SIGUSR1, &handle_signal );
    while (!go) sleep(1);
  }
    
  
  printerror ( "MPI Initialized." );

  char wd[PATH_MAX];
  getcwd( wd, sizeof(wd) );
  printerror ("WorkingDir: %s\n", wd);
  
  for (i = 1; i < fname_count; ++i)
  {
    printerror ( "Reading \"%s\"...", fnames[i] );
    rval = iFace->load_mesh( fnames[i], 0, 0 );
    if (MB_SUCCESS != rval)
    {
      printerror( "Failed to read mesh file: \"%s\"", fnames[i] );
      exit( i+1 );
    }
  }
  
  printerror( "Read %d files.", fname_count - 1 );
  MPI_Barrier( MPI_COMM_WORLD );  
  
  printerror( "Getting/creating tags...");
  rval = iFace->tag_get_handle( MATERIAL_SET_TAG_NAME, blockTag ); assert(!rval);
  rval = iFace->tag_get_handle( GEOM_DIMENSION_TAG_NAME, geomTag ); assert(!rval);
  rval = iFace->tag_get_handle( GLOBAL_ID_TAG_NAME, idTag ); assert(!rval);
  rval = iFace->tag_create( PARALLEL_SHARED_PROC_TAG_NAME, 2*sizeof(int),
                           MB_TAG_SPARSE, ifaceTag, 0 ); assert(!rval);
  rval = iFace->tag_create( PARALLEL_GID_TAG_NAME, sizeof(MBEntityHandle),
                           MB_TAG_SPARSE, gidTag, 0 ); assert(!rval);
  
    // Get the list of geometry volumes this processor is to export
    // (id % numproc == rank).
  printerror( "Getting volumes to export...");
  MBRange myvolumes;
  int dim = 3;
  const void* dimarray[] = {&dim};
  rval = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, dimarray, 1, myvolumes ); assert(!rval);
  if (myvolumes.empty())
  {
    printerror ("Mesh contains no geometry sets.");
    exit(1);
  }
  riter = myvolumes.begin();
  while (riter != myvolumes.end())
  {
    int id;
    rval = iFace->tag_get_data( idTag, &*riter, 1, &id ); assert(!rval);
    if (id % numproc == rank)
      ++riter;
    else
      riter = myvolumes.erase( riter );
  }
  if (myvolumes.empty())
  {
    printerror( "No mesh to export.");
    exit (1);
  }

START_SERIAL;
  
    // Create interface mesh sets.
    // Iterate in descending geometric dimension.
  for (int dimension = 2; dimension >= 0; --dimension)
  {
    MBRange geometry;
    std::vector<MBEntityHandle> children;
    for (riter = myvolumes.begin(); riter != myvolumes.end(); ++riter)
    {
      children.clear();
      rval = iFace->get_child_meshsets( *riter, children, 3 - dimension ); assert(!rval);
      for (viter = children.begin(); viter != children.end(); ++viter)
      {
        int dim;
        rval = iFace->tag_get_data( geomTag, &*viter, 1, &dim );assert(!rval);
        if (dim == dimension)
          geometry.insert(*viter);
      }
    }
    assert(geometry.size());  // assuming no topological spheres...
      
    std::vector<MBEntityHandle> volumes, parents, parent_vols;
    for (riter = geometry.begin(); riter != geometry.end(); ++riter)
    {
      volumes.clear();
      rval = iFace->get_parent_meshsets( *riter, volumes, 3 - dimension );
      assert(MB_SUCCESS == rval);

int id2;
char tmpcstr[32];

rval = iFace->tag_get_data( idTag, &*riter, 1, &id2 ); assert(!rval);
std::string s = dimension == 2 ? "surface" : dimension == 1 ? "curve" : dimension == 0 ? "vertex" : "UNKNOWN";
sprintf(tmpcstr," %d", id2);
s += tmpcstr;
s += " : volumes:";
sprintf(tmpcstr," %d :", volumes.size());
s += tmpcstr;
for (viter = volumes.begin(); viter != volumes.end(); ++viter)
{
  int dim;
  rval = iFace->tag_get_data( geomTag, &*viter, 1, &dim ); assert(!rval);
  if (dim != 3) // not a volume 
    {continue;}
  rval = iFace->tag_get_data( idTag, &*viter, 1, &id2 ); assert(!rval);
  sprintf(tmpcstr," %d", id2);
  s += tmpcstr;
}
printerror("%s", s.c_str());

      for (viter = volumes.begin(); viter != volumes.end(); ++viter)
      {
          // Is the adjacent volume local or remote?
        int id, proc, dim;
        rval = iFace->tag_get_data( geomTag, &*viter, 1, &dim ); assert(!rval);
        if (dim != 3) // not a volume 
          continue;
        rval = iFace->tag_get_data( idTag, &*viter, 1, &id ); assert(!rval);
        proc = id % numproc;
        if (proc == rank) 
          continue;
        
          // Is this entity already part of some parent interface set?
          
          // Get list of immediate parents.
        parents.clear();
        if (2 - dimension > 0)  // empty list for surfaces.
          iFace->get_parent_meshsets( *riter, parents, 1 );
        
          // Check if any parent entity is the remote volume
        bool skip = false;
        for (std::vector<MBEntityHandle>::iterator piter = parents.begin();
             piter != parents.end(); ++piter)
        {
          parent_vols.clear();
          iFace->get_parent_meshsets( *piter, parent_vols, 2 - dimension );
          if (std::find( parent_vols.begin(), parent_vols.end(), *viter ) != parent_vols.end() )
          {
            skip = true;
            break;
          }
        }
          

        if (!skip)
        {
rval = iFace->tag_get_data( idTag, &*riter, 1, &id2 ); assert(!rval);
printerror("Creating interface for %s %d in remote volume %d",
dimension == 2 ? "surface" : dimension == 1 ? "curve" : dimension == 0 ? "vertex" : "UNKNOWN",
id2, id );
          create_interface( *riter, proc );
        }
      }
    }
  }  

END_SERIAL;

    // Find any blocks containing my volumes
  printerror( "Getting element blocks...");
  MBRange blocks, myblocks;
  rval = iFace->get_entities_by_type_and_tag( 0, MBENTITYSET, &blockTag, 0, 1, blocks ); assert(!rval);
  MBRange export_list, block_contents;
  for (riter = myvolumes.begin(); riter != myvolumes.end(); ++riter)
  {
    bool dovol = true;
    for (MBRange::iterator biter = blocks.begin(); biter != blocks.end(); ++biter)
    {
      block_contents.clear();
      rval = iFace->get_entities_by_type( *biter, MBENTITYSET, block_contents ); assert(!rval);
      if (block_contents.find(*riter) != block_contents.end())
      {
        export_list.insert( *biter );
        myblocks.insert(*biter);
        dovol = false;
      }
    }
    if (dovol)
    {
      export_list.insert( *riter );
    }
  }


    // Remove from blocks any volumes not to be exported by this processor
  for (riter = myblocks.begin(); riter != myblocks.end(); ++riter )
  {
    rval = iFace->get_entities_by_type( *riter, MBENTITYSET, block_contents ); assert(!rval);
    for (MBRange::iterator biter = block_contents.begin(); biter != block_contents.end(); ++biter)
    {
      if (myvolumes.find(*biter) == myvolumes.end())
      {
        rval = iFace->remove_entities( *riter, &*biter, 1 );
        assert(MB_SUCCESS == rval);
      }
    }
  }
  
  
    // Print the list of sets to export.
  START_SERIAL;
  printerror( "exporting sets: " );
  for (riter = export_list.begin(); riter != export_list.end(); ++riter)
  {
    int id, dimension;
    if (MB_SUCCESS == iFace->tag_get_data( blockTag, &*riter, 1, &id ))
      printerror( "\tblock %d", id );
    else if(MB_SUCCESS == iFace->tag_get_data( geomTag, &*riter, 1, &dimension ))
    {
      id = -1;
      iFace->tag_get_data( idTag, &*riter, 1, &id );
      printerror( "\t%s %d", dimension == 0 ? "vertex" :
                             dimension == 1 ? "curve" :
                             dimension == 2 ? "surface" :
                             dimension == 3 ? "volume" :
                             "invalid geom dimension", 
                             id );
    }
    else
    {
      printerror( "\tunknown meshset %ul", (unsigned long)*riter );
    }
  }
  END_SERIAL;
        
    // Convert range to vector
  std::vector<MBEntityHandle> list(export_list.size());
  std::copy( export_list.begin(), export_list.end(), list.begin() );
  
    // Assign global ID tag to all entities to make sure
    // all procs are starting out with the same thing.
  MBRange everything, extra;
  iFace->get_entities_by_handle( 0, everything ); // doesn't include meshsets
  iFace->get_entities_by_type( 0, MBENTITYSET, extra );
  everything.merge( extra );
  for (MBRange::iterator all_iter = everything.begin();
       all_iter != everything.end(); ++all_iter) {
    rval = iFace->tag_set_data( gidTag, &*all_iter, 1,  &*all_iter );  assert(!rval); 
  }
     
    // Write all the mesh in a single, serial file to compare with
    // the parallel output.
  if (0 == rank) {
    std::string sname = fnames[0];
    sname += ".serial.h5m";
    rval = iFace->write_mesh( sname.c_str() );
    if (MB_SUCCESS != rval)
    {
      printerror("Failed to write serial file: \"%s\"", sname.c_str());
      exit( 127 );
    }
    else
    {
      printerror("Wrote combined serial file: \"%s\"", sname.c_str());
    }
  }
  
    // Write individual file from each processor
  char str_rank[6];
  sprintf(str_rank, "%02d", rank );
  std::string name( str_rank );
  name += ".";
  name += fnames[0];
  char* ptr = strrchr( fnames[0], '.' );
  if (ptr && strcmp(ptr, ".h5m"))
    name += ".h5m";

  printerror ("Writing individual file: \"%s\"", name.c_str() );
  rval = iFace->write_mesh( name.c_str(), &list[0], list.size() ); 
  //rval = iFace->write_mesh( name.c_str() );
  if (MB_SUCCESS != rval)
  {
    printerror( "Failed to write per-processor file: \"%s\"", name.c_str() );
    exit( 127 );
  }
  printerror( "Wrote per-processor file: \"%s\"", name.c_str() );

  MPI_Barrier( MPI_COMM_WORLD );  

    // Write combined file
  std::vector<std::string> qa;
  qa.push_back( "MOAB Parallel HDF5 Writer." );
  time_t t = time( 0 );
  qa.push_back( ctime( &t ) );
  std::string qa3( "Processor " );
  qa3 += rank;
  qa.push_back( qa3 );
  WriteHDF5Parallel *writer = new WriteHDF5Parallel( iFace );
  
  printerror ("Writing parallel file: \"%s\"", fnames[0] );
  rval = writer->write_file( fnames[0], true, 
                             FileOptions("PARALLEL"),
                             &list[0], list.size(), qa );
  if (MB_SUCCESS != rval)
  {
    printerror( "Failed to write parallel file: \"%s\"", fnames[0] );
    exit( 127 );
  }
  printerror( "Wrote parallel file: \"%s\"", fnames[0] );
  
  H5close();
  MPI_Finalize();

  delete writer;
  delete iFace;
  
  return 0;
}

  
      
void create_interface( MBEntityHandle geom, int proc )
{
   MBEntityHandle iface_set;
   MBErrorCode rval;
   
   int procs[2];
   if (((proc % 2) == (rank % 2)) == (proc < rank))
   {
     procs[0] = proc;
     procs[1] = rank;
   }
   else
   {
     procs[0] = rank;
     procs[1] = proc;
   }

   rval = iFace->create_meshset( MESHSET_SET, iface_set ); assert(!rval);
   rval = iFace->tag_set_data( ifaceTag, &iface_set, 1, procs ); assert(!rval);
   std::vector<MBEntityHandle> children;
   rval = iFace->get_child_meshsets( geom, children, 0 ); assert(!rval);
   children.push_back(geom);
   rval = iFace->add_entities( iface_set, &children[0], children.size() ); assert(!rval);
   
   for (std::vector<MBEntityHandle>::iterator iter = children.begin();
        iter != children.end(); ++iter)
   {
      // Set global ID tag
     std::vector<MBEntityHandle> contents;
     rval = iFace->get_entities_by_handle( *iter, contents ); assert(!rval);
     rval = iFace->tag_set_data( gidTag, &contents[0], contents.size(), &contents[0] );  assert(!rval); 
   }
}
