#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include "MBCore.hpp"
#include "MBTagConventions.hpp"
#include "MBParallelConventions.h"
#include "WriteHDF5Parallel.hpp"

MBCore iFace;
int numproc, rank;
MBTag blockTag, geomTag, pGeomTag, ifaceTag, idTag, gidTag;

static void printerror( const char* format, ... )
{
  fprintf(stderr, "[%d] ", rank );
  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args );
  va_end(args);
  fprintf(stderr, "\n" );
  fflush(stderr);
}

      
void create_interface( MBEntityHandle geom, int proc );


int main( int argc, char* argv[] )
{
  MBErrorCode rval;
  int i;
  MBTag blockTag, geomTag, pGeomTag, IfaceTag, idTag, gidTag;
  MBRange::iterator riter;
  std::vector<MBEntityHandle>::iterator viter;

for (i = 0; i < argc; ++i)
  printf("%s\n", argv[i]);

  MPI_Init( &argc, &argv );
  MPI_Comm_size( MPI_COMM_WORLD, &numproc );
  MPI_Comm_rank( MPI_COMM_WORLD, &rank    );
  
  if (argc < 3)
  {
    printerror( "Usage: %s <output file> <input file list>", argv[0] );
    exit( 1 );
  }
  
  printerror ( "MPI Initialized." );
  
  
  for (i = 2; i < argc; ++i)
  {
    printerror ( "Reading \"%s\"...", argv[i] );
    rval = iFace.load_mesh( argv[i], 0, 0 );
    if (MB_SUCCESS != rval)
    {
      printerror( "Failed to read mesh file: \"%s\"", argv[i] );
      exit( i+1 );
    }
  }
  
  printerror( "Read %d files.", argc-2 );
  MPI_Barrier( MPI_COMM_WORLD );  
  
  printerror( "Getting/creating tags...");
  rval = iFace.tag_get_handle( MATERIAL_SET_TAG_NAME, blockTag ); assert(!rval);
  rval = iFace.tag_get_handle( GEOM_DIMENSION_TAG_NAME, geomTag ); assert(!rval);
  rval = iFace.tag_get_handle( GLOBAL_ID_TAG_NAME, idTag ); assert(!rval);
  rval = iFace.tag_create( PARALLEL_GEOM_TOPO_TAG_NAME, sizeof(int), 
                           MB_TAG_SPARSE, pGeomTag, 0 ); assert(!rval);
  rval = iFace.tag_create( PARALLEL_INTERFACE_TAG_NAME, 2*sizeof(int),
                           MB_TAG_SPARSE, ifaceTag, 0 ); assert(!rval);
  rval = iFace.tag_create( PARALLEL_GLOBAL_ID_TAG_NAME, sizeof(MBEntityHandle),
                           MB_TAG_SPARSE, gidTag, 0 ); assert(!rval);
  
    // Get list meshsets for geometric volumes and surfaces
  printerror( "Getting volumes to export...");
  MBRange geomsets, myvolumes, mysurfaces;
  rval = iFace.get_entities_by_type_and_tag( 0, MBENTITYSET, &geomTag, 0, 1, geomsets ); assert(!rval);
  if (geomsets.empty())
  {
    printerror ("Mesh contains no geometry sets.");
    exit(1);
  }
  for (riter = geomsets.begin(); riter != geomsets.end(); ++riter)
  {
    int dimension;
    rval = iFace.tag_get_data( geomTag, &*riter, 1, &dimension ); assert(!rval);
    if (3 != dimension)
      continue;
    
    int id;
    rval = iFace.tag_get_data( idTag, &*riter, 1, &id ); assert(!rval);
    if (id % numproc != rank)
      continue;
    
    myvolumes.insert( *riter );
    std::vector<MBEntityHandle> surfaces;
    rval = iFace.get_child_meshsets( *riter, surfaces ); assert(!rval);
    for (viter = surfaces.begin(); viter != surfaces.end(); ++viter)
      mysurfaces.insert( *viter );
  }
  if (myvolumes.empty())
  {
    printerror( "No mesh to export.");
    exit (1);
  }
  
    // Create surface interface sets
  printerror( "Getting interface surfaces...");
  MBRange mycurves, finished_vols;
  for (riter = mysurfaces.begin(); riter != mysurfaces.end(); ++riter)
  {
    std::vector<MBEntityHandle> vols, curves;
    rval = iFace.get_parent_meshsets( *riter, vols ); assert(!rval);
    for (viter = vols.begin(); viter != vols.end(); ++viter)
    {
      int id, proc;
      rval = iFace.tag_get_data( idTag, &*riter, 1, &id ); assert(!rval);
      proc = id % numproc;
      if (proc == rank) 
        continue;
      
      create_interface( *riter, proc );
      finished_vols.insert( *viter );
    }
    
    curves.clear();
    rval = iFace.get_child_meshsets( *riter, curves ); assert(!rval);
    for (viter = curves.begin(); viter != curves.end(); ++viter)
      mycurves.insert( *viter );
  }
  
    // Create curve interface sets
  printerror( "Getting interface curves...");
  MBRange myvertices;
  for (riter = mycurves.begin(); riter != mycurves.end(); ++riter)
  {
    std::vector<MBEntityHandle> vols, vertices;
    rval = iFace.get_parent_meshsets( *riter, vols, 2 ); assert(!rval);
    for (viter = vols.begin(); viter != vols.end(); ++viter)
    {
      int id, proc;
      rval = iFace.tag_get_data( idTag, &*riter, 1, &id ); assert(!rval);
      proc = id % numproc;
      if (proc == rank) 
        continue;
      
      if (finished_vols.find( *viter ) != finished_vols.end())
        continue;
      
      create_interface( *riter, proc );
      finished_vols.insert( *viter );
    }
    
    vertices.clear();
    rval = iFace.get_child_meshsets( *riter, vertices ); assert(!rval);
    for (viter = vertices.begin(); viter != vertices.end(); ++viter)
      myvertices.insert( *viter );
  }
  
    // Create vertex interface sets
  printerror( "Getting interface vertices...");
  for (riter = myvertices.begin(); riter != myvertices.end(); ++riter)
  {
    std::vector<MBEntityHandle> vols;
    rval = iFace.get_parent_meshsets( *riter, vols, 2 ); assert(!rval);
    for (viter = vols.begin(); viter != vols.end(); ++viter)
    {
      int id, proc;
      rval = iFace.tag_get_data( idTag, &*riter, 1, &id ); assert(!rval);
      proc = id % numproc;
      if (proc == rank) 
        continue;
      
      if (finished_vols.find( *viter ) != finished_vols.end())
        continue;
      
      create_interface( *riter, proc );
      finished_vols.insert( *viter );
    }
  }

    // Find any blocks containing my volumes
  printerror( "Getting element blocks...");
  MBRange blocks;
  rval = iFace.get_entities_by_type_and_tag( 0, MBENTITYSET, &blockTag, 0, 1, blocks ); assert(!rval);
  MBRange export_list, block_contents;
  for (riter = myvolumes.begin(); riter != myvolumes.end(); ++riter)
  {
    bool dovol = true;
    for (MBRange::iterator biter = blocks.begin(); biter != blocks.end(); ++biter)
    {
      block_contents.clear();
      rval = iFace.get_entities_by_type( *biter, MBENTITYSET, block_contents ); assert(!rval);
      if (block_contents.find(*riter) != block_contents.end())
      {
        export_list.insert( *biter );
        dovol = false;
      }
    }
    if (dovol)
      export_list.insert( *riter );
  }
  
    // Print the list of sets to export.
    // Use loop and MPI_Barrier to serialize writes.
  for (i = 0; i < numproc; ++i)
  {
    MPI_Barrier( MPI_COMM_WORLD ); 
    if (i != rank) 
      continue;
      
    printerror( "exporting sets: " );
    for (riter = export_list.begin(); riter != export_list.end(); ++riter)
    {
      int id, dimension;
      if (MB_SUCCESS == iFace.tag_get_data( blockTag, &*riter, 1, &id ))
        printerror( "\tblock %d", id );
      else if(MB_SUCCESS == iFace.tag_get_data( geomTag, &*riter, 1, &dimension ))
      {
        id = -1;
        iFace.tag_get_data( idTag, &*riter, 1, &id );
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
  }
  MPI_Barrier( MPI_COMM_WORLD );  
        
    // Convert range to vector
  std::vector<MBEntityHandle> list(export_list.size());
  std::copy( export_list.begin(), export_list.end(), list.begin() );
  
    // Write individual file from each processor
  char str_rank[6];
  sprintf(str_rank, ".%02d", rank );
  char* ptr = strrchr( argv[1], '.' );
  if (ptr && strcmp(ptr, ".h5m")) ptr = 0;
  if (ptr)
    *ptr = '\0';
  std::string name = argv[1];
  if (ptr)
    *ptr = '.';
  name += str_rank;
  name += ".h5m";
/*
  printerror ("Writing individual file: \"%s\"", name.c_str() );
  rval = iFace.write_mesh( name.c_str(), &list[0], list.size() ); 
  if (MB_SUCCESS != rval)
  {
    printerror( "Failed to write individual file: \"%s\"", name.c_str() );
    exit( 127 );
  }
  printerror( "Wrote individual file: \"%s\"", name.c_str() );
*/  
  MPI_Barrier( MPI_COMM_WORLD );  
  
    // Write combined file
  std::vector<std::string> qa;
  qa.push_back( "MOAB Parallel HDF5 Writer." );
  time_t t = time( 0 );
  qa.push_back( ctime( &t ) );
  std::string qa3( "Processor " );
  qa3 += rank;
  qa.push_back( qa3 );
  WriteHDF5Parallel writer( &iFace );
  
  printerror ("Writing conbined file: \"%s\"", argv[1] );
  rval = writer.write_file( argv[1], true, &list[0], list.size(), qa );
  if (MB_SUCCESS != rval)
  {
    printerror( "Failed to write combined file: \"%s\"", argv[1] );
    exit( 127 );
  }
  printerror( "Wrote combined file: \"%s\"", argv[1] );
  return 0;
}

  
      
void create_interface( MBEntityHandle geom, int proc )
{
   MBEntityHandle iface_set;
   MBErrorCode rval;
   
   int procs[2];
   if (((proc % 1) == (rank % 1)) == (proc > rank))
   {
     procs[0] = proc;
     procs[1] = rank;
   }
   else
   {
     procs[0] = rank;
     procs[1] = proc;
   }

   rval = iFace.create_meshset( MESHSET_SET, iface_set ); assert(!rval);
   rval = iFace.tag_set_data( ifaceTag, &geom, 1, procs ); assert(!rval);
   std::vector<MBEntityHandle> children;
   rval = iFace.get_child_meshsets( geom, children, 0 ); assert(!rval);
   children.push_back(geom);
   rval = iFace.add_entities( iface_set, &children[0], children.size() ); assert(!rval);
   
   for (std::vector<MBEntityHandle>::iterator iter = children.begin();
        iter != children.end(); ++iter)
   {
     std::vector<MBEntityHandle> contents;
     rval = iFace.get_entities_by_handle( *iter, contents ); assert(!rval);
     rval = iFace.tag_set_data( gidTag, &contents[0], contents.size(), &contents[0] ); assert(!rval); 
     int dim;
     rval = iFace.tag_get_data( geomTag, &*iter, 1, &dim ); assert(!rval);
     rval = iFace.tag_set_data( pGeomTag, &*iter, 1, &dim ); assert(!rval);
   }
}
