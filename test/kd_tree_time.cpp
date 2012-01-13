#include "moab/Core.hpp"
#include "moab/AdaptiveKDTree.hpp"
#include "moab/Range.hpp"
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctime>

using namespace moab;

void usage( const char* argv0 )
{
  fprintf(stderr, "usage: %s [-t] [-d <result_file>] <tree_file> <point_file> [<count>]\n", argv0);
  exit(1);
}

void print_file_stats( Interface& moab )
{
  ErrorCode rval;
  int num_tri;
  Range sets;
  unsigned long set_mem, set_am, tag_mem, tag_am;
  
  rval = moab.get_number_entities_by_type( 0, MBTRI, num_tri );
  if (MB_SUCCESS != rval)
    num_tri = -1;
  rval = moab.get_entities_by_type( 0, MBENTITYSET, sets );
  if (MB_SUCCESS != rval)
    sets.clear();
  
  moab.estimated_memory_use( sets, 0, 0, &set_mem, &set_am, 0, 0, 0, 0, &tag_mem, &tag_am );
  printf( "Triangles:   %d\n", num_tri );
  printf( "Sets:        %lu\n", (unsigned long)sets.size() );
  printf( "Set storage: %lu (%lu)\n", set_mem, set_am );
  printf( "Tag storage: %lu (%lu)\n", tag_mem, tag_am );
}

int main( int argc, char* argv[] ) 
{
  double* values;
  unsigned long length;
  FILE *file, *rfile = 0;
  unsigned long count = 0;
  clock_t t;
  
  const char* tree_file = 0;
  const char* point_file = 0;
  const char* result_file = 0;
  bool query_triangles = false;
  
  if (argc < 3 || argc > 7)
    usage(argv[0]);
  
  for (int i = 1; i < argc; ++i) {
    if (!strcmp("-t", argv[i]))
      query_triangles = true;
    else if (!strcmp("-d", argv[i])) {  
      ++i;
      if (i == argc)
        usage(argv[0]);
      result_file = argv[i];
    }     
    else if (!tree_file)
      tree_file = argv[i];
    else if (!point_file)
      point_file = argv[i];
    else {
      char* endptr;
      count = strtol( argv[i], &endptr, 0 );
      if (*endptr || count < 1)
        usage(argv[0]);
    }
  }
  
  file = fopen( point_file, "rb" );
  if (!file) {
    perror(point_file);
    return 2;
  }
  fseek( file, 0, SEEK_END );
  length = ftell(file) / (3*sizeof(double));
  fseek( file, 0, SEEK_SET );
  values = new double[3*length];
  if (length != fread( values, 3*sizeof(double), length, file )) {
    fprintf( stderr, "Error reading %lu points from file \"%s\"\n", length, argv[2] );
    return 2;
  }
  fclose( file );
  
  if (result_file) {
    rfile = fopen( result_file, "w" );
    if (!rfile) {
      perror(result_file);
      return 2;
    }
  }
  
  if (!count)
    count = length;
 
  printf("Loading tree..."); fflush( stdout );
  t = clock();
  Core moab;
  ErrorCode rval = moab.load_mesh( tree_file, 0, 0 );
  if (MB_SUCCESS != rval) {
    fprintf(stderr,"Failed to read file: %s\n", tree_file );
    return 2;
  }
  printf("%0.2f seconds\n", (clock()-t)/(double)CLOCKS_PER_SEC); fflush( stdout );
  
  Range range;
  AdaptiveKDTree tool(&moab);
  tool.find_all_trees(range);
  if (range.size() != 1) {
    fprintf(stderr,"%s : found %d kd-trees\n", argv[1], (int)range.size() );
    return 3;
  }
  const EntityHandle root = range.front();

  print_file_stats( moab );
  
  printf("Running point queries..."); fflush( stdout );
  t = clock();
  EntityHandle leaf;
  double pt[3];
  for (unsigned long i = 0; i < count; ++i) {
    const double* coords = values + 3 * (i % length);
    if (query_triangles)
      rval = tool.closest_triangle( root, coords, pt, leaf );
    else
      rval = tool.leaf_containing_point( root, coords, leaf );
    if (MB_SUCCESS != rval) {
      fprintf(stderr, "Failure (ErrorCode == %d) for point %d (%f, %f, %f)\n",
        (int)rval, (int)(count % length), coords[0], coords[1], coords[2] );
    }
    else if (rfile) {
      if (query_triangles)
        fprintf( rfile, "%f %f %f %f %f %f %lu\n", 
          coords[0], coords[1], coords[2],
          pt[0], pt[1], pt[2], (unsigned long)leaf );
      else
        fprintf( rfile, "%f %f %f %lu\n", 
          coords[0], coords[1], coords[2], (unsigned long)leaf );
    }
  }
  printf("%0.2f seconds\n", (clock()-t)/(double)CLOCKS_PER_SEC); fflush( stdout );

  return 0;
}

  
