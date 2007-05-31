#include "MBCore.hpp"
#include "MBAdaptiveKDTree.hpp"
#include <errno.h>
#include <stdio.h>

int main( int argc, char* argv[] ) 
{
  double* values;
  unsigned long length;
  FILE* file;
  unsigned long count;
  clock_t t;
  
  if (argc  != 4 && argc != 3) {
    fprintf(stderr,"usage: %s <tree_file> <point_file> <count>\n",argv[0]);
    return 1;
  }
  
  file = fopen( argv[2], "rb" );
  if (!file) {
    perror(argv[1]);
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
  
  if (argc == 4) 
    count = strtol( argv[3], 0, 0 );
  else
    count = length;
 
  printf("Loading tree..."); fflush( stdout );
  t = clock();
  MBCore moab;
  MBErrorCode rval = moab.load_mesh( argv[1], 0, 0 );
  if (MB_SUCCESS != rval) {
    fprintf(stderr,"Failed to read file: %s\n", argv[1] );
    return 2;
  }
  printf("%0.2f seconds\n", (clock()-t)/(double)CLOCKS_PER_SEC); fflush( stdout );
  
  printf("Retrieving tree..."); fflush( stdout );
  t = clock();
  MBRange range;
  MBAdaptiveKDTree tool(&moab);
  tool.find_all_trees(range);
  if (range.size() != 1) {
    fprintf(stderr,"%s : found %d kd-trees\n", argv[1], (int)range.size() );
    return 3;
  }
  const MBEntityHandle root = range.front();
  printf("%0.2f seconds\n", (clock()-t)/(double)CLOCKS_PER_SEC); fflush( stdout );
  
  printf("Running point queries..."); fflush( stdout );
  t = clock();
  MBEntityHandle leaf;
  double pt[3];
  for (unsigned long i = 0; i < count; ++i) {
    const double* coords = values + 3 * (count % length);
    rval = tool.leaf_containing_point( root, coords, leaf );
    //rval = tool.closest_triangle( root, coords, pt, leaf );
    if (MB_SUCCESS != rval) {
      fprintf(stderr, "Failure (MBErrorCode == %d) for point %d (%f, %f, %f)\n",
        (int)rval, (int)(count % length), coords[0], coords[1], coords[2] );
    }
  }
  printf("%0.2f seconds\n", (clock()-t)/(double)CLOCKS_PER_SEC); fflush( stdout );

  return 0;
}

  
