#include "moab_mpi.h"
#include "moab/Core.hpp"
#include <sstream>
#include <iostream>
#include <stdlib.h>

const char usage[] = "[-b|-d|-f] [-p <name>] [-R] [-g <level>] <filename>";

const char DEFAULT_PARTITION_TAG[] = "PARALLEL_PARTITION";

void error( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " " << usage << std::endl;
  exit(1);
}

void help( const char* argv0 ) {
  std::cout << argv0 << " " << usage << std::endl
            << "-p <name>  Tag identifying partition sets (default: \"" << 
                  DEFAULT_PARTITION_TAG << "\")" << std::endl
            << "-a         Assign partitions to processes by matching part number to rank" << std::endl
            << "-R         Do not resolve shared entities" << std::endl
            << "-b         Use broadcast & delete read method" << std::endl
            << "-d         Use read & delete method" << std::endl
            << "-f         Use true parallel read (default)" << std::endl
            << "-g         Set debug output level" << std::endl;
  exit(0);
}

void finalize() { MPI_Finalize(); }

int main( int argc, char* argv[] )
{
  MPI_Init( &argc, &argv );
  atexit( &finalize );
  
  const char BCAST_MODE[] = "BCAST_DELETE";
  const char DELETE_MODE[] = "READ_DELETE";
  const char PART_MODE[] = "READ_PART";
  const char* read_mode = PART_MODE;
  const char* partition_tag_name = DEFAULT_PARTITION_TAG;
  bool assign_by_id = false;
  bool resolve_shared = true;
  int debug_level = 0;
  const char* filename = 0;
  
  bool expect_tag = false;
  bool expect_level = false;
  bool no_more_flags = false;
  for (int i = 1; i < argc; ++i) {
    if (expect_tag) {
      partition_tag_name = argv[i];
      expect_tag = false;
    }
    else if (expect_level) {
      char* endptr;
      debug_level = (int)strtol( argv[i], &endptr, 0 );
      if (*endptr || endptr == argv[i] || debug_level < 0) {
        std::cerr << "Expected positive integer value following '-g' flag" << std::endl;
        error(argv[0]);
      }
      expect_level = false;
    }
    else if (argv[i][0] == '-' && !no_more_flags) {
      for (int j = 1; argv[i][j]; ++j) {
        switch (argv[i][j]) {
          case '-': no_more_flags = true; break;
          case 'p': expect_tag = true; break;
          case 'a': assign_by_id = true; break;
          case 'R': resolve_shared = false; break;
          case 'b': read_mode = BCAST_MODE; break;
          case 'd': read_mode = DELETE_MODE; break;
          case 'f': read_mode = PART_MODE; break;
          case 'g': expect_level = true; break;
          case 'h': help(argv[0]); break;
          default:
            std::cerr << "Unknown flag: -" << argv[i][j] << std::endl;
           error(argv[0]);
        }
      }
    }
    else if (filename) {
      std::cerr << "Unexpected argument: \"" << argv[i] << "\"" << std::endl;
      error(argv[0]);
    }
    else {
      filename = argv[i];
    }
  }
  
  if (expect_tag) {
    std::cerr << "Expected value following -p flag" << std::endl;
    error(argv[0]);
  }
  if (expect_level) {
    std::cerr << "Expected value following -g flag" << std::endl;
    error(argv[0]);
  }
  if (!filename) {
    std::cerr << "No file name specified" << std::endl;
    error(argv[0]);
  }
  
  std::ostringstream options;
  options << "PARTITION=" << partition_tag_name
          << ";PARALLEL=" << read_mode;
  if (resolve_shared)
    options << ";RESOLVE_SHARED_ENTS";
  if (assign_by_id)
    options << ";PARTITION_BY_RANK";
  if (debug_level)
    options << ";DEBUG_IO=" << debug_level;
  
  moab::Core core;
  moab::Interface& mb = core;
  
  std::string opts = options.str();
  std::cout << "Reading \"" << filename << "\" with options=\""
            << opts << "\"." << std::endl;
  double init_time = MPI_Wtime();
  moab::ErrorCode rval = mb.load_file( filename, 0, opts.c_str() );
  double fini_time = MPI_Wtime();
  
  if (moab::MB_SUCCESS == rval) {
    std::cout << "Success!" << std::endl;
  }
  else {
    std::string estr = mb.get_error_string(rval);
    std::string msg;
    mb.get_last_error( msg );
    std::cout << "Read failed with error code " << rval << " ("
              << estr << ")" << std::endl;
    if (!msg.empty()) 
      std::cerr << '"' << msg << '"' << std::endl;
  }
    
  std::cout << "Read returned in " << fini_time - init_time << " seconds" << std::endl;
  return (moab::MB_SUCCESS != rval);
}
