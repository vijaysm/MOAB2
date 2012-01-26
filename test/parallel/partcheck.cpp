#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/Skinner.hpp"
#include "moab/CN.hpp"

#include <iostream>
#include <sstream>
#include <limits>
#include <string>
#include <time.h>

#include "moab_mpi.h"

const char PARTTAG[] = "PARALLEL_PARTITION";
const char READ_DEL_OPT[] = "READ_DELETE";
const char BCAST_DEL_OPT[] = "BCAST_DELETE";
const char READ_PART_OPT[] = "READ_PART";
const char* const DEFAULT_MODE = READ_PART_OPT;

const char USAGE[] = " [-v <n>] [-R] [-p <parttag>[=val]] [-D|-B|-P] <input_file>";

static void usage( const char* argv0 ) {
  std::cerr << "Usage: " << argv0 << USAGE << std::endl
            << "       " << argv0 << " -h" << std::endl;
  exit(1);
}

static const char* defstr( const char* mode )
{
  static const char s[] = " (default)";
  int len = (mode == DEFAULT_MODE) ? 0 : sizeof(s)-1;
  return s + len;
}    

static void help( const char* argv0 ) {
  std::cout << "Usage: " << argv0 << USAGE << std::endl
            << "  -v  : debug output level (DEBUG_IO option)" << std::endl
            << "  -R  : do not resolve shared entities" << std::endl
            << "  -p  : partition tag, with optional '=' value" << std::endl
            << "          (default = \"" << PARTTAG << "\")" << std::endl
            << "  -D  : read mode as \"" << READ_DEL_OPT << "\"" << defstr(READ_DEL_OPT) << std::endl
            << "  -B  : read mode as \"" << BCAST_DEL_OPT << "\"" << defstr(READ_DEL_OPT) << std::endl
            << "  -P  : read mode as \"" << READ_PART_OPT << "\"" << defstr(READ_PART_OPT) << std::endl
            ;
  exit(0);
}

using namespace moab;

int check_parallel_read( Interface& mb, ParallelComm* pcomm, bool shared_ents );

int main( int argc, char* argv[] )
{
  const char* read_mode_opt = DEFAULT_MODE;
  const char* part_tag_opt = PARTTAG;
  bool resolve_shared = true;
  const char* input_file = 0;
  const char* debug_flag_str = 0;
  bool no_more_flags = false;

  MPI_Init( &argc, &argv );
  
    // process command line arguments
  
  for (int i = 1; i < argc; ++i) {
    if (!no_more_flags && argv[i][0] == '-') {
      const char* opts = argv[i]+1;
      for (int j = 0; opts[j]; ++j) {
        switch (opts[j]) {
          case '-':
            no_more_flags = true;
            break;
          case 'v':
            if (++i == argc)
              usage(argv[0]);
            debug_flag_str = argv[i];
            break;
          case 'R':
            resolve_shared = false;
            break;
          case 'p':
            if (++i == argc)
              usage(argv[0]);
            part_tag_opt = argv[i];
            break;
          case 'D':
            read_mode_opt = READ_DEL_OPT;
            break;
          case 'B':
            read_mode_opt = BCAST_DEL_OPT;
            break;
          case 'P':
            read_mode_opt = READ_PART_OPT;
            break;
          case 'h':
            help(argv[0]);
            break;
          default:
            usage(argv[0]);
        }
      }
    }
    else if (!input_file) {
      input_file = argv[i];
    }
    else {
      usage(argv[0]);
    }
  }
  
  if (!input_file) {
    std::cerr << argv[0] << ": no input file specified" << std::endl;
    usage(argv[0]);
  }
  
  
    // build options string
  
  std::ostringstream opts;
  opts << "PARALLEL=" << read_mode_opt;
  
  std::string part_opt( part_tag_opt );
  size_t p = part_opt.find_last_of('=');
  if (p == std::string::npos) {
    opts << ";PARTITION=" << part_opt;
  }
  else {
    char* endptr = 0;
    long n = strtol( part_opt.c_str() + p + 1, &endptr, 0 );
    if (*endptr || p == part_opt.size()-1) {
      std::cerr << "Warning: partition tag option contains an '=' followed "
                   "         by a non-integer value.  Assuming tag name is "
                   "         \"" << part_opt << "\"" << std::endl;
      opts << ";PARTITION=" << part_opt;
    }
    else {
      opts << ";PARTITION=" << part_opt.substr( 0, p );
      opts << ";PARTITION_VAL=" << n;
    }
  }
  
  if (resolve_shared) {
    opts << ";PARALLEL_RESOLVE_SHARED_ENTS";
  }
  
  if (debug_flag_str) {
    char* endptr = 0;
    long n = strtol( debug_flag_str, &endptr, 0 );
    if (*endptr || n < 0 || !*debug_flag_str)
      usage(argv[0]);
    opts << ";DEBUG_IO=" << n;
  }
  
  Core moab;
  Interface& mb = moab;
  ParallelComm* pcomm = new ParallelComm( &mb, MPI_COMM_WORLD );
  if (pcomm->rank() == 0) 
    std::cout << "Loading file: \"" << input_file << "\" with options \"" << opts.str() << '"' << std::endl;
  opts << ";PARALLEL_COMM=" << pcomm->get_id();
  
  clock_t init_time = clock();
  ErrorCode rval = mb.load_file( input_file, 0, opts.str().c_str() );
  if (MB_SUCCESS != rval) {
    std::cerr << input_file << " : file read failed (" << mb.get_error_string(rval) << ")" << std::endl;
    return 1;
  }
  
  clock_t t = clock();
  double sec;
  if (t < init_time) 
    sec = (std::numeric_limits<clock_t>::max() - init_time)/(double)CLOCKS_PER_SEC
      + t/(double)CLOCKS_PER_SEC;
  else
    sec = (t - init_time)/(double)CLOCKS_PER_SEC;
  double allsec;
  MPI_Reduce( &sec, &allsec, 1, MPI_DOUBLE, MPI_MAX, 0, pcomm->comm() );
  if (pcomm->rank() == 0) {
    std::cout << "Read completed in " << allsec << " seconds" << std::endl;
  }
  
  int result = check_parallel_read( mb, pcomm, resolve_shared );
  
  if (MB_SUCCESS != pcomm->check_all_shared_handles(false))
    ++result;
  
  MPI_Finalize();
  return result;
}

int check_parallel_read( Interface& mb, ParallelComm* pcomm, bool /*shared_ents*/ )
{
  int error_count = 0;

  const Range& parts = pcomm->partition_sets();
  if (parts.empty()) 
    std::cout << "No parts for process " << pcomm->rank() << std::endl;
  
    // get all entities from parts
  Range part_ents;
  for (Range::iterator i = parts.begin(); i != parts.end(); ++i) 
    mb.get_entities_by_handle( *i, part_ents );
  
  int dim = 3;
  if (part_ents.empty()) 
    std::cout << "No partitioned entities for process " << pcomm->rank() << std::endl;
  else {
    dim = CN::Dimension( mb.type_from_handle( part_ents.back() ) );
    if (!part_ents.all_of_dimension(dim)) 
      std::cout << "Partitioned ents of mixed dimension for process " << pcomm->rank() << std::endl;
  }
  
  Range all_ents;
  mb.get_entities_by_dimension( 0, dim, all_ents );
  if (!subtract( all_ents, part_ents ).empty()) {
    std::cerr << "Process " << pcomm->rank() << " has entities of dimension "
              << dim << " that are not contained in any part" << std::endl;
    ++error_count;
  }
  
  if (dim == 0) {
    std::cout << "Skipping further tests because mesh is vertex-partitioned" << std::endl;
    return error_count;
  }
  
  Range part_verts;
  mb.get_connectivity( part_ents, part_verts );
  Range all_verts;
  mb.get_entities_by_dimension( 0, 0, all_verts );
  if (!subtract( all_verts, part_verts ).empty()) {
    std::cerr << "Process " << pcomm->rank() << " has vertices "
              << " that are not contained in any partitioned element" << std::endl;
    ++error_count;
  }
    
  //if (!shared_ents) {
  //  std::cout << "Skipping further tests because shared entities were not resolved" << std::endl;
  //  return error_count;
  //}
  
  return error_count;
}
