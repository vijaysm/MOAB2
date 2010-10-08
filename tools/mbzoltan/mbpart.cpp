#include "MBZoltan.hpp"
#include "moab/Core.hpp"

#include <iostream>
#include <stdlib.h>
#include <list>

using namespace moab;

const char DEFAULT_ZOLTAN_METHOD[] = "RCB";
const char ZOLTAN_PARMETIS_METHOD[] = "PARMETIS";
const char ZOLTAN_OCTPART_METHOD[] = "OCTPART";

const char usage[] = " <# parts> <infile> <outfile> [-<n>] [{-z|-p|-o} <method>] [-s|-t|-b]";

void usage_err( const char* argv0 ) {
  std::cerr << argv0 << usage << std::endl;
  std::cerr << argv0 << " -h   # for help" << std::endl;
  exit(1);
}

void help( const char* argv0 )
{
  std::cout << argv0 << usage << std::endl << std::endl
            << "-<n>  Specify dimension of entities to partition.  Default is " << std::endl
            << "       largest dimension entity present in input file" << std::endl
            << "-z    Specify Zoltan partition method.  One of RCB, RIB, HFSC, PHG, " <<std::endl
            << "        or Hypergraph (PHG and Hypergraph are synonymous)." << std::endl
            << "        Default is: " << DEFAULT_ZOLTAN_METHOD << std::endl
            << "-p    Specify Parmetis partition method" << std::endl
            << "-o    Specify OctPart partition method" << std::endl
            << "-s    Write partition as tagged sets (Default)" << std::endl
            << "-t    Write partition by tagging entities" << std::endl
            << "-b    Write partition as both tagged sets and tagged entities" << std::endl
            << std::endl;
  exit(0);
}

int main( int argc, char* argv[] )
{
  int err = MPI_Init( &argc, &argv );
  if (err) {
    std::cerr << "MPI_Init failed.  Aborting." << std::endl;
    return 3;
  }

  Core moab;
  Interface& mb = moab;
  MBZoltan tool(&mb, false, argc, argv);
 
  int num_args = 0;
  const char* args[3] = {0}; // non-flag options (num parts, input_file, output_file)
  const char* zoltan_method = DEFAULT_ZOLTAN_METHOD;
  const char* other_method = 0;
  std::list<const char**> expected;
  int i, j, part_dim = -1;
  long num_parts;
  bool write_sets = true, write_tags = false;
  bool no_more_flags = false;
  
  for (i = 1; i < argc; ++i) {
    if (!expected.empty()) {
      if (!argv[i][0] || (!no_more_flags && argv[i][0] == '-')) {
        std::cerr << "Expected partition method.  Got \"" << argv[i] << '"' << std::endl;
        usage_err(argv[0]);
      }
      *(expected.front()) = argv[i];
      expected.pop_front();
    }
    else if (!no_more_flags && argv[i][0] == '-') {
      for (j = 1; argv[i][j]; ++j) switch (argv[i][j]) {
        case 'h':
          help(argv[0]);
          break;
        case '-':
          no_more_flags = true;
          break;
        case '0': case '1': case '2': case '3':
          part_dim = argv[i][1] - '0';
          break;
        case 'z':
          expected.push_back(&zoltan_method);
          break;
        case 'p':
          zoltan_method = ZOLTAN_PARMETIS_METHOD;
          expected.push_back(&other_method);
          break;
        case 'o':
          zoltan_method = ZOLTAN_OCTPART_METHOD;
          expected.push_back(&other_method);
          break;
        case 's':
          write_sets = true;
          write_tags = false;
          break;
        case 't':
          write_sets = false;
          write_tags = true;
          break;
        case 'b':
          write_sets = true;
          write_tags = true;
          break;
        default:
          std::cerr << "Unknown option: \"" << argv[i] << '"' << std::endl;
          usage_err(argv[0]);
          break;
      }
    }
    else {
      if (num_args >= (int)(sizeof(args)/sizeof(args[0]))) {
        std::cerr << "Unexpected argument: \"" << argv[i] << '"' << std::endl;
        usage_err(argv[0]);
      }
      
      args[num_args++] = argv[i];
    }
  }
  if (!expected.empty()) {
    std::cerr << "Expected partition method before end of argument list" << std::endl;
    usage_err(argv[0]);
  }
  
    // interpret first numerical argument as number of parts
  for (i = 0; i < num_args; ++i) {
    char* endptr = 0;
    num_parts = strtol( args[i], &endptr, 0 );
    if (!*endptr && num_parts > 0) 
      break;
  }
  if (i == num_args) {
    std::cerr << "Number of parts not specified or not understood" << std::endl;
    usage_err(argv[0]);
  }
    // shift down remaining args
  for (++i; i < num_args; ++i)
    args[i-1] = args[i];
  --num_args;
  
  if (num_args < 2) {
    const char* type = num_args ? "input" : "output";
    std::cerr << "No " << type << " file specified." << std::endl;
    usage_err(argv[0]);
  }
  const char* input_file = args[0];
  const char* output_file = args[1];
  
  ErrorCode rval = mb.load_file( input_file );
  if (MB_SUCCESS != rval) {
    std::cerr << input_file << " : failed to read file." << std::endl;
    std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
    std::string errstr;
    mb.get_last_error(errstr);
    if (!errstr.empty())
      std::cerr << "  Error message: " << errstr << std::endl;
    return 2;
  }
  
  if (part_dim < 0) {
    for (int dim = 3; dim >= 0; --dim) {
      int n;
      rval = mb.get_number_entities_by_dimension( 0, dim, n );
      if (MB_SUCCESS == rval && 0 != n) {
        part_dim = dim;
        break;
      }
    }
  }
  if (part_dim < 0) {
    std::cerr << input_file << " : file does not contain any mesh entities" << std::endl;
    return 2;
  }
  
  rval = tool.partition_mesh( num_parts, zoltan_method, other_method,
                              write_sets, write_tags, part_dim );
  if (MB_SUCCESS != rval) {
    std::cerr << "Partitioner failed!" << std::endl;
    std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
    std::string errstr;
    mb.get_last_error(errstr);
    if (!errstr.empty())
      std::cerr << "  Error message: " << errstr << std::endl;
    return 3;
  }
  
  rval = mb.write_file( output_file );
  if (MB_SUCCESS != rval) {
    std::cerr << output_file << " : failed to write file." << std::endl;
    std::cerr << "  Error code: " << mb.get_error_string(rval) << " (" << rval << ")" << std::endl;
    std::string errstr;
    mb.get_last_error(errstr);
    if (!errstr.empty())
      std::cerr << "  Error message: " << errstr << std::endl;
    return 2;
  }
  
  return 0;    
}
