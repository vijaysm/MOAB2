#include "MBRange.hpp"
#include "MBCore.hpp"
#include "MBSkinner.hpp"
#include <iostream>
#include <stdlib.h>

enum {
  NO_ERROR= 0,
  SYNTAX_ERROR = 1,
  FILE_IO_ERROR = 2,
  INTERNAL_ERROR = 3
};

const char* DEFAULT_TAG_NAME = "depth";

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << "[-t <tag name] <input_file> <output_file>" << std::endl
                         << argv0 << "-h" << std::endl;
  exit(SYNTAX_ERROR);
}

void help( const char* argv0 )
{
  std::cout << argv0 << "[-t <tag name] <input_file> <output_file>" << std::endl
            << "Tag elements with integer depth from boundary" << std::endl
            << "-t : tag in which to store depth (default: \"" << DEFAULT_TAG_NAME << "\")" << std::endl;
  exit(NO_ERROR);
}

void check( MBErrorCode rval )
{
  if (MB_SUCCESS != rval) {
    std::cerr << "Internal error.  Aborting." << std::endl;
    exit(INTERNAL_ERROR);
  }
}

void tag_depth( MBInterface& moab, MBTag tag );

int main( int argc, char* argv[] )
{
  const char *input = 0, *output = 0, *tagname = DEFAULT_TAG_NAME;
  bool expect_tag_name = false;
  for (int i = 1; i < argc; ++i) {
    if (expect_tag_name) {
      tagname = argv[i];
      expect_tag_name = false;
    }
    else if (!strcmp("-t", argv[i]))
      expect_tag_name = true;
    else if (input == 0)
      input = argv[i];
    else if (output == 0)
      output = argv[i];
    else {
      std::cerr << "Unexpected argument: '" << argv[i] << "'" << std::endl;
      usage(argv[0]);
    }
  }
  
  if (expect_tag_name) {
    std::cerr << "Expected argument following '-t'" << std::endl;
    usage(argv[0]);
  }
  if (!input) {
    std::cerr << "No input file" << std::endl;
    usage(argv[0]);
  }
  if (!output) {
    std::cerr << "No output file" << std::endl;
    usage(argv[0]);
  }

  MBCore moab;
  MBInterface& mb = moab;
  
  MBEntityHandle file;
  MBErrorCode rval;
  rval = mb.load_file( input, file );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to load file: " << input << std::endl;
    return FILE_IO_ERROR;
  }
  
  int init_val = -1;
  MBTag tag;
  rval = mb.tag_create( tagname, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &init_val );
  if (MB_ALREADY_ALLOCATED == rval) {
    rval = mb.tag_delete( tag ); check(rval);
    rval = mb.tag_create( tagname, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &init_val );
    check(rval);
  }
  
  tag_depth( mb, tag );
  
  rval = mb.write_file( output, 0, 0, &file, 1 );
  if (rval == MB_SUCCESS)
    std::cout << "Wrote file: " << output << std::endl;
  else {
    std::cerr << "Failed to write file: " << output << std::endl;
    return FILE_IO_ERROR;
  }
  
  return NO_ERROR;
}

MBErrorCode get_adjacent_elems( MBInterface& mb, const MBRange& verts, MBRange& elems )
{
  elems.clear();
  MBErrorCode rval;
  for (int dim = 3; dim > 0; --dim) {
    rval = mb.get_adjacencies( verts, dim, false, elems, MBInterface::UNION );
    if (MB_SUCCESS != rval)
      break;
  }
  return rval;
}

void tag_depth( MBInterface& mb, MBTag tag )
{
  MBErrorCode rval;
  int dim;
  
  MBSkinner tool(&mb);
  MBRange verts, elems;
  dim = 3;
  while (elems.empty()) {
    rval = mb.get_entities_by_dimension( 0, dim, elems ); check(rval);
    if (--dim == 0)
      return; // no elements
  }
  rval = tool.find_skin( elems, 0, verts ); check(rval);
  rval = get_adjacent_elems( mb, verts, elems ); check(rval);
  
  std::vector<int> data;
  int val, depth = 0;
  while (!elems.empty()) {
    data.clear();
    data.resize( elems.size(), depth++ );
    rval = mb.tag_set_data( tag, elems, &data[0] ); check(rval);
    
    verts.clear();
    rval = mb.get_adjacencies( elems, 0, false, verts, MBInterface::UNION );
    check(rval);
    
    MBRange tmp;
    rval = get_adjacent_elems( mb, verts, tmp ); check(rval);
    elems.clear();
    for (MBRange::reverse_iterator i = tmp.rbegin(); i != tmp.rend(); ++i) {
      rval = mb.tag_get_data( tag, &*i, 1, &val ); check(rval);
      if (val == -1)
        elems.insert( *i );
    }
  }
}
