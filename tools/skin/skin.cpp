#include <iostream>
#include <vector>
#include "MBInterface.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"


const char FIXED_TAG[] = "fixed"; 

#define CHKERROR( A ) do { if (MB_SUCCESS != (A)) { \
 std::cerr << "Internal error at line " << __LINE__ << std::endl; \
 return 3; } } while(false)

int main( int argc, char* argv[] )
{
  if (argc != 3)
  {
    std::cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << std::endl;
    return 1;
  }
  const char* input_file = argv[1];
  const char* output_file = argv[2];
  
  MBErrorCode result;
  MBCore mbimpl;
  MBInterface* iface = &mbimpl;
  MBSkinner tool( iface );
  
    // read input file
  result = iface->load_mesh( input_file );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to load \"" << input_file << "\"." << std::endl; 
    return 2;
  }
  std::cerr << "Read \"" << input_file << "\"" << std::endl;
  
    // get entities of largest dimension
  int dim = 3;
  MBRange entities;
  while (entities.empty() && dim > 1)
  {
    result = iface->get_entities_by_dimension( 0, dim--, entities );
    CHKERROR(result);
  }
  
    // skin the mesh
  MBRange forward_lower, reverse_lower;
  result = tool.find_skin( entities, forward_lower, reverse_lower );
  MBRange boundary;
  boundary.merge( forward_lower );
  boundary.merge( reverse_lower );
  if (MB_SUCCESS != result || boundary.empty())
  {
    std::cerr << "Mesh skinning failed." << std::endl;
    return 3;
  }
  std::cerr << "Mesh skinning successful." << std::endl;
  
    // get tag handle
  MBTag tag;
  result = iface->tag_get_handle( FIXED_TAG, tag );
  if (result == MB_SUCCESS)
  {
    int size;
    MBDataType type;
    iface->tag_get_size(tag, size);
    iface->tag_get_data_type(tag, type);
    
    if (size != sizeof(int) || type != MB_TYPE_INTEGER)
    {
      std::cerr << '"' << FIXED_TAG << "\" tag defined with incorrect size or type" << std::endl;
      return 3;
    }
  }
  else if (result == MB_TAG_NOT_FOUND)
  {
    int zero = 0;
    result = iface->tag_create( FIXED_TAG, sizeof(int), MB_TAG_DENSE, MB_TYPE_INTEGER, tag, &zero );
    CHKERROR(result);
  }
  else
  {
    CHKERROR(result);
  }
  
    // Set tags
  std::vector<int> ones;
  for (MBRange::iterator i = boundary.begin(); i != boundary.end(); ++i)
  {
    const MBEntityHandle* vertices = 0;
    int num_vertices;
    result = iface->get_connectivity( *i, vertices, num_vertices );
    CHKERROR(result);
    
    if (num_vertices > 0 && ones.size() < (unsigned)num_vertices)
      ones.resize( num_vertices, 1 );
    result = iface->tag_set_data( tag, vertices, num_vertices, &ones[0] );
    CHKERROR(result);
  } 
  std::cerr << "Boundary vertices marked." << std::endl;
  
    // Create a mesh set containing only the highest-dimension elements
  MBEntityHandle set;
  result = iface->create_meshset( MESHSET_SET, set );
  CHKERROR(result);
  result = iface->add_entities( set, entities );
  CHKERROR(result);
    
    // write output file
  result = iface->write_mesh( output_file, &set, 1 );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to write \"" << output_file << "\"." << std::endl; 
    return 2;
  }
  std::cerr << "Wrote \"" << output_file << "\"" << std::endl;
  
  return 0;
}

  
      
  
  
  

