#include <iostream>
#include <vector>
#include "MBInterface.hpp"
#include "MBTagConventions.hpp"
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBSkinner.hpp"


const char FIXED_TAG[] = "fixed"; 

#define CHKERROR( A ) do { if (MB_SUCCESS != (A)) { \
 std::cerr << "Internal error at line " << __LINE__ << std::endl; \
 return 3; } } while(false)

int main( int argc, char* argv[] )
{
  if (argc < 3)
  {
    std::cerr << "Usage: " << argv[0] 
              << " [-b <block_num> [-b ...] ] [-s <sideset_num>] [-t] [-w]"
              << " <input_file> <output_file>" << std::endl;
    std::cerr << "Options: " << std::endl;
    std::cerr << "-b <block_num> : Compute skin only for material set/block <block_num>." << std::endl;
    std::cerr << "-s <sideset_num> : Put skin in neumann set/sideset <sideset_num>." << std::endl;
    std::cerr << "-t : Set 'FIXED' tag on skin vertices." << std::endl;
    std::cerr << "-w : Write out whole mesh (otherwise just writes skin)." << std::endl;
    
    return 1;
  }

  int i = 1;
  std::vector<int> matsets;
  int neuset_num = -1;
  bool write_tag = false, write_whole_mesh = false;
  
  while (i < argc) {
    if (!strcmp(argv[i], "-b")) {
      i++;
      matsets.push_back(atoi(argv[i]));
      i++;
    }
    else if (!strcmp(argv[i], "-s")) {
      i++;
      neuset_num = atoi(argv[i]);
      i++;
    }
    else if (!strcmp(argv[i], "-t")) {
      i++;
      write_tag = true;
    }
    
    else if (!strcmp(argv[i], "-w")) {
      i++;
      write_whole_mesh = true;
    }
    else {
      break;
    }
  }
  
  const char* input_file = argv[i++];
  const char* output_file = argv[i++];
  
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
  int dim = 4;
  MBRange entities;
  while (entities.empty() && dim > 1)
  {
    dim--;
    result = iface->get_entities_by_dimension( 0, dim, entities );
    CHKERROR(result);
  }

  MBRange skin_ents;
  MBTag matset_tag = 0, neuset_tag = 0;
  result = iface->tag_get_handle(MATERIAL_SET_TAG_NAME, matset_tag);
  result = iface->tag_get_handle(NEUMANN_SET_TAG_NAME, neuset_tag);

  if (matsets.empty()) skin_ents = entities;
  else {
      // get all entities in the specified blocks
    if (0 == matset_tag) {
      std::cerr << "Couldn't find any material sets in this mesh." << std::endl;
      return 1;
    }
    
    for (std::vector<int>::iterator vit = matsets.begin(); vit != matsets.end(); vit++) {
      int this_matset = *vit;
      const void *this_matset_ptr = &this_matset;
      MBRange this_range, ent_range;
      result = iface->get_entities_by_type_and_tag(0, MBENTITYSET, &matset_tag,
                                                    &this_matset_ptr, 1, this_range);
      if (MB_SUCCESS != result) {
        std::cerr << "Trouble getting material set #" << *vit << std::endl;
        return 1;
      }
      else if (this_range.empty()) {
        std::cerr << "Warning: couldn't find material set " << *vit << std::endl;
        continue;
      }
      
      result = iface->get_entities_by_dimension(*this_range.begin(), dim, ent_range, true);
      if (MB_SUCCESS != result) continue;
      skin_ents.merge(ent_range);
    }
  }
  
  if (skin_ents.empty()) {
    std::cerr << "No entities for which to compute skin; exiting." << std::endl;
    return 1;
  }
  
    // skin the mesh
  MBRange forward_lower, reverse_lower;
  result = tool.find_skin( skin_ents, forward_lower, reverse_lower );
  MBRange boundary;
  boundary.merge( forward_lower );
  boundary.merge( reverse_lower );
  if (MB_SUCCESS != result || boundary.empty())
  {
    std::cerr << "Mesh skinning failed." << std::endl;
    return 3;
  }

  if (write_tag) {
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
    MBRange bverts;
    result = iface->get_adjacencies(boundary, 0, false, bverts, MBInterface::UNION);
    if (MB_SUCCESS != result) {
      std::cerr << "Trouble getting vertices on boundary." << std::endl;
      return 1;
    }
    ones.resize( bverts.size(), 1 );
    result = iface->tag_set_data( tag, bverts, &ones[0] );
    CHKERROR(result);
  }
  
  if (-1 != neuset_num) {
      // create a neumann set with these entities
    if (0 == neuset_tag) {
      result = iface->tag_create("NEUMANN_SET_TAG_NAME", sizeof(int), MB_TAG_SPARSE,
                                  MB_TYPE_INTEGER, neuset_tag, NULL);
      if (MB_SUCCESS != result || 0 == neuset_tag) return 1;
    }
    

      // always create a forward neumann set, assuming we have something in the set
    MBEntityHandle forward_neuset = 0;
    result = iface->create_meshset(MESHSET_SET, forward_neuset);
    if (MB_SUCCESS != result || 0 == forward_neuset) return 1;
    result = iface->tag_set_data(neuset_tag, &forward_neuset, 1, &neuset_num);
    if (MB_SUCCESS != result) return 1;

    if (!forward_lower.empty()) {
      result = iface->add_entities(forward_neuset, forward_lower);
      if (MB_SUCCESS != result) return 1;
    }
    if (!reverse_lower.empty()) {
      MBEntityHandle reverse_neuset = 1;
      result = iface->create_meshset(MESHSET_SET, reverse_neuset);
      if (MB_SUCCESS != result || 0 == forward_neuset) return 1;

      result = iface->add_entities(reverse_neuset, reverse_lower);
      if (MB_SUCCESS != result) return 1;
      MBTag sense_tag;
      result = iface->tag_get_handle("SENSE", sense_tag);
      if (result == MB_TAG_NOT_FOUND) {
        int dum_sense = 0;
        result = iface->tag_create("SENSE", sizeof(int), MB_TAG_SPARSE, sense_tag, &dum_sense);
      }
      if (result != MB_SUCCESS) return 1;
      int sense_val = -1;
      result = iface->tag_set_data(neuset_tag, &reverse_neuset, 1, &sense_val);
      if (MB_SUCCESS != result) return 0;
      result = iface->add_entities(forward_neuset, &reverse_neuset, 1);
      if (MB_SUCCESS != result) return 0;
    }
  }

  if (write_whole_mesh) {
    
      // write output file
    result = iface->write_mesh( output_file);
    if (MB_SUCCESS != result)
    { 
      std::cerr << "Failed to write \"" << output_file << "\"." << std::endl; 
      return 2;
    }
    std::cerr << "Wrote \"" << output_file << "\"" << std::endl;
  }
  else {
      // write only skin; write them as one set
    MBEntityHandle skin_set;
    result = iface->create_meshset(MESHSET_SET, skin_set);
    if (MB_SUCCESS != result) return 1;
    result = iface->add_entities(skin_set, forward_lower);
    if (MB_SUCCESS != result) return 1;
    result = iface->add_entities(skin_set, reverse_lower);
    if (MB_SUCCESS != result) return 1;

    MBRange this_range, ent_range;
    result = iface->get_entities_by_type_and_tag(0, MBENTITYSET, &matset_tag,
                                                  NULL, 0, this_range);
    if (!this_range.empty()) iface->delete_entities(this_range);

    int dum = 10000;
    result = iface->tag_set_data(matset_tag, &skin_set, 1, &dum);
    

    result = iface->write_mesh( output_file, &skin_set, 1);
    if (MB_SUCCESS != result)
    { 
      std::cerr << "Failed to write \"" << output_file << "\"." << std::endl; 
      return 2;
    }
    std::cerr << "Wrote \"" << output_file << "\"" << std::endl;
  }
  
  
  return 0;
}

  
      
  
  
  

