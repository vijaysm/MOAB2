/* 
   This example takes a .cub mesh file as input and prints out the total area of 
   meshes in each surface of the model. It works for tri and quad elements.
   It makes use of CUBIT's reserved tag - GEOM_DIMENSION and GLOBAL_ID tag.
   Both GLOBAL_ID & GEOM_DIMENSION tag are associated with all the geometric 
   entities in a .cub file. Note: The program would give incorrect result for a 
   non-convex element, since it breaks polygons into triangles for computing the area
*/

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include <iostream>

using namespace moab;

double compute_area(std::vector<EntityHandle>&);

// instantiate
Interface *mb;

int main(int argc, char **argv) {
  if (1 == argc) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }
  
  // declare variables
  Tag gtag, idtag;
  ErrorCode rval;
  const char *tag_geom = "GEOM_DIMENSION";
  const char *tag_gid = "GLOBAL_ID";
  Range sets;
  std::vector<EntityHandle> ents;
  
  // load a file
  mb = new Core();
  rval = mb->load_file(argv[1]);
  if (MB_SUCCESS != rval) return 1;

  // get the tag handle for the tags
  rval = mb->tag_get_handle(tag_geom, 1, MB_TYPE_INTEGER, gtag);
  if (MB_SUCCESS != rval) return 1;
  rval = mb->tag_get_handle(tag_gid, 1, MB_TYPE_INTEGER, idtag);
  if (MB_SUCCESS != rval) return 1;

  // get all the sets with GEOM_DIMESION tag
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &gtag,
					  NULL, 1, sets);
  if (MB_SUCCESS != rval) return 1;
 
  // iterate over each set, getting entities
  Range::iterator set_it;
    
  // loop thru all the geometric entity sets
  for (set_it = sets.begin(); set_it != sets.end(); set_it++)  {

    EntityHandle this_set = *set_it;
    
    // get the id for this set
    int set_id;
    rval = mb->tag_get_data(gtag, &this_set, 1, &set_id);
    if (MB_SUCCESS != rval) return 1;

    // check if it is a surface entities (GEOM_DIMENSION 2) then compute area
    if (set_id ==2){
      
      // area of a surface
      double total_area=0.0;

      //get the global id of this surface
      int gid = 0;
      rval = mb->tag_get_data(idtag, &this_set, 1, &gid);
      if (MB_SUCCESS != rval) return 1;

      // get all entities with dimension 2 in ents
      rval = mb->get_entities_by_dimension(this_set, 2, ents);
      if (MB_SUCCESS != rval) return 1;
 
      // compute the area
      total_area = compute_area(ents);

      ents.clear();
      
      std::cout << "Total area of meshes in surface " << gid << " =  " << total_area << std::endl;
    }
  }
}

// This routine takes all the element entities in a face as input and computes the surface area
// iterating over each element
double compute_area(std::vector<EntityHandle> & entities){
 
  ErrorCode rval= MB_SUCCESS;
  double area = 0.0;
 
  // loop thro' all the elements
  for (int i=0;i<int(entities.size());i++){
    std::vector<EntityHandle> conn;
    EntityHandle handle = entities[i];

    // get the connectivity of this element
    rval = mb->get_connectivity(&handle, 1, conn);
    if (MB_SUCCESS != rval) return -1.0;

    // break polygon into triangles and sum the area - Limitation: Convex polygon
    for (int j = 2; j<=int(conn.size()); ++j){

      EntityHandle vertices[3]={conn[0], conn[j-1],  conn[j-2]};
      CartVect coords[3];
      
      // get 3 coordinates forming the triangle
      rval = mb->get_coords(vertices, 3, coords[0].array());
      if (MB_SUCCESS != rval) return -1.0;
      
      CartVect edge0 = coords[1] - coords [0];
      CartVect edge1 = coords [2] - coords[0];
      
      // using MBCarVect overloaded operators and computing triangle area
      area+=(edge0*edge1).length()/2.0;
    }
  }
  // clear the entities, else old entities remain
  entities.clear();
  return area;
}
