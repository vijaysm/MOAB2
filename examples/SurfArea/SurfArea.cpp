/* 
   This example takes a .cub mesh file as input and prints out the total area of 
   meshes in each surface of the model. It works for tri and quad elements.
   It makes use of CUBIT's reserved tag - GEOM_DIMENSION and GLOBAL_ID tag.
   Both GLOBAL_ID & GEOM_DIMENSION tag are associated with all the geometric 
   entities in a .cub file. Note: The program would give incorrect result for a 
   non-convex element, since it breaks polygons into triangles for computing the area
*/

#include "moab/MBCore.hpp"
#include "moab/MBRange.hpp"
#include "MBCartVect.hpp"
#include <iostream>


double compute_area(std::vector<MBEntityHandle>&);

// instantiate
MBInterface *mb;

int main(int argc, char **argv) {
  if (1 == argc) {
    std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }
  
  // get tag 
  MBTag gtag, idtag;
  MBErrorCode rval;
  const char *tag_geom = "GEOM_DIMENSION";
  const char *tag_gid = "GLOBAL_ID";
  MBRange sets;
  std::vector<MBEntityHandle> ents;
  
  // load a file
  mb = new MBCore();
  rval = mb->load_file(argv[1]);

  // get the tag handle for the tags
  rval = mb->tag_get_handle(tag_geom, gtag);
  rval = mb->tag_get_handle(tag_gid, idtag);

  // get all the sets with GEOM_DIMESION tag
  rval = mb->get_entities_by_type_and_tag(0, MBENTITYSET, &gtag,
					  NULL, 1, sets);
 
  // iterate over each set, getting entities
  MBRange::iterator set_it;
    
  // loop thru all the geometric entity sets
  for (set_it = sets.begin(); set_it != sets.end(); set_it++)  {

    MBEntityHandle this_set = *set_it;
    
    // get the id for this set
    int set_id;
    rval = mb->tag_get_data(gtag, &this_set, 1, &set_id);

    // check if it is a surface entities (GEOM_DIMENSION 2) then compute area
    if (set_id ==2){
      
      // area of a surface
      double total_area=0.0;

      //get the global id of this surface
      int gid = 0;
      rval = mb->tag_get_data(idtag, &this_set, 1, &gid);

      // get all entities with dimension 2 in ents
      rval = mb->get_entities_by_dimension(this_set, 2, ents);
 
      // compute the area
      total_area = compute_area(ents);

      ents.clear();
      
      std::cout << "Total area of meshes in surface " << gid << " =  " << total_area << std::endl;
    }
  }
}

// This routine takes all the element entities in a face as input and computes the surface area
// iterating over each element
double compute_area(std::vector<MBEntityHandle> & entities){
 
  int rval= 0;
  double area = 0.0;
  double coord[9];
 
  // loop thro' all the elements
  for (int i=0;i<entities.size();i++){
    std::vector<MBEntityHandle> conn;
    MBEntityHandle handle = entities[i];

    // get the connectivity of this element
    rval = mb->get_connectivity(&handle, 1, conn);

    // break polygon into triangles and sum the area - Limitation: Convex polygon
    for (int j = 2; j<=conn.size(); ++j){

      MBEntityHandle vertices[3]={conn[0], conn[j-1],  conn[j-2]};
      MBCartVect coords[3];
      
      // get 3 coordinates forming the triangle
      rval = mb->get_coords(vertices, 3, coords[0].array());
      
      MBCartVect edge0 = coords[1] - coords [0];
      MBCartVect edge1 = coords [2] - coords[0];
      
      // using MBCarVect overloaded operators and computing triangle area
      area+=(edge0*edge1).length()/2.0;
    }
  }
  // clear the entities, else old entities remain
  entities.clear();
  return area;
}
