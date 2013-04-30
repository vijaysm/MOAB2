#include <iostream>
#include <stdlib.h>
#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"

// Hold edges in an array of vertex handles.
struct edge {
  MBEntityHandle v0;                                                                                                                      
  MBEntityHandle v1;
};
                                                                                                                                         
// edge structure comparision function for qsort
// If the first vertex handle is the same, compare the second.
int compare_edge(const void *a, const void *b) {
  struct edge *ia = (struct edge *)a;
  struct edge *ib = (struct edge *)b;
  if(ia->v0 == ib->v0) {
    return (int)(100.f*ia->v1 - 100.f*ib->v1);
  } else {
    return (int)(100.f*ia->v0 - 100.f*ib->v0);
  }
}

// This skinner is fast partly because it assumes that no edges exist in the MOAB 
// instance. Checking to see if an edge exists before creating a new one is slow. 
MBErrorCode skin_tris( MBInterface *mb, MBRange tris, MBRange &skin_edges ) {    
  
  // Empty the output range and make sure that the input range is only tris
  skin_edges.clear(); 
  if(tris.empty()) return MB_ENTITY_NOT_FOUND;
  if(!tris.all_of_type(MBTRI)) return MB_FAILURE;

  // Remove edges from the instance.
  int n_edges;
  MBErrorCode rval = mb->get_number_entities_by_type( 0, MBEDGE, n_edges );
  if(MB_SUCCESS != rval) return rval;
  if(0 != n_edges) {
    std::cerr << "skin_tris: failed because " << n_edges 
              << " edges exist in the MOAB instance" << std::endl;
    return MB_FAILURE;
  }      

  // Get connectivity. Do not create MBEdges.
  edge *edges = new edge[3*tris.size()];
  int n_verts;
  int ii = 0;
  for(MBRange::iterator i=tris.begin(); i!=tris.end(); i++) {
    const MBEntityHandle *conn;
    rval = mb->get_connectivity( *i, conn, n_verts );
    if(MB_SUCCESS != rval) return rval;
    if(3 != n_verts) return MB_FAILURE;
    // points should not be degenerate
    if(conn[0]==conn[1] || conn[1]==conn[2] || conn[2]==conn[0]) {
      std::cerr << "skin_tris: degenerate triangle" << std::endl;
      return MB_FAILURE;
    }

    // make edges
    edges[3*ii+0].v0 = conn[0];
    edges[3*ii+0].v1 = conn[1];
    edges[3*ii+1].v0 = conn[1];
    edges[3*ii+1].v1 = conn[2];
    edges[3*ii+2].v0 = conn[2];
    edges[3*ii+2].v1 = conn[0];
    ii++;
  }

  // Ensure that the first vertex handle is the lowest
  for(unsigned int i=0; i<3*tris.size(); ++i) {
    if(edges[i].v0 > edges[i].v1) {
      MBEntityHandle temp = edges[i].v0;
      edges[i].v0 = edges[i].v1;
      edges[i].v1 = temp;
    }
  }

  // Sort by first handle, then second handle.
  qsort(edges, 3*tris.size(), sizeof(struct edge), compare_edge);    

  // Go through array, saving edges that are not paired.
  for(unsigned int i=0; i<3*tris.size(); i++) {
    // If the last edge has not been paired, create it. This avoids overrunning
    // the edges array with i+1.
    if(3*tris.size()-1 == i) {
      const MBEntityHandle conn[2] = {edges[i].v0, edges[i].v1};
      MBEntityHandle edge;
      rval = mb->create_element( MBEDGE, conn, 2, edge );
      if(MB_SUCCESS != rval) return rval;
      skin_edges.insert(edge);
    
    // If a match exists, skip ahead
    } else if(edges[i].v0==edges[i+1].v0 && edges[i].v1==edges[i+1].v1) {
      i++;
      // test to make sure surface is manifold
      while( edges[i].v0==edges[i+1].v0 && edges[i].v1==edges[i+1].v1 ) {
        std::cout << "find_skin WARNING: non-manifold edge" << std::endl;
        mb->list_entity( edges[i].v0 );
        mb->list_entity( edges[i].v1 );
        ++i;
      }
    // otherwise a skin edge has been found
    } else {
      const MBEntityHandle conn[2] = {edges[i].v0, edges[i].v1};
      MBEntityHandle edge;
      rval = mb->create_element( MBEDGE, conn, 2, edge );
      if(MB_SUCCESS != rval) return rval;
      skin_edges.insert( edge );
    } 
  }
  delete[] edges;
  return MB_SUCCESS;
}


// Skin triangles to recover edges.
// Triangles are contained in surface sets.
int main(int argc, char **argv) {
  if (1 == argc) {
  std::cout << "Usage: " << argv[0] << " <filename>" << std::endl;
    return 0;
  }     

  // get MOAB instance and read the file                                                                                                  
  MBCore *mb = new MBCore();
  MBErrorCode rval = mb->load_file(argv[1]);
  if(MB_SUCCESS != rval) return 0;

  // this optimized skinner requires removing all MBEdges from the MOAB instance
  MBRange edges;
  rval = mb->get_entities_by_type( 0, MBEDGE, edges );
  if(MB_SUCCESS != rval) return 0;
  if( !edges.empty() ) std::cout << "Warning: deleting all MBEdges" << std::endl;
  rval = mb->delete_entities( edges );
  if(MB_SUCCESS != rval) return 0;

  // get surface sets
  MBTag geom_tag;
  rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, geom_tag );
  if(MB_SUCCESS != rval) return 0;
  MBRange surf_sets;
  int two = 2;
  void *dim[] = {&two};                                                                                   
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,                                                                
                                           dim, 1, surf_sets );
  if(MB_SUCCESS != rval) return 0;

  // skin each surface
  for(MBRange::iterator i=surf_sets.begin(); i!=surf_sets.end(); ++i) {

    // get triangles in the surface set
    MBRange tris;
    rval = mb->get_entities_by_type( *i, MBTRI, tris );
    if(MB_SUCCESS != rval) return 0;

    // call the skinning function
    MBRange skin_edges;
    rval = skin_tris( mb, tris, skin_edges );
    if(MB_SUCCESS != rval) return 0;

    // do something with the result
    std::cout << "surface has " << skin_edges.size() << " skin edges" << std::endl;

    // remove the edges for the optimized skinner
    rval = mb->delete_entities( skin_edges );
    if(MB_SUCCESS != rval) return 0;    
  }
}  
