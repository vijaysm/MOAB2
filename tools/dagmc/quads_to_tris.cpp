// Brandon Smith
// January 20, 2009

/*
   loop over all surface meshsets
     get all quads of surface meshset
     loop over all quads
       make tris from quad
       add tris to the surface meshset
       remove quad from surface meshset
       delete quad
     end loop
   end loop
*/

#include "quads_to_tris.hpp"

// Generic function to create two tris from a quad. This can be improved later.
MBErrorCode make_tris_from_quad( MBInterface *MBI,
                                 MBEntityHandle quad,  /* intput */
                                 MBEntityHandle &tri0, /* output */
				 MBEntityHandle &tri1  /* output */) {
  
  // get connectivity (ordered counterclockwise for 2D elements in MOAB)
  MBErrorCode result;    
  const MBEntityHandle *quad_conn;
  int n_verts=0;
  result = MBI->get_connectivity( quad, quad_conn, n_verts );
    assert( 4 == n_verts );
    assert( MB_SUCCESS == result);
   
  // make tris from quad
  MBEntityHandle tri0_conn[] = {quad_conn[0], quad_conn[1], quad_conn[3]};
  MBEntityHandle tri1_conn[] = {quad_conn[1], quad_conn[2], quad_conn[3]};
  result = MBI->create_element(MBTRI, tri0_conn, 3, tri0);
    assert( MB_SUCCESS == result);
  result = MBI->create_element(MBTRI, tri1_conn, 3, tri1);
    assert( MB_SUCCESS == result);

  return MB_SUCCESS;
}

MBErrorCode make_tris_from_quads( MBInterface *MBI,
                                  const MBRange quads,
                                  MBRange &tris ) {
  MBErrorCode result;
  tris.clear();
  for(MBRange::const_iterator i=quads.begin(); i!=quads.end(); ++i) {
    MBEntityHandle tri0, tri1;
    result = make_tris_from_quad( MBI, *i, tri0, tri1 );
    assert(MB_SUCCESS == result);
    tris.insert( tri0 );
    tris.insert( tri1 );
  }
  return MB_SUCCESS;
}  

MBErrorCode quads_to_tris( MBInterface *MBI, MBEntityHandle input_meshset ) {

  // create a geometry tag to find the surfaces with
  MBErrorCode result;
  MBTag geom_tag, id_tag;
  result = MBI->tag_create( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                            MB_TYPE_INTEGER, geom_tag, 0, true );
    assert( MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);

  // create an id tag to find the surface id with
  result = MBI->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                            MB_TYPE_INTEGER, id_tag, 0, true );
    assert( MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);

  // get all surface meshsets
  MBRange surface_meshsets;
  int dim = 2;
  void* input_dim[] = {&dim};
  result = MBI->get_entities_by_type_and_tag( input_meshset, MBENTITYSET, &geom_tag,
                                              input_dim, 1, surface_meshsets);
    assert( MB_SUCCESS == result );
  std::cout << surface_meshsets.size() << " surfaces found." << std::endl;

  // ******************************************************************
  // Loop over every surface meshset and grab each surface's quads.
  // ******************************************************************
  for( MBRange::iterator i=surface_meshsets.begin(); i!=surface_meshsets.end(); i++ ) {

    // get the surface id of the surface meshset
    int surf_id=0;
    result = MBI->tag_get_data( id_tag, &(*i), 1, &surf_id );
      assert(MB_SUCCESS == result);
    std::cout << "  Surface " << surf_id << " has ";

    // get all quads of the surface
    MBRange quads;
    result = MBI->get_entities_by_type( *i, MBQUAD, quads );
      assert( MB_SUCCESS == result );
    std::cout << quads.size() << " quads." << std::endl;

    // ******************************************************************
    // For each quad, make two triangles then delete the quad.
    // ******************************************************************
    for(MBRange::iterator j=quads.begin(); j!=quads.end(); j++ ) {
    
      // make the tris
      MBEntityHandle tri0 = 0, tri1 = 0;
      result = make_tris_from_quad( MBI, *j, tri0, tri1 );
        assert( MB_SUCCESS == result );
      
      // add all the tris to the same surface meshset as the quads were inside.
      result = MBI->add_entities( *i, &tri0, 1 );
        if ( MB_SUCCESS != result ) std::cout << "result=" << result << std::endl;
	assert( MB_SUCCESS == result);
      result = MBI->add_entities( *i, &tri1, 1 );
        assert( MB_SUCCESS == result);

      // remove the quad from the surface meshset
      result = MBI->remove_entities( *i, &(*j), 1 );
        assert( MB_SUCCESS == result);

      // delete the quad
      result = MBI->delete_entities( &(*j), 1 );
        assert( MB_SUCCESS == result);

    } // end quad loop
  }   // end surface meshset loop
  return MB_SUCCESS; 
}
