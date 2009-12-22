#include "GeometryQueryTool.hpp"
#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "MBCore.hpp"
#include "MBCartVect.hpp"
#include "cubfile.h"
#include "Tqdcfr.hpp"
#include "FileOptions.hpp"
#include "ReadNCDF.hpp"
#include "MBSkinner.hpp"
#include "quads_to_tris.hpp"
#include <limits>

#define GF_CUBIT_FILE_TYPE    "CUBIT"
#define GF_STEP_FILE_TYPE     "STEP"
#define GF_IGES_FILE_TYPE     "IGES"
#define GF_ACIS_TXT_FILE_TYPE "ACIS_SAT"
#define GF_ACIS_BIN_FILE_TYPE "ACIS_SAB"
#define GF_OCC_BREP_FILE_TYPE "OCC"

// Given parent volume senses, an id, and a set handle, this function creates a
// new surface set with dimension, geometry category, id, and sense tags.
MBErrorCode build_new_surface( MBInterface *MBI,
                               MBEntityHandle &new_surf,
                               const MBEntityHandle forward_parent_vol,
                               const MBEntityHandle reverse_parent_vol,
                               const int new_surf_id,
                               const MBTag dimTag, const MBTag idTag, 
                               const MBTag categoryTag, const MBTag senseTag ) {
 
  MBErrorCode result;
  result = MBI->create_meshset( 0, new_surf );
  assert(MB_SUCCESS == result);
  if(0 != forward_parent_vol) {
    result = MBI->add_parent_child( forward_parent_vol, new_surf );
    assert(MB_SUCCESS == result);
  }
  if(0 != reverse_parent_vol) {
    result = MBI->add_parent_child( reverse_parent_vol, new_surf );
    assert(MB_SUCCESS == result);
  }
  const int two = 2;
  result = MBI->tag_set_data( dimTag, &new_surf, 1, &two );
  assert(MB_SUCCESS == result);
  result = MBI->tag_set_data( idTag, &new_surf, 1, &new_surf_id );
  assert(MB_SUCCESS == result);  
  const char geom_category[CATEGORY_TAG_SIZE] = {"Surface\0"};
  result = MBI->tag_set_data( categoryTag, &new_surf, 1, &geom_category );
  assert(MB_SUCCESS == result);
  MBEntityHandle vols[2] = { forward_parent_vol, reverse_parent_vol };
  result = MBI->tag_set_data( senseTag, &new_surf, 1, vols );
  assert(MB_SUCCESS == result);
  return MB_SUCCESS;
}


// Given a face, orient it outward wrt its adjacent mesh element.
// Each face must be adjacent to exactly one mesh element.
MBErrorCode orient_faces_outward( MBInterface *MBI, const MBRange faces,
                                  const bool debug ) {

  MBErrorCode result;
  for(MBRange::const_iterator i=faces.begin(); i!=faces.end(); ++i) {
    MBRange adj_elem;
    result = MBI->get_adjacencies( &(*i), 1, 3, false, adj_elem );
    assert(MB_SUCCESS == result);
    assert(1 == adj_elem.size());
        
    // get center of element
    const MBEntityHandle *elem_conn;
    int n_nodes;   
    result = MBI->get_connectivity( adj_elem.front(), elem_conn, n_nodes );
    assert(MB_SUCCESS == result);
    MBCartVect elem_coords[n_nodes];
    result = MBI->get_coords( elem_conn, n_nodes, elem_coords[0].array() );
    assert(MB_SUCCESS == result);
    MBCartVect elem_center(0.0);
    for(int j=0; j<n_nodes; ++j) elem_center += elem_coords[j];
    elem_center /= n_nodes;

    // get the center of the face
    const MBEntityHandle *face_conn;
    result = MBI->get_connectivity( *i, face_conn, n_nodes );
    assert(MB_SUCCESS == result);
    MBCartVect face_coords[n_nodes];
    result = MBI->get_coords( face_conn, n_nodes, face_coords[0].array() );
    assert(MB_SUCCESS == result);
    MBCartVect face_center(0.0);
    for(int j=0; j<n_nodes; ++j) face_center += face_coords[j];
    face_center /= n_nodes;
    if(debug) std::cout << "      center of mesh face exposed by dead element=" 
                        << face_center << std::endl;

    // get the normal of the face
    MBCartVect face_normal, a, b;
    a = face_coords[2] - face_coords[1];
    b = face_coords[0] - face_coords[1];
    face_normal = a*b;
    face_normal.normalize();

    // get the direction into the element
    MBCartVect elem_dir = elem_center - face_center;
    elem_dir.normalize();

    // the dot product of the face_normal and elem_dir should be ~-1 if the face
    // is oriented outward wrt the elem.
    double dot_prod = face_normal % elem_dir;
    
    // If the face is not oriented outward wrt the element, reverse it
    if(0 < dot_prod) {
      MBEntityHandle new_face_conn[4] = { face_conn[3], face_conn[2], 
                                          face_conn[1], face_conn[0] };
      result = MBI->set_connectivity( *i, new_face_conn, 4 );
      assert(MB_SUCCESS == result);
    }
  }
  return MB_SUCCESS;
}

/* qsort int comparison function */
int handle_compare(const void *a, const void *b)
{
  const MBEntityHandle *ia = (const MBEntityHandle *)a; // casting pointer types
  const MBEntityHandle *ib = (const MBEntityHandle *)b;
  return *ia  - *ib; 
  /* integer comparison: returns negative if b > a 
     and positive if a > b */
}

// qsort face comparison function. assume each face has 4 nodes
int compare_face(const void *a, const void *b) {                                       
  MBEntityHandle *ia = (MBEntityHandle *)a;                                                  
  MBEntityHandle *ib = (MBEntityHandle *)b;                                                  
  if(*ia == *ib) {
    if(*(ia+1) == *(ib+1)) {
      if(*(ia+2) == *(ib+2)) {
        return (int)(*(ia+3) - *(ib+3));                                         
      } else {
        return (int)(*(ia+2) - *(ib+2));                                         
      }
    } else {
      return (int)(*(ia+1) - *(ib+1));                                         
    }
  } else {                                                                             
    return (int)(*ia - *ib);                                         
  }                                                                                    
}  

// Use this to get quad faces from hex elems. MBSkinner does not produce the 
// expected result.
MBErrorCode skin_hex_elems(MBInterface *MBI, MBRange elems, const int dim, 
                           MBRange &faces ) {
  // get faces of each hex
  const int nodes_per_face = 4;
  const unsigned int faces_per_elem = 6;
  unsigned int n_faces = faces_per_elem*elems.size();
  MBEntityHandle f[n_faces][nodes_per_face];
  MBErrorCode result;
  int counter = 0;
  for(MBRange::iterator i=elems.begin(); i!=elems.end(); ++i) {
    MBRange elem_faces;
    result = MBI->get_adjacencies( &(*i), 1, 2, true, elem_faces );
    assert(MB_SUCCESS == result);
    assert(faces_per_elem == elem_faces.size());
    for(MBRange::iterator j=elem_faces.begin(); j!=elem_faces.end(); ++j) {
      const MBEntityHandle *conn;
      int n_nodes;
      MBErrorCode result = MBI->get_connectivity( *j, conn, n_nodes );
      assert(MB_SUCCESS == result);
      assert(nodes_per_face == n_nodes);
      // Sort the node handles of the face
      for(int k=0; k<nodes_per_face; ++k) f[counter][k] = conn[k];
      qsort( &f[counter][0], nodes_per_face, sizeof(MBEntityHandle), handle_compare );
      ++counter;
    }
  }
  
  // Sort the faces by the first node handle, then second node, then third node...
  qsort( &f[0][0], n_faces, nodes_per_face*sizeof(MBEntityHandle), compare_face );

  // if a face has 1 or more duplicates, it is not on the skin
  faces.clear();
  for(unsigned int i=0; i<n_faces; ++i) {
    // if the last face is tested, it must be on the skin
    if(n_faces-1 == i) {
      MBRange face_handle;
      result = MBI->get_adjacencies( &(f[i][0]), nodes_per_face, 2, false, face_handle );
      assert(MB_SUCCESS == result);
      assert(1 == face_handle.size());
      faces.insert( face_handle.front() );
      // Due to sort, if a duplicate exists it must be next
    } else if( f[i][0]==f[i+1][0] && f[i][1]==f[i+1][1] && 
               f[i][2]==f[i+1][2] && f[i][3]==f[i+1][3] ) {
      ++i;
      while( f[i][0]==f[i+1][0] && f[i][1]==f[i+1][1] && 
             f[i][2]==f[i+1][2] && f[i][3]==f[i+1][3] ) {
	std::cout << "    skin WARNING: non-manifold face" << std::endl;
        ++i;
      }
      // otherwise it is on the skin
    } else {
      MBRange face_handle;
      result = MBI->get_adjacencies( &(f[i][0]), nodes_per_face, 2, false, face_handle );
      assert(MB_SUCCESS == result);
      assert(1 == face_handle.size());
      faces.insert( face_handle.front() );
    }
  }

  return MB_SUCCESS;
}

// Given a 1D array of data, axis labels, title, and number of bins, create a
// histogram.
void plot_histogram( const std::string title, 
                     const std::string x_axis_label, const std::string y_axis_label,
                     const int n_bins, const double data[], const int n_data ) {
  // find max and min
  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::max();
  for(int i=0; i<n_data; ++i) {
    if(min > data[i]) min = data[i];
    if(max < data[i]) max = data[i];

  }

  // make bins for histogram
  double bin_width = (max-min)/n_bins;
  int bins[n_bins];
  for(int i=0; i<n_bins; ++i) bins[i] = 0;
  
  // fill the bins
  for(int i=0; i<n_data; ++i) {
    double diff = data[i] - min;
    int bin = diff/bin_width;
    if(9<bin) bin = 9; // cheap fix for numerical precision error
    if(0>bin) bin = 0; // cheap fix for numerical precision error
    ++bins[bin];
  }
 
  // create bars
  int max_bin = 0;
  for(int i=0; i<n_bins; ++i) if(max_bin < bins[i]) max_bin = bins[i];
  int bar_height;
  int max_bar_chars = 72;
  std::string bars[n_bins];
  for(int i=0; i<n_bins; ++i) {
    bar_height = (max_bar_chars*bins[i])/max_bin;
    for(int j=0; j<bar_height; ++j) bars[i] += "*";
  }

  // print histogram header
  std::cout << std::endl; 
  std::cout << "                                 " << title << std::endl;

  // print results
  std::cout.width(15);
  std::cout << min << std::endl;
  for(int i=0; i<n_bins; ++i) {
    std::cout.width(15);
    std::cout << min+((i+1)*bin_width);
    std::cout.width(max_bar_chars);
    std::cout << bars[i] << bins[i] << std::endl;
  }

  // print histogram footer
  std::cout.width(15);
  std::cout << y_axis_label;
  std::cout.width(max_bar_chars);
  std::cout << " " << x_axis_label << std::endl;
  std::cout << std::endl; 
}

// This is a helper function that creates data and labels for histograms.
void generate_plots( const double orig[], const double defo[], 
                     const int n_elems, const std::string time_step ) {

  // find volume ratio then max and min
  double ratio[n_elems];
  for(int i=0; i<n_elems; ++i) ratio[i] = (defo[i]-orig[i])/orig[i];

  plot_histogram( "Predeformed Element Volume", "Num_Elems", "Volume", 10, orig, n_elems );
  std::string title = "Element Volume Change Ratio at Time Step " + time_step;
  plot_histogram( title, "Num_Elems", "Volume Ratio", 10, ratio, n_elems );

}

// Given four nodes, calculate the tet volume.
inline static double tet_volume( const MBCartVect& v0,
                                 const MBCartVect& v1,
                                 const MBCartVect& v2, 
                                 const MBCartVect& v3 )
{
  return 1./6. * ( ((v1 - v0) * (v2 - v0)) % (v3 - v0) );
}

// Measure and tet volume are taken from measure.cpp
double measure(MBInterface *MBI, const MBEntityHandle element) {
  MBEntityType type = MBI->type_from_handle( element );

  const MBEntityHandle *conn;
  int num_vertices;
  MBErrorCode result = MBI->get_connectivity( element, conn, num_vertices );
  assert(MB_SUCCESS == result);

  MBCartVect coords[num_vertices];
  result = MBI->get_coords( conn, num_vertices, coords[0].array() );
  assert(MB_SUCCESS == result);

  switch( type )
    {
    case MBEDGE:
      return (coords[0] - coords[1]).length();
    case MBTRI:
      return 0.5 * ((coords[1] - coords[0]) * (coords[2] - coords[0])).length();
    case MBQUAD:
      num_vertices = 4;
    case MBPOLYGON:
      {
	MBCartVect mid(0,0,0);
	for (int i = 0; i < num_vertices; ++i)
	  mid += coords[i];
	mid /= num_vertices;
      
	double sum = 0.0;
	for (int i = 0; i < num_vertices; ++i)
	  {
	    int j = (i+1)%num_vertices;
	    sum += ((mid - coords[i]) * (mid - coords[j])).length();
	  }
	return 0.5 * sum;
      }
    case MBTET:
      return tet_volume( coords[0], coords[1], coords[2], coords[3] ) ;
    case MBPYRAMID:
      return tet_volume( coords[0], coords[1], coords[2], coords[4] ) +
	tet_volume( coords[0], coords[2], coords[3], coords[4] ) ;
    case MBPRISM:
      return tet_volume( coords[0], coords[1], coords[2], coords[5] ) +
	tet_volume( coords[3], coords[5], coords[4], coords[0] ) +
	tet_volume( coords[1], coords[4], coords[5], coords[0] ) ;
    case MBHEX:
      return tet_volume( coords[0], coords[1], coords[3], coords[4] ) +
	tet_volume( coords[7], coords[3], coords[6], coords[4] ) +
	tet_volume( coords[4], coords[5], coords[1], coords[6] ) +
	tet_volume( coords[1], coords[6], coords[3], coords[4] ) +
	tet_volume( coords[2], coords[6], coords[3], coords[1] ) ;
    default:
      return 0.0;
    }
}

/* Calculate the signed volumes beneath the surface (x 6.0). Use the triangle's
   cannonical sense. Do not take sense tags into account. Code taken from
   DagMC::measure_volume. */
MBErrorCode get_signed_volume( MBInterface *MBI, 
                               const MBEntityHandle surf_set, 
                               double &signed_volume) {
  MBErrorCode rval;
  MBRange tris;
  rval = MBI->get_entities_by_type( surf_set, MBTRI, tris ); 
  if(MB_SUCCESS != rval) return rval;
  signed_volume = 0.0;                                                                 
  const MBEntityHandle *conn;                                                            
  int len;                                                                               
  MBCartVect coords[3];                                                                  
  for (MBRange::iterator j = tris.begin(); j != tris.end(); ++j) {             
    rval = MBI->get_connectivity( *j, conn, len, true );                                 
    if (MB_SUCCESS != rval) return rval;                                                 
    assert(3 == len);                                                                    
    rval = MBI->get_coords( conn, 3, coords[0].array() );                                
    if (MB_SUCCESS != rval) return rval;                                                 
                                                                                           
    coords[1] -= coords[0];                                                              
    coords[2] -= coords[0];                                                              
    signed_volume += (coords[0] % (coords[1] * coords[2]));                                   
  }                                                                                      
  return MB_SUCCESS;
} 

// The cgm and cub surfaces may not have the same sense. Create triangles that
// represent the quads in the cub surface. Calculate the signed volume of both
// the cgm and cub surface. If they are different, change the cgm sense so that
// it matches the sense of the cub surface.
MBErrorCode fix_surface_senses( MBInterface *MBI, 
                                const MBEntityHandle cgm_file_set, 
                                const MBEntityHandle cub_file_set,
                                const MBTag idTag, const MBTag dimTag, const MBTag senseTag,
                                const bool debug ) {
  MBErrorCode result;
  const int two = 2;
  const void* const two_val[] = {&two};
  MBRange cgm_surfs;
  result = MBI->get_entities_by_type_and_tag(cgm_file_set, MBENTITYSET, &dimTag,
						 two_val, 1, cgm_surfs );
  if(MB_SUCCESS != result) return result;    
  for(MBRange::iterator i=cgm_surfs.begin(); i!=cgm_surfs.end(); i++) {
    int surf_id;
    result = MBI->tag_get_data(idTag, &(*i), 1, &surf_id);
    if(MB_SUCCESS != result) return result;    
            
    // Find the meshed surface set with the same id
    MBRange cub_surf;
    const MBTag tags[] = {idTag, dimTag};
    const void* const tag_vals[] = { &surf_id, &two };
    result = MBI->get_entities_by_type_and_tag(cub_file_set, MBENTITYSET, tags,
						   tag_vals, 2, cub_surf );
    assert(MB_SUCCESS == result);
    if(1 != cub_surf.size()) {
      std::cout << "  surf_id=" << surf_id << " no meshed surface found" 
                << std::endl;
      continue;
    }

    // Get tris that represent the quads of the cub surf
    MBRange quads;
    result = MBI->get_entities_by_type( cub_surf.front(), MBQUAD, quads ); 
    assert(MB_SUCCESS == result);  
    MBRange cub_tris; 
    result = make_tris_from_quads( MBI, quads, cub_tris );

    // Add the tris to the same surface meshset as the quads are inside.            
    result = MBI->add_entities( cub_surf.front(), cub_tris );                       
    assert( MB_SUCCESS == result);                                                     
    
    // get the signed volume for each surface representation
    double cgm_signed_vol, cub_signed_vol;
    result = get_signed_volume( MBI,               *i, cgm_signed_vol ); 
    assert( MB_SUCCESS == result);                                      
    result = get_signed_volume( MBI, cub_surf.front(), cub_signed_vol ); 
    assert( MB_SUCCESS == result);                                      
    if(debug) std::cout << "  surf_id=" << surf_id << " cgm_signed_vol=" 
      << cgm_signed_vol << " cub_signed_vol=" << cub_signed_vol << std::endl;

    // If the sign is different, reverse the cgm senses so that both
    // representations have the same signed volume.
    if( (cgm_signed_vol<0 && cub_signed_vol>0) ||
	(cgm_signed_vol>0 && cub_signed_vol<0) ) {
      MBEntityHandle cgm_surf_volumes[2], reversed_cgm_surf_volumes[2];
      result = MBI->tag_get_data( senseTag, &(*i), 1, cgm_surf_volumes );
      assert(MB_SUCCESS == result);  
      if(MB_SUCCESS != result) return result;    

      reversed_cgm_surf_volumes[0] = cgm_surf_volumes[1];
      reversed_cgm_surf_volumes[1] = cgm_surf_volumes[0];

      result = MBI->tag_set_data( senseTag, &(*i), 1, reversed_cgm_surf_volumes );
      assert(MB_SUCCESS == result);
    }
  }

  return MB_SUCCESS;
}

  // The quads in the cub_file_set have been updated for dead elements. For each
  // cgm_surf, if there exists a cub_surf with the same id, replace the cgm tris
  // with cub_tris (created from the quads). Note the a surface that is not 
  // meshed (in cub file) will not be effected.
MBErrorCode replace_faceted_cgm_surfs( MBInterface *MBI,
                                       const MBEntityHandle cgm_file_set, 
                                       const MBEntityHandle cub_file_set, 
                                       const MBTag idTag, const MBTag dimTag, 
                                       const bool debug ) {
  MBErrorCode result;
  const int two = 2;
  const void* const two_val[] = {&two};
  MBRange cgm_surfs;
  result = MBI->get_entities_by_type_and_tag(cgm_file_set, MBENTITYSET, &dimTag,
						 two_val, 1, cgm_surfs );
  if(MB_SUCCESS != result) return result;    

  for(MBRange::iterator i=cgm_surfs.begin(); i!=cgm_surfs.end(); ++i) {
    int surf_id;
    result = MBI->tag_get_data(idTag, &(*i), 1, &surf_id);
    if(MB_SUCCESS != result) return result;    
    if(debug) std::cout << "surf_id=" << surf_id << std::endl;
            
    // Find the meshed surface set with the same id
    MBRange cub_surf;
    const MBTag tags[] = {idTag, dimTag};
    const void* const tag_vals[] = { &surf_id, &two };
    result = MBI->get_entities_by_type_and_tag(cub_file_set, MBENTITYSET, tags,
						   tag_vals, 2, cub_surf );
    assert(MB_SUCCESS == result);
    if(1 != cub_surf.size()) {
      std::cout << " Surface " << surf_id << ": no meshed representation found" << std::endl;
      continue;
    }

    // Get tris that represent the quads of the cub surf
    MBRange quads;
    result = MBI->get_entities_by_type( cub_surf.front(), MBQUAD, quads ); 
    assert(MB_SUCCESS == result);  
    MBRange cub_tris; 
    result = make_tris_from_quads( MBI, quads, cub_tris );
	
    // Remove the tris from the cgm surf. Don't forget to remove them from the
    // cgm_file_set because it is not TRACKING.
    MBRange cgm_tris;
    result = MBI->get_entities_by_type( *i, MBTRI, cgm_tris );
    assert( MB_SUCCESS == result);                              
    result = MBI->remove_entities( *i, cgm_tris );                     
    assert( MB_SUCCESS == result);                              
    result = MBI->remove_entities( cgm_file_set, cgm_tris );                     
    assert( MB_SUCCESS == result);                              
    result = MBI->delete_entities( cgm_tris );
    assert( MB_SUCCESS == result);                                      

    // Add the cub_tris to the cgm_surf
    result = MBI->add_entities( *i, cub_tris );
    assert(MB_SUCCESS == result);                                                                                          
  }
  
  return MB_SUCCESS;
}

// Dead elements need removed from the simulation. This is done by removing them
// from their volume set and adding them to the implicit complement. New surfaces
// must be created for this. 
// IF MODIFYING THIS CODE, BE AWARE THAT DEAD ELEMENTS CAN BE ADJACENT TO MORE
// THAN ONE SURFACE, MAKING THE ASSOCIATION BETWEEN NEWLY EXPOSED AND EXISTING
// SURFACES AMBIGUOUS.
MBErrorCode add_dead_elems_to_impl_compl( MBInterface *MBI, 
                                          const MBEntityHandle cgm_file_set,
                                          const MBEntityHandle cub_file_set,
                                          const MBTag idTag, const MBTag dimTag,
                                          const MBTag categoryTag, const MBTag senseTag,
                                          const bool debug ) {

  // Get the cgm surfaces
  MBErrorCode result;
  const int two = 2;
  const void* const two_val[] = {&two};
  MBRange cgm_surfs;
  result = MBI->get_entities_by_type_and_tag(cgm_file_set, MBENTITYSET, &dimTag,
					     two_val, 1, cgm_surfs );
  if(MB_SUCCESS != result) return result;    

  // Get the maximum surface id. This is so that new surfaces to do have
  // duplicate ids.
  int max_surf_id = -1;
  for(MBRange::const_iterator i=cgm_surfs.begin(); i!=cgm_surfs.end(); ++i) {
    int surf_id;
    result = MBI->tag_get_data( idTag, &(*i), 1, &surf_id );
    assert(MB_SUCCESS == result);
    if(max_surf_id < surf_id) max_surf_id = surf_id;
  }
  std::cout << "  Maximum surface id=" << max_surf_id << std::endl;

  // For each cgm volume, does a cub volume with the same id exist?
  const int three = 3;
  const void* const three_val[] = {&three};
  MBRange cgm_vols;
  result = MBI->get_entities_by_type_and_tag(cgm_file_set, MBENTITYSET, &dimTag,
					     three_val, 1, cgm_vols );
  assert(MB_SUCCESS == result);
  // get the corresponding cub volume
  for(MBRange::iterator i=cgm_vols.begin(); i!=cgm_vols.end(); i++) {
    int vol_id;
    result = MBI->tag_get_data(idTag, &(*i), 1, &vol_id);
    if(MB_SUCCESS != result) return result;    
    std::cout << "  Volume " << vol_id << std::endl;
            
    // Find the meshed vol set with the same id
    MBRange cub_vol;
    const MBTag tags[] = {idTag, dimTag};
    const void* const tag_vals[] = { &vol_id, &three };
    result = MBI->get_entities_by_type_and_tag(cub_file_set, MBENTITYSET, tags,
					       tag_vals, 2, cub_vol );
    assert(MB_SUCCESS == result);
    if(1 != cub_vol.size()) {
      std::cout << "    no meshed volume found" << std::endl;
      continue;
    }

    // get the mesh elements of the volume.
    MBRange elems;
    result = MBI->get_entities_by_type( cub_vol.front(), MBHEX, elems ); 
    assert(MB_SUCCESS == result);  
    if(debug) std::cout << "    found " << elems.size() << " hex elems" << std::endl;

    // skin the volumes
    // BEWARE: THE MBSKINNER DOES NOT PRODUCE THE EXPECTED RESULT. CONFIRMED IN VISIT.
    MBSkinner tool(MBI);
    MBRange moab_faces, my_faces;
    result = tool.find_skin( elems, 2, moab_faces, true );
    assert(MB_SUCCESS == result);
    skin_hex_elems( MBI, elems, 2, my_faces );
    assert(MB_SUCCESS == result);
    if(debug) std::cout << "    moab_faces=" << moab_faces.size() << " my_faces=" 
	      << my_faces.size() << std::endl;
   
    // get cub child surfaces.
    MBRange cub_surfs;
    result = MBI->get_child_meshsets( cub_vol.front(), cub_surfs );
    assert(MB_SUCCESS == result);
    for(MBRange::iterator j=cub_surfs.begin(); j!=cub_surfs.end(); ++j) {
      // get the quads on each surface
      MBRange cub_faces;
      result = MBI->get_entities_by_type( *j, MBQUAD, cub_faces );
      assert(MB_SUCCESS == result);
      // Meshed volumes must have meshed surfaces
      assert(!cub_faces.empty());
      // get the faces common to both the skin and this surface
      MBRange common_faces = intersect( cub_faces, my_faces );
      // find the surfaces faces not on the skin - these are orphaned and need removed
      MBRange orphaned_faces = subtract( cub_faces, common_faces );
      result = MBI->remove_entities( *j, orphaned_faces );
      assert(MB_SUCCESS == result);
      int surf_id;
      result = MBI->tag_get_data( idTag, &(*j), 1, &surf_id );
      assert(MB_SUCCESS == result);
      // remove the common faces from the skin faces
      my_faces = subtract( my_faces, common_faces );
      // If no orphaned faces exist we are done
      if(orphaned_faces.empty()) continue;
      std::cout << "    Surface " << surf_id << " had " << orphaned_faces.size() 
		<< " orphaned faces removed" << std::endl;
      // place the orphaned faces in a new surface. Get the parent vols of the
      // surf.
      MBRange cgm_surf;
      const MBTag tags[] = {idTag, dimTag};
      const void* const tag_vals[] = { &surf_id, &two };
      result = MBI->get_entities_by_type_and_tag(cgm_file_set, MBENTITYSET, tags,
						 tag_vals, 2, cgm_surf );
      assert(MB_SUCCESS == result);
      assert(1 == cgm_surf.size());
      MBEntityHandle cgm_vols[2], cub_vols[2];
      result = MBI->tag_get_data( senseTag, cgm_surf, &cgm_vols );
      assert(MB_SUCCESS == result);
      result = MBI->tag_get_data( senseTag, &(*j), 1, &cub_vols );
      assert(MB_SUCCESS == result);
      // for the new surf, replace the current volume with the impl compl vol.
      // This is because the faces that no longer exist will become adjacent to
      // the impl compl
      if(*i == cgm_vols[0]) {
	cgm_vols[0] = 0;
	cub_vols[0] = 0;
      }
      if(*i == cgm_vols[1]) {
	cgm_vols[1] = 0;
	cub_vols[1] = 0;
      }
      // If both sides of the surface are the impl comp, do not create the surface.
      if(0==cgm_vols[0] && 0==cgm_vols[1]) {
	std::cout << "    New surface was not created for orphaned faces because both parents are impl_compl volume " << std::endl;
	continue;
      }
      // build the new surface, convert quads to tris, and add the faces.
      MBEntityHandle new_cgm_surf, new_cub_surf;
      ++max_surf_id;
      result = build_new_surface( MBI, new_cgm_surf,
				  cgm_vols[0], cgm_vols[1], max_surf_id,
				  dimTag, idTag, 
				  categoryTag, senseTag );
      assert(MB_SUCCESS == result);
      result = build_new_surface( MBI, new_cub_surf,
				  cub_vols[0], cub_vols[1], max_surf_id,
				  dimTag, idTag, 
				  categoryTag, senseTag );
      assert(MB_SUCCESS == result);
      // add the new surface to the file set and populate it with faces
      result = MBI->add_entities( cgm_file_set, &new_cgm_surf, 1 );
      assert(MB_SUCCESS == result);
      result = MBI->add_entities( cub_file_set, &new_cub_surf, 1 );
      assert(MB_SUCCESS == result);       
      result = MBI->add_entities( new_cub_surf, orphaned_faces ); 
      assert(MB_SUCCESS == result);
      std::cout << "    Surface " << max_surf_id << " was created for " 
		<< orphaned_faces.size() << " orphaned faces" << std::endl;
    }

    // remaining skin faces must be assigned to a surface
    if(my_faces.empty()) continue;
    std::cout << "    Surface " << max_surf_id+1 << " will be created for " 
	      << my_faces.size() << " dead skin faces" << std::endl;

    // Ensure that faces are oriented outwards
    result = orient_faces_outward( MBI, my_faces, debug );
    assert(MB_SUCCESS == result);

    // Create the new surface.
    MBEntityHandle new_cgm_surf, new_cub_surf;
    ++max_surf_id;
    result = build_new_surface( MBI, new_cgm_surf,
				*i, 0, max_surf_id,
				dimTag, idTag, 
				categoryTag, senseTag );
    assert(MB_SUCCESS == result);
    result = build_new_surface( MBI, new_cub_surf,
				cub_vol.front(), 0, max_surf_id,
				dimTag, idTag, 
				categoryTag, senseTag );
    assert(MB_SUCCESS == result);
    // Insert the new surf into file sets and populate it with faces.
    result = MBI->add_entities( cgm_file_set, &new_cgm_surf, 1 );
    assert(MB_SUCCESS == result);
    result = MBI->add_entities( cub_file_set, &new_cub_surf, 1 );
    assert(MB_SUCCESS == result);       
    result = MBI->add_entities( new_cub_surf, my_faces );
    assert(MB_SUCCESS == result);
  }

  return MB_SUCCESS;
}
  
 
/* Get the type of a file.
   Return value is one of the above constants
*/
const char* get_geom_file_type( const char* filename );
const char* get_geom_fptr_type( FILE* file );

int is_cubit_file( FILE* file );
int is_step_file( FILE* file );
int is_iges_file( FILE* file );
int is_acis_txt_file( FILE* file );
int is_acis_bin_file( FILE* file );
int is_occ_brep_file( FILE* file );

double DEFAULT_DISTANCE = 0.001;
double DEFAULT_LEN = 0.0;
int DEFAULT_NORM = 5;

// load cub file 
// load cgm file
// for each surface
//   convert cub surf quads to tris
//   get signed volume from cgm and cub surf         MUST COME BEFORE COORD UPDATE, NEEDS TRIS
//   reverse cgm surface sense if needed
//   replace cgm surf tris with cub surf tris
//   X remove quads from cub surf
// measure volume of predeformed cub elements
// convert cub volumes sets to tracking so that dead elems are removed from vol sets
// update coordinates and delete dead elems
// measure volume of deformed cub elems
// print histogram of volume change
// for each cub volume
//   skin volume elems to get faces
//   for each child cub surface
//     remove surface faces that are not in skin
//   orient each skin face outward
//   assign new (leftover) skin faces to a surface
// for each surface
//   remove existing tris (from before the update)
//   convert quads to tris
int main( int argc, char* argv[] )
{
  const bool debug = false;
  const char *file_type = NULL;
  
  const char* cub_name = 0;
  const char* exo_name = 0;
  const char* out_name = 0;
  const char* time_step = 0;
  const char* sat_name = 0;
  double dist_tol = 0.001, len_tol = 0.0;
  int norm_tol = 5;

  if(4!=argc && 6!=argc && 7!=argc)
    {
      std::cerr << "To read meshed geometry for DagMC:" << std::endl;
      std::cerr << "$> cub_file acis_file output_file" << std::endl;
      std::cerr << "To read meshed geometry for DagMC and update node coordinates:" << std::endl;
      std::cerr << "$> <cub_file.cub> <acis_file.sat> <output_file.h5m> <deformed_exo_file.e> time_step<int> check_vol_change<bool>" 
		<< std::endl;
      exit(4);
    }

  // check filenames for proper suffix
  std::string temp;
  cub_name = argv[1];
  temp.assign(cub_name);
  if(std::string::npos == temp.find(".cub")) {
    std::cerr << "cub_file does not have correct suffix" << std::endl;
    return 1;
  }
  sat_name = argv[2]; // needed because the cub file's embedded sat file does not have groups
  temp.assign(sat_name);
  if(std::string::npos == temp.find(".sat")) {
    std::cerr << "sat_file does not have correct suffix" << std::endl;
    return 1;
  }
  out_name = argv[3];
  temp.assign(out_name);
  if(std::string::npos == temp.find(".h5m")) {
    std::cerr << "out_file does not have correct suffix" << std::endl;
    return 1;
  }

  // Should the nodes be updated?
  bool update_coords = false;
  if(6 <= argc) {
    exo_name = argv[4];
    temp.assign(exo_name);
    if(std::string::npos == temp.find(".e")) {
      std::cerr << "e_file does not have correct suffix" << std::endl;
      return 1;
    }
    time_step= argv[5];
    update_coords = true;
  }

  // Should the volume change be determined?
  bool determine_volume_change = false;
  if(7 == argc) {
    temp.assign(argv[6]);
    if(std::string::npos != temp.find("true")) determine_volume_change = true;
  }
  
  // Get CGM file type
  if (!file_type) {
    file_type = get_geom_file_type( cub_name );
    if (!file_type) {
      std::cerr << cub_name << " : unknown file type, try '-t'" << std::endl;
      exit(1);
    }
  }

  // Read the mesh from the cub file with Tqcdfr 
  MBCore *MBI = new MBCore();
  MBErrorCode result;
  MBEntityHandle cub_file_set;
  result = MBI->create_meshset( 0, cub_file_set );
  assert(MB_SUCCESS == result);
  //char cub_options[256] = "120";
  //result = MBI->load_file(cub_name, cub_file_set, cub_options, NULL, 0, 0);
  result = MBI->load_file(cub_name, &cub_file_set, 0, NULL, 0, 0);
  if(MB_SUCCESS != result) return result;
  std::cout << "Mesh file read." << std::endl;

  // Read the ACIS file with ReadCGM
  char cgm_options[256];
  sprintf(cgm_options,
	  "CGM_ATTRIBS=yes;FACET_DISTANCE_TOLERANCE=%g;FACET_NORMAL_TOLERANCE=%d;MAX_FACET_EDGE_LENGTH=%g;",
	  dist_tol,norm_tol,len_tol); 
  MBEntityHandle cgm_file_set;
  result = MBI->create_meshset( 0, cgm_file_set );
  assert(MB_SUCCESS == result);
  result = MBI->load_file(sat_name, &cgm_file_set, cgm_options,NULL,0,0);
  if(MB_SUCCESS != result) return result;
  std::cout << "Geometry file read." << std::endl;
    
  // Create tags
  MBTag dimTag, idTag, categoryTag, senseTag;
  result = MBI->tag_create(GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
			       MB_TYPE_INTEGER, dimTag, NULL, true );
  if(MB_SUCCESS != result) return result;
  result = MBI->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
			       MB_TYPE_INTEGER, idTag, NULL, true );
  if(MB_SUCCESS != result) return result;
  result = MBI->tag_create(CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, MB_TAG_SPARSE, 
			       MB_TYPE_OPAQUE, categoryTag, NULL, true );
  if(MB_SUCCESS != result) return result;
  result = MBI->tag_create("GEOM_SENSE_2", 2*sizeof(MBEntityHandle), MB_TAG_DENSE, 
			       MB_TYPE_HANDLE, senseTag, NULL, true );
  if(MB_SUCCESS != result) return result;
    
  // Create triangles from the quads of the cub surface sets and add them to the
  // cub surface sets. Get the signed volume of each surface for both cgm and 
  // cub representations. Change the sense of the cgm representation to match 
  // the cub representation.
  result = fix_surface_senses( MBI, cgm_file_set, cub_file_set, idTag, dimTag, senseTag, debug );
  assert(MB_SUCCESS == result);
  std::cout << "Fixed geometry surface senses to match meshed surface senses." << std::endl;

  // Get the 3D elements in the cub file and measure their volume.
  MBRange orig_elems;
  std::vector<double> orig_size;
  if(determine_volume_change) {
    result = MBI->get_entities_by_dimension( 0, 3, orig_elems );   
    assert(MB_SUCCESS == result);
    orig_size.resize(orig_elems.size());
    for(unsigned int i=0; i<orig_elems.size(); ++i) {
      orig_size[i] = measure( MBI, orig_elems[i] );
    }
  }

  // Before updating the nodes and removing dead elements, force the cub volume
  // sets to track ownership so that dead elements will be deleted from the sets.
  const int three = 3;
  const void* const three_val[] = {&three};
  MBRange cub_vols;
  result = MBI->get_entities_by_type_and_tag(cub_file_set, MBENTITYSET, &dimTag,
						 three_val, 1, cub_vols );
  assert(MB_SUCCESS == result);
  for(MBRange::const_iterator i=cub_vols.begin(); i!=cub_vols.end(); ++i) {
    result = MBI->set_meshset_options( *i, MESHSET_TRACK_OWNER );
    assert(MB_SUCCESS == result);
  }

  // Update the coordinates if needed. Do not do this before checking surface
  // sense, because the coordinate update could deform the surfaces too much
  // to make an accurate comparison.
  // The cub node ids are unique because cgm vertex ids are tagged on the vertex
  // meshset and not the vertex itself.
  //result = MBI->delete_entities( &cub_file_set, 1 );
  //assert(MB_SUCCESS == result);
  // Assume dead elements exist until I think of something better.
  bool dead_elements_exist = true;
  if(update_coords) {
    ReadNCDF my_ex_reader(MBI);
    char exo_options[120] = "tdata=coord,";
    strcat(exo_options, time_step);
    strcat(exo_options,",set");
    FileOptions exo_opts(exo_options)  ;
    //opts = "tdata=coord, 100, sum, temp.exo";
    //result =  my_ex_reader.load_file(exo_name, cgm_file_set, exo_opts, NULL, 0 , 0);
    //result =  my_ex_reader.load_file(exo_name, cub_file_set, exo_opts, NULL, 0 , 0);
    result = my_ex_reader.load_file(exo_name, &cub_file_set, exo_opts, NULL, 0 , 0);
    if(MB_SUCCESS != result) {
      std::cout << "coordinate update failed" << std::endl;
      return result;
    }
    std::cout << "Updated mesh nodes with deformed coordinates from exodus file." << std::endl;
  }

  if(determine_volume_change) {
    // Dead elements have been removed by the deformation. Get the elements that 
    // still exist.
    MBRange defo_elems;
    result = MBI->get_entities_by_dimension( 0, 3, defo_elems );   
    assert(MB_SUCCESS == result);
    
    // Determine the volume of the elements now that a deformation has been
    // applied. Condense the original array by removing dead elements.
    double orig_size_condensed[defo_elems.size()];
    double defo_size_condensed[defo_elems.size()];
    int j=0;
    for(unsigned int i=0; i<orig_elems.size(); ++i) {
      if(orig_elems[i] == defo_elems[j]) {
	orig_size_condensed[j] = orig_size[i];
	defo_size_condensed[j] = measure( MBI, defo_elems[j] );
	++j;
      }
    }
    generate_plots( orig_size_condensed, defo_size_condensed, 
		    defo_elems.size(), std::string(time_step) );
  }

  // Deal with dead elements. For now, add them to the impl_compl volume.
  // Extra surfaces are created to do this.
  if(update_coords && dead_elements_exist) {
    result = add_dead_elems_to_impl_compl( MBI, cgm_file_set, cub_file_set,
                                           idTag, dimTag, categoryTag, senseTag, debug );
    assert(MB_SUCCESS == result);
    std::cout << "Placed dead elements to implicit complement volume and added required surfaces." << std::endl;
  }

  // The quads in the cub_file_set have been updated for dead elements. For each
  // cgm_surf, if there exists a cub_surf with the same id, replace the cgm tris
  // with cub_tris (created from the quads). Note the a surface that is not 
  // meshed (in cub file) will not be effected.
  result = replace_faceted_cgm_surfs( MBI, cgm_file_set, cub_file_set, idTag, dimTag, debug );
  assert(MB_SUCCESS == result);
  std::cout << "Replaced faceted geometry surfaces with meshed surfaces of triangles." << std::endl;

  result = MBI->write_mesh( out_name, &cgm_file_set, 1 );
  if(MB_SUCCESS != result) {
    std::cout << "write mesh failed" << std::endl;
    return result;
  }
  std::cout << "Saved output file for mesh-based analysis." << std::endl;
  std::cout << std::endl;
  
  return 0;
}

const char* get_geom_file_type( const char* name )
{
  FILE* file;
  const char* result = 0;
  
  file = fopen( name, "r" );
  if (file) {
    result = get_geom_fptr_type( file );
    fclose( file );
  }
  
  return result;
}

const char* get_geom_fptr_type( FILE* file )
{
  static const char* CUBIT_NAME = GF_CUBIT_FILE_TYPE;
  static const char*  STEP_NAME = GF_STEP_FILE_TYPE;
  static const char*  IGES_NAME = GF_IGES_FILE_TYPE;
  static const char*   SAT_NAME = GF_ACIS_TXT_FILE_TYPE;
  static const char*   SAB_NAME = GF_ACIS_BIN_FILE_TYPE;
  static const char*  BREP_NAME = GF_OCC_BREP_FILE_TYPE;
  
  if (is_cubit_file(file))
    return CUBIT_NAME;
  else if (is_step_file(file))
    return STEP_NAME;
  else if (is_iges_file(file))
    return IGES_NAME;
  else if (is_acis_bin_file(file))
    return SAB_NAME;
  else if (is_acis_txt_file(file))
    return SAT_NAME;
  else if (is_occ_brep_file(file))
    return BREP_NAME;
  else
    return 0;
}

int is_cubit_file( FILE* file )
{
  unsigned char buffer[4];
  return !fseek(file, 0, SEEK_SET) &&
    fread(buffer, 4, 1, file) &&
    !memcmp(buffer, "CUBE", 4);
}

int is_step_file( FILE* file )
{
  unsigned char buffer[9];
  return !fseek(file, 0, SEEK_SET) &&
    fread(buffer, 9, 1, file) &&
    !memcmp(buffer, "ISO-10303", 9);
}

int is_iges_file( FILE* file )
{
  unsigned char buffer[10];
  return !fseek(file, 72, SEEK_SET) &&
    fread(buffer, 10, 1, file) &&
    !memcmp(buffer, "S      1\r\n", 10);
}

int is_acis_bin_file( FILE* file )
{
  char buffer[15];
  return !fseek(file, 0, SEEK_SET) &&
    fread(buffer, 15, 1, file) &&
    !memcmp(buffer, "ACIS BinaryFile", 9);
}

int is_acis_txt_file( FILE* file )
{
  char buffer[5];
  int version, length;
  
  if (fseek(file,0,SEEK_SET) || 
      2 != fscanf( file, "%d %*d %*d %*d %d ", &version, &length ))
    return 0;
    
  if (version < 1 || version >0xFFFF)
    return 0;
  
  // Skip appliation name
  if (fseek(file, length, SEEK_CUR))
    return 0;
    
  // Read length of version string followed by first 5 characters
  if (2 != fscanf(file, "%d %4s", &length, buffer))
    return 0;
    
  return !strcmp( buffer, "ACIS" );
}

int is_occ_brep_file( FILE* file )
{
  unsigned char buffer[6];
  return !fseek(file, 0, SEEK_SET) &&
    fread(buffer, 6, 1, file) &&
    !memcmp(buffer, "DBRep_", 6);
}
