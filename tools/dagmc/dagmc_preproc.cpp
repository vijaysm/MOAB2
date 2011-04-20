#include "ProgOptions.hpp"

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include "MBTagConventions.hpp"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <cstdlib>
#include <algorithm>

using namespace moab;

static bool verbose = false;

/**
 * Estimate the volume of the surface (actually multiplied by a factor of six).
 * See DagMC::measure_volume, from which this code is borrowed, for algorithmic details.
 * For our purpose, all that matters is the signedness.
 *
 * @param offset Offset to apply to surface to avoid a zero result.
 */
static ErrorCode get_signed_volume( Interface *MBI, 
                                    const EntityHandle surf_set,
                                    const CartVect offset, 
                                    double &signed_volume) 
{
  ErrorCode rval;
  Range tris;
  rval = MBI->get_entities_by_type( surf_set, MBTRI, tris ); 
  if(MB_SUCCESS != rval) return rval;

  signed_volume = 0.0;                                                                 
  const EntityHandle *conn;                                                            
  int len;                                                                               
  CartVect coords[3];                                                                  
  for (Range::iterator j = tris.begin(); j != tris.end(); ++j) {             
    rval = MBI->get_connectivity( *j, conn, len, true );                                 
    if(MB_SUCCESS != rval) return rval;                                                 
    if(3 != len) return MB_INVALID_SIZE;                                                                    
    rval = MBI->get_coords( conn, 3, coords[0].array() );                                
    if(MB_SUCCESS != rval) return rval;                                                 

    // Apply offset to avoid calculating 0 for cases when the origin is in the
    // plane of the surface.
    for(unsigned int k=0; k<3; ++k) {
      coords[k] += offset;        
    }
                                                                             
    coords[1] -= coords[0];                                                              
    coords[2] -= coords[0];                                                              
    signed_volume += (coords[0] % (coords[1] * coords[2]));                                   
  }                                                                                      
  return MB_SUCCESS;
}

/**
 * Replace the triangles in an old surface with those in a new surface, ensuring
 * that their surfaces senses match before the replacement
 */
static ErrorCode replace_surface( Interface *mbi, EntityHandle old_surf, EntityHandle old_file_set, 
                                  EntityHandle new_surf, const Tag& senseTag )
{

  ErrorCode rval;

  // Get the signed volume for each surface representation. If a volume comes 
  // back as zero, it's probably because a planar surface passes through the
  // origin.  In that case, try applying random offsets until reasonable 
  // values are returned.

  CartVect offset; // start with no offset
  const double min_vol = 0.1; // try again if abs(vol) < this value

  double old_vol = 0, new_vol = 0;

  bool success = false;
  int num_attempts = 100;

  while( num_attempts-- > 0 ){
    
    rval = get_signed_volume( mbi, old_surf, offset, old_vol );
    if( MB_SUCCESS != rval ) return rval;
    
    rval = get_signed_volume( mbi, new_surf, offset, new_vol );
    if( MB_SUCCESS != rval ) return rval;

    if( std::fabs(old_vol) >= min_vol && std::fabs(new_vol) >= min_vol ){
      success = true; break;
    }

    // haven't succeeded yet: randomize the offset vector 
    const int max_random = 10;
    for( int i = 0; i < 3; ++i ){ offset[i] = std::rand() % max_random; }
    
  } 

  if( !success ){
    std::cerr << "Error: could not calculate a surface volume" << std::endl;
    return MB_FAILURE;
  }
  
  // If the sign is different, reverse the old surf senses so that both
  // representations have the same signed volume.
  if( (old_vol<0 && new_vol>0) || (old_vol>0 && new_vol<0) ) {

    EntityHandle old_surf_volumes[2];
    rval = mbi->tag_get_data( senseTag, &old_surf, 1, old_surf_volumes );
    if(MB_SUCCESS != rval) return rval;  
    
    std::swap( old_surf_volumes[0], old_surf_volumes[1] );
    
    rval = mbi->tag_set_data( senseTag, &old_surf, 1, old_surf_volumes );
    if(MB_SUCCESS != rval) return rval;
  }
  
  int num_old_tris, num_new_tris;

  // Remove the tris from the old surf. Also remove them from the
  // old_file_set because it is not TRACKING.
  Range old_tris;
  rval = mbi->get_entities_by_type( old_surf, MBTRI, old_tris );
  if(MB_SUCCESS != rval) return rval;                              
  num_old_tris = old_tris.size();
  rval = mbi->remove_entities( old_surf, old_tris );                     
  if(MB_SUCCESS != rval) return rval;                              
  rval = mbi->remove_entities( old_file_set, old_tris );                     
  if(MB_SUCCESS != rval) return rval;                              
  rval = mbi->delete_entities( old_tris );
  if(MB_SUCCESS != rval) return rval;                                      

  // Add the new_surf's triangles to the old_surf  
  Range new_tris;
  rval = mbi->get_entities_by_type( new_surf, MBTRI, new_tris );
  if( MB_SUCCESS != rval) return rval;
  num_new_tris = new_tris.size();
  rval = mbi->add_entities( old_surf, new_tris );
  if(MB_SUCCESS != rval) return rval;     

  if( verbose ){ std::cout << num_old_tris << " tris -> " << num_new_tris << " tris"  << std::endl; }

  return MB_SUCCESS;

}

/**
 * Given an "old" file and a "new" file, replace the facets in any surface of the old
 * file with facets from the new file.  
 */ 
static ErrorCode merge_input_surfs( Interface *mbi, 
                                    const EntityHandle old_file_set, const EntityHandle new_file_set,
                                    const Tag& idTag, const Tag& dimTag, const Tag& senseTag )
{
  ErrorCode rval;

  const int two = 2;
  const Tag tags[2] = {dimTag, idTag};
  const void* tag_vals[2] = {&two, NULL};  

  Range old_surfs;
  rval = mbi->get_entities_by_type_and_tag( old_file_set, MBENTITYSET, &dimTag, 
                                            tag_vals, 1, old_surfs );
  if( MB_SUCCESS != rval ) return rval;

  int count = 0; 

  for( Range::iterator i = old_surfs.begin(); i!=old_surfs.end(); ++i ){
    EntityHandle old_surf = *i;

    int surf_id;
    rval = mbi->tag_get_data( idTag, &old_surf, 1, &surf_id );
    if( MB_SUCCESS != rval ) return rval;

    Range new_surf_range;
    tag_vals[1] = &surf_id;
    rval = mbi->get_entities_by_type_and_tag( new_file_set, MBENTITYSET, 
                                              tags, tag_vals, 2, new_surf_range );
    if( MB_SUCCESS != rval ) return rval;

    if( new_surf_range.size() != 1 ){
      if( new_surf_range.size() > 1 ){
        std::cerr << "Warning: surface " << surf_id 
                  << " has more than one representation in new file" << std::endl;
      }
      continue;
    }
    
    // Now we have found a surf in new_file_set to replace an old surf
    EntityHandle new_surf = new_surf_range[0];
    

    // TODO: check for quads and convert to triangles

    if( verbose ){ std::cout << "Surface " << surf_id << ": " << std::flush; } 
    rval = replace_surface( mbi, old_surf, old_file_set, new_surf, senseTag );
    if( MB_SUCCESS != rval ) return rval;
    count++;
  }

  std::cout << "Replaced " << count << " surface" << (count==1?".":"s.") << std::endl;

  return MB_SUCCESS;
}

/**
 * Kill the program informatively if code != MB_SUCCESS 
 */
static void CHECKERR( Interface& mbi, ErrorCode code ){
  if( code != MB_SUCCESS ){
    std::cerr << "Error: " << mbi.get_error_string(code) << " (" << code << ")" << std::endl;
    std::string message;
      if (MB_SUCCESS == mbi.get_last_error(message) && !message.empty())
        std::cerr << "Error message: " << message << std::endl;
    std::exit( EXIT_FAILURE );
  }
}


int main( int argc, char* argv[] ){

  ProgOptions po("dagmc_preproc: a tool for preprocessing CAD and mesh files for DAGMC analysis");

  std::string input_file;
  std::string output_file = "dagmc_preproc_out.h5m";

  po.addOpt<void>( ",v", "Verbose output", &verbose );
  po.addOpt<std::string>( ",o", "Specify output file name (default "+output_file+")", &output_file );
  po.addOpt<std::string>( ",m", "Specify alternate input mesh to override surfaces in input_file" );
  po.addOptionHelpHeading("Options for loading CAD files");
  po.addOpt<double>( "ftol,f", "Faceting distance tolerance", po.add_cancel_opt );
  po.addOpt<double>( "ltol,l", "Faceting edge length tolerance", po.add_cancel_opt );
  po.addOpt<int>( "atol,a", "Faceting normal angle tolerance (degrees)", po.add_cancel_opt );
  po.addOpt<void>( "all-curve-warnings", "Verbose warnings about curve tolerances" );
  po.addOpt<void>( "no-attribs", "Do not actuate CGM attributes" );

  po.addRequiredArg<std::string>( "input_file", "Path to input file for preprocessing", &input_file );

  po.parseCommandLine( argc, argv );

  /* Load input file, with CAD processing options, if specified */
  std::string options;
#define OPTION_APPEND(X) { if( options.length() ) options += ";"; options += (X); }

  if( po.numOptSet("no-attribs") ){
    OPTION_APPEND( "CGM_ATTRIBS=no" );
  }

  if( po.numOptSet("all-curve-warnings" ) ){
    OPTION_APPEND( "VERBOSE_CGM_CURVE_WARNINGS" );
  }

  // This is more roundabout than necessary, but we don't want *any* of the CGM-specific options
  //   to appear in the option string unless they were explicitly requested 
  double tol;
  static const int tol_prec = 12;
  if( po.getOpt( "ftol", &tol ) ){
    std::stringstream s;
    s << "FACET_DISTANCE_TOLERANCE=" << std::setprecision(tol_prec) << tol;
    OPTION_APPEND( s.str() );
  }

  if( po.getOpt( "ltol", &tol ) ){
    std::stringstream s;
    s << "MAX_FACET_EDGE_LENGTH=" << std::setprecision(tol_prec) << tol;
    OPTION_APPEND( s.str() );
  }

  int atol;
  if( po.getOpt( "atol", &atol ) ){
    std::stringstream s;
    s << "FACET_NORMAL_TOLERANCE=" << atol;
    OPTION_APPEND( s.str() );
  }


#undef OPTION_APPEND

  /* Load main input file */
  if( verbose ){
    std::cout << "Loading file " << input_file << std::endl;
    if( options.length() ) std::cout << "Option string: " << options << std::endl;
  }

  EntityHandle input_file_set;
  ErrorCode ret;
  Core mbi;

  ret = mbi.create_meshset( 0, input_file_set );
  CHECKERR( mbi, ret );

  ret = mbi.load_file( input_file.c_str(), &input_file_set, options.c_str() );
  if( ret == MB_UNHANDLED_OPTION ){
    std::cerr << "Warning: unhandled option while loading input_file, will proceed anyway" << std::endl;
  }
  else{ 
    CHECKERR( mbi, ret );
  }

  /* Iterate through any -m alternate mesh files and replace surfaces */

  std::vector<std::string> m_list;
  po.getOptAllArgs( ",m", m_list );

  if( m_list.size() > 0 ){
    // Create tags
    Tag dimTag, idTag, senseTag;
    ret = mbi.tag_create(GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
			 MB_TYPE_INTEGER, dimTag, NULL, true );
    CHECKERR( mbi, ret );

    ret = mbi.tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
			 MB_TYPE_INTEGER, idTag, NULL, true );
    CHECKERR( mbi, ret );

    ret = mbi.tag_create("GEOM_SENSE_2", 2*sizeof(EntityHandle), MB_TAG_SPARSE, 
			 MB_TYPE_HANDLE, senseTag, NULL, true );
    CHECKERR( mbi, ret );
  
    for( std::vector<std::string>::iterator i = m_list.begin(); i!=m_list.end(); ++i){
      std::cout << "Loading alternate mesh file " << *i << std::endl;
      
      EntityHandle alt_file_set;
      ret = mbi.create_meshset( 0, alt_file_set );
      CHECKERR( mbi, ret );
      
      ret = mbi.load_file( (*i).c_str(), &alt_file_set );
      CHECKERR( mbi, ret );
      
      if( verbose ) std::cout << "Merging input surfaces..." << std::flush;

      ret = merge_input_surfs( &mbi, input_file_set, alt_file_set, idTag, dimTag, senseTag );
      CHECKERR( mbi, ret );
      
      if( verbose ) std::cout << "done." << std::endl;

    }
  }

  /* Write output file */

  if( verbose ){ std::cout << "Writing " << output_file << std::endl; } 
  ret= mbi.write_file( output_file.c_str(), NULL, NULL, &input_file_set, 1 );
  CHECKERR( mbi, ret );

  return 0; 

}
