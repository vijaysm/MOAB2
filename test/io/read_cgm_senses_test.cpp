
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "moab/GeomTopoTool.hpp"
#include "InitCGMA.hpp"
#include "GeometryQueryTool.hpp"
#include "cgm_version.h"

#define SENSE_FORWARD 1
#define SENSE_REVERSE -1 
#define SENSE_UNKNOWN 0
using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)


#ifdef MESHDIR
#ifdef HAVE_OCC_STEP
static const char input_cylcube[] = STRINGIFY(MESHDIR) "/io/cylcube.stp";
#else
static const char input_cylcube[] = STRINGIFY(MESHDIR) "/io/cylcube.sat";
#endif
#else
#ifdef HAVE_OCC_STEP
static const char input_cylcube[] = "cylcube.stp";
#else
static const char input_cylcube[] = "cylcube.sat";
#endif
#endif


// Function used to load the test file
void read_file( Interface* moab, const char* input_file );

// Functions containing known sense data
ErrorCode load_sat_curve_sense_data( Interface* moab, EntityHandle curve,  std::vector<int>& surf_ids_out, std::vector<int>& senses_out );
ErrorCode load_stp_curve_sense_data( Interface* moab, EntityHandle curve,  std::vector<int>& surf_ids_out, std::vector<int>& senses_out );
ErrorCode load_precgm14_stp_curve_sense_data( Interface* moab, EntityHandle curve,  std::vector<int>& surf_ids_out, std::vector<int>& senses_out );
ErrorCode load_sat_surf_sense_data( Interface* moab, EntityHandle surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out );
ErrorCode load_stp_surf_sense_data( Interface* moab, EntityHandle surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out );
ErrorCode load_precgm14_stp_surf_sense_data( Interface* moab, EntityHandle surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out );

// Functions used to compare sense information found in 
// the model to reference information
void check_sense_data( Interface* moab, std::vector<EntityHandle> wrt_ents, std::vector<int> senses, 
                             std::vector<int> known_wrt_ids, std::vector<int> known_senses );

//Function used to get id's from entity handles
int geom_id_by_handle( Interface* moab, const EntityHandle set );

// List of tests in this file
void read_cylcube_curve_senses_test();
void read_cylcube_surf_senses_test();
void delete_mesh_test();


int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST(read_cylcube_surf_senses_test);  
  result += RUN_TEST(read_cylcube_curve_senses_test);  
 
  return result;
}


void read_file( Interface* moab, const char* input_file )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();

  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

// Gets the sense data for each curve from a file 
// containing a cube and a cylinder. It then checks
// that this sense data matches the reference
// sense data from prior file reads.
void read_cylcube_curve_senses_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cylcube );
  
  //Get all curve handles
  Tag geom_tag;
  rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1,
			MB_TYPE_INTEGER, geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);
  
  // Check that the proper number of curves exist
  int dim = 1;
  void *val[] = {&dim};
  int number_of_curves;
  rval = mb->get_number_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
    					          val, 1, number_of_curves );
  CHECK_ERR(rval);
  //Step format adds a surface on the barrel of the cylinder.
  //This created 4 extra surfaces in comparison to the .sat format from Cubit. 
  //(New surface breaks the barrel of the cylinder into two half-pipes)
#ifdef HAVE_OCC_STEP
  CHECK_EQUAL( 18, number_of_curves );
#else
  CHECK_EQUAL( 14, number_of_curves );
#endif

  //Get curve handles
  Range curves;
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, curves );
  CHECK_ERR(rval);

  // Establish GeomTopoTool instance needed to get curve data 
  moab::GeomTopoTool gt( mb, false );  
  // Initialize vectors for sense checking
  std::vector<EntityHandle> surfs;
  std::vector<int> senses;  
  std::vector<int> known_surf_ids;
  std::vector<int> known_senses;

for(unsigned int i = 0; i < curves.size() ; i++)
  {

   //Clean data from previous curve
   surfs.clear();
   senses.clear();
   //Get sense info for the current curve
   rval = gt.get_senses( curves[i], surfs, senses );
   CHECK_ERR(rval);

   //Clear reference data from previous curve
   known_surf_ids.clear();
   known_senses.clear();
   //Load known curve-sense ID data
#ifdef HAVE_OCC_STEP
   if (CGM_MAJOR_VERSION >= 14)
   {
     rval = load_stp_curve_sense_data( mb, curves[i], known_surf_ids, known_senses );
   }
   else
   {
     rval = load_precgm14_stp_curve_sense_data( mb, curves[i], known_surf_ids, known_senses );
   }
   CHECK_ERR(rval);
#else
   rval = load_sat_curve_sense_data( mb, curves[i], known_surf_ids, known_senses );
   CHECK_ERR(rval);
#endif

   //Check that each surf and sense has a match in the references
   check_sense_data( mb, surfs, senses, known_surf_ids, known_senses);
  }
}


int geom_id_by_handle( Interface* moab, const EntityHandle set )
{
    
    ErrorCode rval;
    //Get the id_tag handle
    Tag id_tag;
    rval = moab->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE );
    CHECK_ERR(rval);
    //Load the ID for the EntHandle given to the function                  
    int id;
    rval = moab->tag_get_data( id_tag, &set, 1, &id );                  
    CHECK_ERR(rval);                        
    return id;
 }

void check_sense_data( Interface* moab, std::vector<EntityHandle> wrt_ents, std::vector<int> senses, 
                             std::vector<int> known_wrt_ids, std::vector<int> known_senses )
{
  
  //Get ID's of the wrt entities
  std::vector<int> wrt_ent_ids;
  for(unsigned int i=0 ; i<wrt_ents.size() ; i++)
  {
      wrt_ent_ids.push_back( geom_id_by_handle( moab, wrt_ents[i] ) );
  }

  for(unsigned int i=0; i< wrt_ent_ids.size() ; i++)
  {
     for(unsigned int j=0; j< known_wrt_ids.size(); j++)
     {
       if( wrt_ent_ids[i] == known_wrt_ids [j] )
         {
          // Make sure the senses of the matching wrt entities
          // are correct
          CHECK_EQUAL( senses[i], known_senses[j] );
          //Once a wrt entity is matched with a known entity,
          // remove it from the list
          wrt_ent_ids.erase( wrt_ent_ids.begin()+i );
          senses.erase( senses.begin()+i );
          --i;
          break;
         }
     }
  }

  // After both loops are complete, known_wrt_ents should be empty 
  int leftovers = wrt_ent_ids.size();
  CHECK_EQUAL( leftovers, 0 );

}

//Loads two vectors with reference curve and curve_sense data
ErrorCode load_sat_curve_sense_data( Interface* moab, EntityHandle curve, std::vector<int>& surf_ids_out, std::vector<int>& senses_out )
{

  int curve_id = geom_id_by_handle( moab, curve );
  switch(curve_id)
  {
    case 1:
          surf_ids_out.push_back(1); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 2:
          surf_ids_out.push_back(1); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 3:
          surf_ids_out.push_back(1); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 4:
          surf_ids_out.push_back(1); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 5:
          surf_ids_out.push_back(2); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 6:
          surf_ids_out.push_back(2); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 7:
          surf_ids_out.push_back(2); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 8:
          surf_ids_out.push_back(2); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 9:
          surf_ids_out.push_back(3); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 10:
          surf_ids_out.push_back(3); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 11:
          surf_ids_out.push_back(4); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 12:
          surf_ids_out.push_back(5); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 13:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;

    case 14:
      surf_ids_out.push_back(7); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    default:
          return MB_FAILURE;

  } 
  return MB_SUCCESS;
}

//Loads two vectors with reference curve and curve_sense data
ErrorCode load_stp_curve_sense_data( Interface* moab, EntityHandle curve, std::vector<int>& surf_ids_out, std::vector<int>& senses_out )
{

  int curve_id = geom_id_by_handle( moab, curve );
  switch(curve_id)
  {
    case 1:
          surf_ids_out.push_back(1); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_FORWARD);
          break;

    case 2:
          surf_ids_out.push_back(1); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_FORWARD);
          break;

    case 3:
          surf_ids_out.push_back(1); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_FORWARD);
          break;

    case 4:
          surf_ids_out.push_back(1); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_FORWARD);
          break;

    case 5:
          surf_ids_out.push_back(2); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 6:
          surf_ids_out.push_back(2); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 7:
          surf_ids_out.push_back(2); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 8:
          surf_ids_out.push_back(2); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 9:
          surf_ids_out.push_back(3); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 10:
          surf_ids_out.push_back(3); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 11:
          surf_ids_out.push_back(4); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 12:
          surf_ids_out.push_back(5); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 13:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;

    case 14:
      surf_ids_out.push_back(7); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 15:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 16:
      surf_ids_out.push_back(7); surf_ids_out.push_back(10);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 17:
      surf_ids_out.push_back(8); surf_ids_out.push_back(10);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 18:
      surf_ids_out.push_back(8); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    default:
          return MB_FAILURE;

  } 
  return MB_SUCCESS;
}

ErrorCode load_precgm14_stp_curve_sense_data( Interface* moab, EntityHandle curve, std::vector<int>& surf_ids_out, std::vector<int>& senses_out )
{

  int curve_id = geom_id_by_handle( moab, curve );
  switch(curve_id)
  {
    case 1:
          surf_ids_out.push_back(1); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 2:
          surf_ids_out.push_back(1); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 3:
          surf_ids_out.push_back(1); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 4:
          surf_ids_out.push_back(1); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 5:
          surf_ids_out.push_back(2); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 6:
          surf_ids_out.push_back(2); surf_ids_out.push_back(3);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 7:
          surf_ids_out.push_back(2); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 8:
          surf_ids_out.push_back(2); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 9:
          surf_ids_out.push_back(3); surf_ids_out.push_back(4);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 10:
          surf_ids_out.push_back(3); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
          break;

    case 11:
          surf_ids_out.push_back(4); surf_ids_out.push_back(5);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 12:
          surf_ids_out.push_back(5); surf_ids_out.push_back(6);
          senses_out.push_back(SENSE_FORWARD); senses_out.push_back(SENSE_REVERSE);
          break;

    case 13:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;

    case 14:
      surf_ids_out.push_back(7); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 15:
      surf_ids_out.push_back(7); surf_ids_out.push_back(8);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 16:
      surf_ids_out.push_back(7); surf_ids_out.push_back(10);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 17:
      surf_ids_out.push_back(8); surf_ids_out.push_back(10);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    case 18:
      surf_ids_out.push_back(8); surf_ids_out.push_back(9);
      senses_out.push_back(SENSE_REVERSE); senses_out.push_back(SENSE_FORWARD);
      break;
    default:
          return MB_FAILURE;

  } 
  return MB_SUCCESS;
}

// SURFACE SENSE CHECKING
// Gets the sense data for each surface from a file 
// containing a cube and a cylinder. It then checks
// that this sense data matches the reference
// sense data from prior file reads.
void read_cylcube_surf_senses_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cylcube );
  
  //Get geometry tag for gathering surface information from the mesh
  Tag geom_tag;
  rval = mb->tag_get_handle( GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, 
                             geom_tag, moab::MB_TAG_DENSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);
  
  // Check that the proper number of surfaces exist
  int dim = 2;
  void *val[] = {&dim};
  int number_of_surfs;
  rval = mb->get_number_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, number_of_surfs );
  CHECK_ERR(rval);
  //Step format adds a surface on barrel of the cylinder.
  // (Breaks it into two half-pipes)
#ifdef HAVE_OCC_STEP
  CHECK_EQUAL( 10, number_of_surfs );
#else
  CHECK_EQUAL( 9, number_of_surfs );
#endif
  // Get surface handles
  Range surfs;
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, surfs );
  CHECK_ERR(rval);

  // Establish GeomTopoTool instance needed to get surf data 
  moab::GeomTopoTool gt( mb, false );  
  std::vector<EntityHandle> vols;
  std::vector<int> senses;  
  std::vector<int> known_vol_ids;
  std::vector<int> known_senses;

for(unsigned int i = 0; i < surfs.size(); i++)
  {
   //Clean data from previous surface
   vols.clear();
   senses.clear();
   // Get sense information for the current
   // surface from the mesh
   rval = gt.get_senses( surfs[i], vols, senses );
   CHECK_ERR(rval);
   //Clear previous reverence data
   known_vol_ids.clear();
   known_senses.clear();
   // Load known surface-volume data 
   // for this surface and check that it's correct
#ifdef HAVE_OCC_STEP
   if (CGM_MAJOR_VERSION >= 14)
   {
     rval = load_stp_surf_sense_data( mb, surfs[i], known_vol_ids, known_senses );
   }
   else
   {
     rval = load_precgm14_stp_surf_sense_data( mb, surfs[i], known_vol_ids, known_senses );
   }
   CHECK_ERR(rval);
#else
   rval = load_sat_surf_sense_data( mb, surfs[i], known_vol_ids, known_senses );
   CHECK_ERR(rval);
#endif
   // Check sense information from the loaded mesh against 
   // reference sense information
   check_sense_data( mb, vols, senses, known_vol_ids, known_senses );

  }

}

//Loads reference surface to volume sense data into the reference vectors
ErrorCode load_sat_surf_sense_data( Interface* moab, EntityHandle surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out ){

  int surf_id = geom_id_by_handle( moab, surf );
  switch(surf_id)
  {
    case 1:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 2:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 3:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 4:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 5:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 6:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 7:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
  
    case 8:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;

    case 9:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
    default:
          return MB_FAILURE;

   }
  return MB_SUCCESS;
}

//Loads reference surface to volume sense data into the reference vectors
ErrorCode load_stp_surf_sense_data( Interface* moab, EntityHandle surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out ){

  int surf_id = geom_id_by_handle( moab, surf );
  switch(surf_id)
  {
    case 1:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 2:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_REVERSE); 
          break;

    case 3:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_REVERSE); 
          break;

    case 4:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_REVERSE); 
          break;

    case 5:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_REVERSE); 
          break;

    case 6:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_REVERSE); 
          break;

    case 7:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
  
    case 8:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;

    case 9:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;

    case 10:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
    default:
      std::cout << "Failure to find surface sense reference data. Returning failure..." << std::endl;
          return MB_FAILURE;
   }
  return MB_SUCCESS;
}

//Loads reference surface to volume sense data into the reference vectors
ErrorCode load_precgm14_stp_surf_sense_data( Interface* moab, EntityHandle surf, std::vector<int>& vol_ids_out, std::vector<int>& senses_out ){

  int surf_id = geom_id_by_handle( moab, surf );
  switch(surf_id)
  {
    case 1:
    case 2:
    case 3:
    case 4:
    case 5:
    case 6:
          vol_ids_out.push_back(1);
          senses_out.push_back(SENSE_FORWARD); 
          break;

    case 7:
    case 8:
    case 9:
    case 10:
          vol_ids_out.push_back(2);
          senses_out.push_back(SENSE_FORWARD);
          break;
  
    default:
      std::cout << "Failure to find surface sense reference data. Returning failure..." << std::endl;
          return MB_FAILURE;
   }
  return MB_SUCCESS;
}
