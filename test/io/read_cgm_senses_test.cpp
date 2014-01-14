
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
using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)


#ifdef MESHDIR
static const char input_cube[] = STRINGIFY(MESHDIR) "/io/cube.sat";
#else
static const char input_cube[] = "/io/cube.sat";
#endif

// Function used to load the test file
void read_file( Interface* moab, const char* input_file );

// Functions containing known sense data
void load_curve_sense_data( Interface* moab, EntityHandle curve,  std::vector<int>& surf_ids_out, std::vector<int>& senses_out );

// Functions used to compare sense information found in 
// the model to reference information

void check_curve_sense_data( Interface* moab, std::vector<EntityHandle> wrt_ents, std::vector<int> senses, 
                             std::vector<int> known_wrt_ids, std::vector<int> known_senses);

//Function used to get id's from entity handles
int geom_id_by_handle( Interface* moab, const EntityHandle set );

// List of tests in this file
void read_cube_curve_senses_test();
void delete_mesh_test();


int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result += RUN_TEST( read_cube_curve_senses_test );  

  return result;
}


void read_file( Interface* moab, const char* input_file )
{
  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

void read_cube_curve_senses_test()
{
  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cube );
  
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
  CHECK_EQUAL(12 , number_of_curves);
  
  // Get curve handles
  Range curves;
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &geom_tag,
	  					    val, 1, curves );
  CHECK_ERR(rval);

  // Establish GeomTopoTool instance needed to get curve data 
  moab::GeomTopoTool gt(mb, false);  
  std::vector<EntityHandle> surfs;
  std::vector<int> senses;  
  std::vector<int> known_surf_ids;
  std::vector<int> known_senses;

for(unsigned int i = 0; i < curves.size() ; i++)
  {
  
  surfs.clear();
  senses.clear();
  //Curve 1
  gt.get_senses( curves[i], surfs, senses);
  CHECK_ERR(rval);

  //Load known curve-sense data

   known_surf_ids.clear();
   known_senses.clear();
   load_curve_sense_data( mb, curves[i], known_surf_ids, known_senses );

   check_curve_sense_data( mb, surfs, senses, known_surf_ids, known_senses);
  }

}


int geom_id_by_handle( Interface* moab, const EntityHandle set ) {
    
ErrorCode rval;
    Tag id_tag;
    rval = moab->tag_get_handle( GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, id_tag, moab::MB_TAG_DENSE);
    assert(MB_SUCCESS==result || MB_ALREADY_ALLOCATED==result);                       
    int id;
    rval = moab->tag_get_data( id_tag, &set, 1, &id );                  
    assert(MB_SUCCESS == result);                           
    return id;
 }

void check_curve_sense_data( Interface* moab, std::vector<EntityHandle> wrt_ents, std::vector<int> senses, 
                             std::vector<int> known_wrt_ids, std::vector<int> known_senses)
{
  
  std::vector<int> wrt_ent_ids;
  for( unsigned int i=0 ; i<wrt_ents.size() ; i++)
  {
      wrt_ent_ids.push_back(geom_id_by_handle(moab, wrt_ents[i]));
  }

  for (unsigned int i=0; i< wrt_ent_ids.size() ; i++)
  {
     for(unsigned int j=0; j< known_wrt_ids.size(); j++)
     {
       if(wrt_ent_ids[i] == known_wrt_ids [j])
         {
          CHECK_EQUAL(senses[i],known_senses[j]);
          known_wrt_ids.erase(known_wrt_ids.begin()+j);
          known_senses.erase(known_senses.begin()+j);
         }
     }
  }

  int leftovers = known_wrt_ids.size();

  CHECK_EQUAL(leftovers, 0 );

}
void load_curve_sense_data( Interface* moab, EntityHandle curve, std::vector<int>& surf_ids_out, std::vector<int>& senses_out ){

  int curve_id = geom_id_by_handle( moab, curve);

  std::cout << curve_id << std::endl;


  switch(curve_id)
  {
    case 1:
          surf_ids_out.push_back(1); surf_ids_out.push_back(6);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 2:
          surf_ids_out.push_back(1); surf_ids_out.push_back(5);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 3:
          surf_ids_out.push_back(1); surf_ids_out.push_back(4);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 4:
          surf_ids_out.push_back(1); surf_ids_out.push_back(3);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 5:
          surf_ids_out.push_back(2); surf_ids_out.push_back(6);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 6:
          surf_ids_out.push_back(2); surf_ids_out.push_back(3);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 7:
          surf_ids_out.push_back(2); surf_ids_out.push_back(4);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 8:
          surf_ids_out.push_back(2); surf_ids_out.push_back(5);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 9:
          surf_ids_out.push_back(3); surf_ids_out.push_back(4);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 10:
          surf_ids_out.push_back(3); surf_ids_out.push_back(6);
          senses_out.push_back(-1); senses_out.push_back(1);

    break;
    case 11:
          surf_ids_out.push_back(4); surf_ids_out.push_back(5);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
    case 12:
          surf_ids_out.push_back(5); surf_ids_out.push_back(6);
          senses_out.push_back(1); senses_out.push_back(-1);

    break;
  } 

}

