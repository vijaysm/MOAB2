
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "TestUtil.hpp"
#include "Internals.hpp"
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"
#include "InitCGMA.hpp"
#include "GeometryQueryTool.hpp"

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

//Function for getting entity ids
int geom_id_by_handle( Interface* moab, const EntityHandle set );

//Function for checking retrieved group data
void check_group_data( std::vector<int> & group_ids, std::vector<std::string> & group_names,
#ifdef HAVE_OCC_STEP
                      std::vector<int> & );
#else
                      std::vector<int> & group_ent_ids );
#endif

//Function for loading all reference data
void load_group_references( std::vector<int>& ids, std::vector<std::string>& names, std::vector<int>& ent_ids);

// List of tests in this file
void read_cylcube_groups_test();


int main(int /* argc */, char** /* argv */)
{
  int result = 0;

  result+=RUN_TEST(read_cylcube_groups_test);
 
  return result;
}



void read_file( Interface* moab, const char* input_file )
{
  InitCGMA::initialize_cgma();
  GeometryQueryTool::instance()->delete_geometry();

  ErrorCode rval = moab->load_file( input_file );
  CHECK_ERR(rval);
}

// Checks that group information is being read correctly if MOAB is 
// being build with CGM. If MOAB is built with OCC, it makes sure
// no erroneous group data is loaded as STEP files do not hold 
// information about groups.
void read_cylcube_groups_test()
{

  ErrorCode rval;
  //Open the test file
  Core moab;
  Interface* mb = &moab;
  read_file( mb, input_cylcube);  

  //Get (or create) the name and category tags
  Tag name_tag, category_tag;

  rval = mb->tag_get_handle( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TYPE_OPAQUE,
                             name_tag, moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);

  rval = mb->tag_get_handle( CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE,
                             MB_TYPE_OPAQUE, category_tag,
                             moab::MB_TAG_SPARSE|moab::MB_TAG_CREAT );
  CHECK_ERR(rval);

  //Get the group entity handles
  Range group_sets;
  char query[CATEGORY_TAG_SIZE] = "Group\0";   
  //Has to be this way because of the way tags are created
  void* val[] = {&query};
  rval = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &category_tag, val, 1, group_sets);
  CHECK_ERR(rval);
  //Get group names and IDs
  std::vector<int> g_ids;
  std::vector<std::string> g_names;
  std::vector<int> g_ent_ids;

  for(Range::iterator i=group_sets.begin(); i!=group_sets.end(); ++i)
    {
      int group_id = geom_id_by_handle( mb, *i );
      g_ids.push_back(group_id);
      //Get the group name
      char group_name[NAME_TAG_SIZE+1];
      rval = mb->tag_get_data( name_tag, &(*i), 1, &group_name);
      CHECK_ERR(rval);
      //Store group name
      std::string temp(group_name);
      g_names.push_back(temp);
      //Get all entities in the group
      Range group_ents;
      rval = mb->get_entities_by_type( *i, MBENTITYSET, group_ents, false );
      CHECK_ERR(rval);
      if( group_ents.size() != 1) CHECK(false);
      int grp_ent_id = geom_id_by_handle( mb, group_ents[0] );
      g_ent_ids.push_back(grp_ent_id);
  
    }
  check_group_data( g_ids, g_names, g_ent_ids );
}

void check_group_data(std::vector<int>& group_ids, std::vector<std::string>& group_names,
#ifdef HAVE_OCC_STEP
                      std::vector<int>& /*group_ent_ids*/ )
#else
                      std::vector<int>& group_ent_ids )
#endif
{
  // Step files do not contain group data, MOAB shouldn't return errors when trying to access
  // this data but there shouldn't be any found.
#ifdef HAVE_OCC_STEP
  int num_g_ids = group_ids.size();
  int num_g_names = group_names.size();

  CHECK_EQUAL( 0, num_g_ids );
  CHECK_EQUAL( 0, num_g_names );
#else


  //Initialize reference data
  std::vector<int> group_ref_ids;
  std::vector<std::string> group_ref_names;
  std::vector<int> group_ref_ent_ids;
  load_group_references( group_ref_ids, group_ref_names, group_ref_ent_ids );

  // check that the correct number of entities were found
  CHECK_EQUAL ( group_ref_ids.size(), group_ids.size() );
  CHECK_EQUAL ( group_ref_names.size(), group_names.size() );  
  CHECK_EQUAL ( group_ref_ent_ids.size(), group_ent_ids.size() );

  //now make sure that each group has a matching group
  for(unsigned int i=0 ; i<group_ids.size(); i++)
    {
      for(unsigned int j=0; j<group_ref_ids.size(); j++)
	{
          if( group_ids[i]==group_ref_ids[j]
	      && group_names[i] ==  group_ref_names[j]
              && group_ent_ids[i]==group_ref_ent_ids[j])
	    {
              group_ref_ids.erase(group_ref_ids.begin()+j);
              group_ref_names.erase(group_ref_names.begin()+j);
              group_ref_ent_ids.erase(group_ref_ent_ids.begin()+j);
              continue;
	    }
	}
    }

  // Check sizes of reference vectors after matching
  // (all should be zero)
   int leftovers = group_ref_ids.size();
  CHECK_EQUAL( 0, leftovers );
  leftovers = group_ref_names.size();
  CHECK_EQUAL( 0, leftovers );
  leftovers = group_ref_ent_ids.size();
  CHECK_EQUAL( 0, leftovers );
#endif
}

void load_group_references( std::vector<int>& ids, std::vector<std::string>& names, std::vector<int>& ent_ids )
{
  //First set of group info
  names.push_back("Group 3"); ids.push_back(3); ent_ids.push_back(2);
  
  //Second set of group info
  names.push_back("Group 2"); ids.push_back(2); ent_ids.push_back(1);
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
