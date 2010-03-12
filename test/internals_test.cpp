#include "Internals.hpp"
#include <iostream>
using namespace std;

using namespace moab;

HandleUtils handleUtils(0,1);

bool internal_assert( bool c ) { return !c; }

#define handle_test_assert( C ) \
  if ( internal_assert(C) ) {\
    cout << "Test: " #C " failed ." << endl; \
    return false; \
  }

bool handle_test( EntityType type, EntityID id, int proc, bool should_fail )
{
  int err = 0;
  EntityHandle handle = CREATE_HANDLE( type, handleUtils.create_id(id, proc), err );
  if (should_fail) {
    handle_test_assert( err )
    return true;
  }
  handle_test_assert( !err )
  
  EntityType type_from_handle = TYPE_FROM_HANDLE(handle);
  handle_test_assert( type_from_handle == type )
  
  EntityID id_from_handle = handleUtils.id_from_handle(handle);
  handle_test_assert( id_from_handle == id )
  
  int proc_from_handle = handleUtils.rank_from_handle(handle);
  handle_test_assert( proc_from_handle == proc )
  
  return true;
}

#define tag_test_assert( C ) \
  if ( internal_assert(C) ) {\
    cout << "Test: " #C " failed." << endl; \
    return false; \
  }

bool tag_test( TagId id, TagType prop )
{
  Tag tag = TAG_HANDLE_FROM_ID( id, prop );
  
  unsigned long id_from_handle = ID_FROM_TAG_HANDLE(tag);
  tag_test_assert( id_from_handle == id )
  
  TagType prop_from_handle = PROP_FROM_TAG_HANDLE(tag);
  tag_test_assert( prop_from_handle == prop )
  
  return true;
}

  
//! Sample code on a 32-bit system.
int main()
{
  const unsigned cpus[] = { 1, 4, 16, 5, 20 };
  const int num_cpus = sizeof(cpus)/sizeof(cpus[0]);
  unsigned errors = 0, tests = 0;
  const int num_prop = MB_TAG_LAST+1;
  
  ++tests;
  if (MB_TAG_LAST > num_prop) {
    cout << "MB_TAG_PROP_WIDTH insufficient for size of TagType" << endl;
    ++errors;
  }
  
  ++tests;
  if (MBMAXTYPE > 1<<MB_TYPE_WIDTH) {
    cout << "MB_TYPE_WIDTH insufficient for size of EntityType" << endl;
    ++errors;
  }
  
    // if any errors so far, abort because everything else will
    // probably fail.
  if (errors)
    return errors;
  
  for (int num_cpu = 0; num_cpu < num_cpus; ++num_cpu) {
    
    handleUtils = HandleUtils( 0, cpus[num_cpu] );
    
    // init these after setting num_cpu, because max id depends on num cpu.
    const EntityID ids[] = {0, 1, handleUtils.max_id()/2, handleUtils.max_id()-1, handleUtils.max_id()};
    const TagId tids[] = {0, 1, MB_TAG_PROP_MASK/2, MB_TAG_PROP_MASK-1, MB_TAG_PROP_MASK};
    const int num_ids = sizeof(ids)/sizeof(ids[0]);
    const int num_tids = sizeof(tids)/sizeof(tids[0]);
    
    for (unsigned cpu = 0; cpu < cpus[num_cpu]; ++cpu) {
      for (EntityType type = MBVERTEX; type < MBMAXTYPE; ++type) 
        for (int id = 0; id < num_ids; ++id) {
        ++tests;
         if (!handle_test( type, ids[id], cpu, false )) {
            cout << "Test of handle with type=" << type << ", id=" << ids[id]
                 << ", proc=" << cpu << ", and numproc=" << cpus[num_cpu] << endl;
            ++errors;
          }
        }
          
      for (int prop = 0; prop < num_prop; ++prop)
        for (int id = 0; id < num_tids; ++id) {
          ++tests;
          if (!tag_test( tids[id], (TagType)prop )) {
            cout << "Test of tag handle with prop=" << prop << ", id=" << tids[id]
                 << ", proc=" << cpu << ", and numproc=" << cpus[num_cpu] << endl;
            ++errors;
          }
        }
    }
  }
  
    // test some stuff that should fail
  handleUtils = HandleUtils(0, 16);
  ++tests;
  if (!handle_test( MBVERTEX, MB_END_ID+1, 0, true)) {
    cout << "Failed to catch ID overflow" << endl;
    ++errors;
  }
  ++tests;
  if (!handle_test( (EntityType)(MBMAXTYPE+1), 1, 0, true)) {
    cout << "Failed to catch type overflow" << endl;
    ++errors;
  }
//  ++tests;
//  if (!handle_test( MBHEX, 1, 17, true)) {
//    cout << "Failed to catch Proc# overflow" << endl;
//    ++errors;
//  }
  
  if (errors) 
    cout << endl << errors << " of " << tests << " tests failed." << endl << endl;
  else
    cout << endl << tests << " tests passed." << endl << endl;
  return errors;
}
