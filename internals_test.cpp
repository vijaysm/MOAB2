#include "MBInternals.hpp"
#include <iostream>
using namespace std;

bool internal_assert( bool c ) { return !c; }

#define handle_test_assert( C ) \
  if ( internal_assert(C) ) {\
    cout << "Test: " #C " failed ." << endl; \
    return false; \
  }

bool handle_test( MBEntityType type, MBEntityHandle id, int proc, bool should_fail )
{
  int err = 0;
  MBEntityHandle handle = CREATE_HANDLE( type, id, proc, err );
  if (should_fail) {
    handle_test_assert( err )
    return true;
  }
  handle_test_assert( !err )
  
  MBEntityType type_from_handle = TYPE_FROM_HANDLE(handle);
  handle_test_assert( type_from_handle == type )
  
  unsigned long id_from_handle = ID_FROM_HANDLE(handle);
  handle_test_assert( id_from_handle == id )
  
  int proc_from_handle = PROC_FROM_HANDLE(handle);
  handle_test_assert( proc_from_handle == proc )
  
  return true;
}

#define tag_test_assert( C ) \
  if ( internal_assert(C) ) {\
    cout << "Test: " #C " failed." << endl; \
    return false; \
  }

bool tag_test( MBTagId id, MBTagType prop )
{
  MBTag tag = TAG_HANDLE_FROM_ID( id, prop );
  
  unsigned long id_from_handle = ID_FROM_TAG_HANDLE(tag);
  tag_test_assert( id_from_handle == id )
  
  MBTagType prop_from_handle = PROP_FROM_TAG_HANDLE(tag);
  tag_test_assert( prop_from_handle == prop )
  
  return true;
}

  
//! Sample code on a 32-bit system.
int main()
{
  const unsigned cpus[] = { 1, 4, 16, 5, 20 };
  const int num_cpus = sizeof(cpus)/sizeof(cpus[0]);
  unsigned errors = 0, tests = 0;
  const int num_prop = MB_TAG_LAST;
  
  ++tests;
  if (MB_TAG_LAST > num_prop) {
    cout << "MB_TAG_PROP_WIDTH insufficient for size of MBTagType" << endl;
    ++errors;
  }
  
  ++tests;
  if (MBMAXTYPE > 1<<MB_TYPE_WIDTH) {
    cout << "MB_TYPE_WIDTH insufficient for size of MBEntityType" << endl;
    ++errors;
  }
  
    // if any errors so far, abort because everything else will
    // probably fail.
  if (errors)
    return errors;
  
  for (int num_cpu = 0; num_cpu < num_cpus; ++num_cpu) {
    
    MB_SET_PROC(cpus[num_cpu], 0);
    
    // init these after setting num_cpu, because max id depends on num cpu.
    const MBEntityHandle ids[] = {0, 1, MB_END_ID/2, MB_END_ID-1, MB_END_ID};
    const MBTagId tids[] = {0, 1, MB_TAG_PROP_MASK/2, MB_TAG_PROP_MASK-1, MB_TAG_PROP_MASK};
    const int num_ids = sizeof(ids)/sizeof(ids[0]);
    const int num_tids = sizeof(tids)/sizeof(tids[0]);
    
    for (unsigned cpu = 0; cpu < cpus[num_cpu]; ++cpu) {
      for (MBEntityType type = MBVERTEX; type < MBMAXTYPE; ++type) 
        for (int id = 0; id < num_ids; ++id) {
        ++tests;
         if (!handle_test( type, ids[id], cpu, false )) {
            cout << "Test of handle with type=" << type << ", id=" << ids[id]
                 << ", proc=" << cpu << ", and numproc=" << num_cpu << endl;
            ++errors;
          }
        }
          
      for (int prop = 0; prop < num_prop; ++prop)
        for (int id = 0; id < num_tids; ++id) {
          ++tests;
          if (!tag_test( tids[id], (MBTagType)prop )) {
            cout << "Test of tag handle with prop=" << prop << ", id=" << tids[id]
                 << ", proc=" << cpu << ", and numproc=" << num_cpu << endl;
            ++errors;
          }
        }
    }
  }
  
    // test some stuff that should fail
  MB_SET_PROC(16, 0);
  ++tests;
  if (!handle_test( MBVERTEX, MB_END_ID+1, 0, true)) {
    cout << "Failed to catch ID overflow" << endl;
    ++errors;
  }
  ++tests;
  if (!handle_test( (MBEntityType)(MBMAXTYPE+1), 1, 0, true)) {
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
