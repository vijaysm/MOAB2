
#include <iostream>
#include "moab/Interface.hpp"
#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif
#include "Internals.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/Core.hpp"
#include "moab/Range.hpp"

using namespace moab;

#define CHKERR(A) do { if (MB_SUCCESS != (A)) { \
  std::cerr << "Failure (error code " << (A) << ") at " __FILE__ ":" \
            << __LINE__ << std::endl; \
  return A; } } while(false)

ErrorCode gather_related_test() 
{
    // create an entityset structure and test related entities function
    // sets: A
    //       |
    //       B  C
    //       |/ |
    //       D  E
    // if D is passed in to gather_related_ents, A-D should be returned, and E should not be
    //
  EntityHandle A, B, C, D, E;
  Core mb;
  ErrorCode rval = mb.create_meshset( MESHSET_SET, A ); CHKERR(rval);
  rval = mb.create_meshset( MESHSET_SET, B ); CHKERR(rval);
  rval = mb.create_meshset( MESHSET_SET, C ); CHKERR(rval);
  rval = mb.create_meshset( MESHSET_SET, D ); CHKERR(rval);
  rval = mb.create_meshset( MESHSET_SET, E ); CHKERR(rval);

  rval = mb.add_parent_child(A, B); CHKERR(rval);
  rval = mb.add_parent_child(B, D); CHKERR(rval);
  rval = mb.add_parent_child(C, D); CHKERR(rval);
  rval = mb.add_parent_child(C, E); CHKERR(rval);
  
    // now test it
  ReadUtilIface* readMeshIface;
  mb.Interface::query_interface(readMeshIface);

  Range init_range, set_range, all_sets(A, E);
  init_range.insert(D);
  rval = readMeshIface->gather_related_ents(init_range, set_range); CHKERR(rval);
  
  if (set_range.size() != 4) return MB_FAILURE;
  all_sets -= set_range;
  if (1 != all_sets.size() || *all_sets.begin() != E) return MB_FAILURE;
  
  return MB_SUCCESS;
}
  
int number_tests = 0;
int number_tests_failed = 0;
#define RUN_TEST( A ) _run_test( (A), #A )


typedef ErrorCode (*TestFunc)();
static void _run_test( TestFunc func, const char* func_str ) 
{
  ++number_tests;
  std::cout << "   " << func_str << ": ";
  std::cout.flush();
  ErrorCode error = func( );
  
  if (MB_SUCCESS == error)
    std::cout << "Success" << std::endl;
  else if (MB_FAILURE == error)
    std::cout << "Failure" << std::endl;
  else {
    Core moab;
    std::cout << "Failed: " << moab.get_error_string( error ) << std::endl;
  }
  
  if (MB_SUCCESS != error) {
    ++number_tests_failed;
  }
}

int main(int argc, char* argv[])
{
  RUN_TEST( gather_related_test );

  std::cout << "\nMB TEST SUMMARY: \n"
       << "   Number Tests:           " << number_tests << "\n"
       << "   Number Successful:      " << number_tests - number_tests_failed << "\n"
       << "   Number Failed:          " << number_tests_failed 
       << "\n\n";

  return number_tests_failed;
}
