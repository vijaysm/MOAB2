#ifndef TEST_RUNNER_HPP
#define TEST_RUNNER_HPP

#include "TestUtil.hpp"



/***************************************************************************************
 * Begin test runner implememtation.
 * This is a higher-level API that can be used to register tests,
 * test dependencies, and to run-time select a subset of tests to 
 * run.
 ***************************************************************************************/

/* Register a test to be run */
#define REGISTER_TEST( TEST_FUNC ) \
  runner_register_test( __FILE__, __LINE__, #TEST_FUNC, (TEST_FUNC), NULL )

/* Mark a dependency between tests.  The second argument must be
 * an alredy-regsitered test.  The first argument will be registered
 * as a test if it has not already been registered.  The test specified
 * by the first argument will be run only if the test specified by
 * the second argument is run and succeeds.
 */
#define REGISTER_DEP_TEST( TEST_FUNC, REQUIRED_FUNC ) \
  runner_register_test( __FILE__, __LINE__, #TEST_FUNC, (TEST_FUNC), (REQUIRED_FUNC) )

/* Run registered tests.  
 * Arguments should be argc and argv passed to main.
 * If ARGC is less than or equal to 1 then all tests are run.
 * Otherwse only tests specified in the argument list are run.
 * Returns number of failed tests.
 */
#define RUN_TESTS( ARGC, ARGV ) \
  runner_run_tests( (ARGC), (ARGV) )



/***************************************************************************************
 * NOTE: The remainder of this file contains the implementation of the above macros.
 *       The above macros constitute the entire intended API.
 ***************************************************************************************/



static void runner_register_test( const char* filename, int line_number,
                                  const char* name, test_func function, 
                                  test_func requisite = 0 );
static int runner_run_tests( int argc, char* argv[] );

static void runner_usage( FILE* str, int argc, char* argv[] );
static void runner_list_tests(int long_format);
static void runner_help( int argc, char* argv[] );

enum RunnerStatus { PASSED, FAILED, DESELECTED, SELECTED };
struct RunnerTest {
  test_func testFunc;
  char* testName;
  enum RunnerStatus testStatus;
  int* testRequisites;
  size_t numRequisites;
};


struct RunnerTest* RunnerTestList = 0;
size_t RunnerTestCount = 0;
const size_t RUNNER_NOT_FOUND = ~(size_t)0;
static size_t runner_find_test_func( test_func f );
static size_t runner_find_test_name( const char* name );
static size_t runner_add_test( test_func f, const char* name );
static void runner_add_requisite( size_t idx, size_t req );
static void free_test_list();

static size_t runner_find_test_func( test_func f ) {
  for (size_t i = 0; i < RunnerTestCount; ++i)
    if (RunnerTestList[i].testFunc == f)
      return i;
  return RUNNER_NOT_FOUND;
}
static size_t runner_find_test_name( const char* name ) {
  for (size_t i = 0; i < RunnerTestCount; ++i)
    if (!strcmp(RunnerTestList[i].testName, name))
      return i;
  return RUNNER_NOT_FOUND;
}
static size_t runner_add_test( test_func f, const char* name ) {
  size_t idx = runner_find_test_func( f );
  if (idx == RUNNER_NOT_FOUND) {
    if (!RunnerTestCount)
      atexit( &free_test_list );
    idx = RunnerTestCount++;
    RunnerTestList = (RunnerTest*)realloc( RunnerTestList, RunnerTestCount * sizeof(RunnerTest) );
    RunnerTestList[idx].testFunc = f;
    RunnerTestList[idx].testName = strdup(name);
    RunnerTestList[idx].testStatus = SELECTED;
    RunnerTestList[idx].testRequisites = 0;
    RunnerTestList[idx].numRequisites = 0;
  }
  return idx;
}
static void runner_add_requisite( size_t idx, size_t req )
{
  size_t i;
  for (i = 0; i < RunnerTestList[idx].numRequisites; ++i)
    if (RunnerTestList[idx].testRequisites[i] == (int)req)
      return;
  ++RunnerTestList[idx].numRequisites;
  RunnerTestList[idx].testRequisites = (int*)realloc( RunnerTestList[idx].testRequisites,
                            RunnerTestList[idx].numRequisites * sizeof(*RunnerTestList[idx].testRequisites) );
  RunnerTestList[idx].testRequisites[RunnerTestList[idx].numRequisites-1] = req;
}
static void free_test_list()
{
  for (size_t i = 0; i < RunnerTestCount; ++i) {
    free( RunnerTestList[i].testName );
    free( RunnerTestList[i].testRequisites );
  }
  free( RunnerTestList );
}

void runner_register_test( const char* filename,
                           int line_number,
                           const char* name, 
                           test_func test,
                           test_func req )
{
  size_t i = runner_add_test( test, name );
  size_t req_idx;
  if (req) {
    req_idx = runner_find_test_func( req );
    if (RUNNER_NOT_FOUND == req_idx) {
      fprintf( stderr, "Error registering requisite for test: \"%s\"\n"
                       "\tat %s:%d\n"
                       "\tRequisite test funciton not registered.\n",
                       name, filename, line_number );
      abort();
    }
    runner_add_requisite( i, req_idx );
  }
}

void runner_usage( FILE* str, int argc, char* argv[] ) 
{
  fprintf( str, "%s [-l|-L] [-h] [-r] [<test_name> [<test_name> ...]]\n", argv[0] );
}

void runner_help( int argc, char* argv[] ) 
{
  runner_usage( stdout, argc, argv );
  fprintf( stdout, "-l : List test names and exit\n"
                   "-L : List test names and requisites and exit\n"
                   "-h : This help text\n"
                   "-r : Recursively run requisite tests for any specified test\n"
                   "\n");
}

void runner_list_tests( int long_format ) 
{
  size_t i, j;
  printf("Test List:\n");
  for (i = 0; i < RunnerTestCount; ++i) {
    if (RunnerTestList[i].testStatus == DESELECTED)
      continue;
    printf( " o %s\n", RunnerTestList[i].testName );
    if (!long_format || ! RunnerTestList[i].numRequisites) 
      continue;
    if (RunnerTestList[i].numRequisites == 1)
      printf( "  Requires : %s\n", RunnerTestList[RunnerTestList[i].testRequisites[0]].testName );
    else {
      printf( "  Requires : \n" );
      for (j = 0; j < RunnerTestList[i].numRequisites; ++j)
        printf( "    - %s\n", RunnerTestList[ RunnerTestList[i].testRequisites[j] ].testName );
    }
  }
}

int runner_run_tests( int argc, char* argv[] )
{
    /* Counters */
  int error_count = 0;
  int fail_count = 0;
  int num_selected = 0;
  int num_run = 0;

    /* Flags from parsed arguments */
  int run_requisites = 0;
  int list_tests = 0;
  int first_selected = 1;
  
    /* Misc iterator vars and such */
  int changed_some, ran_some, can_run, fail;
  int k;
  const char* c;
  size_t i, j;
  
    /* Pricess command line arguments */
  for (k = 1; k < argc; ++k) {
    if (argv[k][0] == '-') {
      for (c = argv[k] + 1; *c; ++c) {
        switch (*c) {
          case 'l': list_tests = 1; break;
          case 'L': list_tests = 2; break;
          case 'r': run_requisites = true; break;
          case 'h': runner_help( argc, argv ); return 0;
          default:
            runner_usage( stderr, argc, argv );
            fprintf( stderr, "Unknown flag: '%c'\n", *c );
            return 1;
        }
      }
    }
    else {
        // If user has specified some tests to run, begin
        // by marking all tests as de-selected.
      if (first_selected) {
        for (i = 0; i < RunnerTestCount; ++i)
          RunnerTestList[i].testStatus = DESELECTED;
        first_selected = 0;
      }
        // Mark specified test as selected.
      i = runner_find_test_name( argv[k] );
      if (RUNNER_NOT_FOUND == i) {
        fprintf( stderr, "Unknown test name: \"%s\"\n", argv[k] );
        ++error_count;
      }
      else {
        RunnerTestList[i].testStatus = SELECTED;
      }
    }
  }
  
    /* If recursively running requisite tests, select those also. */
  if (run_requisites) {
    do {
      changed_some = 0;
      for (i = 0; i < RunnerTestCount; ++i) {
        if (RunnerTestList[i].testStatus == DESELECTED)
          continue;
        
        for (j = 0; j < RunnerTestList[i].numRequisites; ++j) {
          if (RunnerTestList[ RunnerTestList[i].testRequisites[j] ].testStatus == DESELECTED) {
            RunnerTestList[ RunnerTestList[i].testRequisites[j] ].testStatus = SELECTED;
            changed_some = 1;
          }
        }
      }
    } while(changed_some);
  }
  
    // Count number of selected tests
  num_selected = 0;
  for (i = 0; i < RunnerTestCount; ++i) 
    if (RunnerTestList[i].testStatus == SELECTED)
      ++num_selected;
  
  if (list_tests) {
    runner_list_tests( list_tests - 1 );
    return error_count;
  }
  
    // Now run the tests
  num_run = 0;
  do {
    ran_some = 0;
    for (i = 0; i < RunnerTestCount; ++i) {
      if (RunnerTestList[i].testStatus != SELECTED)
        continue;
      can_run = 1;
      for (j = 0; j < RunnerTestList[i].numRequisites; ++j) {
        k = RunnerTestList[i].testRequisites[j];
        if (RunnerTestList[k].testStatus != PASSED &&
            RunnerTestList[k].testStatus != DESELECTED) {
          can_run = 0;
          break;
        }
      }
      if (!can_run)
        continue;
      
      ran_some = 1;
      ++num_run;
      fail = run_test( RunnerTestList[i].testFunc, RunnerTestList[i].testName );
      if (fail) {
        error_count++;
        fail_count++;
        RunnerTestList[i].testStatus = FAILED;
      }
      else {
        RunnerTestList[i].testStatus = PASSED;
      }
    }
  } while (ran_some);
  
    // Print brief summary
  if (num_run == (int)RunnerTestCount && !fail_count) {
    printf("All %d tests passed.\n", num_run);
  }
  else if (num_run == num_selected && !fail_count) {
    printf("All %d selected tests passed.\n", num_run );
    printf("Skipped %d non-selected tests\n", (int)(RunnerTestCount - num_selected));
  }
  else {
    printf( "%2d tests registered\n", (int)RunnerTestCount );
    if (num_selected == num_run) 
      printf( "%2d tests selected\n", num_selected );
    else
      printf( "%2d of %2d tests ran\n", num_run, num_selected );
    if (num_run < (int)RunnerTestCount)
      printf( "%2d of %2d tests skipped\n", 
        (int)RunnerTestCount - num_run, (int)RunnerTestCount );
    printf( "%2d of %2d tests passed\n", num_run - fail_count, num_run );
    if (fail_count) 
      printf( "%2d of %2d tests FAILED\n", fail_count, num_run );
  }
  
  return error_count;
}

#endif
