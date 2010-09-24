#ifndef TEST_UTIL_HPP
#define TEST_UTIL_HPP

/* Define these here because they are used by many tests
 * to find the add directory for input files */
#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

/* How to use this test suite utility:
 * 1) Write tests that use the CHECK and CHECK_* macros defined below to assert test conditions.
 * 2) Write a main routine that invokes each test through the RUN_TEST macro
 * 3) RUN_TEST evaluates to 1 if test failed, zero otherwize.  Count failures and print summary.
 */

/** Check that A is MB_SUCCESS */
#define CHECK_ERR( A )                    check_equal( MB_SUCCESS, (A), "MB_SUCCESS", #A, __LINE__, __FILE__ )
/**  Ensure that A is true */
#define CHECK( A )                        check_true( (A), #A, __LINE__, __FILE__ )
/** Check that two values are equal */
#define CHECK_EQUAL( EXP, ACT )           check_equal( (EXP), (ACT), #EXP, #ACT, __LINE__, __FILE__ )
/** Check that two real (float or double) values are equal within EPS */
#define CHECK_REAL_EQUAL( EXP, ACT, EPS ) check_equal( (EXP), (ACT), (EPS), #EXP, #ACT, __LINE__, __FILE__ )
/** Check that two arrays contain the same values in the same order */
#define CHECK_ARRAYS_EQUAL( EXP, EXP_LEN, ACT, ACT_LEN ) check_array_equal( (EXP), (EXP_LEN), (ACT), (ACT_LEN), #EXP, #ACT, __LINE__, __FILE__ )
/** Run a test
 *  Argument should be a function with the signature:  void func(void)
 *  Evaluates to zero if test is successful, one otherwise.
 */
#define RUN_TEST( FUNC )           run_test( &FUNC, #FUNC )


// Use C++ exceptions to return error state to test runner
// Portable, but whole test suite stops if any test segfaults, asserts, etc.
#define EXCEPTION_MODE 1   

// Test runner forks separate process for each test.
// Difficult to debug tests (with debugger).  Not portable to Windows.  
// Very robust (no test can distrub test running code)
#define FORK_MODE 2

// Use signal handler and long jumps to return error state to test runner.
// Might be portable to Windows (not sure).  Possibly undefined behavior (e.g. continuing 
// with next test after catching segfault is technically undefined behavior.)
// Also, tests can corrupt heap memory management, interferring with later tests.
// Leaks memory on test failure (no stack unwind).  This is actually a feature, as
// we don't care too much about tests leaking memory and trying to reconver memory
// might make things worse, depending on why the test failed.
#define LONGJMP_MODE 3      

// If test application hasn't set MODE, set to default
#ifndef MODE
#if defined(_MSC_VER) || defined(__MINGW32__)
#    define MODE EXCEPTION_MODE
#  else
#    define MODE LONGJMP_MODE
#  endif
#endif


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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#ifdef __cplusplus
#include <iostream>
#include <vector>
#endif

/***************************************************************************************
 *                     Define What to do when a test fails.
 ***************************************************************************************/

// For EXCEPTION_MODE, throw an exception when a test fails.
// This will unwind stack, recover memory, etc. 
#if MODE == EXCEPTION_MODE
   struct ErrorExcept{};
#  define FLAG_ERROR throw ErrorExcept()
// For FORK_MODE, the test is running in its own processs.  Just
// terminate the process with a non-zero exit code when the test
// fails.
#elif MODE == FORK_MODE
#  include <sys/types.h>
#  include <sys/wait.h>
#  include <unistd.h>
#  include <errno.h>
#  define FLAG_ERROR exit(1)
// For LONGJMP_MODE, we do a long jump to just before the test is
// run, with a return value of -1 to indicate failures (positive
// return codes are used if the test caused a segfault or other
// signal.)
#elif MODE == LONGJMP_MODE
#  include <signal.h>
#  include <setjmp.h>
#  define FLAG_ERROR siglongjmp( jmpenv, -1 )
#else
#  error "MODE not set"
#endif

/***************************************************************************************
 *                              Setup for LONGJMP_MODE
 ***************************************************************************************/

#if MODE == LONGJMP_MODE

// Variable to hold stack state for longjmp
sigjmp_buf jmpenv;

// Define signal handler used to catch errors such as segfaults.
// Signal handler does longjmp with the signal number as the 
// return value.
extern "C" {
  void sighandler( int sig ) {
    signal( sig, sighandler );
    siglongjmp(jmpenv, sig);
    // should never return from longjmp
    exit(1);
  }
  typedef void (*sigfunc_t)(int);
} // extern "C"

// Helper function to register signal handlers.  
int sethandler( int sig ) {
  sigfunc_t h = signal( sig, &sighandler );
  if (h == SIG_ERR)
    return  1;
   // If user-defined signal handler (or signal is ignored),
   // than unregister our handler.
  else if (h != SIG_DFL)
    signal( sig, h );
  return 0;
}

// Register signal handlers for all defined signals that typicall result
// in process termination.
int init_signal_handlers()
{
  int result = 0;
#ifdef SIGHUP
  result += sethandler( SIGHUP );
#endif
#ifdef SIGINT
  result += sethandler( SIGINT );
#endif
#ifdef SIGQUIT
  result += sethandler( SIGQUIT );
#endif
#ifdef SIGILL
  result += sethandler( SIGILL );
#endif
#ifdef SIGTRAP
  result += sethandler( SIGTRAP );
#endif
#ifdef SIGABRT
  result += sethandler( SIGABRT );
#endif
#ifdef SIGBUS
  result += sethandler( SIGBUS );
#endif
#ifdef SIGFPE
  result += sethandler( SIGFPE );
#endif
#ifdef SIGUSR1
  result += sethandler( SIGUSR1 );
#endif
#ifdef SIGSEGV
  result += sethandler( SIGSEGV );
#endif
#ifdef SIGUSR2
  result += sethandler( SIGUSR2 );
#endif
#ifdef SIGPIPE
  result += sethandler( SIGPIPE );
#endif
#ifdef SIGTERM
  result += sethandler( SIGTERM );
#endif
#ifdef SIGCHLD
  result += sethandler( SIGCHLD );
#endif
#ifdef SIGIO
  result += sethandler( SIGIO );
#endif
#ifdef SIGSYS
  result += sethandler( SIGSYS );
#endif
  return result;
}

// Declare a garbage global variable.  Use variable initialization to
// force call to init_signal_handlers().  
int junk_init_var = init_signal_handlers();

#endif // LONGJMP_MODE


/***************************************************************************************
 *                            Function to handle failed tests
 ***************************************************************************************/

// use a function rather than substituting FLAG_ERROR directly
// so we have a convenient place to set a break point
inline void flag_error() 
  { FLAG_ERROR; }


/***************************************************************************************
 *                            The Code to Run Tests
 ***************************************************************************************/


/* Make sure IS_BUILDING_MB is defined so we can include MBInternals.hpp */
#include "moab/Types.hpp"
#ifndef IS_BUILDING_MB
#  define IS_BUILDING_MB
#  include "Internals.hpp"
#  undef IS_BUILDING_MB
#else
#  include "Internals.hpp"
#endif

typedef void (*test_func)(void);
int run_test( test_func test, const char* func_name )
{
  printf("Running %s ...\n", func_name );
  
#if MODE == EXCEPTION_MODE
  /* On Windows, run all tests in same process.
     Flag errors by throwing an exception.
   */
  try {
    (*test)();
    return 0;
  }
  catch (ErrorExcept) {
    printf( "  %s: FAILED\n", func_name );
    return 1;
  }
  catch (...) {
    printf( "  %s: UNCAUGHT EXCEPTION\n", func_name );
    return 1;
  }
    
#elif MODE == FORK_MODE
    /* For non-Windows OSs, fork() and run test in child process. */
  pid_t pid = fork();
  int status;
  
    /* Fork failed? */
  if (pid == -1) {  
    perror( "fork()" );
    abort(); /* abort all tests (can't fork child processes) */
  }
  
    /* If child process*/
  if (pid == 0) {
    (*test)();  /* call test function */
    exit(0);    /* if function returned, then it succeeded */
  }
  
    /* If here, then parent process */
    
    /* Wait until child process exits */
  waitpid( pid, &status, 0 );
  
    /* Check child exit status */
  if (WIFSIGNALED(status)) {
    if (WTERMSIG(status))
      printf("  %s: TERMINATED (signal %d)\n", func_name, (int)WTERMSIG(status) );
    if (WCOREDUMP(status))
      printf("  %s: CORE DUMP\n", func_name);
    return 1;
  }
  else if(WEXITSTATUS(status)) {
    printf( "  %s: FAILED\n", func_name );
    return 1;
  }
  else {
    return 0;
  }
  
#elif MODE == LONGJMP_MODE
    // Save stack state at this location.
  int rval = sigsetjmp( jmpenv, 1 );
    // If rval is zero, then we haven't run the test yet. 
    // If rval is non-zero then
    // a) we ran the test
    // b) the test failed
    // c) we did a longjmp back to the location where we called setsigjmp.
    
    // run test
  if (!rval) {
    (*test)();
    return 0;
  }
    // some check failed
  else if (rval == -1) {
    printf( "  %s: FAILED\n", func_name );
    return 1;
  }
    // a signal was raised (e.g. segfault)
  else {
    printf( "  %s: TERMINATED (signal %d)\n", func_name, rval );
    return 1;
  }
#else
  #error "MODE not set"
#endif // MODE
}



/***************************************************************************************
 *                            CHECK_EQUAL implementations
 ***************************************************************************************/

// Common implementatation for most types
#define EQUAL_TEST_IMPL( TEST, TYPE ) if( !(TEST) ) { \
  printf( "Equality Test Failed: %s == %s\n", sA, sB ); \
  printf( "  at line %d of '%s'\n", line, file ); \
  printf( "  Expected value: %" #TYPE "\n", A ); \
  printf( "  Actual value:   %" #TYPE "\n", B ); \
  printf( "\n" ); \
  flag_error(); \
}

void check_equal( int A, int B, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( A == B, d ) }

void check_equal( unsigned A, unsigned B, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( A == B, u ) }

void check_equal( long A, long B, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( A == B, ld ) }

void check_equal( unsigned long A, unsigned long B, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( A == B, lu ) }

void check_equal( void* A, void* B, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( A == B, p ) }

void check_equal( const char* A, const char* B, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( !strcmp((A),(B)), s ) }

void check_equal( const std::string& A, const std::string& B, const char* sA, const char* sB, int line, const char* file )
  {  check_equal( A.c_str(), B.c_str(), sA, sB, line, file); }

void check_equal( float A, float B, float eps, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( fabsf(A - B) <= eps, f ) }

void check_equal( double A, double B, float eps, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( fabs(A - B) <= eps, f ) }

const char* mb_error_str( moab::ErrorCode err )
{
  switch (err) {
    case moab::MB_SUCCESS                 : return "Success";
    case moab::MB_INDEX_OUT_OF_RANGE      : return "Index Out of Range";
    case moab::MB_TYPE_OUT_OF_RANGE       : return "Type Out of Range";
    case moab::MB_MEMORY_ALLOCATION_FAILED: return "Memory Alloc. Failed";
    case moab::MB_ENTITY_NOT_FOUND        : return "Entity Not Found";
    case moab::MB_MULTIPLE_ENTITIES_FOUND : return "Multiple Entities Found";
    case moab::MB_TAG_NOT_FOUND           : return "Tag Not Found";
    case moab::MB_FILE_DOES_NOT_EXIST     : return "File Not Found";
    case moab::MB_FILE_WRITE_ERROR        : return "File Write Error";
    case moab::MB_NOT_IMPLEMENTED         : return "Not Implemented";
    case moab::MB_ALREADY_ALLOCATED       : return "Already Allocated";
    case moab::MB_VARIABLE_DATA_LENGTH    : return "Variable Data Length";
    case moab::MB_INVALID_SIZE            : return "Invalid Size";
    case moab::MB_UNSUPPORTED_OPERATION   : return "Unsupported Operation";
    case moab::MB_UNHANDLED_OPTION        : return "Unhandled Option";
    case moab::MB_FAILURE                 : return "Failure";
    default                         : return "(unknown)";
  }
}


// Special case for MBErrorCode, use mb_error_str() to print the 
// string name of the error code.
void check_equal( moab::ErrorCode A, moab::ErrorCode B, const char* sA, const char* sB, int line, const char* file )
{
  if (A == B)
    return;
  
  printf( "MBErrorCode Test Failed: %s == %s\n", sA, sB ); 
  printf( "  at line %d of '%s'\n", line, file ); 
  printf( "  Expected value: %s (%d)\n", mb_error_str(A), (int)A ); 
  printf( "  Actual value:   %s (%d)\n", mb_error_str(B), (int)B ); 
  printf( "\n" ); 
  flag_error(); 
}

const char* mb_type_str( moab::EntityType type )
{
  switch(type) {
    case moab::MBVERTEX    : return "Vertex";
    case moab::MBEDGE      : return "Edge";
    case moab::MBTRI       : return "Triangle";
    case moab::MBQUAD      : return "Quadrilateral";
    case moab::MBPOLYGON   : return "Polygon";
    case moab::MBTET       : return "Tetrahedron";
    case moab::MBPYRAMID   : return "Pyramid";
    case moab::MBPRISM     : return "Prism (wedge)";
    case moab::MBKNIFE     : return "Knife";
    case moab::MBHEX       : return "Hexahedron";
    case moab::MBPOLYHEDRON: return "Polyhedron";
    case moab::MBENTITYSET : return "Entity (Mesh) Set";
    case moab::MBMAXTYPE   : return "(max type)";
    default          : return "(unknown)";
  }
}

const char* mb_type_str( moab::EntityHandle a )
  { return mb_type_str( moab::TYPE_FROM_HANDLE(a) ); }
/*
void check_equal( moab::EntityHandle A, moab::EntityHandle B, const char* sA, const char* sB, int line, const char* file )
{
  if (A == B)
    return;
  
  printf( "Entity handles not equal: %s == %s\n", sA, sB );
  printf( "  at line %d of '%s'\n", line, file ); 
  if (A) 
    printf( "  Expected value: %lx (%s %ld)\n", (unsigned long)A, mb_type_str( A ), (long)ID_FROM_HANDLE(A) ); 
  else 
    printf( "  Expected value: 0\n" ); 
  if (B)
    printf( "  Actual value:   %lx (%s %ld)\n", (unsigned long)B, mb_type_str( B ), (long)ID_FROM_HANDLE(B) ); 
  else 
    printf( "  Actual value: 0\n" ); 
  printf( "\n" ); 
  flag_error(); 
}  
*/

void check_true( bool cond, const char* str, int line, const char* file )
{
  if( !cond ) { 
    printf( "Test Failed: %s\n", str ); 
    printf( "  at line %d of '%s'\n", line, file ); 
    printf( "\n" ); 
    flag_error(); 
  }
}

#ifdef __cplusplus

template <typename T>
void check_array_equal( const T* A, size_t A_size,
                        const T* B, size_t B_size, 
                        const char* sA, const char* sB, 
                        int line, const char* file )
{
  size_t i = 0;
  for (;;) {
    if (i == A_size && i == B_size)
      return; // equal
    else if (i == A_size || i == B_size)
      break; // differene lengths
    else if (A[i] != B[i])
      break;
    ++i;
  }
  
  std::cout << "Equality Test Failed: " << sA << " == " << sB << std::endl;
  std::cout << "  at line " << line << " of '" << file << "'" << std::endl;
  std::cout << "  Vectors differ at position " << i << std::endl;
  
    // print at most 10 values, roughly centered on the unequal one
  size_t count = 10, num_front_values = std::min(count/2,i);
  size_t max_len = std::max(A_size,B_size);
  if (i + count - num_front_values > max_len) {
    if (count > max_len) {
      num_front_values = i;
      count = max_len;
    }
    else {
      num_front_values = count - (max_len - i);
    }
  }
  
  std::cout << "  Expected: ";
  if (!A_size) {
    std::cout << "(empty)" << std::endl;
  }
  else {
    size_t j = i - num_front_values;
    size_t end = std::min(j + count, A_size);
    if (j) 
      std::cout << "... ";
    for (; j < end; ++j) {
      if (j == i)
        std::cout << '>' << A[j] << "< ";
      else
        std::cout << A[j] << " ";
    }
    if (end != A_size)
      std::cout << "...";
    std::cout << std::endl;
  }
  
  std::cout << "  Actual:   ";
  if (!B_size) {
    std::cout << "(empty)" << std::endl;
  }
  else {
    size_t j = i - num_front_values;
    size_t end = std::min(j + count, B_size);
    if (j) 
      std::cout << "... ";
    for (; j < end; ++j) {
      if (j == i)
        std::cout << '>' << B[j] << "< ";
      else
        std::cout << B[j] << " ";
    }
    if (end != B_size)
      std::cout << ", ...";
    std::cout << std::endl;
  }
  
  flag_error(); 
}
  
 
template <typename T>
void check_equal( const std::vector<T>& A, const std::vector<T>& B, 
                  const char* sA, const char* sB, 
                  int line, const char* file )
{
  check_array_equal( &A[0], A.size(), &B[0], B.size(), sA, sB, line, file );
}

#ifdef MOAB_RANGE_HPP

void check_equal( const moab::Range& A, const moab::Range& B, const char* sA, const char* sB, int line, const char* file )
{
  if (A == B)
    return;
    
  std::cout << "moab::ErrorCode Test Failed: " << sA << " == " << sB << std::endl;
  std::cout << "  at line " << line << " of '" << file << "'" << std::endl;
  std::cout << "   Expected: " << A << std::endl;
  std::cout << "   Actual  : " << B << std::endl;
  std::cout << std::endl;
  flag_error();
}

#endif  /* ifdef MOAB_RANGE_HPP */
    
#endif /* ifdef __cplusplus */

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
    printf( "%2d registered tests\n", (int)RunnerTestCount );
    if (num_selected == num_run) 
      printf( "%2d tests selected and ran\n", num_selected );
    else
      printf( "%2d of %2d selected tests were ran\n", num_run, num_selected );
    if (num_run < (int)RunnerTestCount)
      printf( "%2d of %2d registered tests skipped\n", 
        (int)RunnerTestCount - num_run, (int)RunnerTestCount );
    printf( "%2d of %2d tests passed\n", num_run - fail_count, num_run );
    if (fail_count) 
      printf( "%2d of %2d tests FAILED\n", fail_count, num_run );
  }
  
  return error_count;
}

#endif
