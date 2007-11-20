#ifndef TEST_UTIL_HPP
#define TEST_UTIL_HPP

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
/** Run a test
 *  Argument should be a function with the signature:  void func(void)
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
#define LONGJMP_MODE 3      

// If test application hasn't set MODE, set to default
#ifndef MODE
#  ifdef _MSC_VER
#    define MODE EXCEPTION_MODE
#  else
#    define MODE LONGJMP_MODE
#  endif
#endif


/***************************************************************************************
 * NOTE: The remainder of this file contains the implementation of the above macros.
 *       The above macros constitute the entire intended API.
 ***************************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#if MODE == EXCEPTION_MODE
   struct ErrorExcept{};
#  define FLAG_ERROR throw ErrorExcept()
#elif MODE == FORK_MODE
#  include <sys/types.h>
#  include <sys/wait.h>
#  include <unistd.h>
#  include <errno.h>
#  define FLAG_ERROR exit(1)
#elif MODE == LONGJMP_MODE
#  include <signal.h>
#  include <setjmp.h>
#  ifdef SIGUSR2
#    define FLAG_ERROR raise(SIGUSR2)
#  else
#    define FLAG_ERROR abort()
#  endif
#else
#  error "MODE not set"
#endif

#if MODE == LONGJMP_MODE

jmp_buf jmpenv;
extern "C" {
  void sighandler( int sig ) {
    signal( sig, sighandler );
    longjmp(jmpenv, sig);
    // should never return from longjmp
    exit(1);
  }
  typedef void (*sigfunc_t)(int);
} // extern "C"
int sethandler( int sig ) {
  sigfunc_t h = signal( sig, &sighandler );
  if (h == SIG_ERR)
    return  1;
  else if (h != SIG_DFL)
    signal( sig, h );
  return 0;
}

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

// initialize global to force call to init_signal_handlers
int junk_init_var = init_signal_handlers();

#endif // LONGJMP_MODE



/* Make sure IS_BUILDING_MB is defined so we can include MBInternals.hpp */
#include "MBTypes.h"
#ifndef IS_BUILDING_MB
#  define IS_BUILDING_MB
#  include "MBInternals.hpp"
#  undef IS_BUILDING_MB
#else
#  include "MBInternals.hpp"
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
    
    /* Wait until child process exists */
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
  int sig = setjmp( jmpenv );
  if (!sig) {
    (*test)();
    return 0;
  }
#ifdef SIGUSR2
  else if(sig == SIGUSR2) {
    printf( "  %s: FAILED\n", func_name );
    return 1;
  }
#endif
  else {
    printf( "  %s: TERMINATED (signal %d)\n", func_name, sig );
    return 1;
  }
#else
  #error "MODE not set"
#endif // MODE
}


#define EQUAL_TEST_IMPL( TEST, TYPE ) if( !(TEST) ) { \
  printf( "Equality Test Failed: %s == %s\n", sA, sB ); \
  printf( "  at line %d of '%s'\n", line, file ); \
  printf( "  Expected value: %" #TYPE "\n", A ); \
  printf( "  Actual value:   %" #TYPE "\n", B ); \
  printf( "\n" ); \
  FLAG_ERROR; \
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

void check_equal( float A, float B, float eps, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( fabsf(A - B) <= eps, f ) }

void check_equal( double A, double B, float eps, const char* sA, const char* sB, int line, const char* file )
  {  EQUAL_TEST_IMPL( fabs(A - B) <= eps, f ) }

const char* mb_error_str( MBErrorCode err )
{
  switch (err) {
    case MB_SUCCESS                 : return "Success";
    case MB_INDEX_OUT_OF_RANGE      : return "Index Out of Range";
    case MB_TYPE_OUT_OF_RANGE       : return "Type Out of Range";
    case MB_MEMORY_ALLOCATION_FAILED: return "Memory Alloc. Failed";
    case MB_ENTITY_NOT_FOUND        : return "Entity Not Found";
    case MB_MULTIPLE_ENTITIES_FOUND : return "Multiple Entities Found";
    case MB_TAG_NOT_FOUND           : return "Tag Not Found";
    case MB_FILE_DOES_NOT_EXIST     : return "File Not Found";
    case MB_FILE_WRITE_ERROR        : return "File Write Error";
    case MB_NOT_IMPLEMENTED         : return "Not Implemented";
    case MB_ALREADY_ALLOCATED       : return "Already Allocated";
    case MB_FAILURE                 : return "Failure";
    default                         : return "(unknown)";
  }
}


void check_equal( MBErrorCode A, MBErrorCode B, const char* sA, const char* sB, int line, const char* file )
{
  if (A == B)
    return;
  
  printf( "MBErrorCode Test Failed: %s == %s\n", sA, sB ); 
  printf( "  at line %d of '%s'\n", line, file ); 
  printf( "  Expected value: %s (%d)\n", mb_error_str(A), (int)A ); 
  printf( "  Actual value:   %s (%d)\n", mb_error_str(B), (int)B ); 
  printf( "\n" ); 
  FLAG_ERROR; 
}

const char* mb_type_str( MBEntityType type )
{
  switch(type) {
    case MBVERTEX    : return "Vertex";
    case MBEDGE      : return "Edge";
    case MBTRI       : return "Triangle";
    case MBQUAD      : return "Quadrilateral";
    case MBPOLYGON   : return "Polygon";
    case MBTET       : return "Tetrahedron";
    case MBPYRAMID   : return "Pyramid";
    case MBPRISM     : return "Prism (wedge)";
    case MBKNIFE     : return "Knife";
    case MBHEX       : return "Hexahedron";
    case MBPOLYHEDRON: return "Polyhedron";
    case MBENTITYSET : return "Entity (Mesh) Set";
    case MBMAXTYPE   : return "(max type)";
    default          : return "(unknown)";
  }
}

const char* mb_type_str( MBEntityHandle a )
  { return mb_type_str( TYPE_FROM_HANDLE(a) ); }
/*
void check_equal( MBEntityHandle A, MBEntityHandle B, const char* sA, const char* sB, int line, const char* file )
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
  FLAG_ERROR; 
}  
*/

void check_true( bool cond, const char* str, int line, const char* file )
{
  if( !cond ) { 
    printf( "Test Failed: %s\n", str ); 
    printf( "  at line %d of '%s'\n", line, file ); 
    printf( "\n" ); 
    FLAG_ERROR; 
  }
}

#endif
