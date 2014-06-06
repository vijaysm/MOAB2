#include "TestUtil.hpp"
#include "TestRunner.hpp"
#include "moab/ProgOptions.hpp"
#include <limits>
#include <stdlib.h>
#ifdef USE_MPI
# include "moab_mpi.h"
#endif

void test_flag_opt_short();
void test_flag_opt_long_short();
void test_flag_opt_long();
void test_flag_cancel();
void test_flag_store_false();

void test_int_opt();
void test_int_arg();

void test_real_opt();
void test_real_arg();

void test_string_opt();
void test_string_arg();
void test_string_rank_subst();

void test_int_vect_opt();
void test_int_vect_arg();
void test_optional_args();
void test_optional_arg();
void test_squashed_short();

#define ARGCV(A)  (sizeof(A)/sizeof(A[0])), const_cast<char**>(A)

int main( int argc, char* argv[] )
{
  // make ProgOptions abort() rather than exiting with an
  // error code for invalid options so that we can catch
  // the signal and continue with later tests
#ifndef WIN32
  setenv("MOAB_PROG_OPT_ABORT","1",0);
#endif
  
  REGISTER_TEST( test_flag_opt_short );
  REGISTER_TEST( test_flag_opt_long_short );
  REGISTER_TEST( test_flag_opt_long );
  REGISTER_TEST( test_flag_cancel );
  REGISTER_TEST( test_flag_store_false );

  REGISTER_TEST( test_int_opt );
  REGISTER_TEST( test_int_arg );
  REGISTER_TEST( test_real_opt );
  REGISTER_TEST( test_real_arg );

  REGISTER_TEST( test_string_opt );
  REGISTER_TEST( test_string_arg );
  REGISTER_TEST( test_string_rank_subst );

  REGISTER_TEST( test_int_vect_opt );
  REGISTER_TEST( test_int_vect_arg );
  REGISTER_TEST( test_optional_args );
  REGISTER_TEST( test_optional_arg );
  REGISTER_TEST( test_squashed_short );

#ifdef USE_MPI
  MPI_Init( &argc, &argv );
#endif
  int result = RUN_TESTS( argc, argv );
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return result;
}
  
  
void test_flag_opt_short()
{
  ProgOptions opts1;
  bool value1 = false;
  opts1.addOpt<void>( ",s", "short opt", &value1 );
  const char* argv1[] = { "prog", "-s" };
  opts1.parseCommandLine( ARGCV(argv1) );
  CHECK(value1);
  CHECK_EQUAL(1,opts1.numOptSet( ",s" ));
}

void test_flag_opt_long_short()
{
  ProgOptions opts2;
  bool value = false;
  opts2.addOpt<void>( "long,l", "long opt", &value );
  const char* argv2[] = { "prog", "-l", "-l" };
  opts2.parseCommandLine( ARGCV(argv2) );
  CHECK(value);
  CHECK_EQUAL(2,opts2.numOptSet( ",l" ));
  CHECK_EQUAL(2,opts2.numOptSet( "long,l" ));
  CHECK_EQUAL(2,opts2.numOptSet( "long" ));
}

void test_flag_opt_long()
{
  ProgOptions opts3;
  bool value = false;
  opts3.addOpt<void>( "long,l", "long opt", &value );
  const char* argv3[] = { "prog", "--long", "--long" };
  opts3.parseCommandLine( ARGCV(argv3) );
  CHECK(value);
  CHECK_EQUAL(2,opts3.numOptSet( ",l" ));
  CHECK_EQUAL(2,opts3.numOptSet( "long,l" ));
  CHECK_EQUAL(2,opts3.numOptSet( "long" ));
}

void test_flag_cancel() 
{
  ProgOptions opts1;
  bool value = false;
  opts1.addOpt<void>( "flag", "my flag", &value, ProgOptions::add_cancel_opt );
  const char* argv1[] = { "test", "--flag" };
  opts1.parseCommandLine( ARGCV(argv1) );
  CHECK(value);
  
  ProgOptions opts2;
  value = true;
  opts2.addOpt<void>( "flag", "my flag", &value, ProgOptions::add_cancel_opt );
  const char* argv2[] = { "test", "--flag", "--no-flag" };
  opts2.parseCommandLine( ARGCV(argv2) );
  CHECK(!value);
  
  
  ProgOptions opts3;
  value = false;
  opts3.addOpt<void>( "flag", "my flag", &value, ProgOptions::add_cancel_opt );
  const char* argv3[] = { "test", "--flag", "--no-flag", "--flag" };
  opts3.parseCommandLine( ARGCV(argv3) );
  CHECK(value);
}

  
void test_flag_store_false()
{
  ProgOptions opts1;
  bool value1 = true;
  opts1.addOpt<void>( ",s", "short opt", &value1, ProgOptions::store_false );
  const char* argv1[] = { "prog", "-s" };
  opts1.parseCommandLine( ARGCV(argv1) );
  CHECK(!value1);
}


void test_int_opt()
{
  ProgOptions opts;
  int val1 = -1;
  opts.addOpt( "long,s", "my int opt", &val1 );
  const char* argv[] = { "test",
                   "-s","2",
                   "--long","-0xA",
                   "--long=5" };
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_EQUAL( 5, val1 );
  int val2 = -1;
  CHECK(opts.getOpt(",s",&val2));
  CHECK_EQUAL( 5, val2 );
  val2 = -1;
  CHECK(opts.getOpt("long,s",&val2));
  CHECK_EQUAL( 5, val2 );
  val2 = -1;
  CHECK(opts.getOpt("long",&val2));
  CHECK_EQUAL( 5, val2 );
  std::vector<int> list;
  opts.getOptAllArgs( ",s", list );
  CHECK_EQUAL( (size_t)3, list.size() );
  CHECK_EQUAL(   2, list[0] );
  CHECK_EQUAL( -10, list[1] );
  CHECK_EQUAL(   5, list[2] );
}

void test_int_arg()
{
  ProgOptions opts;
  int val1 = 5;
  opts.addRequiredArg( "arg", "my test arg", &val1 );
  opts.addRequiredArg<int>( "arg2", "my other test arg" );
  const char* argv[] = { "test",
                   "--",
                   "-1",
                   "-010" };
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_EQUAL( -1, val1 );
  CHECK_EQUAL( -010, opts.getReqArg<int>("arg2") );  // octal -10 == decimal -8
  CHECK_EQUAL( -1, opts.getReqArg<int>("arg") );
}


void test_real_opt()
{
  const double EPS = std::numeric_limits<double>::epsilon();
  ProgOptions opts;
  double val1 = -1;
  opts.addOpt( "long,s", "my real opt", &val1 );
  const char* argv[] = { "test",
                   "-s","2",
                   "--long","2e5",
                   "--long=-0.01" };
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_REAL_EQUAL( -0.01, val1, EPS );
  double val2 = -1;
  CHECK(opts.getOpt(",s",&val2));
  CHECK_REAL_EQUAL( -0.01, val2, EPS );
  val2 = -1;
  CHECK(opts.getOpt("long,s",&val2));
  CHECK_REAL_EQUAL( -0.01, val2, EPS );
  val2 = -1;
  CHECK(opts.getOpt("long",&val2));
  CHECK_REAL_EQUAL( -0.01, val2, EPS );
  std::vector<double> list;
  opts.getOptAllArgs( ",s", list );
  CHECK_EQUAL( (size_t)3, list.size() );
  CHECK_REAL_EQUAL(  2   , list[0], EPS );
  CHECK_REAL_EQUAL(  2e5 , list[1], EPS );
  CHECK_REAL_EQUAL( -0.01, list[2], EPS );
}

void test_real_arg()
{
  const double EPS = std::numeric_limits<double>::epsilon();
  ProgOptions opts;
  double val1 = 5;
  opts.addRequiredArg( "arg", "my test arg", &val1 );
  opts.addRequiredArg<double>( "arg2", "my other test arg" );
  const char* argv[] = { "test",
                   "--",
                   "-1.2",
                   "1.01e-3" };
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_REAL_EQUAL( -1.2, val1, EPS );
  CHECK_REAL_EQUAL( 1.01e-3, opts.getReqArg<double>("arg2"), EPS );
  CHECK_REAL_EQUAL( -1.2, opts.getReqArg<double>("arg"), EPS );
}

void test_string_opt( )
{
  ProgOptions opts;
  std::string val1;
  opts.addOpt( "long,s", "my first opt", &val1 );
  opts.addOpt<std::string>( "second,2", "my second opt" );

  const char* argv[] = { "test",
                   "--long","2",
                   "-s","foobar",
                   "-2","two",
                   "--second=testval" };
  opts.parseCommandLine( ARGCV(argv) );

  CHECK_EQUAL( "foobar", val1 );
  std::string val2;
  CHECK(opts.getOpt(",s",&val2));
  CHECK_EQUAL( "foobar", val2 );
  val2.clear();
  CHECK(opts.getOpt("long,s",&val2));
  CHECK_EQUAL( "foobar", val2 );
  val2.clear();
  CHECK(opts.getOpt("long",&val2));
  CHECK_EQUAL( "foobar", val2 );
  
  val2.clear();
  CHECK(opts.getOpt(",2",&val2));
  CHECK_EQUAL( "testval", val2 );
  val2.clear();
  CHECK(opts.getOpt("second,2",&val2));
  CHECK_EQUAL( "testval", val2 );
  val2.clear();
  CHECK(opts.getOpt("second",&val2));
  CHECK_EQUAL( "testval", val2 );
  
  std::vector<std::string> list;
  opts.getOptAllArgs( "long", list );
  CHECK_EQUAL( (size_t)2, list.size() );
  CHECK_EQUAL( "2", list[0] );
  CHECK_EQUAL( "foobar", list[1] );
  
  list.clear();
  opts.getOptAllArgs( ",2", list );
  CHECK_EQUAL( (size_t)2, list.size() );
  CHECK_EQUAL( "two", list[0] );
  CHECK_EQUAL( "testval", list[1] );
}

void test_string_arg( )
{
  ProgOptions opts;
  std::string val1, val2;
  opts.addRequiredArg<std::string>( "arg", "my test arg" );
  opts.addRequiredArg( "arg2", "my other test arg",&val2 );
  const char* argv[] = { "test",
                   "my_string",
                   "with spaces" };
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_EQUAL( "with spaces", val2 );
  CHECK_EQUAL( "my_string", opts.getReqArg<std::string>("arg") );
  CHECK_EQUAL( "with spaces", opts.getReqArg<std::string>("arg2") );
}

void test_string_rank_subst( )
{
  std::string exp1 = "ddd%";
  std::string exp2 = "%";
  std::string exp2b = "foo%bar";
  std::string exp3 = "%string%";
  std::string exp4 = "file.%";

  ProgOptions opts;
  std::string val1, val2, val3, val4;
  opts.addOpt( ",n", "no subst flag", &val1, 0 );
  opts.addOpt( "dosub,s",    "subst flag", &val2, ProgOptions::rank_subst );
  opts.addRequiredArg( "nos", "no subst arg", &val3, 0 );
  opts.addRequiredArg( "sub",    "subst arg", &val4, ProgOptions::rank_subst );
  std::string eqflg("--dosub="); eqflg += exp2;
  const char* argv[] = { "test",
                   exp3.c_str(),
                   exp4.c_str(),
                   "-s",exp2b.c_str(),
                   "-n",exp1.c_str(),
                   eqflg.c_str() };
  opts.parseCommandLine( ARGCV(argv) );
  
#ifdef USE_MPI
  int rank, size;
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  char buffer[64];
  int width = 1;
  while (size > 10) {
    width++;
    size/=10;
  }
  sprintf(buffer,"%0*d",width,rank);
  exp2 = buffer;
  exp2b = std::string("foo") + buffer + "bar";
  exp4 = std::string("file.") + buffer;
#endif
  
  CHECK_EQUAL( exp1, val1 );
  CHECK_EQUAL( exp2, val2 );
  CHECK_EQUAL( exp3, val3 );
  CHECK_EQUAL( exp4, val4 );
  
  val1.clear();
  val2.clear();
  CHECK( opts.getOpt(",n",    &val1 ) );
  CHECK( opts.getOpt("dosub", &val2 ) );
  
  CHECK_EQUAL( exp1, val1 );
  CHECK_EQUAL( exp2, val2 );
  CHECK_EQUAL( exp3, opts.getReqArg<std::string>("nos") );
  CHECK_EQUAL( exp4, opts.getReqArg<std::string>("sub") );
  
  std::vector<std::string> list;
  opts.getOptAllArgs( ",s", list );
  CHECK_EQUAL( (size_t)2, list.size() );
  CHECK_EQUAL( exp2b, list[0] );
  CHECK_EQUAL( exp2,  list[1] );
}

void test_int_vect_opt()
{
  ProgOptions opts;
  std::vector<int> list1, list2;
  opts.addOpt( "ids,d", "id list", &list1 );
  const char* argv[] = { "test",
                   "--ids=1,2,3",
                   "--ids","4-10",
                   "-d","4-5,2" };
  
  const int exp1[] = { 1, 2, 3 };
  const int exp1_len = sizeof(exp1)/sizeof(exp1[0]);
  const int exp2[] = { 4, 5, 6, 7, 8, 9, 10 };
  const int exp2_len = sizeof(exp2)/sizeof(exp2[0]);
  const int exp3[] = { 4, 5, 2 };
  const int exp3_len = sizeof(exp3)/sizeof(exp3[0]);
  std::vector<int> all;
  std::copy( exp1, exp1+exp1_len, std::back_inserter(all) );
  std::copy( exp2, exp2+exp2_len, std::back_inserter(all) );
  std::copy( exp3, exp3+exp3_len, std::back_inserter(all) );
 
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_EQUAL( all, list1 );
  CHECK(opts.getOpt(",d", &list2 ));
  CHECK_ARRAYS_EQUAL( exp3, exp3_len, &list2[0], list2.size() );
  
  std::vector< std::vector<int> > lists;
  opts.getOptAllArgs( "ids", lists );
  CHECK_EQUAL( (size_t)3, lists.size() );
  CHECK_ARRAYS_EQUAL( exp1, exp1_len, &lists[0][0], lists[0].size() );
  CHECK_ARRAYS_EQUAL( exp2, exp2_len, &lists[1][0], lists[1].size() );
  CHECK_ARRAYS_EQUAL( exp3, exp3_len, &lists[2][0], lists[2].size() );
  
  list2.clear();
  opts.getOptAllArgs( "ids", list2 );
  CHECK_EQUAL( all, list2 );
}

void test_int_vect_arg()
{
  ProgOptions opts;
  std::vector<int> list1, list2;
  opts.addRequiredArg( "ints", "int list", &list1 );
  const char* argv[] = { "test", "5,6,-3--1,10" };
  
  const int exp1[] = { 5, 6, -3, -2, -1, 10 };
  const int exp1_len = sizeof(exp1)/sizeof(exp1[0]);
  
  opts.parseCommandLine( ARGCV(argv) );
  CHECK_ARRAYS_EQUAL( exp1, exp1_len, &list1[0], list1.size() );
  CHECK_EQUAL(list1, opts.getReqArg< std::vector<int> >("ints"));
}

void test_optional_args()
{
  ProgOptions opts1;
  std::string arg;
  opts1.addOptionalArgs<std::string>( 0, "opts", "optional arguments" );
  opts1.addRequiredArg( "req", "required final argument", &arg );
  const char* argv1[] = { "test", "arg" }; 
  opts1.parseCommandLine( ARGCV(argv1) );
  std::vector<std::string> list;
  CHECK_EQUAL( "arg", arg );
  opts1.getArgs( "opts", list );
  CHECK(list.empty());
  
  ProgOptions opts2;
  arg.clear();
  opts2.addOptionalArgs<std::string>( 0, "opts", "optional arguments" );
  opts2.addRequiredArg( "req", "required final argument", &arg );
  const char* argv2[] = { "test", "arg1", "arg2", "arg" }; 
  opts2.parseCommandLine( ARGCV(argv2) );
  CHECK_EQUAL( "arg", arg );
  list.clear();
  opts2.getArgs( "opts", list );
  CHECK_EQUAL((size_t)2, list.size());
  CHECK_EQUAL("arg1", list[0] );
  CHECK_EQUAL("arg2", list[1] );
}

void test_optional_arg()
{
  ProgOptions opts1;
  std::string init, fini;
  opts1.addRequiredArg( "init", "required initial argument", &init );
  opts1.addOptionalArgs<std::string>( 1, "mid", "optional arguments" );
  opts1.addRequiredArg( "fini", "required final argument", &fini );
  const char* argv1[] = { "test", "arg1", "arg2" }; 
  opts1.parseCommandLine( ARGCV(argv1) );
  std::vector<std::string> list;
  CHECK_EQUAL( "arg1", init );
  CHECK_EQUAL( "arg2", fini );
  opts1.getArgs( "mid", list );
  CHECK(list.empty());
  
  ProgOptions opts2;
  init.clear();
  fini.clear();
  opts2.addRequiredArg( "init", "required initial argument", &init );
  opts2.addOptionalArgs<std::string>( 1, "mid", "optional arguments" );
  opts2.addRequiredArg( "fini", "required final argument", &fini );
  const char* argv2[] = { "test", "arg1", "arg2", "arg3" }; 
  opts2.parseCommandLine( ARGCV(argv2) );
  CHECK_EQUAL( "arg1", init );
  CHECK_EQUAL( "arg3", fini );
  list.clear();
  opts2.getArgs( "mid", list );
  CHECK_EQUAL((size_t)1, list.size());
  CHECK_EQUAL("arg2", list[0] );
}

void test_squashed_short()
{
  ProgOptions opts;
  int intval = 0;
  double realval = 0;
  bool flagval = false;
  std::string strval;
  opts.addOpt<void>( ",f", "flag", &flagval );
  opts.addOpt( ",i", "int",  & intval );
  opts.addOpt( ",r", "real", &realval );
  opts.addOpt( ",s", "str",  & strval );
  
  const char* argv[] = { "test",
                   "-ifsrff",
                   "-0xBEEF",
                   "string",
                   "-1.0e55" };
  opts.parseCommandLine( ARGCV(argv) );
  
  CHECK_EQUAL( -0xBEEF, intval );
  CHECK_REAL_EQUAL( -1.0e55, realval, 1e-6 );
  CHECK_EQUAL( "string", strval );
  CHECK( flagval );
  CHECK_EQUAL( 3, opts.numOptSet( ",f" ) );
  CHECK_EQUAL( 1, opts.numOptSet( ",i" ) );
  CHECK_EQUAL( 1, opts.numOptSet( ",r" ) );
  CHECK_EQUAL( 1, opts.numOptSet( ",s" ) );
}


