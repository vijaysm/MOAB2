#include "moab/Range.hpp"
#include "Internals.hpp"
#include "TestUtil.hpp"

using namespace moab;

void insert_hint_test();
void intersect_test();
void merge_test();
void const_pair_iterator_test();
void swap_test();
void subtract_test();
void subset_by_dimension_test();
void erase_test();
void contains_test();

int main()
{
  int rval = 0;
  rval += RUN_TEST(insert_hint_test);
  rval += RUN_TEST(intersect_test);
  rval += RUN_TEST(merge_test);
  rval += RUN_TEST(const_pair_iterator_test);
  rval += RUN_TEST(swap_test);
  rval += RUN_TEST(subtract_test);
  rval += RUN_TEST(subset_by_dimension_test);
  rval += RUN_TEST(erase_test);
  rval += RUN_TEST(contains_test);
  return rval;
}

void insert_hint_test()
{
  const EntityHandle pairs[][2] = {
    { 4980, 4981 },
    { 4985, 4986 },
    { 4990, 4990 },
    { 5886, 5886 },
    { 5888, 5890 },
    { 5890, 5890 },
    { 5886, 5888 },
    { 5890, 5890 },
    { 5894, 5894 },
    { 5899, 5899 } 
  };
  const int num_pairs = sizeof(pairs)/sizeof(pairs[0]);
  const EntityHandle exp_pairs[][2] = {
    { 4980, 4981 },
    { 4985, 4986 },
    { 4990, 4990 },
    { 5886, 5890 },
    { 5894, 5894 },
    { 5899, 5899 } 
  };
  const int num_exp = sizeof(exp_pairs)/sizeof(exp_pairs[0]);
  
  Range range;
  Range::iterator it = range.begin();
  for (int i = 0; i < num_pairs; ++i) 
    it = range.insert( it, pairs[i][0], pairs[i][1] );
  
  Range::const_pair_iterator pit = range.const_pair_begin();
  for (int i = 0; i < num_exp; ++i) {
    CHECK( pit != range.const_pair_end() );
    CHECK_EQUAL( exp_pairs[i][0], pit->first );
    CHECK_EQUAL( exp_pairs[i][1], pit->second );
    ++pit;
  }
  CHECK( pit == range.const_pair_end() );
}

// common constants used for a bunch of tests below
const EntityHandle h1 = CREATE_HANDLE(MBVERTEX, 1);
const EntityHandle h4 = CREATE_HANDLE(MBVERTEX, 4);
const EntityHandle h5 = CREATE_HANDLE(MBVERTEX, 5);
const EntityHandle h9 = CREATE_HANDLE(MBVERTEX, 9);
const EntityHandle h10 = CREATE_HANDLE(MBVERTEX, 10);
const EntityHandle h15 = CREATE_HANDLE(MBVERTEX, 15);
const EntityHandle h16 = CREATE_HANDLE(MBVERTEX, 16);
const EntityHandle h20 = CREATE_HANDLE(MBVERTEX, 20);
const EntityHandle hh1 = CREATE_HANDLE(MBHEX, 1);
const EntityHandle hh3 = CREATE_HANDLE(MBHEX, 3);

void intersect_test()
{
  Range r1, r2, rhs;

    // equal start/end
  r1.insert(h1, h5);
  r1.insert(h10, h15);
  r2.insert(h5, h10);
  r2.insert(h15, h20);
  rhs = intersect( r1, r2);
  CHECK_EQUAL( (size_t)3, rhs.size() );
  
    // identical ranges test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h5); r1.insert(h10, h20); 
  r2.insert(h1, h5); r2.insert(h10, h20); 
  rhs = intersect( r1, r2);
  CHECK_EQUAL(r1.size(), rhs.size());
  CHECK_EQUAL(r2.size(), rhs.size());
  
    // off by one test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h4); r1.insert(h10, h15); 
  r2.insert(h5, h9); r2.insert(h16, h20);
  rhs = intersect( r1, r2);
  CHECK( rhs.empty() );
  
    // interior test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h20); 
  r2.insert(h5, h10); 
  rhs = intersect( r1, r2);
  CHECK_EQUAL(r2.size(), rhs.size());
  
    // half-above test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h1, h10);
  r2.insert(h5, h20);
  rhs = intersect( r1, r2);
  CHECK_EQUAL( (size_t)(h10-h5+1), rhs.size() );
  
    // half-below test
  r1.clear(); r2.clear(); rhs.clear();
  r1.insert(h5, h20);
  r2.insert(h1, h10);
  rhs = intersect( r1, r2);
  CHECK_EQUAL( (size_t)(h10-h5+1), rhs.size() );
}

void merge_test()
{
  Range r1, r2;
   
    // merge all test
  r1.clear();
  r2.clear();
  r1.insert(h1, h5);
  r1.insert(h10, h15);
  r2 = r1;
  r1.merge( r2 );
  CHECK_EQUAL( r2.size(), r1.size() );
  
    // merge subset test
  r1.clear();
  r2.clear();
  r1.insert(h1, h5);
  r1.insert(h10, h15);
  Range::const_iterator i1 = r1.begin();
  Range::const_iterator i2 = r1.end();
  i1 += 2;
  i2 -= 2;
  r2.merge( i1, i2 );
  CHECK_EQUAL( r1.size() - 4, r2.size() );
}

void const_pair_iterator_test()
{  
    // const_pair_iterator test
  Range r1;
  r1.insert( h1 );
  r1.insert( h4 );
  r1.insert( h5 );
  Range::const_pair_iterator pair_iter = r1.const_pair_begin();
  EntityHandle cpi_h1 = (*pair_iter).first;
  EntityHandle cpi_h2 = (*pair_iter).second;
  ++pair_iter;
  EntityHandle cpi_h3 = (*pair_iter).first;
  EntityHandle cpi_h4 = (*pair_iter).second;
  
  ++pair_iter;
  CHECK_EQUAL( cpi_h1, h1 );
  CHECK_EQUAL( cpi_h2, h1 );
  CHECK_EQUAL( cpi_h3, h4 );
  CHECK_EQUAL( cpi_h4, h5 );
}

void swap_test()
{
  Range r1, r2;

    // both non-empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r2.insert(h9);
  r1.swap(r2);
  CHECK_EQUAL( (size_t)1, r1.size() );
  CHECK_EQUAL( (size_t)3, r2.size() );
  CHECK( r1.find(h9) != r1.end() );
  CHECK( r2.find(h1) != r2.end() );
  CHECK( r2.find(h4) != r2.end() );
  CHECK( r2.find(h5) != r2.end() );
  
    // destination empty
  r1.clear();
  r2.clear();
  r2.insert(h1);
  r2.insert(h4);
  r2.insert(h5);
  r1.swap(r2);
  CHECK_EQUAL( (size_t)3, r1.size() );
  CHECK_EQUAL( (size_t)0, r2.size() );
  CHECK( r1.find(h1) != r1.end() );
  CHECK( r1.find(h4) != r1.end() );
  CHECK( r1.find(h5) != r1.end() );
  
    // source empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r1.swap(r2);
  CHECK_EQUAL( (size_t)0, r1.size() );
  CHECK_EQUAL( (size_t)3, r2.size() );
  CHECK( r2.find(h1) != r2.end() );
  CHECK( r2.find(h4) != r2.end() );
  CHECK( r2.find(h5) != r2.end() );
  
    // both empty
  r1.clear();
  r2.clear();
  r1.swap(r2);
  CHECK( r1.empty() );
  CHECK( r2.empty() );
}

void subtract_test()
{
  Range r1, r2, r3;
    
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r2.insert(h4);
  r3 = subtract( r1, r2);
  CHECK_EQUAL( (size_t)2, r3.size() );
  CHECK( r3.find(h1) != r3.end() );
  CHECK( r3.find(h5) != r3.end() );
  CHECK( r3.find(h4) == r3.end() );
  
    // destination empty
  r1.clear();
  r2.clear();
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(h5);
  r3 = subtract( r1, r2);
  CHECK_EQUAL( (size_t)3, r3.size() );
  CHECK( r3.find(h1) != r3.end() );
  CHECK( r3.find(h4) != r3.end() );
  CHECK( r3.find(h5) != r3.end() );
}

void subset_by_dimension_test()
{
  Range r1, r2;

    // subset_by_dimension test
  r1.insert(h1);
  r1.insert(h4);
  r1.insert(hh1);
  r1.insert(hh3);
  r2 = r1.subset_by_dimension(3);
  CHECK_EQUAL( (size_t)2, r2.size() );
  CHECK_EQUAL( hh1, r2.front() );
  CHECK_EQUAL( hh3, r2.back() );
}

void erase_test() 
{
  Range range;
  Range::iterator result;
  
    // test erase from first node
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin(), range.begin()+2 );
  CHECK_EQUAL( (EntityHandle)7, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)13, range.size() );
  CHECK_EQUAL( (EntityHandle)7, *result );
  
    // test erase first node
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin(), range.begin()+6 );
  CHECK_EQUAL( (EntityHandle)12, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)9, range.size() );
  CHECK_EQUAL( (EntityHandle)12, *result );
  
    // test erase from back of first node
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+2, range.begin()+6 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)11, range.size() );
  CHECK_EQUAL( (EntityHandle)12, *result );
  
    // test erase from middle of first node
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+2, range.begin()+5 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)12, range.size() );
  CHECK_EQUAL( (EntityHandle)10, *result );
  
    // test erase spanning two nodes
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+3, range.begin()+7 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)11, range.size() );
  CHECK_EQUAL( (EntityHandle)13, *result );
  
    // test erase of first node and part of second
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin(), range.begin()+7 );
  CHECK_EQUAL( (EntityHandle)13, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)8, range.size() );
  CHECK_EQUAL( (EntityHandle)13, *result );
  
    // test erase spanning three nodes
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  range.insert( 100, 101 );
  result = range.erase( range.begin()+3, range.begin()+16 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)101, range.back() );
  CHECK_EQUAL( (size_t)4, range.size() );
  CHECK_EQUAL( (EntityHandle)101, *result );
  
    // test erase from start of second node
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+6, range.begin()+8 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)13, range.size() );
  CHECK_EQUAL( (EntityHandle)14, *result );
  
    // test erase from back of last node
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+13, range.end() );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)18, range.back() );
  CHECK_EQUAL( (size_t)13, range.size() );
  CHECK( result == range.end() );
  
    // test erase part of first node through end
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+4, range.end() );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)8, range.back() );
  CHECK_EQUAL( (size_t)4, range.size() );
  CHECK( result == range.end() );
  
    // test erase of single node
  range.clear();
  range.insert( 5, 10 );
  result = range.erase( range.begin(), range.end() );
  CHECK_EQUAL( (size_t)0, range.size() );
  CHECK( result == range.end() );
  
    // test erase of multi-node range
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  range.insert( 100, 101 );
  result = range.erase( range.begin(), range.end() );
  CHECK_EQUAL( (size_t)0, range.size() );
  CHECK( result == range.end() );
  
    // test erase nothing
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  result = range.erase( range.begin()+3, range.begin()+3 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)15, range.size() );
  CHECK_EQUAL( (EntityHandle)8, *result );
  
    // test iterators before erase remain valid
  Range::iterator a, b, c;
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  a = range.begin();
  b = range.begin() + 6;  
  c = range.begin() + 8;
  result = range.erase( range.begin() + 9, range.end() );
  CHECK( a == range.begin() );
  CHECK( b == range.begin() + 6 );
  CHECK( c == range.begin() + 8 );
  CHECK( result == range.end() );
  
    // test iterators before erase remain valid, single value case
  range.clear();
  range.insert( 5, 10 );
  range.insert( 12, 20 );
  a = range.begin();
  b = range.begin() + 6;  
  c = range.begin() + 8;
  result = range.erase( range.begin() + 9 );
  CHECK_EQUAL( (EntityHandle)5, range.front() );
  CHECK_EQUAL( (EntityHandle)20, range.back() );
  CHECK_EQUAL( (size_t)14, range.size() );
  CHECK( a == range.begin() );
  CHECK( b == range.begin() + 6 );
  CHECK( c == range.begin() + 8 );
  CHECK_EQUAL( (EntityHandle)16, *result );
}

void contains_test() 
{
  Range r1, r2;

    // simple test cases: one range each
    
  r1.clear();
  r2.clear();
  r1.insert( 1, 20 );
  r2.insert( 1, 21 );
  CHECK( !r1.contains(r2) );
  CHECK(  r2.contains(r1) );
  
  r1.clear();
  r2.clear();
  r1.insert( 2, 20 );
  r2.insert( 1, 20 );
  CHECK( !r1.contains(r2) );
  CHECK(  r2.contains(r1) );
  
  r1.clear();
  r2.clear();
  r1.insert( 5 );
  r2.insert( 1, 6 );
  CHECK( !r1.contains(r2) );
  CHECK(  r2.contains(r1) );
  
  r1.clear();
  r2.clear();
  r1.insert( 5 );
  r2.insert( 6 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r1.clear();
  r2.clear();
  r1.insert( 18 );
  r2.insert( 18 );
  CHECK( r1.contains(r2) );
  CHECK( r2.contains(r1) );
  
    // empty range test cases
  
  r1.clear();
  r2.clear();
  CHECK( r1.contains(r2) );
  CHECK( r2.contains(r1) );

  r1.clear();
  r2.clear();
  r1.insert( 18 );
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
    // slightly more complex tests: one range in one container
    // and multiple in the other

  r1.clear();
  r1.insert( 10, 100 );
  r2.clear();
  r2.insert( 20, 30 );
  r2.insert( 40, 50 );
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r2.insert( 10, 12 );
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );

  r2.insert( 90, 100 );  
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );

  r2.insert( 9 );  
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );

  r2.erase(9);
  r2.insert( 101 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r2.erase( 101 );
  r2.insert( 103, 110 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r2.insert( 1, 5 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
    // most complex case: both containers have several ranges

  r1.clear();
  r1.insert( 10, 30 );
  r1.insert( 40, 50 );
  r1.insert( 90, 100 );
  r2.clear();
  r2.insert( 20, 30 );
  r2.insert( 40, 50 );
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r2.insert( 10, 12 );
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );

  r2.insert( 90, 100 );  
  CHECK(  r1.contains(r2) );
  CHECK( !r2.contains(r1) );

  r2.insert( 9 );  
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );

  r2.erase(9);
  r2.insert( 101 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r2.erase( 101 );
  r2.insert( 103, 110 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
  
  r2.insert( 1, 5 );
  CHECK( !r1.contains(r2) );
  CHECK( !r2.contains(r1) );
}
