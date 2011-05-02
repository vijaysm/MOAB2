#include "TestUtil.hpp"
#include "VarLenTag.hpp"

using namespace moab;

#include <iostream>

void test_valid_struct();
void test_inline();
void test_non_inline();
void test_resize_ii();
void test_resize_in();
void test_resize_ni();
void test_resize_nn();

int main()
{
  int count = 0;
  count += RUN_TEST( test_valid_struct );
  if (count) {
    std::cerr << "ABORTING VarLenTag TEST" << std::endl
              << "Structure is not valid" << std::endl;
    return count;
  }
  
  count += RUN_TEST( test_inline );
  count += RUN_TEST( test_non_inline );
  count += RUN_TEST( test_resize_ii );
  count += RUN_TEST( test_resize_in );
  count += RUN_TEST( test_resize_ni );
  count += RUN_TEST( test_resize_nn );
  
  return count;
}
  

#define OFFSET( A ) ((char*)(&(A)) - (char*)this)
class GetOffsets : public VarLenTag
{
public:
  unsigned pointer_array_offset() { return OFFSET(mData.mData.mPointer.array); }
  unsigned pointer_size_offset() { return OFFSET(mData.mData.mPointer.size); }
#ifdef VAR_LEN_TAG_ELIDE_DATA
  unsigned inline_array_offset() { return OFFSET(mData.mData.mInline.array); }
  unsigned inline_size_offset() { return OFFSET(mData.mData.mInline.size); }
#endif
};
struct ExpectedSize
{
  unsigned char* pointer;
  unsigned size;
};


void test_valid_struct()
{
  GetOffsets off;
  CHECK_EQUAL( 0u, off.pointer_array_offset() );
#ifdef VAR_LEN_TAG_ELIDE_DATA
  CHECK_EQUAL( 0u, off.inline_array_offset() );
  CHECK_EQUAL( off.pointer_size_offset(), off.inline_size_offset() );
#endif
  CHECK_EQUAL( sizeof(ExpectedSize), sizeof(VarLenTag) );
}

void test_inline()
{
  VarLenTag tag( sizeof(void*) );
  CHECK_EQUAL( (unsigned char*)&tag, tag.data() );
}

void test_non_inline()
{
  VarLenTag tag( 2*sizeof(void*) );
  CHECK( (unsigned char*)&tag != tag.data() );
}

void test_resize_ii()
{
  VarLenTag tag( 1 );
  tag.data()[0] = 'X';
  unsigned char* ptr = tag.resize( 3 );
  CHECK_EQUAL( tag.data(), ptr );
  CHECK_EQUAL( (unsigned char*)&tag, tag.data() );
  CHECK_EQUAL( tag.data()[0], 'X' );
}

void test_resize_in()
{
  VarLenTag tag( sizeof(void*) );
  memcpy( tag.data(), "ABCDEFGHIJKLMNOPQRST", sizeof(void*) );
  unsigned char* ptr = tag.resize( 2*sizeof(void*) );
  CHECK_EQUAL( tag.data(), ptr );
  CHECK( (unsigned char*)&tag != tag.data() );
  CHECK( !memcmp( tag.data(), "ABCDEFGHIJKLMNOPQRST", sizeof(void*) ) );
}

void test_resize_ni()
{
  VarLenTag tag( 2*sizeof(void*) );
  memcpy( tag.data(), "12345678901234567890", sizeof(void*) );
  unsigned char* ptr = tag.resize( sizeof(void*) );
  CHECK_EQUAL( tag.data(), ptr );
  CHECK_EQUAL( (unsigned char*)&tag, tag.data() );
  CHECK( !memcmp( tag.data(), "12345678901234567890", sizeof(void*) ) );
}
  
void test_resize_nn()
{
  VarLenTag tag( 2*sizeof(void*) );
  memcpy( tag.data(), "TSRQPONMLKJIHGFEDCBA", 2*sizeof(void*) );
  unsigned char* ptr = tag.resize( 4*sizeof(void*) );
  CHECK_EQUAL( tag.data(), ptr );
  CHECK( (unsigned char*)&tag != tag.data() );
  CHECK( !memcmp( tag.data(), "TSRQPONMLKJIHGFEDCBA", sizeof(void*) ) );
}
 

