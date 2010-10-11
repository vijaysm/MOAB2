#include "TestUtil.hpp"
#include "ElemUtil.hpp"
#include <iostream>

using namespace moab;

void test_tet();
void test_hex();

int main()
{
  int rval = 0;
  rval += RUN_TEST(test_tet);
  rval += RUN_TEST(test_hex);
  return rval;
}

void test_tet() {
  moab::Element::LinearTet tet;
}// test_tet()

void test_hex() {
  moab::Element::LinearHex hex;
}// test_hex()
