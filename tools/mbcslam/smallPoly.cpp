/*
 * smallPoly.cpp
 *
 */

#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "TestUtil.hpp"

using namespace moab;

int main(int argc, char **argv)
{

  // read an MPAS file, extract few polys and write them out
  std::string fileN="mpas_p8.h5m";
  const char *filename_mesh1 = fileN.c_str();
  Core moab;
  Interface & mb = moab;
  EntityHandle euler_set;
  ErrorCode rval;
  rval = mb.create_meshset(MESHSET_SET, euler_set);
  CHECK_ERR(rval);

  rval = mb.load_file(filename_mesh1);
  CHECK_ERR(rval);

  // get one / two vertices , and polys adjacent to them;
  Range verts;
  rval = mb.get_entities_by_type(0, MBVERTEX, verts);
  CHECK_ERR(rval);

  EntityHandle v[3]={verts[0], verts[1], verts[2]};
  Range  elems;
  rval= mb.get_adjacencies(v, 3, 2, false, elems, Interface::UNION);
  CHECK_ERR(rval);

  rval = mb.add_entities(euler_set, elems);
  CHECK_ERR(rval);

  rval = mb.write_mesh("fewPolys.h5m", &euler_set, 1);
  CHECK_ERR(rval);
    return 0;

}


