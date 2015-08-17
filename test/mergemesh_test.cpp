#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MergeMesh.hpp"
#include <iostream>
#include "TestUtil.hpp"

#ifndef MESHDIR
#error Specify MESHDIR to compile test
#endif

#ifdef MOAB_HAVE_MPI
#include "moab_mpi.h"
#endif

using namespace moab;

const char* meshfile = STRINGIFY(MESHDIR) "/16_unmerged_hex.h5m";
const char* meshfile2 = STRINGIFY(MESHDIR) "/merge_with_tag.h5m";
const char *outfile = "mm_out.h5m";

void mergesimple_test();
void merge_with_tag_test();

#ifdef MOAB_HAVE_MPI
int main(int argc, char** argv)
#else
int main()
#endif
{
#ifdef MOAB_HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  int result = 0;

  result += RUN_TEST(mergesimple_test);
  result += RUN_TEST(merge_with_tag_test);

#ifdef MOAB_HAVE_MPI
  MPI_Finalize();
#endif

  return result;
}

void mergesimple_test()
{
  ErrorCode rval;
  Interface* iface = new Core();
  // can be generalized to load user defined input/output file

  rval = iface->load_mesh(meshfile);
  CHECK_ERR(rval);
  int dim = 3;
  moab::Range ents;
  iface->get_entities_by_dimension(0, dim, ents);

  // Make sure that mm is destroyed before deleting iface
  {
    MergeMesh mm(iface);
    double merge_tol = 1e-3;

    rval = mm.merge_entities(ents, merge_tol);
    CHECK_ERR(rval);
  }

  // Fixed for now

  rval = iface->write_file( outfile);
  CHECK_ERR(rval);

  delete iface;

  return;
}

void merge_with_tag_test()
{
  ErrorCode rval;
  Interface* iface = new Core();
  // can be generalized to load user defined input/output file

  rval = iface->load_mesh(meshfile2);
  CHECK_ERR(rval);
  int dim = 0;
  moab::Range verts;
  iface->get_entities_by_dimension(0, dim, verts);
  Tag  tag_for_merge;
  rval = iface->tag_get_handle("IDFTAG", tag_for_merge);
  CHECK_ERR(rval);

  // Make sure that mm is destroyed before deleting iface
  {
    MergeMesh mm(iface);
    rval = mm.merge_using_integer_tag(verts, tag_for_merge);
    CHECK_ERR(rval);
  }

  rval = iface->write_file( outfile);
  CHECK_ERR(rval);

  verts.clear();
  iface->get_entities_by_dimension(0, dim, verts);
  CHECK_EQUAL( 405, (int)verts.size()) ;

  delete iface;

  return;
}
