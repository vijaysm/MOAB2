/*This function tests the AHF datastructures on CST meshes*/
#include <iostream>
#include <vector>
#include <algorithm>
#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/MeshTopoUtil.hpp"
#include "moab/HalfFacetRep.hpp"
#include "TestUtil.hpp"

using namespace moab;

#ifdef MESHDIR
static const char example[] = STRINGIFY(MESHDIR) "/hexes_mixed.vtk";
#else
static const char example[] = "/hexes_mixed.vtk";
#endif

void ahf_mbintf_test()
{

  Core moab;
  Interface* mbImpl = &moab;
  MeshTopoUtil mtu(mbImpl);

  ErrorCode error = mbImpl->load_file(example);
  CHECK_ERR(error);

  /*Create ranges for handles of explicit elements of the mixed mesh*/
  Range verts, edges, faces, cells;
  error = mbImpl->get_entities_by_dimension( 0, 0, verts);
  error = mbImpl->get_entities_by_dimension( 0, 1, edges);
  error = mbImpl->get_entities_by_dimension( 0, 2, faces);
  error = mbImpl->get_entities_by_dimension( 0, 3, cells);

  //Perform queries
  std::vector<EntityHandle> adjents;
  Range mbents, ahfents;

  //1D Queries //
  //IQ1: For every vertex, obtain incident edges
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, mbents );
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(),mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

  //NQ1:  For every edge, obtain neighbor edges
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 1, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mtu.get_bridge_adjacencies( *i, 0, 1, mbents);
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(), mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

  // 2D Queries
  // IQ21: For every vertex, obtain incident faces
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
      adjents.clear();
      error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
      CHECK_ERR(error);
      mbents.clear();
      error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
      CHECK_ERR(error);

      CHECK_EQUAL(adjents.size(), mbents.size());

      std::sort(adjents.begin(), adjents.end());
      std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
      mbents = subtract(mbents, ahfents);
      CHECK(!mbents.size());
  }

  //IQ22: For every edge, obtain incident faces
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, mbents);
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(), mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

  //NQ2: For every face, obtain neighbor faces
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 2, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mtu.get_bridge_adjacencies( *i, 1, 2, mbents);
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(), mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

  // 3D Queries
  //IQ 31: For every vertex, obtain incident cells
  for (Range::iterator i = verts.begin(); i != verts.end(); ++i) {
      adjents.clear();
      error = mbImpl->get_adjacencies( &*i, 1, 3, false, adjents);
      CHECK_ERR(error);
      mbents.clear();
      error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
      CHECK_ERR(error);

      CHECK_EQUAL(adjents.size(), mbents.size());

      std::sort(adjents.begin(), adjents.end());
      std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
      mbents = subtract(mbents, ahfents);
      CHECK(!mbents.size());
  }

  // IQ 32: For every edge, obtain incident cells
  for (Range::iterator i = edges.begin(); i != edges.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 3, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(), mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

  //IQ32: For every face, obtain incident cells
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 3, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mbImpl->get_adjacencies(&*i, 1, 3, false, mbents);
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(), mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

  //NQ3: For every cell, obtain neighbor cells
  for (Range::iterator i = cells.begin(); i != cells.end(); ++i) {
    adjents.clear();
    error = mbImpl->get_adjacencies( &*i, 1, 3, false, adjents);
    CHECK_ERR(error);
    mbents.clear();
    error = mtu.get_bridge_adjacencies( *i, 2, 3, mbents);
    CHECK_ERR(error);

    CHECK_EQUAL(adjents.size(), mbents.size());

    std::sort(adjents.begin(), adjents.end());
    std::copy(adjents.begin(), adjents.end(), range_inserter(ahfents));
    mbents = subtract(mbents, ahfents);
    CHECK(!mbents.size());
  }

}

int main(int argc, char *argv[])
{
  int result = 0;

  argv[0] = argv[argc - argc]; // Followed read_mpas_nc.cpp test for removing warnings in serial mode about unused variables.

  result += RUN_TEST(ahf_mbintf_test);

  return result;
}

