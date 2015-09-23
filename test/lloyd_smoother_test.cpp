#include "moab/Core.hpp"
#include "moab/LloydSmoother.hpp"
#include "moab/CartVect.hpp"
#include "TestUtil.hpp"

std::string filename = TestDir + "/surfrandomtris-4part.h5m";

using namespace moab;

int main(int argc, char**argv) 
{
  if (argc > 1) filename = std::string(argv[1]);
  Core mb;
  ErrorCode rval = mb.load_file(filename.c_str());
  CHECK_ERR(rval);
  
  Range elems;
  rval = mb.get_entities_by_dimension(0, 3, elems);
  CHECK_ERR(rval);
  if (elems.empty()) {
    rval = mb.get_entities_by_dimension(0, 2, elems);
    CHECK_ERR(rval);
  }
  if (elems.empty()) {
    std::cout << "Mesh must have faces or regions for this test." << std::endl;
    CHECK_ERR(MB_FAILURE);
  }

    // get the vertex positions and set on an intermediate tag, to test smoothing for
    // a tag instead of for coords
  Tag ctag;
  rval = mb.tag_get_handle("vcentroid", 3, MB_TYPE_DOUBLE, ctag, MB_TAG_CREAT|MB_TAG_DENSE);
  CHECK_ERR(rval);
  Range verts;
  rval = mb.get_entities_by_dimension(0, 0, verts);
  CHECK_ERR(rval);
  std::vector<double> coords(3*verts.size());
  rval = mb.get_coords(verts, &coords[0]);
  CHECK_ERR(rval);
  rval = mb.tag_set_data(ctag, verts, &coords[0]);
  CHECK_ERR(rval);

  LloydSmoother ll(&mb, NULL, elems, ctag);
  ll.report_its(10);
  rval = ll.perform_smooth();
  CHECK_ERR(rval);
  std::cout << "Mesh smoothed in " << ll.num_its() << " iterations." << std::endl;

    // now, set vertex coords to almost their converged positions, then re-smooth; should take fewer
    // iterations
  std::vector<double> new_coords(3*verts.size());
  rval = mb.tag_get_data(ctag, verts, &new_coords[0]);
  CHECK_ERR(rval);
  unsigned int i;
  Range::iterator vit;
  for (vit = verts.begin(), i = 0; vit != verts.end(); ++vit, i+=3) {
    CartVect old_pos(&coords[i]), new_pos(&new_coords[i]);
    CartVect almost_pos = old_pos + .99 * (new_pos - old_pos);
    almost_pos.get(&new_coords[i]);
  }
  rval = mb.set_coords(verts, &new_coords[0]);
  CHECK_ERR(rval);

  LloydSmoother ll2(&mb, NULL, elems);
  ll2.report_its(10);
  rval = ll2.perform_smooth();
  CHECK_ERR(rval);
  std::cout << "Mesh smoothed in " << ll2.num_its() << " iterations." << std::endl;
}

