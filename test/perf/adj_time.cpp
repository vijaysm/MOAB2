#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include <iostream>
#include <assert.h>
#include <time.h>

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)
#ifndef SRCDIR
# define SRCDIR .
#endif

const char* default_input_file = "../mb_big_test.g";

int main()
{
  ErrorCode rval;
  Core moab;
  Interface& mb = moab;
  
    // load test file
  rval = mb.load_file( default_input_file );
  if (MB_SUCCESS != rval) {
    std::cerr << default_input_file <<": failed to load file." << std::endl;
    return 1;
  }
  
    // get all region elements
  Range vols;
  rval = mb.get_entities_by_dimension( 0, 3, vols );
  if (MB_SUCCESS != rval)
    return 2;
  if (vols.empty())
    return 1;
  
    // create internal face elements
  Range faces;
  rval = mb.get_adjacencies( vols, 2, true, faces, Interface::UNION );
  if (MB_SUCCESS != rval)
    return 2;
  assert(faces.size() > vols.size());
  
    // time query of all adjacent volumes
  std::vector<EntityHandle> adj;
  clock_t t_0 = clock();
  for (Range::iterator i = faces.begin(); i != faces.end(); ++i) {
    adj.clear();
    rval = mb.get_adjacencies( &*i, 1, 3, false, adj );
    if (MB_SUCCESS != rval)
      return 2;
    assert( adj.size() == 1 || adj.size() == 2 );
  }
  clock_t t_up = clock() - t_0;
  std::cout << "Querying of volumes for " << faces.size() << " faces: "
            << t_up/(double)CLOCKS_PER_SEC << " seconds" << std::endl;
  
    // time downward adjacency query from volumes to faces
  t_0 = clock();
  for (Range::iterator i = vols.begin(); i != vols.end(); ++i) {
    adj.clear();
    rval = mb.get_adjacencies( &*i, 1, 1, false, adj );
    if (MB_SUCCESS != rval)
      return 2;
    assert( adj.size() > 3 );
  }
  clock_t t_down = clock() - t_0;
  std::cout << "Querying of faces for " << vols.size() << " volumes: "
            << t_down/(double)CLOCKS_PER_SEC << " seconds" << std::endl;
  
  
  return 0;
}

