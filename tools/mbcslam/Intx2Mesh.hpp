/*
 * Intx2Mesh.hpp
 *
 *  Created on: Oct 2, 2012
 *
 */

#ifndef INTX2MESH_HPP_
#define INTX2MESH_HPP_

#include <iostream>
#include <sstream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"

// these are intersection utils
#include "IntxUtils.hpp"

namespace moab {


class Intx2Mesh
{
public:
  Intx2Mesh(Interface * mbimpl);
  virtual ~Intx2Mesh();

  virtual ErrorCode intersect_meshes(EntityHandle mbs1, EntityHandle mbs2,
        EntityHandle & outputSet)=0; // pure abstract method, needs to be implemented in
  // the derived classes

  virtual void createTags();
  ErrorCode GetOrderedNeighbors(EntityHandle set, EntityHandle quad,
      EntityHandle neighbors[4]);

  void SetErrorTolerance(double eps) { epsilon_1=eps;}
  // clean some memory allocated
  void clean();
// private: everything public?
  Interface * mb;

  EntityHandle mbs1;
  EntityHandle mbs2;

  EntityHandle outSet; // will contain intersection

  // tags used in computation, advanced front
  Tag BlueFlagTag; // to mark blue quads already considered

  Tag RedFlagTag; // to mark red quads already considered

  Range RedEdges; //

  // red parent and blue parent tags
  // these will be on the out mesh
  Tag redParentTag;
  Tag blueParentTag;
  Tag countTag;

  const EntityHandle * redConn;
  const EntityHandle * blueConn;
  CartVect redCoords[4];
  CartVect blueCoords[4];
  double redQuad[8]; // these are in plane
  double blueQuad[8]; // these are in plane

  std::ofstream mout_1[6]; // some debug files
  int dbg_1;
  // for each red edge, we keep a vector of extra nodes, coming from intersections
  // use the index in RedEdges range, instead of a map, as before
  // std::map<EntityHandle, std::vector<EntityHandle> *> extraNodesMap;
  std::vector<std::vector<EntityHandle> *> extraNodesVec;

  double epsilon_1;

};

} /* namespace moab */
#endif /* INTX2MESH_HPP_ */
