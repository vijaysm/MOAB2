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
#include "CslamUtils.hpp"

namespace moab {

// forward declarations
class ParallelComm;
class AdaptiveKDTree;
class TupleList;

class Intx2Mesh
{
public:
  Intx2Mesh(Interface * mbimpl);
  virtual ~Intx2Mesh();

  ErrorCode intersect_meshes(EntityHandle mbs1, EntityHandle mbs2,
       EntityHandle & outputSet);

  // mark could be 3 or 4, depending on type
  // this is pure abstract, this need s to be implemented by
  // all derivations
  virtual int computeIntersectionBetweenRedAndBlue(EntityHandle red,
      EntityHandle blue, double * P, int & nP, double & area,
      int markb[4], int markr[4])=0;

  // this is also abstract
  virtual int findNodes(EntityHandle red, EntityHandle blue,
      double * iP, int nP)=0;

  virtual void createTags();
  ErrorCode GetOrderedNeighbors(EntityHandle set, EntityHandle quad,
      EntityHandle neighbors[4]);

  void SetErrorTolerance(double eps) { epsilon_1=eps;}

  void SetEntityType (EntityType tp) { type=tp;}

  // clean some memory allocated
  void clean();

  ErrorCode initialize_local_kdtree(EntityHandle euler_set);

  // this will work in parallel
  ErrorCode locate_departure_points(EntityHandle euler_set); // get the points and elements from the local set
  ErrorCode locate_departure_points(Range & local_verts);
  ErrorCode test_local_box(double *xyz, int from_proc, int remote_index, TupleList *tl);
  ErrorCode inside_entities(double xyz[3], std::vector<EntityHandle> &entities);

  // this will depend on the problem and element type; return true if on the border edge too
  virtual bool is_inside_element(double xyz[3], EntityHandle eh) = 0;
  void set_box_error(double berror)
   {box_error = berror;}

protected: // so it can be accessed in derived classes, InPlane and OnSphere
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

  EntityType type;

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

  ParallelComm * parcomm;

  AdaptiveKDTree *myTree;
  std::vector<double> allBoxes;
  double box_error;
  /* \brief Local root of the kdtree */
  EntityHandle localRoot;
  Range localEnts;// this range is for local elements of interest

};

} /* namespace moab */
#endif /* INTX2MESH_HPP_ */
