/*
 * Intx2MeshOnSphere.hpp
 *
 *  Created on: Oct 3, 2012
 *      Author: iulian
 */

#ifndef INTX2MESHONSPHERE_HPP_
#define INTX2MESHONSPHERE_HPP_

#include "Intx2Mesh.hpp"

namespace moab {

class Intx2MeshOnSphere: public moab::Intx2Mesh
{
public:
  Intx2MeshOnSphere(Interface * mbimpl);
  virtual ~Intx2MeshOnSphere();

  void SetRadius(double radius) { R=radius ;}
  // main method to intersect meshes on a sphere
  ErrorCode intersect_meshes(EntityHandle mbset1, EntityHandle mbset2,
      EntityHandle & outputSet);


private:
  int plane; // current gnomonic plane
  double R; // radius of the sphere

  // these methods are really private, we should not expose them yet
  // maybe later?
  int computeIntersectionBetweenRedAndBlue(EntityHandle red, EntityHandle blue,
        double * P, int & nP, double & area, int markb[4], int markr[4]);

  int findNodes(EntityHandle red, EntityHandle blue, double * iP, int nP);
};

} /* namespace moab */
#endif /* INTX2MESHONSPHERE_HPP_ */
