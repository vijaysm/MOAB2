/*
 * Intx2MeshInPlane.hpp
 *
 *  Created on: Oct 24, 2012
 *      Author: iulian
 */

#ifndef INTX2MESHINPLANE_HPP_
#define INTX2MESHINPLANE_HPP_

#include "Intx2Mesh.hpp"
namespace moab {

class Intx2MeshInPlane: public moab::Intx2Mesh {
public:
  Intx2MeshInPlane(Interface * mbimpl);
  virtual ~Intx2MeshInPlane();

  int computeIntersectionBetweenRedAndBlue(EntityHandle red, EntityHandle blue,
          double * P, int & nP, double & area, int markb[4], int markr[4]);

  int findNodes(EntityHandle red, EntityHandle blue, double * iP, int nP);
};
} // end namespace moab
#endif /* INTX2MESHINPLANE_HPP_ */
