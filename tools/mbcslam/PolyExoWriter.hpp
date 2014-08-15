/*
 * PolyExoWriter.hpp
 *
 *  Created on: Aug 15, 2014
 */

#ifndef POLYEXOWRITER_HPP_
#define POLYEXOWRITER_HPP_

// use netcdf cpp interface, maybe it is easier to copy code
// later on, maybe use c interface
#include "moab/Interface.hpp"

namespace moab {

class PolyExoWriter {
public:
  PolyExoWriter(Interface * moab);
  virtual ~PolyExoWriter();

  ErrorCode write_poly_set(EntityHandle set, const char * filename);
private:
  Interface* mb;
};



} /* namespace moab */
#endif /* POLYEXOWRITER_HPP_ */
