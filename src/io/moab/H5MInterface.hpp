/*
 * H5MInterface.hpp
 *
 *  Created on: May 5, 2015
 *      Author: iulian
 */

#ifndef H5MINTERFACE_HPP
#define H5MINTERFACE_HPP

namespace moab {

class H5MInterface {
public:
  H5MInterface(const char * filena);
  virtual ~H5MInterface();

  int get_mesh_info( int* verts, int *edges, int*faces, int* regions, int *numdim,  int* parts);

private:
  const char * filename;
  /*mhdf_FileHandle fileHandle;
  unsigned long max_id;*/
};

} /* namespace moab */

#endif /* H5MINTERFACE_HPP */
