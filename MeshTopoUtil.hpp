/*!
 *  \class   MeshTopoUtil
 *  \authors Tim Tautges
 *  \date    2/04
 *  \brief   MeshTopoUtil contains general mesh utility functions 
 *          
 */ 

#ifndef MESH_TOPO_UTIL_HPP
#define MESH_TOPO_UTIL_HPP

#include "MBInterface.hpp"

class MeshTopoUtil
{
public:
  MeshTopoUtil(MBInterface *impl) : mbImpl(impl) {}
  
  ~MeshTopoUtil() {}

    //! generate all the AEntities bounding the vertices
  MBErrorCode construct_aentities(const MBRange &vertices);

    //! given an entity, get its average position (avg vertex locations)
  MBErrorCode get_average_position(const MBEntityHandle entity,
                                   double *avg_position);
  
  
private:
  MBInterface *mbImpl;
  
};


#endif

