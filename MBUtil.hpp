//-------------------------------------------------------------------------
// Filename      : MBUtil.hpp
//
// Purpose       : This file contains utility functions that can be used
//                 with MB
//
// Special Notes : This is a pure virtual class, to prevent instantiation.
//                 All functions are static, called like this:
//                 MBUtil::function_name();
//
// Creator       : Ray J. Meyers
//
// Date          : 09/01/02
//
// Owner         : Ray J. meyers
//-------------------------------------------------------------------------


#ifndef MB_UTIL_HPP
#define MB_UTIL_HPP

#include "MBInterface.hpp"

struct  MBCoord
{
  double x;
  double y;
  double z;
};


//! utility function to be use throughout VERDE
class MBUtil
{
public:

   
  static void normal(MBInterface* MB, MBEntityHandle handle, double& x, double& y, double& z);

  static void centroid(MBInterface *MB, MBEntityHandle handle,MBCoord &coord);
 // static void edge_centers(MBInterface *MB, MBEntityHandle handle, std::vector<MBCoord> &coords_list);

  //static void face_centers(MBInterface *MB, MBEntityHandle handle, std::vector<MBCoord> &coords_list);

private:

  MBUtil(){};

};

#endif
