#ifndef VTK_UTIL
#define VTK_UTIL

//
// ExoIIUtil class: utility class for functions used by both reader
// and writer

#ifndef IS_BUILDING_MB
#error "ExoIIUtil.hpp isn't supposed to be included into an application"
#endif

#include "MBInterface.hpp"


class VtkUtil 
{

public:

    //! vtk element type numbers
  static const int vtkElemType[];

    //! vtk element type names
  static const char *vtkElemNames[];
};

#endif
