#include "VtkUtil.hpp"
#include "MBInternals.hpp"
#include "MBInterface.hpp"

const int VtkUtil::vtkElemType[] =
{
    1,  //  MBVERTEX
    3,  //  MBEDGE,     
    5,  //  MBTRI,      
    9,  //  MBQUAD,     
    7,  //  MBPOLYGON,  
    10, //  MBTET,      
    14, //  MBPYRAMID,  
    13, //  MBPRISM,    
    -1, //  MBKNIFE,    
    12, //  MBHEX,      
    -1, //  MBPOLYHEDRON
    -1, //  MBENTITYSET,
    -1  //  MBMAXTYPE    
};

const char *VtkUtil::vtkElemNames[] = 
{
  "VTK_VERTEX", //  MBVERTEX     
  "VTK_LINE", //  MBEDGE,      
  "VTK_TRIANGLE", //  MBTRI,       
  "VTK_QUAD", //  MBQUAD,      
  "VTK_POLYGON", //  MBPOLYGON,   
  "VTK_TETRA", //  MBTET,       
  "VTK_PYRAMID", //  MBPYRAMID,   
  "VTK_WEDGE", //  MBPRISM,     
  "VTK_NO_EQUIVALENT", //  MBKNIFE,     
  "VTK_HEXAHEDRON", //  MBHEX,       
  "VTK_NO_EQUIVALENT", //  MBPOLYHEDRON 
  "VTK_NO_EQUIVALENT", //  MBENTITYSET, 
  "VTK_NO_EQUIVALENT"  //  MBMAXTYPE    
};
