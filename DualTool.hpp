/*!
 *  \class   DualTool
 *  \authors Tim Tautges
 *  \date    2/04
 *  \brief   Tools for constructing and working with mesh duals (both tet- and hex-based,
 *           though some functions may not make sense for tet duals)
 *          
 */ 

#ifndef DUAL_TOOL_HPP
#define DUAL_TOOL_HPP

#include "MBInterface.hpp"

class DualTool
{
public:
  DualTool(MBInterface *impl);
  
  ~DualTool();

    //! construct the dual entities for the entire mesh
  MBErrorCode construct_dual();
  
    //! construct the dual entities for a hex mesh, including dual surfaces & curves
  MBErrorCode construct_hex_dual();
  
    //! get the sets representing dual surfaces (sheets) in hex dual
  MBErrorCode get_dual_surfaces(MBRange &dual_surfs);
  
    //! get the sets representing dual curves (chords) in hex dual
  MBErrorCode get_dual_curves(MBRange &dual_curves);
  
    //! get the cells of the dual
  MBErrorCode get_dual_cells(MBRange &dual_cells);

    //! get/set the tag for dual surfaces
  MBTag dualSurface_tag();
  MBErrorCode dualSurface_tag(const MBTag tag);
  
    //! get/set the tag for dual curves
  MBTag dualCurve_tag();
  MBErrorCode dualCurve_tag(const MBTag tag);

    //! get/set the tag for dual cells
  MBTag isDualCell_tag();
  MBErrorCode isDualCell_tag(const MBTag tag);

    //! get/set the tag for dual entities
  MBTag dualEntity_tag();
  MBErrorCode dualEntity_tag(const MBTag tag);

    //! tag name for dual surfaces
  static char *DUAL_SURFACE_TAG_NAME;

    //! tag name for dual curves
  static char *DUAL_CURVE_TAG_NAME;

    //! tag name for dual cells
  static char *IS_DUAL_CELL_TAG_NAME;

    //! tag name for dual entitys
  static char *DUAL_ENTITY_TAG_NAME;

  
private:

    //! construct dual vertices for specified regions
  MBErrorCode construct_dual_vertices(const MBRange &all_regions,
                                      MBRange &new_dual_ents);
  
    //! construct dual edges for specified faces
  MBErrorCode construct_dual_edges(const MBRange &all_faces,
                                      MBRange &new_dual_ents);
  
    //! construct dual faces for specified edges
  MBErrorCode construct_dual_faces(const MBRange &all_edges,
                                      MBRange &new_dual_ents);
  
    //! construct dual cells for specified vertices
  MBErrorCode construct_dual_cells(const MBRange &all_verts,
                                   MBRange &new_dual_ents);
  
  //! given an edge handle, return a list of dual vertices in radial order 
  //! around the edge
  MBErrorCode get_radial_dverts(const MBEntityHandle edge,
                                std::vector<MBEntityHandle> &rad_verts);
  
    //! private copy of interface *
  MBInterface *mbImpl;
  
    //! tags used for dual surfaces, curves, cells, entities
  MBTag dualCurveTag;
  MBTag dualSurfaceTag;
  MBTag isDualCellTag;
  MBTag dualEntityTag;
};

inline MBTag DualTool::dualSurface_tag()
{
  return dualSurfaceTag;
}

inline MBTag DualTool::dualCurve_tag()
{
  return dualCurveTag;
}

inline MBTag DualTool::isDualCell_tag()
{
  return isDualCellTag;
}

inline MBTag DualTool::dualEntity_tag()
{
  return dualEntityTag;
}

  //! get/set the tag for dual surfaces
inline MBErrorCode DualTool::dualSurface_tag(const MBTag tag) 
{
  MBErrorCode result = MB_FAILURE;
  if (0 == dualSurfaceTag && tag || dualSurfaceTag != tag) {
    dualSurfaceTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual curves
inline MBErrorCode DualTool::dualCurve_tag(const MBTag tag)
{
  MBErrorCode result = MB_FAILURE;
  if (0 == dualCurveTag && tag || dualCurveTag != tag) {
    dualCurveTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual cells
inline MBErrorCode DualTool::isDualCell_tag(const MBTag tag)
{
  MBErrorCode result = MB_FAILURE;
  if (0 == isDualCellTag && tag || isDualCellTag != tag) {
    isDualCellTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual entities
inline MBErrorCode DualTool::dualEntity_tag(const MBTag tag)
{
  MBErrorCode result = MB_FAILURE;
  if (0 == dualEntityTag && tag || dualEntityTag != tag) {
    dualEntityTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
#endif

