/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

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
    //! tag name for dual surfaces
  static char *DUAL_SURFACE_TAG_NAME;

    //! tag name for dual curves
  static char *DUAL_CURVE_TAG_NAME;

    //! tag name for dual cells
  static char *IS_DUAL_CELL_TAG_NAME;

    //! tag name for dual entitys
  static char *DUAL_ENTITY_TAG_NAME;

    //! tag name for dual entitys
  static char *EXTRA_DUAL_ENTITY_TAG_NAME;

    //! tag name for dual entitys
  static char *DUAL_GRAPHICS_POINT_TAG_NAME;

    //! struct for storing a graphics pt
  class GraphicsPoint 
  {
  public:
    GraphicsPoint() 
      {xyz[0] = 0.0; xyz[1] = 0.0; xyz[2] = 0.0; id = -1;}

    GraphicsPoint(float xi, float yi, float zi, int idi) 
      {xyz[0] = xi; xyz[1] = yi; xyz[2] = zi; id = idi;}

    GraphicsPoint(float xyzi[3], int idi)
      {xyz[0] = xyzi[0]; xyz[1] = xyzi[1]; xyz[2] = xyzi[2]; id = idi;}

    GraphicsPoint(double xyzi[3], int idi)
      {xyz[0] = xyzi[0]; xyz[1] = xyzi[1]; xyz[2] = xyzi[2]; id = idi;}

    GraphicsPoint(const GraphicsPoint &gp) 
      {xyz[0] = gp.xyz[0]; xyz[1] = gp.xyz[1]; xyz[2] = gp.xyz[2]; id = gp.id;}
    
    float xyz[3];
    int id;
  };
  
  DualTool(MBInterface *impl);
  
  ~DualTool();

    //! construct the dual entities for the entire mesh
  MBErrorCode construct_dual();
  
    //! construct the dual entities for a hex mesh, including dual surfaces & curves
  MBErrorCode construct_hex_dual();
  
  //! get the faces of the dual
  MBErrorCode get_dual_entities(const int dim, MBRange &dual_ents);
  
  //! get the faces of the dual
  MBErrorCode get_dual_entities(const int dim, 
                                std::vector<MBEntityHandle> &dual_ents);
  
    //! get the d-dimensional hyperplane sets; static 'cuz it's easy to do without an active
    //! dualtool
  static MBErrorCode get_dual_hyperplanes(const MBInterface *impl, const int dim, 
                                          MBRange &dual_ents);

    //! get the graphics points for single entity (dual_ent CAN'T be a set);
    //! returns multiple facets, each with npts[i] points
  MBErrorCode get_graphics_points(MBEntityHandle dual_ent,
                                  std::vector<int> &npts,
                                  std::vector<GraphicsPoint> &gpoints);
  
    //! get the graphics points for a range of entities or sets (if set, the
    //! entities in those sets); optionally reset ids on points
  MBErrorCode get_graphics_points(const MBRange &in_range,
                                  std::vector<GraphicsPoint> &gpoints,
                                  const bool assign_ids = false,
                                  const int start_id = 0);
  
    //! given a last_v (possibly zero) and this_v, find the next loop vertex on 
    //! this dual surface
  static MBEntityHandle next_loop_vertex(const MBEntityHandle last_v,
                                         const MBEntityHandle this_v,
                                         const MBEntityHandle dual_surf);
  
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

    //! get/set the tag for dual entities
  MBTag extraDualEntity_tag();
  MBErrorCode extraDualEntity_tag(const MBTag tag);

    //! get/set the tag for dual entities
  MBTag dualGraphicsPoint_tag();
  MBErrorCode dualGraphicsPoint_tag(const MBTag tag);

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
  
    //! traverse dual faces of input dimension, constructing
    //! dual hyperplanes of them in sets as it goes
  MBErrorCode construct_dual_hyperplanes(const int dim);
  
    //! connect dual surfaces with dual curves using parent/child connections
  MBErrorCode construct_hp_parent_child();
  
  //! given an edge handle, return a list of dual vertices in radial order 
  //! around the edge; also returns whether this edge is on the boundary
  MBErrorCode get_radial_dverts(const MBEntityHandle edge,
                                std::vector<MBEntityHandle> &rad_verts,
                                bool &bdy_edge);
  
  MBErrorCode construct_dual_vertex(MBEntityHandle entity, 
                                    MBEntityHandle &dual_ent, 
                                    const bool extra = false,
                                    const bool add_graphics_pt = true);

    //! add a graphics point to an entity (on a tag)
  MBErrorCode add_graphics_point(MBEntityHandle entity,
                                 double *avg_pos = NULL);
  
    //! get points defining facets of a 2cell
  MBErrorCode get_cell_points(MBEntityHandle dual_ent,
                              std::vector<int> &npts,
                              std::vector<GraphicsPoint> &points);

    //! if this_ent is an edge, is a dual entity, and has quads as
    //! its vertices' dual entities, return true, otherwise false
  bool check_1d_loop_edge(MBEntityHandle this_ent);
  
    //! private copy of interface *
  MBInterface *mbImpl;

    //! static constant number of points bounding any cell
  enum
  {
    GP_SIZE=20
  };
  
    //! tags used for dual surfaces, curves, cells, entities
  MBTag dualCurveTag;
  MBTag dualSurfaceTag;
  MBTag isDualCellTag;
  MBTag dualEntityTag;
  MBTag extraDualEntityTag;
  MBTag dualGraphicsPointTag;
  MBTag categoryTag;
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

inline MBTag DualTool::extraDualEntity_tag()
{
  return extraDualEntityTag;
}

inline MBTag DualTool::dualGraphicsPoint_tag()
{
  return dualGraphicsPointTag;
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
  
  //! get/set the tag for dual entities
inline MBErrorCode DualTool::extraDualEntity_tag(const MBTag tag)
{
  MBErrorCode result = MB_FAILURE;
  if (0 == extraDualEntityTag && tag || extraDualEntityTag != tag) {
    extraDualEntityTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
  //! get/set the tag for dual entities
inline MBErrorCode DualTool::dualGraphicsPoint_tag(const MBTag tag)
{
  MBErrorCode result = MB_FAILURE;
  if (0 == dualGraphicsPointTag && tag || dualGraphicsPointTag != tag) {
    dualGraphicsPointTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
#endif

