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
  static const char *DUAL_SURFACE_TAG_NAME;

    //! tag name for dual curves
  static const char *DUAL_CURVE_TAG_NAME;

    //! tag name for dual cells
  static const char *IS_DUAL_CELL_TAG_NAME;

    //! tag name for dual entitys
  static const char *DUAL_ENTITY_TAG_NAME;

    //! tag name for dual entitys
  static const char *EXTRA_DUAL_ENTITY_TAG_NAME;

    //! tag name for dual entitys
  static const char *DUAL_GRAPHICS_POINT_TAG_NAME;

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
  MBErrorCode construct_dual(MBEntityHandle *entities, 
                             const int num_entities);
  
    //! construct the dual entities for a hex mesh, including dual surfaces & curves
  MBErrorCode construct_hex_dual(MBEntityHandle *entities,
                                 const int num_entities);
  
    //! construct the dual entities for a hex mesh, including dual surfaces & curves
  MBErrorCode construct_hex_dual(MBRange &entities);
  
  //! get the dual entities; if non-null, only dual of entities passed in are returned
  MBErrorCode get_dual_entities(const int dim, 
                                MBEntityHandle *entities, 
                                const int num_entities, 
                                MBRange &dual_ents);
  
  //! get the dual entities; if non-null, only dual of entities passed in are returned
  MBErrorCode get_dual_entities(const int dim, 
                                MBEntityHandle *entities, 
                                const int num_entities, 
                                std::vector<MBEntityHandle> &dual_ents);

    //! return the corresponding dual entity
  MBEntityHandle get_dual_entity(const MBEntityHandle this_ent) const;
  
    //! return the corresponding extra dual entity
  MBEntityHandle get_extra_dual_entity(const MBEntityHandle this_ent);
  
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
  MBEntityHandle next_loop_vertex(const MBEntityHandle last_v,
                                  const MBEntityHandle this_v,
                                  const MBEntityHandle dual_surf);
  
    //! get/set the tag for dual surfaces
  MBTag dualSurface_tag() const;
  MBErrorCode dualSurface_tag(const MBTag tag);
  
    //! get/set the tag for dual curves
  MBTag dualCurve_tag() const;
  MBErrorCode dualCurve_tag(const MBTag tag);

    //! get/set the tag for dual cells
  MBTag isDualCell_tag() const;
  MBErrorCode isDualCell_tag(const MBTag tag);

    //! get/set the tag for dual entities
  MBTag dualEntity_tag() const;
  MBErrorCode dualEntity_tag(const MBTag tag);

    //! get/set the tag for dual entities
  MBTag extraDualEntity_tag() const;
  MBErrorCode extraDualEntity_tag(const MBTag tag);

    //! get/set the tag for dual entities
  MBTag dualGraphicsPoint_tag() const;
  MBErrorCode dualGraphicsPoint_tag(const MBTag tag);

    //! get/set the global id tag
  MBTag globalId_tag() const;
  MBErrorCode globalId_tag(const MBTag tag);

    //! given an entity, return any dual surface or curve it's in
  MBEntityHandle get_dual_hyperplane(const MBEntityHandle ncell);

    //! returns true if first & last vertices are dual to hexes (not faces)
  bool is_blind(const MBEntityHandle chord);
  
    //! set the dual surface or curve for an entity
  MBErrorCode set_dual_surface_or_curve(MBEntityHandle entity, 
                                        const MBEntityHandle dual_hyperplane,
                                        const int dimension);
  
    //! effect atomic pillow operation
  MBErrorCode atomic_pillow(MBEntityHandle odedge, MBEntityHandle &new_hp);

    //! effect reverse atomic pillow operation
  MBErrorCode rev_atomic_pillow(MBEntityHandle pillow, MBRange &chords);

    //! effect face shrink operation
  MBErrorCode face_shrink(MBEntityHandle odedge);
  
    //! effect a face open-collapse operation
  MBErrorCode face_open_collapse(MBEntityHandle ocl, MBEntityHandle ocr,
                                 MBEntityHandle tcm);
  
    //! given a 1-cell and a chord, return the neighboring vertices on the
    //! chord, in the same order as the 1-cell's vertices
  MBErrorCode get_opposite_verts(const MBEntityHandle middle_edge, 
                                 const MBEntityHandle chord, 
                                 MBEntityHandle *verts);

    //! given a dual surface or curve, return the 2-cells, 1-cells, 0-cells, and
    //! loop 0/1-cells, if requested; any of those range pointers can be NULL,
    //! in which case that range isn't returned
  MBErrorCode get_dual_entities(const MBEntityHandle dual_ent,
                                MBRange *dcells,
                                MBRange *dedges,
                                MBRange *dverts,
                                MBRange *dverts_loop,
                                MBRange *dedges_loop);
  
  MBErrorCode list_entities(const MBRange &entities) const;
  MBErrorCode list_entities(const MBEntityHandle *entities,
                            const int num_entities) const;
  
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
  MBErrorCode construct_dual_hyperplanes(const int dim, 
                                         MBEntityHandle *entities, 
                                         const int num_entities);

    //! order 1cells on a chord 
  MBErrorCode order_chord(MBEntityHandle chord_set);
  
    //! make a new dual hyperplane with the specified id; if the id specified is -1,
    //! set the new one's id to the max found
  MBErrorCode construct_new_hyperplane(const int dim, MBEntityHandle &new_hyperplane,
                                       int &id);
  
    //! traverse the cells of a dual hyperplane, starting with this_ent (dimension
    //! of this_ent determines hyperplane dimension)
    //! simpler method for traversing hyperplane, using same basic algorithm but
    //! using MeshTopoUtil::get_bridge_adjacencies
  MBErrorCode traverse_hyperplane(const MBTag hp_tag, 
                                  MBEntityHandle &this_hp, 
                                  MBEntityHandle this_ent);
  
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
  
    //! go through potential dual equivalent edges (edges whose nodes define
    //! multiple edges), and add explicit adjacencies to corrent 2cells
  MBErrorCode check_dual_equiv_edges(MBRange &dual_edges);
  
    //! delete a dual entity; updates primal to no longer point to it
  MBErrorCode delete_dual_entities(MBEntityHandle *entities, const int num_entities);

    //! delete a range of dual entities; updates primal to no longer point to them
  MBErrorCode delete_dual_entities(MBRange &entities);
  
  MBErrorCode foc_gather_data(const MBEntityHandle ocl, const MBEntityHandle ocr, 
                              const MBEntityHandle tcm, 
                                // 0-cells, left & right
                              MBEntityHandle &zclf, MBEntityHandle &zclb, 
                              MBEntityHandle &zcrf, MBEntityHandle &zcrb,
                                // 2-cells, left & right
                              MBEntityHandle &tclu, MBEntityHandle &tclm, MBEntityHandle &tcll, 
                              MBEntityHandle &tcru, MBEntityHandle &tcrm, MBEntityHandle &tcrl,
                                // 3-cells, left & right
                              MBEntityHandle &thclu, MBEntityHandle &thcll, MBEntityHandle &thcmu, 
                              MBEntityHandle &thcml, MBEntityHandle &thcru, MBEntityHandle &thcrl,
                                // sheets
                              MBEntityHandle &sl, MBEntityHandle &sm, MBEntityHandle &sr,
                                // chords     
                              MBEntityHandle &cl, MBEntityHandle &cr);
  
  MBErrorCode foc_get_neighbor_23cells(const MBEntityHandle oc,
                                       const MBEntityHandle tcm,
                                       const MBEntityHandle thcmu,
                                       const MBEntityHandle thcml,
                                       MBEntityHandle &tcu, 
                                       MBEntityHandle &tcm, 
                                       MBEntityHandle &tcl, 
                                       MBEntityHandle &thcu, 
                                       MBEntityHandle &thcl);
  
  MBErrorCode foc_1cells(MBEntityHandle zclf, MBEntityHandle zclb, 
                         MBEntityHandle ocl, MBEntityHandle cl,
                         MBEntityHandle zcrf, MBEntityHandle zcrb, 
                         MBEntityHandle ocr, MBEntityHandle cr,
                         MBEntityHandle sm, MBEntityHandle sr,
                         MBEntityHandle &new_ocb, MBEntityHandle &new_ocf,
                         MBEntityHandle &new_cb, MBEntityHandle &new_cf);
  
    //! break a chord such that first_1cell is the first 1cell, the next 1cell
    //! shares next_0cell, and the new chord starts with the 1cell in the other
    //! direction; if chord is blind, new_chord is not created
  MBErrorCode foc_break_chord(MBEntityHandle chord,
                              MBEntityHandle first_1cell,
                              MBEntityHandle next_0cell,
                              std::vector<MBEntityHandle> &chord_1cells,
                              std::vector<MBEntityHandle> &new_chord_1cells);

    //! replace zc_old with zc_new in 1cell oc
  MBErrorCode foc_replace_0cell(MBEntityHandle oc,
                                MBEntityHandle zc_old,
                                MBEntityHandle zc_new);
  
    //! check sense of connect arrays, and reverse/rotate if necessary
  MBErrorCode fs_check_quad_sense(MBEntityHandle hex0,
                                  MBEntityHandle quad0,
                                  std::vector<MBEntityHandle> *connects);
  
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
  MBTag globalIdTag;
};

inline MBTag DualTool::dualSurface_tag() const
{
  return dualSurfaceTag;
}

inline MBTag DualTool::dualCurve_tag() const
{
  return dualCurveTag;
}

inline MBTag DualTool::isDualCell_tag() const
{
  return isDualCellTag;
}

inline MBTag DualTool::dualEntity_tag() const
{
  return dualEntityTag;
}

inline MBTag DualTool::extraDualEntity_tag() const
{
  return extraDualEntityTag;
}

inline MBTag DualTool::dualGraphicsPoint_tag() const
{
  return dualGraphicsPointTag;
}

inline MBTag DualTool::globalId_tag() const
{
  return globalIdTag;
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
  
  //! get/set the tag for dual entities
inline MBErrorCode DualTool::globalId_tag(const MBTag tag)
{
  MBErrorCode result = MB_FAILURE;
  if (0 == globalIdTag && tag || globalIdTag != tag) {
    globalIdTag = tag;
    result = MB_SUCCESS;
  }
  
  return result;
}
  
#endif

