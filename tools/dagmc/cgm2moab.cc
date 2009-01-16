#include "GeometryQueryTool.hpp"
#include "GeometryQueryEngine.hpp"
#include "ModelQueryEngine.hpp"
#include "RefEntityName.hpp"
#include "GMem.hpp"

#include "RefGroup.hpp"
#include "RefVolume.hpp"
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVertex.hpp"

#include "SenseEntity.hpp"
#include "Surface.hpp"
#include "Curve.hpp"

#include "MBCore.hpp"
#include "MBInterface.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"

#include "GeomTopoTool.hpp"

#include <stdio.h>
#include <algorithm>

#include "cgm2moab.hpp"

// default settings
const char DEFAULT_NAME[] = "cgm2moab";

// settings
const char* name = DEFAULT_NAME;
const char* file_type = 0;

// copy geometry into mesh database
bool cgm2moab(MBInterface* iface,
              double faceting_tol,
              int norm_tol,
              double len_tol,
              int actuate_attribs)
{
  MBErrorCode rval;
  GeomTopoTool geomTool(iface);

    // get some tag handles
  MBTag geom_tag, id_tag, name_tag, category_tag;
  rval = iface->tag_create( GEOM_DIMENSION_TAG_NAME, sizeof(int), MB_TAG_SPARSE, 
                            MB_TYPE_INTEGER, geom_tag, 0, true );
  assert(!rval);
  rval = iface->tag_create( GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE, 
                            MB_TYPE_INTEGER, id_tag, 0, true );
  assert(!rval);
  rval = iface->tag_create( NAME_TAG_NAME, NAME_TAG_SIZE, MB_TAG_SPARSE, 
                            MB_TYPE_OPAQUE, name_tag, 0, true );
  assert(!rval);
  std::vector<MBTag> extra_name_tags;
  
  rval = iface->tag_create( CATEGORY_TAG_NAME, CATEGORY_TAG_SIZE, 
                            MB_TAG_SPARSE, MB_TYPE_OPAQUE, category_tag, 0, true );
  assert(!rval);
  
    // CGM data
  std::map<RefEntity*,MBEntityHandle> entmap[5]; // one for each dim, and one for groups
  std::map<RefEntity*,MBEntityHandle>::iterator ci;
  const char geom_categories[][CATEGORY_TAG_SIZE] = 
    {"Vertex\0", "Curve\0", "Surface\0", "Volume\0", "Group\0"};
  const char* const names[] = { "Vertex", "Curve", "Surface", "Volume"};
  DLIList<RefEntity*> entlist;
  DLIList<ModelEntity*> me_list;

    // create entity sets for all geometric entities
  for (int dim = 0; dim < 4; ++dim) {
    entlist.clean_out();
    GeometryQueryTool::instance()->ref_entity_list( names[dim], entlist, true );
    
    entlist.reset();
    for (int i = entlist.size(); i--; ) {
      RefEntity* ent = entlist.get_and_step();
      MBEntityHandle handle;
      rval = iface->create_meshset( dim == 1 ? MESHSET_ORDERED : MESHSET_SET, handle );
      if (MB_SUCCESS != rval)
        return false;
    
      entmap[dim][ent] = handle;
      
      rval = iface->tag_set_data( geom_tag, &handle, 1, &dim );
      if (MB_SUCCESS != rval)
        return rval;
      int id = ent->id();
      rval = iface->tag_set_data( id_tag, &handle, 1, &id );
      if (MB_SUCCESS != rval)
        return rval;

      rval = iface->tag_set_data( category_tag, &handle, 1, &geom_categories[dim] );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  
    // create topology for all geometric entities
  for (int dim = 1; dim < 4; ++dim) {
    for (ci = entmap[dim].begin(); ci != entmap[dim].end(); ++ci) {
      entlist.clean_out();
      ci->first->get_child_ref_entities( entlist );
    
      entlist.reset();
      for (int i = entlist.size(); i--; ) {
        RefEntity* ent = entlist.get_and_step();
        MBEntityHandle h = entmap[dim-1][ent];
        rval = iface->add_parent_child( ci->second, h );
        if (MB_SUCCESS != rval)
          return false;
      }
    }
  }
  
    // store CoFace senses
  for (ci = entmap[2].begin(); ci != entmap[2].end(); ++ci) {
    RefFace* face = (RefFace*)(ci->first);
    BasicTopologyEntity *forward = 0, *reverse = 0;
    for (SenseEntity* cf = face->get_first_sense_entity_ptr();
         cf; cf = cf->next_on_bte()) {
      BasicTopologyEntity* vol = cf->get_parent_basic_topology_entity_ptr();
      if (cf->get_sense() == CUBIT_UNKNOWN || 
          cf->get_sense() != face->get_surface_ptr()->bridge_sense()) {
        if (reverse) {
          std::cout << "Surface " << face->id() << " has reverse senes " <<
                       "with multiple volume " << reverse->id() << " and " <<
                       "volume " << vol->id() << std::endl;
          return false;
        }
        reverse = vol;
      }
      if (cf->get_sense() == CUBIT_UNKNOWN || 
          cf->get_sense() == face->get_surface_ptr()->bridge_sense()) {
        if (forward) {
          std::cout << "Surface " << face->id() << " has forward senes " <<
                       "with multiple volume " << forward->id() << " and " <<
                       "volume " << vol->id() << std::endl;
          return false;
        }
        forward = vol;
      }
    }
    
    if (forward) {
      rval = geomTool.set_sense( ci->second, entmap[3][forward], true );
      if (MB_SUCCESS != rval)
        return rval;
    }
    if (reverse) {
      rval = geomTool.set_sense( ci->second, entmap[3][reverse], false );
      if (MB_SUCCESS != rval)
        return rval;
    }
  }
  
    // create entity sets for all ref groups
  DLIList<CubitString*> name_list;
  entlist.clean_out();
  GeometryQueryTool::instance()->ref_entity_list( "group", entlist );
  entlist.reset();
  for (int i = entlist.size(); i--; ) {
    RefEntity* grp = entlist.get_and_step();
    name_list.clean_out();
    RefEntityName::instance()->get_refentity_name( grp, name_list, true );
    if (name_list.size() == 0)
      continue;
    name_list.reset();
    CubitString name = *name_list.get();
    
    MBEntityHandle h;
    rval = iface->create_meshset( MESHSET_SET, h );
    if (MB_SUCCESS != rval)
      return false;
    
    char namebuf[NAME_TAG_SIZE];
    memset( namebuf, '\0', NAME_TAG_SIZE );
    strncpy( namebuf, name.c_str(), NAME_TAG_SIZE - 1 );
    if (name.length() >= (unsigned)NAME_TAG_SIZE) 
      std::cout << "WARNING: group name '" << name.c_str() 
                << "' truncated to '" << namebuf << "'" << std::endl;
    rval = iface->tag_set_data( name_tag, &h, 1, namebuf );
    if (MB_SUCCESS != rval)
      return false;
      
    int id = grp->id();
    rval = iface->tag_set_data( id_tag, &h, 1, &id );
    if (MB_SUCCESS != rval)
      return false;
      
    rval = iface->tag_set_data( category_tag, &h, 1, &geom_categories[4] );
    if (MB_SUCCESS != rval)
      return false;
      
    if (name_list.size() > 1) {
      for (int j = extra_name_tags.size(); j < name_list.size(); ++j) {
        sprintf( namebuf, "EXTRA_%s%d", NAME_TAG_NAME, j );
        MBTag t;
        rval = iface->tag_create( namebuf, NAME_TAG_SIZE, MB_TAG_SPARSE, MB_TYPE_OPAQUE, t, 0, true );
        assert(!rval);
        extra_name_tags.push_back(t);
      }
        
      for (int j = 0; j < name_list.size(); ++j) {
        name = *name_list.get_and_step();
        memset( namebuf, '\0', NAME_TAG_SIZE );
        strncpy( namebuf, name.c_str(), NAME_TAG_SIZE - 1 );
        if (name.length() >= (unsigned)NAME_TAG_SIZE) 
          std::cout << "WARNING: group name '" << name.c_str() 
                    << "' truncated to '" << namebuf << "'" << std::endl;
        rval = iface->tag_set_data( extra_name_tags[j], &h, 1, namebuf );
        if (MB_SUCCESS != rval)
          return false;
      }
    }
      
    entmap[4][grp] = h;
  }
  
    // store contents for each group
  entlist.reset();
  for (ci = entmap[4].begin(); ci != entmap[4].end(); ++ci) {
    RefGroup* grp = (RefGroup*)(ci->first);
    entlist.clean_out();
    grp->get_child_ref_entities( entlist );
    
    MBRange entities;
    while (entlist.size()) {
      RefEntity* ent = entlist.pop();
      int dim = ent->dimension();
      if (dim < 0) {
        if (entmap[4].find(ent) != entmap[4].end())
          entities.insert( entmap[4][ent] );
      }
      else if (dim < 4) {
        if (entmap[dim].find(ent) != entmap[dim].end())
          entities.insert( entmap[dim][ent] );
      }
    }
    
    if (!entities.empty()) {
      rval = iface->add_entities( ci->second, entities );
      if (MB_SUCCESS != rval)
        return false;
    }
  }
  
    // done with volumes and groups
  entmap[3].clear();
  entmap[4].clear();
  
    // create geometry for all vertices and replace 
    // vertex set handles with vertex handles in map
  for (ci = entmap[0].begin(); ci != entmap[0].end(); ++ci) {
    CubitVector pos = dynamic_cast<RefVertex*>(ci->first)->coordinates();
    double coords[3] = {pos.x(), pos.y(), pos.z()};
    MBEntityHandle vh;
    rval = iface->create_vertex( coords, vh );
    if (MB_SUCCESS != rval)
      return false;
    
    rval = iface->add_entities( ci->second, &vh, 1 );
    if (MB_SUCCESS != rval)
      return false;
    
    ci->second = vh;
  }
  
    // create geometry for all curves
  GMem data;
  for (ci = entmap[1].begin(); ci != entmap[1].end(); ++ci) {
    RefEdge* edge = dynamic_cast<RefEdge*>(ci->first);
    Curve* curve = edge->get_curve_ptr();
    GeometryQueryEngine* gqe = curve->get_geometry_query_engine();
    data.clean_out();
    int count;
    CubitStatus s = gqe->get_graphics( curve, count, &data, faceting_tol);
    if (CUBIT_SUCCESS != s)
      return false;
      
    std::vector<CubitVector> points;
    for (int i = 0; i < data.pointListCount; ++i)
      points.push_back( CubitVector( data.point_list()[i].x,
                                     data.point_list()[i].y,
                                     data.point_list()[i].z ) );

      // need to reverse data?
    if (curve->bridge_sense() == CUBIT_REVERSED) 
      std::reverse( points.begin(), points.end() );
    
       // check for closed curve
    RefVertex *start_vtx, *end_vtx;
    start_vtx = edge->start_vertex();
    end_vtx = edge->end_vertex();
    
      // Special case for point curve
    if (points.size() < 2) {
      if (start_vtx != end_vtx || curve->measure() > gqe->get_sme_resabs_tolerance()) {
        std::cerr << "Warning: No facetting for curve " << edge->id() << std::endl;
        continue;
      }
      MBEntityHandle h = entmap[0][start_vtx];
      rval = iface->add_entities( ci->second, &h, 1 );
      if (MB_SUCCESS != rval)
        return false;
      continue;
    }
    
    const bool closed = (points.front() - points.back()).length() < gqe->get_sme_resabs_tolerance();
    if (closed != (start_vtx == end_vtx)) {
      std::cerr << "Warning: topology and geometry inconsistant for possibly closed curve "
                << edge->id() << std::endl;
    }
    
      // check proximity of vertices to end coordinates
    if ((start_vtx->coordinates() - points.front()).length() > gqe->get_sme_resabs_tolerance()
     || (  end_vtx->coordinates() - points.back() ).length() > gqe->get_sme_resabs_tolerance()) {
      std::cerr << "Warning: vertices not at ends of curve " << edge->id() << std::endl;
    }
    
      // create interior points
    std::vector<MBEntityHandle> verts, edges;
    verts.push_back( entmap[0][start_vtx] );
    for (size_t i = 1; i < points.size() - 1; ++i) {
      double coords[] = { points[i].x(), points[i].y(), points[i].z() };
      MBEntityHandle h;
      rval = iface->create_vertex( coords, h );
      if (MB_SUCCESS != rval)
        return false;
      verts.push_back( h );
    }
    verts.push_back( entmap[0][end_vtx] );
    
      // create edges
    for (size_t i = 0; i < verts.size()-1; ++i) {
      MBEntityHandle h;
      rval = iface->create_element( MBEDGE, &verts[i], 2, h );
      if (MB_SUCCESS != rval)
        return false;
      edges.push_back( h );
    }
    
      // if closed, remove duplicate
    if (verts.front() == verts.back())
      verts.pop_back();
    
    rval = iface->add_entities( ci->second, &verts[0], verts.size() );
    if (MB_SUCCESS != rval)
      return false;
    rval = iface->add_entities( ci->second, &edges[0], edges.size() );
    if (MB_SUCCESS != rval)
      return false;
  }
  
    // create geometry for all surfaces
  for (ci = entmap[2].begin(); ci != entmap[2].end(); ++ci) {
    RefFace* face = dynamic_cast<RefFace*>(ci->first);
    Surface* surf = face->get_surface_ptr();
    GeometryQueryEngine* gqe = surf->get_geometry_query_engine();
    data.clean_out();
    int nt, np, nf;
    CubitStatus s = gqe->get_graphics( surf, nt, np, nf, &data, 
                                       norm_tol, faceting_tol, len_tol);
    if (CUBIT_SUCCESS != s)
      return false;

      // declare array of all vertex handles
    std::vector<MBEntityHandle> verts( data.pointListCount, 0 );
    
      // get list of vertices in surface
    me_list.clean_out();
    ModelQueryEngine::instance()->query_model( *face, DagType::ref_vertex_type(), me_list );

      // for each vertex, find coincident point in facets
    for (int i = me_list.size(); i--; ) {
      RefVertex* vtx = dynamic_cast<RefVertex*>(me_list.get_and_step());
      CubitVector pos = vtx->coordinates();

      for (int j = 0; j < data.pointListCount; ++j) {
        CubitVector vpos( data.point_list()[j].x,
                          data.point_list()[j].y,
                          data.point_list()[j].z );
        if ((pos - vpos).length_squared() < gqe->get_sme_resabs_tolerance()*gqe->get_sme_resabs_tolerance()) {
          if (verts[j])
            std::cerr << "Warning: Coincident vertices in surface " << face->id() << std::endl;
          verts[j] = entmap[0][vtx];
          break;
        }
      }
    }
    
      // now create vertices for the remaining points in the facetting
    for (int i = 0; i < data.pointListCount; ++i) {
      if (verts[i]) // if a geometric vertex
        continue;
      double coords[] = { data.point_list()[i].x,
                          data.point_list()[i].y,
                          data.point_list()[i].z };
      rval = iface->create_vertex( coords, verts[i] );
      if (MB_SUCCESS != rval)
        return rval;
    }
    
      // now create facets
    MBRange facets;
    std::vector<MBEntityHandle> corners;
    for (int i = 0; i < data.fListCount; i += data.facet_list()[i]+1) {
      int* facet = data.facet_list() + i;
      corners.resize( *facet );
      for (int j = 1; j <= *facet; ++j) {
        if (facet[j] >= (int)verts.size()) {
          std::cerr << "ERROR: Invalid facet data for surface " << face->id() << std::endl;
          return false;
        }
        corners[j-1] = verts[facet[j]];
      }
      MBEntityType type;
      if (*facet == 3)
        type = MBTRI;
      else {
        std::cerr << "Warning: non-triangle facet in surface " << face->id() << std::endl;
        if (*facet == 4)
          type = MBQUAD;
        else
          type = MBPOLYGON;
      }
      
      // if (surf->bridge_sense() == CUBIT_REVERSED)
      //   std::reverse( corners.begin(), corners.end() );
      
      MBEntityHandle h;
      rval = iface->create_element( type, &corners[0], corners.size(), h );
      if (MB_SUCCESS != rval)
        return false;
        
      facets.insert( h );
    }
    
      // add vertices and facets to surface set
    rval = iface->add_entities( ci->second, &verts[0], verts.size() );
    if (MB_SUCCESS != rval)
      return false;
    rval = iface->add_entities( ci->second, facets );
    if (MB_SUCCESS != rval)
      return false;
  }
  
  return true;
}

    
