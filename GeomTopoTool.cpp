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

#include "GeomTopoTool.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include <assert.h>

    //! Restore parent/child links between GEOM_TOPO mesh sets
MBErrorCode GeomTopoTool::restore_topology() 
{
  
    // look for geometric topology sets and restore parent/child links between them
    // algorithm:
    // - for each entity
    //   . get nodes inclusive in range and store as tag on entity
    // - for each entity of dimension, starting low & working upward:
    //   . for each vertex: if part of boundary:
    //     - get all connected entities of dim d-1
    //     - look for intersections in ranges of (d-1)-entities with d-entity (whole intersection)
    //     - if whole intersection, make parent/child link between (d-1)- and d-entities

    // get the geom topology tag
  MBTag geom_tag, verts_tag;
  MBErrorCode result = mdbImpl->tag_create(GEOM_DIMENSION_TAG_NAME, 4, 
                                            MB_TAG_SPARSE, geom_tag, NULL);
  if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result)
    return result;
  
    // get all sets with this tag
  MBRange geom_sets;
  result = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag, NULL, 1,
                                                geom_sets);
  if (MB_SUCCESS != result || geom_sets.empty()) 
    return result;

    // create the tag holding the vertex ranges
  MBRange *dum = NULL;
  result = mdbImpl->tag_create("__vertex_range", sizeof(MBRange*), MB_TAG_SPARSE, verts_tag, &dum);
  if (MB_SUCCESS != result)
    return result;

  result = construct_vertex_ranges(geom_sets, verts_tag);
  if (MB_SUCCESS != result)
    return result;

  MBRange entities[4];
  result = separate_by_dimension(geom_sets, entities, geom_tag);
  if (MB_SUCCESS != result)
    return result;
  
    // loop over dimensions
  std::vector<MBRange*> dm1_ranges, d_ranges, d0_ranges;
  MBRange::iterator d_it, v_it;
  std::vector<MBRange*>::iterator d_rit, v_rit;
  MBRange inters;
  std::vector<MBEntityHandle> dm1_ents;
  MBRange real_dm1_ents;

    // pre-load d_ranges with vertex ranges
  d_ranges.resize(entities[0].size());
  result = mdbImpl->tag_get_data(verts_tag, entities[0], &d_ranges[0]);
  if (MB_SUCCESS != result)
    return result;
  d0_ranges = d_ranges;
  
  for (int dim = 1; dim <= 3; dim++) {

        // get the range * for these entities
    dm1_ranges.swap(d_ranges);
    d_ranges.resize(entities[dim].size());
    result = mdbImpl->tag_get_data(verts_tag, entities[dim], &d_ranges[0]);
    if (MB_SUCCESS != result)
      return result;
    
    for (d_it = entities[dim].begin(), d_rit = d_ranges.begin(); 
         d_it != entities[dim].end(); d_it++, d_rit++) {

        // iterate over vertices, finding any with intersected ranges
        //dm1_ents.clear();
      real_dm1_ents.clear();
      for (v_it = entities[dim-1].begin(), v_rit = dm1_ranges.begin(); 
           v_it != entities[dim-1].end(); v_it++, v_rit++) {
        inters = (*v_rit)->intersect(*(*d_rit));
        if (!inters.empty() && inters.size() == (*v_rit)->size()) {
            // non-zero intersection; get possible parent sets
/*
          if (dim == 1) dm1_ents.push_back(*v_it);
          else {
            result = mdbImpl->get_parent_meshsets(*v_it, dm1_ents);
            if (MB_SUCCESS != result)
              return result;
          }
*/
          real_dm1_ents.insert(*v_it);
        }
      }

/*
        // ok, we have possible children; check for real overlap, but only
        // if we're not doing edges
      if (dim == 1) std::copy(dm1_ents.begin(), dm1_ents.end(),
                              mb_range_inserter(real_dm1_ents));
      else {
          // reuse v_it, v_rit here
        for (v_it = entities[dim-1].begin(), v_rit = dm1_ranges.begin();
             v_it != entities[dim-1].end(); v_it++, v_rit++) {

            // if this isn't one of the possible dm1 entities, go on
          if (std::find(dm1_ents.begin(), dm1_ents.end(), *v_it) == dm1_ents.end())
            continue;

            // check for real overlap
          inters = (*v_rit)->intersect(*(*d_rit));
          if (!inters.empty()) real_dm1_ents.insert(*v_it);
        }
      }
  
*/    
        // ok, we have the real children; add parent/child links; reuse v_it again
      for (v_it = real_dm1_ents.begin(); v_it != real_dm1_ents.end(); v_it++) {
        result = mdbImpl->add_parent_child(*d_it, *v_it);
        if (MB_SUCCESS != result)
          return result;
      }
      
        
    } // d_it
  } // dim

  return result;
}

MBErrorCode GeomTopoTool::separate_by_dimension(const MBRange &geom_sets,
                                                 MBRange *entities, MBTag geom_tag) 
{
  MBErrorCode result;
  
  if (0 == geom_tag) {
    
    result = mdbImpl->tag_get_handle(GEOM_DIMENSION_TAG_NAME, geom_tag);
    if (MB_SUCCESS != result)
      return result;
  }

    // get the data for those tags
  std::vector<int> tag_vals(geom_sets.size());
  result = mdbImpl->tag_get_data(geom_tag, geom_sets, &tag_vals[0]);
  if (MB_SUCCESS != result)
    return result;

  MBRange::const_iterator git;
  std::vector<int>::iterator iit;
  
  for (git = geom_sets.begin(), iit = tag_vals.begin(); git != geom_sets.end(); 
       git++, iit++) {
    if (0 <= *iit && 3 >= *iit) 
      entities[*iit].insert(*git);
    else
      assert(false);
  }
  
  return MB_SUCCESS;
}

MBErrorCode GeomTopoTool::construct_vertex_ranges(const MBRange &geom_sets,
                                                   const MBTag verts_tag) 
{
    // construct the vertex range for each entity and put on that tag
  MBRange *temp_verts, temp_elems;
  MBErrorCode result = MB_SUCCESS;
  for (MBRange::const_iterator it = geom_sets.begin(); it != geom_sets.end(); it++) {
      // make the new range
    temp_verts = new MBRange();
    assert(NULL != temp_verts);
    temp_elems.clear();
    
      // get all the elements in the set, recursively
    result = mdbImpl->get_entities_by_handle(*it, temp_elems, true);
    if (MB_SUCCESS != result) 
      return result;
    
      // get all the verts of those elements; use get_adjacencies 'cuz it handles ranges better
    result = mdbImpl->get_adjacencies(temp_elems, 0, false, *temp_verts,
                                      MBInterface::UNION);
    if (MB_SUCCESS != result) 
      return result;

      // store this range as a tag on the entity
    result = mdbImpl->tag_set_data(verts_tag, &(*it), 1, &temp_verts);
    if (MB_SUCCESS != result) 
      return result;
    
  }
  
  return result;
}

  
  
