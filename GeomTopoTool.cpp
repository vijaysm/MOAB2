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
#include "MBInterface.hpp"
#include <assert.h>

    //! Restore parent/child links between GEOM_TOPO mesh sets
MBErrorCode GeomTopoTool::restore_topology() 
{
  
    // look for geometric topology sets and restore parent/child links between them
    // algorithm:
    // - for each entity of dimension d=D-1..0:
    //   . get d-dimensional entity in entity
    //   . get all (d+1)-dim adjs to that entity
    //   . for each geom entity if dim d+1, if it contains any of the ents,
    //     add it to list of parents
    //   . make parent/child links with parents

    // get the geom topology tag
  MBTag geom_tag;
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

  MBRange entities[4];
  result = separate_by_dimension(geom_sets, entities, geom_tag);
  if (MB_SUCCESS != result)
    return result;

  std::vector<MBEntityHandle> parents;
  MBRange tmp_parents;
  
    // loop over dimensions
  for (int dim = 2; dim >= 0; dim--) {
      // mark entities of next higher dimension with their owners; regenerate tag
      // each dimension so prev dim's tag data goes away
    MBTag owner_tag;
    MBEntityHandle dum_val = 0;
    result = mdbImpl->tag_create("__owner_tag", sizeof(MBEntityHandle), MB_TAG_DENSE,
                                 MB_TYPE_HANDLE, owner_tag, &dum_val);
    if (MB_SUCCESS != result) continue;
    MBRange dp1ents;
    std::vector<MBEntityHandle> owners;
    for (MBRange::iterator rit = entities[dim+1].begin(); rit != entities[dim+1].end(); rit++) {
      dp1ents.clear();
      result = mdbImpl->get_entities_by_dimension(*rit, dim+1, dp1ents);
      if (MB_SUCCESS != result) continue;
      owners.resize(dp1ents.size());
      std::fill(owners.begin(), owners.end(), *rit);
      result = mdbImpl->tag_set_data(owner_tag, dp1ents, &owners[0]);
      if (MB_SUCCESS != result) continue;
    }
    
    for (MBRange::iterator d_it = entities[dim].begin(); 
         d_it != entities[dim].end(); d_it++) {
      MBRange dents;
      result = mdbImpl->get_entities_by_dimension(*d_it, dim, dents);
      if (MB_SUCCESS != result) continue;
      if (dents.empty()) continue;
      
        // get (d+1)-dimensional adjs
      dp1ents.clear();
      result = mdbImpl->get_adjacencies(&(*dents.begin()), 1, dim+1, 
                                        false, dp1ents);
      if (MB_SUCCESS != result || dp1ents.empty()) continue;

        // get owner tags
      parents.resize(dp1ents.size());
      result = mdbImpl->tag_get_data(owner_tag, dp1ents, &parents[0]);
      assert(MB_TAG_NOT_FOUND != result);
      if (MB_SUCCESS != result) continue;
      
        // compress to a range to remove duplicates
      tmp_parents.clear();
      std::copy(parents.begin(), parents.end(), mb_range_inserter(tmp_parents));
      for (MBRange::iterator pit = tmp_parents.begin(); pit != tmp_parents.end(); pit++) {
        result = mdbImpl->add_parent_child(*pit, *d_it);
        if (MB_SUCCESS != result) return result;
      }
    }
    
      // now delete owner tag on this dimension, automatically removes tag data
    result = mdbImpl->tag_delete(owner_tag);
    if (MB_SUCCESS != result) return result;
    
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
    else {
        // assert(false);
        // do nothing for now
    }
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

  
  
