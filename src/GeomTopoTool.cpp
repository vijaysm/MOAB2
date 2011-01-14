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

#include "moab/GeomTopoTool.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/Interface.hpp"
#include "moab/CN.hpp"
#include "Internals.hpp"
#include <assert.h>
#include <iostream>

namespace moab {

// Tag name used for saving sense of faces in volumes.
// We assume that the surface occurs in at most two volumes.
// Code will error out if more than two volumes per surface.
// The tag data is a pair of tag handles, representing the
// forward and reverse volumes, respectively.  If a surface
// is non-manifold in a single volume, the same volume will
// be listed for both the forward and reverse slots.
const char GEOM_SENSE_2_TAG_NAME[] = "GEOM_SENSE_2";

const char GEOM_SENSE_N_ENTS_TAG_NAME[] = "GEOM_SENSE_N_ENTS";
const char GEOM_SENSE_N_SENSES_TAG_NAME[] = "GEOM_SENSE_N_SENSES";

GeomTopoTool::GeomTopoTool(Interface *impl, bool find_geoments)
 : mdbImpl(impl), 
   sense2Tag(0), 
   senseNEntsTag(0), 
   senseNSensesTag(0), 
   obbTree(impl), 
   contiguous(true),
   oneVolRootSet(0)
{
   ErrorCode result = mdbImpl->tag_create(GEOM_DIMENSION_TAG_NAME, 4,
         MB_TAG_SPARSE, geomTag, NULL);
   if (MB_SUCCESS != result && MB_ALREADY_ALLOCATED != result) {
      std::cerr << "Error: Failed to create geometry dimension tag."
            << std::endl;
   }

   if (find_geoments)
      find_geomsets();
}

ErrorCode GeomTopoTool::set_sense(EntityHandle surface, EntityHandle volume,
      bool forward) {
   ErrorCode rval;
   if (!sense2Tag) {
      rval = mdbImpl->tag_create(GEOM_SENSE_2_TAG_NAME, 2
            * sizeof(EntityHandle), MB_TAG_SPARSE, MB_TYPE_HANDLE, sense2Tag,
            0, true);
      if (MB_SUCCESS != rval && (MB_ALREADY_ALLOCATED != rval || !sense2Tag))
         return rval;
   }

   EntityHandle sense_data[2] = { 0, 0 };
   rval = mdbImpl->tag_get_data(sense2Tag, &surface, 1, sense_data);
   if (MB_TAG_NOT_FOUND != rval && MB_SUCCESS != rval)
      return MB_FAILURE;

   if (sense_data[!forward] == volume)
      return MB_SUCCESS;
   else if (sense_data[!forward])
      return MB_MULTIPLE_ENTITIES_FOUND;

   sense_data[!forward] = volume;
   return mdbImpl->tag_set_data(sense2Tag, &surface, 1, sense_data);
}

ErrorCode GeomTopoTool::set_senses(EntityHandle edge,
      std::vector<EntityHandle> &faces, std::vector<int> &senses) {
   ErrorCode rval;
   if (!senseNEntsTag) {
      rval = mdbImpl->tag_create_variable_length(GEOM_SENSE_N_ENTS_TAG_NAME,
            MB_TAG_SPARSE, MB_TYPE_HANDLE, senseNEntsTag);
      if (MB_SUCCESS != rval && (MB_ALREADY_ALLOCATED != rval || !senseNEntsTag))
         return rval;

      rval = mdbImpl->tag_create_variable_length(GEOM_SENSE_N_SENSES_TAG_NAME,
            MB_TAG_SPARSE, MB_TYPE_INTEGER, senseNSensesTag);
      if (MB_SUCCESS != rval && (MB_ALREADY_ALLOCATED != rval || !senseNSensesTag))
         return rval;
   }

   int dum_size = faces.size() * sizeof(EntityHandle);
   void *dum_ptr = &faces[0];
   rval = mdbImpl->tag_set_data(senseNEntsTag, &edge, 1, &dum_ptr, &dum_size);
   if (MB_SUCCESS != rval)
      return rval;

   dum_ptr = &senses[0];
   dum_size = faces.size() * sizeof(int);
   rval = mdbImpl->tag_set_data(senseNSensesTag, &edge, 1, &dum_ptr, &dum_size);
   if (MB_SUCCESS != rval)
      return rval;

   return rval;
}

ErrorCode GeomTopoTool::get_sense(EntityHandle surface, EntityHandle volume,
      bool& forward) {
   ErrorCode rval;
   if (!sense2Tag) {
      rval = mdbImpl->tag_get_handle(GEOM_SENSE_2_TAG_NAME, sense2Tag);
      if (MB_SUCCESS != rval) {
         sense2Tag = 0;
         return MB_FAILURE;
      }
   }

   EntityHandle sense_data[2] = { 0, 0 };
   rval = mdbImpl->tag_get_data(sense2Tag, &surface, 1, sense_data);
   if (MB_SUCCESS != rval)
      return rval;

   if (sense_data[0] == volume)
      forward = true;
   else if (sense_data[1] == volume)
      forward = false;
   else
      return MB_ENTITY_NOT_FOUND;

   return MB_SUCCESS;
}

ErrorCode GeomTopoTool::get_senses(EntityHandle edge,
      std::vector<EntityHandle> &faces, std::vector<int> &senses) {
   ErrorCode rval;
   if (!senseNEntsTag) {
      rval = mdbImpl->tag_get_handle(GEOM_SENSE_N_ENTS_TAG_NAME, senseNEntsTag);
      if (MB_SUCCESS != rval) {
         senseNEntsTag = 0;
         return MB_FAILURE;
      }
   }

   if (!senseNSensesTag) {
      rval = mdbImpl->tag_get_handle(GEOM_SENSE_N_SENSES_TAG_NAME,
            senseNSensesTag);
      if (MB_SUCCESS != rval) {
         senseNSensesTag = 0;
         return MB_FAILURE;
      }
   }

   const void *dum_ptr;
   int num_ents;
   rval = mdbImpl->tag_get_data(senseNEntsTag, &edge, 1, &dum_ptr, &num_ents);
   if (MB_SUCCESS != rval)
      return rval;

   faces.clear();
   num_ents /= sizeof(EntityHandle);
   const EntityHandle *ents_data = static_cast<const EntityHandle*> (dum_ptr);
   std::copy(ents_data, ents_data + num_ents, std::back_inserter(faces));

   rval = mdbImpl->tag_get_data(senseNSensesTag, &edge, 1, &dum_ptr, &num_ents);
   if (MB_SUCCESS != rval)
      return rval;

   senses.clear();
   num_ents /= sizeof(int);
   const int *senses_data = static_cast<const int*> (dum_ptr);
   std::copy(senses_data, senses_data + num_ents, std::back_inserter(senses));

   return MB_SUCCESS;
}

ErrorCode GeomTopoTool::find_geomsets(Range *ranges) {
   // get all sets with this tag
   Range geom_sets;
   ErrorCode result = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET,
         &geomTag, NULL, 1, geom_sets);
   if (MB_SUCCESS != result || geom_sets.empty())
      return result;

   result = separate_by_dimension(geom_sets, geomRanges, geomTag);
   if (MB_SUCCESS != result)
      return result;

   if (ranges) {
      for (int i = 0; i < 4; i++)
         ranges[i] = geomRanges[i];
   }

   return MB_SUCCESS;
}

ErrorCode GeomTopoTool::construct_obb_trees(bool make_one_vol) {
   ErrorCode rval;

   // get all surfaces and volumes
   Range surfs, vols, vol_trees;
   const int three = 3;
   const void* const three_val[] = { &three };
   rval = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &geomTag,
         three_val, 1, vols);
   if (MB_SUCCESS != rval)
      return rval;

   const int two = 2;
   const void* const two_val[] = { &two };
   rval = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &geomTag,
         two_val, 1, surfs);
   if (MB_SUCCESS != rval)
      return rval;

   if (vols.empty() && !surfs.empty()) {
      setOffset = surfs.front();
   } else if (!vols.empty() && surfs.empty()) {
      setOffset = vols.front();
   } else {
      setOffset = (surfs.front() < vols.front() ? surfs.front() : vols.front());
   }
   EntityHandle minSet = setOffset;
   EntityHandle maxSet = setOffset;
   Range::iterator it;
   for (it = surfs.begin(); it != surfs.end(); ++it) {
      EntityHandle sf = *it;
      if (sf > maxSet)
         maxSet = sf;
      if (sf < minSet)
         minSet = sf;
   }
   for (it = vols.begin(); it != vols.end(); ++it) {
      EntityHandle sf = *it;
      if (sf > maxSet)
         maxSet = sf;
      if (sf < minSet)
         minSet = sf;
   }
   if (surfs.size() + vols.size() == maxSet - minSet + 1)
      contiguous = true;
   else
      contiguous = false; // need map arrangements
   // for surface
   EntityHandle root;
   if (contiguous)
      rootSets.resize(surfs.size() + vols.size());
   for (Range::iterator i = surfs.begin(); i != surfs.end(); ++i) {
      Range tris;
      rval = mdbImpl->get_entities_by_dimension(*i, 2, tris);
      if (MB_SUCCESS != rval)
         return rval;

      if (tris.empty()) {
         std::cerr << "WARNING: Surface has no facets." << std::endl;
      }

      rval = obbTree.build(tris, root);
      if (MB_SUCCESS != rval)
         return rval;

      rval = mdbImpl->add_entities(root, &*i, 1);
      if (MB_SUCCESS != rval)
         return rval;

      //surfRootSets[*i - surfOffset] = root;
      if (contiguous)
         rootSets[*i - setOffset] = root;
      else
         mapRootSets[*i] = root;
   }

   // for volumes
   Range trees;
   for (Range::iterator i = vols.begin(); i != vols.end(); ++i) {
      // get all surfaces in volume
      Range tmp_surfs;
      rval = mdbImpl->get_child_meshsets(*i, tmp_surfs);
      if (MB_SUCCESS != rval)
         return rval;

      // get OBB trees for each surface
      if (!make_one_vol)
         trees.clear();
      for (Range::iterator j = tmp_surfs.begin(); j != tmp_surfs.end(); ++j) {
         rval = get_root(*j, root);
         if (MB_SUCCESS != rval || !root)
            return MB_FAILURE;
         trees.insert(root);
      }

      // build OBB tree for volume
      if (!make_one_vol) {
         rval = obbTree.join_trees(trees, root);
         if (MB_SUCCESS != rval)
            return rval;
         if (contiguous)
            rootSets[*i - setOffset] = root;
         else
            mapRootSets[*i] = root;
      }
   }

   // build OBB tree for volume
   if (make_one_vol) {
      rval = obbTree.join_trees(trees, root);
      if (MB_SUCCESS != rval)
         return rval;
      oneVolRootSet = root;
   }

   return rval;
}

//! Restore parent/child links between GEOM_TOPO mesh sets
ErrorCode GeomTopoTool::restore_topology() {

   // look for geometric topology sets and restore parent/child links between them
   // algorithm:
   // - for each entity of dimension d=D-1..0:
   //   . get d-dimensional entity in entity
   //   . get all (d+1)-dim adjs to that entity
   //   . for each geom entity if dim d+1, if it contains any of the ents,
   //     add it to list of parents
   //   . make parent/child links with parents

   // get the geom topology tag
   Tag geom_tag;
   ErrorCode result = mdbImpl->tag_create(GEOM_DIMENSION_TAG_NAME, 4,
         MB_TAG_SPARSE, geom_tag, NULL);
   if (MB_SUCCESS != result && (MB_ALREADY_ALLOCATED != result || !geom_tag))
      return result;

   // get all sets with this tag
   Range geom_sets;
   result = mdbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &geom_tag,
         NULL, 1, geom_sets);
   if (MB_SUCCESS != result || geom_sets.empty())
      return result;

   Range entities[4];
   result = separate_by_dimension(geom_sets, entities, geom_tag);
   if (MB_SUCCESS != result)
      return result;

   std::vector<EntityHandle> parents;
   Range tmp_parents;

   // loop over dimensions
   for (int dim = 2; dim >= 0; dim--) {
      // mark entities of next higher dimension with their owners; regenerate tag
      // each dimension so prev dim's tag data goes away
      Tag owner_tag;
      EntityHandle dum_val = 0;
      result = mdbImpl->tag_create("__owner_tag", sizeof(EntityHandle),
            MB_TAG_DENSE, MB_TYPE_HANDLE, owner_tag, &dum_val);
      if (MB_SUCCESS != result && (MB_ALREADY_ALLOCATED != result || !owner_tag))
         continue;
      Range dp1ents;
      std::vector<EntityHandle> owners;
      for (Range::iterator rit = entities[dim + 1].begin(); rit != entities[dim
            + 1].end(); rit++) {
         dp1ents.clear();
         result = mdbImpl->get_entities_by_dimension(*rit, dim + 1, dp1ents);
         if (MB_SUCCESS != result)
            continue;
         owners.resize(dp1ents.size());
         std::fill(owners.begin(), owners.end(), *rit);
         result = mdbImpl->tag_set_data(owner_tag, dp1ents, &owners[0]);
         if (MB_SUCCESS != result)
            continue;
      }

      for (Range::iterator d_it = entities[dim].begin(); d_it
            != entities[dim].end(); d_it++) {
         Range dents;
         result = mdbImpl->get_entities_by_dimension(*d_it, dim, dents);
         if (MB_SUCCESS != result)
            continue;
         if (dents.empty())
            continue;

         // get (d+1)-dimensional adjs
         dp1ents.clear();
         result = mdbImpl->get_adjacencies(&(*dents.begin()), 1, dim + 1,
               false, dp1ents);
         if (MB_SUCCESS != result || dp1ents.empty())
            continue;

         // get owner tags
         parents.resize(dp1ents.size());
         result = mdbImpl->tag_get_data(owner_tag, dp1ents, &parents[0]);
         assert(MB_TAG_NOT_FOUND != result);
         if (MB_SUCCESS != result)
            continue;

         // compress to a range to remove duplicates
         tmp_parents.clear();
         std::copy(parents.begin(), parents.end(), range_inserter(tmp_parents));
         for (Range::iterator pit = tmp_parents.begin(); pit
               != tmp_parents.end(); pit++) {
            result = mdbImpl->add_parent_child(*pit, *d_it);
            if (MB_SUCCESS != result)
               return result;
         }

         // store surface senses
         if (dim != 2)
            continue;
         const EntityHandle *conn3, *conn2;
         int len3, len2, err, num, sense, offset;
         for (size_t i = 0; i < parents.size(); ++i) {
            result = mdbImpl->get_connectivity(dp1ents[i], conn3, len3, true);
            if (MB_SUCCESS != result)
               return result;
            result
                  = mdbImpl->get_connectivity(dents.front(), conn2, len2, true);
            if (MB_SUCCESS != result)
               return result;
            assert(len2 <= 4);
            err = CN::SideNumber(TYPE_FROM_HANDLE(dp1ents[i]), conn3, conn2,
                  len2, dim, num, sense, offset);
            if (err)
               return MB_FAILURE;

            result = set_sense(*d_it, parents[i], sense == 1);
            if (MB_MULTIPLE_ENTITIES_FOUND == result) {
               std::cerr
                     << "Warning: Multiple volumes use surface with same sense."
                     << std::endl << "         Some geometric sense data lost."
                     << std::endl;
            } else if (MB_SUCCESS != result) {
               return result;
            }
         }
      }

      // now delete owner tag on this dimension, automatically removes tag data
      result = mdbImpl->tag_delete(owner_tag);
      if (MB_SUCCESS != result)
         return result;

   } // dim

   return result;
}

ErrorCode GeomTopoTool::separate_by_dimension(const Range &geom_sets,
      Range *entities, Tag geom_tag) {
   ErrorCode result;

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

   Range::const_iterator git;
   std::vector<int>::iterator iit;

   for (git = geom_sets.begin(), iit = tag_vals.begin(); git != geom_sets.end(); git++, iit++) {
      if (0 <= *iit && 3 >= *iit)
         entities[*iit].insert(*git);
      else {
         // assert(false);
         // do nothing for now
      }
   }

   return MB_SUCCESS;
}

ErrorCode GeomTopoTool::construct_vertex_ranges(const Range &geom_sets,
      const Tag verts_tag) {
   // construct the vertex range for each entity and put on that tag
   Range *temp_verts, temp_elems;
   ErrorCode result = MB_SUCCESS;
   for (Range::const_iterator it = geom_sets.begin(); it != geom_sets.end(); it++) {
      // make the new range
      temp_verts = new Range();
      assert(NULL != temp_verts);
      temp_elems.clear();

      // get all the elements in the set, recursively
      result = mdbImpl->get_entities_by_handle(*it, temp_elems, true);
      if (MB_SUCCESS != result)
         return result;

      // get all the verts of those elements; use get_adjacencies 'cuz it handles ranges better
      result = mdbImpl->get_adjacencies(temp_elems, 0, false, *temp_verts,
            Interface::UNION);
      if (MB_SUCCESS != result)
         return result;

      // store this range as a tag on the entity
      result = mdbImpl->tag_set_data(verts_tag, &(*it), 1, &temp_verts);
      if (MB_SUCCESS != result)
         return result;

   }

   return result;
}

} // namespace moab


