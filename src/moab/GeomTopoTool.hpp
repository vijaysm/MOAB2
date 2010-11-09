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



#ifndef MOAB_GEOM_TOPO_TOOL_HPP
#define MOAB_GEOM_TOPO_TOOL_HPP

#include "moab/Forward.hpp"
#include "moab/Range.hpp"
#include "moab/OrientedBoxTreeTool.hpp"

#include <map>

namespace moab {

class GeomTopoTool
{
public:
  GeomTopoTool(Interface *impl, bool find_geoments = false);
  ~GeomTopoTool() {}
  
    //! Restore parent/child links between GEOM_TOPO mesh sets
  ErrorCode restore_topology();
  
    //! Store sense of surface relative to volume.
    //!\return MB_MULTIPLE_ENTITIES_FOUND if surface already has a forward volume.
    //!        MB_SUCCESS if successful
    //!        otherwise whatever internal error code occured.
  ErrorCode set_sense( EntityHandle surface,
                         EntityHandle volume,
                         bool forward );

  ErrorCode get_sense( EntityHandle surface,
                         EntityHandle volume,
                         bool& forward );

    //! Store senses and ents for edges
  ErrorCode get_senses(EntityHandle edge,
                       std::vector<EntityHandle> &faces,
                       std::vector<int> &senses);
  
  ErrorCode set_senses(EntityHandle edge,
                       std::vector<EntityHandle> &faces,
                       std::vector<int> &senses);
  
  ErrorCode find_geomsets(Range *ranges = NULL);

  ErrorCode construct_obb_trees(bool make_one_vol = false);

  ErrorCode get_root(EntityHandle vol_or_surf, EntityHandle &root);

  EntityHandle get_one_vol_root();

  OrientedBoxTreeTool *obb_tree() {return &obbTree;}

private:
  Interface *mdbImpl;
  Tag sense2Tag;
  Tag senseNEntsTag, senseNSensesTag;
  Tag geomTag;
  Range geomRanges[4];

  OrientedBoxTreeTool obbTree;
  EntityHandle setOffset;
  std::vector<EntityHandle> rootSets;

  bool contiguous;
  std::map<EntityHandle, EntityHandle>  mapRootSets;
  EntityHandle oneVolRootSet;

    //! compute vertices inclusive and put on tag on sets in geom_sets
  ErrorCode construct_vertex_ranges(const Range &geom_sets,
				      const Tag verts_tag);
  
    //! given a range of geom topology sets, separate by dimension
  ErrorCode separate_by_dimension(const Range &geom_sets,
				    Range *entities, Tag geom_tag = 0);
};

// get the root of the obbtree for a given entity
inline ErrorCode GeomTopoTool::get_root(EntityHandle vol_or_surf, EntityHandle &root) 
{
   if(contiguous)
   {
     unsigned int index = vol_or_surf - setOffset;
     root = (index < rootSets.size() ? rootSets[index] : 0);
   }
   else
      root = mapRootSets[vol_or_surf];
   return (root ? MB_SUCCESS : MB_INDEX_OUT_OF_RANGE);
}

inline EntityHandle GeomTopoTool::get_one_vol_root()
{
  return oneVolRootSet;
}

}
 // namespace moab 

#endif

