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



#ifndef GEOM_TOPO_UTIL_HPP
#define GEOM_TOPO_UTIL_HPP

#include "moab/Forward.hpp"
#include "moab/Range.hpp"

using namespace moab;

class GeomTopoUtil
{
public:
  GeomTopoUtil(Interface *impl, bool find_geoments = false);
  ~GeomTopoUtil() {}
  
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

  ErrorCode find_geomsets(Range *ranges = NULL);

private:
  Interface *mdbImpl;
  Tag sense2Tag;
  Range geomRanges[4];
  
    //! compute vertices inclusive and put on tag on sets in geom_sets
  ErrorCode construct_vertex_ranges(const Range &geom_sets,
                                       const Tag verts_tag);
  
    //! given a range of geom topology sets, separate by dimension
  ErrorCode separate_by_dimension(const Range &geom_sets,
                                     Range *entities, Tag geom_tag = 0);
};

inline GeomTopoUtil::GeomTopoUtil(Interface *impl, 
                                  bool find_geoments) 
        : mdbImpl(impl), sense2Tag(0) 
{
  if (find_geoments) find_geomsets();
}

#endif

