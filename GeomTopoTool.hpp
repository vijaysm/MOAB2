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



#ifndef GEOM_TOPO_TOOL_HPP
#define GEOM_TOPO_TOOL_HPP

#include "MBForward.hpp"

class GeomTopoTool
{
public:
  GeomTopoTool(MBInterface *impl) : mdbImpl(impl), sense2Tag(0) {}
  
  ~GeomTopoTool() {}
  
    //! Restore parent/child links between GEOM_TOPO mesh sets
  MBErrorCode restore_topology();
  
    //! Store sense of surface relative to volume.
    //!\return MB_MULTIPLE_ENTITIES_FOUND if surface already has a forward volume.
    //!        MB_SUCCESS if successful
    //!        otherwise whatever internal error code occured.
  MBErrorCode set_sense( MBEntityHandle surface,
                         MBEntityHandle volume,
                         bool forward );

  MBErrorCode get_sense( MBEntityHandle surface,
                         MBEntityHandle volume,
                         bool& forward );

private:
  MBInterface *mdbImpl;
  MBTag sense2Tag;
  
    //! compute vertices inclusive and put on tag on sets in geom_sets
  MBErrorCode construct_vertex_ranges(const MBRange &geom_sets,
                                       const MBTag verts_tag);
  
    //! given a range of geom topology sets, separate by dimension
  MBErrorCode separate_by_dimension(const MBRange &geom_sets,
                                     MBRange *entities, MBTag geom_tag = 0);
  
};


#endif

