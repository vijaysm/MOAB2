

#ifndef GEOM_TOPO_TOOL_HPP
#define GEOM_TOPO_TOOL_HPP

#include "MBInterface.hpp"

class GeomTopoTool
{
public:
  GeomTopoTool(MBInterface *impl) : mdbImpl(impl) {}
  
  ~GeomTopoTool() {}
  
    //! Restore parent/child links between GEOM_TOPO mesh sets
  MBErrorCode restore_topology();

private:
  MBInterface *mdbImpl;
  
    //! compute vertices inclusive and put on tag on sets in geom_sets
  MBErrorCode construct_vertex_ranges(const MBRange &geom_sets,
                                       const MBTag verts_tag);
  
    //! given a range of geom topology sets, separate by dimension
  MBErrorCode separate_by_dimension(const MBRange &geom_sets,
                                     MBRange *entities, MBTag geom_tag = 0);
  
};


#endif

