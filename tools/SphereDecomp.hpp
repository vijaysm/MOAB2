#ifndef SPHERE_DECOMP_HPP
#define SPHERE_DECOMP_HPP

#include "MBInterface.hpp"

class SphereDecomp 
{
public:
  SphereDecomp(MBInterface *impl);

  MBErrorCode build_sphere_mesh(const char *sphere_radii_tag_name, 
                                MBEntityHandle *hex_set = NULL);
  
private:

    //! compute subdivision vertices on entities of specified dimension
  MBErrorCode compute_nodes(const int dim);

    //! subdivide tets based on subdiv vertices, returning in lists according
    //! to whether they're inside or outside spheres
  MBErrorCode build_hexes(std::vector<MBEntityHandle> &sphere_hexes,
                          std::vector<MBEntityHandle> &interstic_hexes);
  
    //! subdivide an individual tet
  MBErrorCode subdivide_tet(MBEntityHandle tet, 
                            std::vector<MBEntityHandle> &sphere_hexes,
                            std::vector<MBEntityHandle> &interstic_hexes);
  
    //! retrieve the subdivision vertices for a given entity in a given tet,
    //! placing them in the array oriented wrt the tet
  MBErrorCode retrieve_subdiv_verts(MBEntityHandle tet, MBEntityHandle this_ent,
                                    const MBEntityHandle *tet_conn,
                                    const int dim, MBEntityHandle *subdiv_verts);
  
    //! tag used to hold sphere radii (assigned to vertices)
  MBTag sphereRadiiTag;

    //! used to store subdiv vertices for a given d>0 entity
  MBTag subdivVerticesTag;
  
    //! MOAB interface ptr
  MBInterface *mbImpl;
  
};
#endif
