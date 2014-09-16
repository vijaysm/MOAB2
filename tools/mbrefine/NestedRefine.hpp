#ifndef NESTED_REFINE_HPP
#define NESTED_REFINE_HPP

#include "moab/Range.hpp"
#include "moab/HalfFacetRep.hpp"
#include "RefTemplates.hpp"

namespace moab
{
  
  enum METHOD{
    MULTI_LEVEL = 0,
    ALL_LEVELS
  };
  
  class Core;
  class HalfFacetRep;

  class NestedRefine{
  protected: 
    Core * mb;
    HalfFacetRep * ahf;
    
  public:
    NestedRefine(Core *mesh_in);
    
    ~NestedRefine() {}
    
    //User interface functions
    //1st class: Basic functionalities

    /* Generate a hierarchy of meshes from the input mesh.
     * level is a sequence of degrees used for each level in the hierarchy. There is no upper bound on the number of levels.
     * However, in practice, maximum of 5 or 6 levels are used.
     */
    ErrorCode generate_mesh_hierarchy(int *level_seq);
    //ErrorCode get_coordinates();
    //ErrorCode get_connectivities();
    //ErrorCode get_adjacencies();
    //ErrorCode get_tag();
    //ErrorCode set_tag();

    //2nd class: Interlevel
    //ErrorCode interpolate_data();
    //ErrorCode restriction_operator();
    //ErrorCode prolongation_operator();

  protected:

    ErrorCode estimate_hm_storage(int *level_seq, int num_level, int hmest[][4]);
    ErrorCode create_hm_storage(int num_level, int hmest[][4]);
    ErrorCode generate_hm(int *level_seq);

    
    struct level_memory_pointer{
      EntityHandle start_vertex, start_edge, start_face, start_cell;
      int num_verts, num_edges, num_faces, num_cells;
      std::vector<double *> coordinates;
      EntityHandle *edge_conn, *face_conn, *cell_conn;
    };

    level_memory_pointer access_lmem[20];


  };
} //name space moab
#endif
