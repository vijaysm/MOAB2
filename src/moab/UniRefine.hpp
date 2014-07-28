#ifndef UNI_REFINE_HPP
#define UNI_REFINE_HPP

#include "moab/Range.hpp"
#include "moab/HalfFacetRep.hpp"

namespace moab
{

#define MAX_LEVEL 3
#define MAX_VERTS 100
#define MAX_CHILDRENS 100
#define MAX_LOC 19
#define MAX_CONN 8
  
  enum METHOD{
    MULTI_LEVEL = 0,
    ALL_LEVELS
  };
  
  class UniRefine: public  HalfFacetRep{
  protected: 
    //Interface * mb;
    //ArrayHalfFacetMDS * ahf;
    
  public:
    UniRefine(Interface *mesh_in) : ArrayHalfFacetMDS(mesh_in) {}
    
    ~UniRefine() {}
    
    
    //ArrayHalfFacetMDS ahf(Interface * mb);
    
    //ErrorCode initialize();
    //  ErrorCode deinitialize();
    
    struct LevelTemplates{
      short int num_new_verts_per_edge; // Number of new vertices on edge 
      short int num_new_verts_per_face; // Number of new vertices on face
      short int num_new_verts_per_cell; // Number of new vertices in cell
      short int total_new_verts;
      
      short int num_new_ents;
      
      int vert_indices[2];
      double vert_params[MAX_VERTS][2];
      
      int new_entsConn[MAX_CHILDRENS][MAX_CONN];
      
      int vert_identify[MAX_LOC][MAX_VERTS];
    };
    
    static const LevelTemplates UnirefPatterns[3][MAX_LEVEL];
    
    int get_index_from_type(EntityHandle ent);

    ErrorCode get_total_new_verts_and_ents(EntityHandle edg, int nedges, int level, METHOD method, int *newverts, int *newsubents);
    
    ErrorCode get_total_new_verts_and_ents(EntityHandle face, int nedges, int nfaces, int level, METHOD method, int *newverts, int *newsubents);
    
    ErrorCode get_total_new_verts_and_ents(EntityHandle cell, int nedges, int nfaces, int ncells, int level, METHOD method, int *newverts, int *newsubents);
    
    ErrorCode uniform_refinement( Range &ents, int level, METHOD method);

    ErrorCode uniform_refinement_mixed_2d(Range &verts, Range &edges, Range &faces, int level, METHOD method);

    ErrorCode uniform_refinement_mixed_single_level_2d(Range &edges, Range &faces, int level, EntityHandle vert_bnds[2], std::vector<double *> coords, EntityHandle edge_bnds[2], EntityHandle *econnect, EntityHandle face_bnds[2], EntityHandle *fconnect );
    
    ErrorCode uniform_refinement_mixed_3d(Range &verts, Range &edges, Range &faces, Range &cells, int level, METHOD method);

    ErrorCode uniform_refinement_mixed_single_level_3d(Range &edges, Range &faces, Range &cells, int level, EntityHandle vert_bnds[2], std::vector<double *> coords, EntityHandle edge_bnds[2], EntityHandle *econnect, EntityHandle face_bnds[2], EntityHandle *fconnect , EntityHandle cell_bnds[2], EntityHandle *cconnect);

  };
} //name space moab
#endif
