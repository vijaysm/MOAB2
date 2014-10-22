#ifndef NESTED_REFINE_HPP
#define NESTED_REFINE_HPP

#include "moab/Range.hpp"
#include "moab/HalfFacetRep.hpp"
//#include "Templates.hpp"

namespace moab
{
  
#define MAX_DEGREE 3
#define MAX_VERTS 100
#define MAX_CHILDRENS 100
#define MAX_HF 12
#define MAX_CONN 8
#define MAX_VHF 20
#define MAX_LEVELS 10


  //class Core;
  //class HalfFacetRep;

  class NestedRefine: public HalfFacetRep
  {
    
  public:
    NestedRefine(Core *impl);
    
    ~NestedRefine() {}
    
    //User interface functions
    //1st class: Basic functionalities

    /* Generate a hierarchy of meshes from the input mesh.
     * level is a sequence of degrees used for each level in the hierarchy. There is no upper bound on the number of levels.
     * However, in practice, maximum of 5 or 6 levels are used.
     */
    ErrorCode generate_mesh_hierarchy(int *level_degrees, int num_level, EntityHandle *hm_set);
    ErrorCode get_connectivity(EntityHandle ent, int level, std::vector<EntityHandle> &conn);
    ErrorCode get_coordinates(std::vector<EntityHandle> verts, int num_verts,  int cur_level, double *coords);

    //ErrorCode get_adjacencies(); // Called directly from the AHF class
    //ErrorCode tag_get_data(); // Get meta data for the new levels
    //ErrorCode tag_set_data(); // Set meta data for the new levels

  protected:

    // Refinement Patterns
    struct refPatterns{
      short int nv_edge; // Number of new vertices on edge
      short int nv_face; // Number of new vertices on face, does not include those on edge
      short int nv_cell; // Number of new vertices in cell
      short int total_new_verts; // Total number of new vertices per entity
      short int total_new_ents; // Total number of new child entities

      int vert_index_bnds[2]; //Lower and upper indices of the new vertices
      double vert_nat_coord[MAX_VERTS][3]; //Natural coordinates of the new vertices
      int ents_conn[MAX_CHILDRENS][MAX_CONN]; //Connectivity of the new entities

      int v2hf[MAX_VERTS][2]; //Vertex to half-facet map of the new vertices
      int ents_opphfs[MAX_CHILDRENS][2*MAX_CONN]; // Opposite half-facet map of the new entities

      int vert_on_edges[MAX_HF][MAX_VHF]; //Helper: storing the local ids of vertices on each local edge
      int vert_on_faces[MAX_HF][MAX_VHF]; // Helper: storing local ids of verts on each local face, doesnt include those on edges of the face
      int ents_on_pent[MAX_HF][MAX_CHILDRENS]; //Helper: stores child half-facets incident on parent half-facet. First column contain the number of such children
    };

    static const refPatterns refTemplates[9][MAX_DEGREE];

    int get_index_from_degree(int degree);

    // Octahedron Tessellation
    struct tessellate_octahedron{
      int combination;
      int diag_conn[2]; //Connectivity of the diagonal
      int tet_conn[4][4]; //Connectivity of the tets
      int tet_opphfs[4][8]; //Opposite half-face map between the children tets
      int olfid_to_tlfid[8][2]; // Map each local face of the parent octahedron to the corresponding local face of the incident child tet
      int tlfid_to_olfid[4][4]; // Map each local face of each child tet to the corresponding local face of the oct
    };

    static const tessellate_octahedron oct_tessellation[3];

    //Permutation matrix
    struct pmat{
      short int num_comb;
      int mat[MAX_HF][MAX_HF];
    };

    static const pmat permute_matrix[2];

    // HM Storage Helper
    struct level_memory{
      int num_verts, num_edges, num_faces, num_cells;
      EntityHandle start_vertex, start_edge, start_face, start_cell;
      std::vector<double *> coordinates;
      EntityHandle *edge_conn, *face_conn, *cell_conn;
      Range verts, edges, faces, cells;
    };

    level_memory level_mesh[MAX_LEVELS];

    //Basic Functions
    //Generate HM
    ErrorCode generate_hm(int *level_degrees, int num_level, int *hmest, EntityHandle *hm_set);

    //Estimate and create storage for the levels
    ErrorCode estimate_hm_storage(int *level_degrees, int num_level, int *hmest);
    ErrorCode create_hm_storage_single_level(EntityHandle *set, int cur_level, int *estL);
    ErrorCode add_entities(EntityHandle *hm_set, int num_level, int *estL);

    //Construct the hierarchical mesh: 1D, 2D, 3D
    ErrorCode construct_hm_entities(int cur_level, int deg);
    ErrorCode construct_hm_1D(int cur_level, int deg);
    ErrorCode construct_hm_2D(int cur_level, int deg);
    ErrorCode construct_hm_3D(int cur_level, int deg);

    ErrorCode subdivide_cells(EntityType type, std::vector<EntityHandle> conn, int cur_level, int deg, EntityHandle  *vbuffer, int *count_ents);
    ErrorCode subdivide_tets(std::vector<EntityHandle> conn, int cur_level, int deg, EntityHandle *vbuffer, int *count_ents);

    // Helper functions
    ErrorCode copy_vertices_from_prev_level(int cur_level);
    ErrorCode update_tracking_verts(EntityHandle cidl, int cur_level, int deg, std::vector<EntityHandle> trackvertsC_edg, std::vector<EntityHandle> trackvertsC_face, EntityHandle *vbuffer);
    ErrorCode match_and_reorder_vertices(EntityType type, int cur_level, int deg, EntityHandle cell, int lfid, EntityHandle sib_cell, int sib_lfid, int *id_sib);
    int find_shortest_diagonal_octahedron( double *coords);
    int get_local_vid(EntityHandle vid, EntityHandle ent, int level);
    ErrorCode print_tags_1D(int level);

    // Coordinates
    ErrorCode compute_coordinates(int cur_level, int deg, EntityType type, EntityHandle *vbuffer, int vtotal, double *corner_coords);

    // Update the ahf maps

    ErrorCode update_local_ahf(int deg, EntityType type, EntityHandle *vbuffer, EntityHandle *ent_buffer, int etotal);

    ErrorCode update_local_ahf(int cur_level, int deg, std::vector<int> nents_flag, std::vector<int> idx_buffer);

    ErrorCode update_global_ahf(EntityType type, int cur_level, int deg);

    ErrorCode update_global_ahf_1D(int cur_level, int deg);

    ErrorCode update_global_ahf_2D(int cur_level, int deg);

    ErrorCode update_global_ahf_3D(int cur_level, int deg);

    ErrorCode get_sibling_tag(EntityType type, EntityHandle ent, EntityHandle *sib_entids, int *sib_lids);

    ErrorCode set_sibling_tag(EntityType type, EntityHandle ent, EntityHandle *set_entids, int *set_lids);

    ErrorCode get_incident_tag(EntityType type, EntityHandle vid, EntityHandle *inci_entid, int *inci_lid);

    ErrorCode set_incident_tag(EntityType type, EntityHandle vid, EntityHandle *set_entid, int *set_lid);


  };
} //name space moab
#endif
