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


#ifndef MOAB_HALF_FACET_REP_HPP
#define MOAB_HALF_FACET_REP_HPP

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"

namespace moab {

/*! 
 *  \brief   HalfFacetRep class implements the Array-Based Half-Facet(AHF) Mesh data structure on top of MOAB.
 *  \        It is based on a generalized notion of a half-facet derived from that a half-edge/half-face data structure for 2D/3D.
 *  \        The core construct of AHF are its two maps: sibling-half-facets(SIBHFS) and vertex-to-half-facet(V2HF)
 *  \        1. SIBHFS: Maps every half-facet in the mesh to its sibling half-facets,
 *  \        2. V2HF: Maps each vertex to an incident half-facet
 *  \        Using these two maps, a range of adjacency queries is performed. The maps are stored in dense tags over entities and vertices.
 *  \
 *  \        Current support for adjacency functions:
 *  \        1. upward-incidence queries: vertex -> edge, edge -> faces, edge -> cells, face ->cells
 *  \        2. neighborhood (same-dimensional) adjacency queries: edge -> edges, face -> faces, cell -> cells, etc.
 *  \ 
 *  \        Mesh types supported: 
 *  \        1D(edges), 2D(triangles, quads), 3D(tet, pyramid, prism, hex), Mixed dimensional meshes
 *  \        NOT SUPPORTED: Meshes with mixed entity types of same dimension. Ex. a volume mesh with both tets and prisms.
 *  \
 */ 

  //! ENUM for the type of input mesh.
  enum MESHTYPE{
   CURVE = 0, //Homogeneous curve mesh
   SURFACE, // Homogeneous surface mesh
   SURFACE_MIXED, // Mixed surface with embedded curves
   VOLUME, // Homogeneous volume mesh
   VOLUME_MIXED_1, // Volume mesh with embedded curves
   VOLUME_MIXED_2, // Volume mesh with embedded surface
   VOLUME_MIXED //Volume mesh with embedded curves and surfaces
  };

  class Core;

  class HalfFacetRep{
       
  public:
    
    HalfFacetRep(Core *impl);
    
    ~HalfFacetRep();

    // User interface functions

    //! Creates all the necessary tags to store the maps. Constructs the sibling-half-facet and vertex-to-incident-half-facet maps for each dimension present in the input.

    ErrorCode initialize();

    //! Deletes all the created tags.

    ErrorCode deinitialize();

    //! Prints the tag values.
    ErrorCode print_tags();
    

    ErrorCode get_adjacencies(const EntityHandle source_entity,
                              const unsigned int target_dimension,
                              std::vector<EntityHandle> &target_entities);


    //! Get the upward incidences associated with an entity.
    /** Given an entity of dimension <em>d</em>, gather all the incident <em>D(>d)</em> dimensional entities.
     * Parameters:
     * \param ent EntityHandle to which incident entities have to be found.
     * \param out_dim Dimension of the desired incidence information.
     * \param adjents Vector in which the incident entities are returned.
     * \param local_id Set to false by default. If true, returns the local id's of the half-facets
     * \param lids Vector in which the local id's are returned.
     */

    ErrorCode get_up_adjacencies(EntityHandle ent,
                                 int out_dim,
                                 std::vector<EntityHandle> &adjents,
                                 bool local_id = false,
                                 std::vector<int> * lids = NULL );

    //! Get the same-dimensional entities connected with an entity.
    /** Given an entity of dimension <em>d</em>, gather all the entities connected via <em>d-1</em> dimensional entities.
     *  Same as bridge_adjacencies in MOAB.
     * Parameters:
     * \param ent EntityHandle to which neighbor entities have to be found.
     * \param adjents Vector in which the neighbor entities are returned.
     */

    ErrorCode get_neighbor_adjacencies(EntityHandle ent,
                                       std::vector<EntityHandle> &adjents);

    // 1D Maps and queries

    //! Given a range of edges, determines the map for sibling half-verts and stores them into SIBHVS_EID, SIBHVS_LVID tags.
    /** Compute all sibling half-vertices for all half-vertices in the given curve. The sibling half-verts is
     *  defined in terms of the containing edge and the local id of the vertex w.r.t that edge.
     *  That is, the map consists of two pieces of information: <EntityHandle eid, int lvid>
     *  Parameters:
     * \param edges Range of edges.
    */
    
    ErrorCode determine_sibling_halfverts(Range &edges);

    //! Given a range of edges, determines the map for incident half-verts and stores them into V2HV_EID, V2HV_LVID tags.
    /** Compute a map between a vertex and an incident half-vertex. This map is not always required, but is
     * essential for local neighborhood searching as it acts like an anchor to start the search.
     * Parameters:
     * \param edges Range of edges
    */

    ErrorCode determine_incident_halfverts(Range &edges);

    //! Given a vertex, finds the edges incident on it.
    /** Given a vertex handle, it starts by first finding an incident half-vert by using the incident
     * half-vert map, and then obtaining all the sibling half-verts of the corresponding half-vertex.
     * Parameters:
     * \param vid EntityHandle of the query vertex
     * \param adjents Vector returning the incident edges
     * \param local_id False by default. If true, returns the local vertex id's corresponding to vid
     * \param lvids Vector returning the local vertex id's
    */

    ErrorCode get_up_adjacencies_1d(EntityHandle vid,
                                    std::vector<EntityHandle> &adjents,
                                    bool local_id = false,
                                    std::vector<int> * lvids = NULL);

    //! Given an edge, finds vertex-connected neighbor edges
    /** Given an edge, it gathers all the incident edges of each vertex of the edge.
     * Parameters:
     * \param eid EntityHandle of the query edge
     * \param adjents Vector returning neighbor edges
    */
    
    ErrorCode get_neighbor_adjacencies_1d(EntityHandle eid,
                                          std::vector<EntityHandle> &adjents);
    
    
    // 2D Maps and queries    

    //! Given a range of faces, determines the sibling half-edges and stores them into SIBHES_FID, SIBHES_LEID tags.
    /** Compute all sibling half-edges for all half-edges in the given surface.
     * The sibling half-edges is defined in terms of the containing face and the local id of the edge w.r.t that entity.
     * That is, the map consists of two pieces of information: <EntityHandle fid, int leid>
     * Parameters:
     * \param faces Range of faces
    */
    
    ErrorCode determine_sibling_halfedges(Range &faces);

    //! Given a range of faces, determines the incident half-edges and stores them into V2HE_FID, V2HE_LEID tags.
    /** Compute a map between a vertex and an incident half-edge.
     * This map is not always required, but is essential for local neighborhood searching as it acts
     * like an anchor to start the search.
     * 	Parameters:
     * \param faces Range of faces
    */

    ErrorCode determine_incident_halfedges(Range &faces);

   /* ErrorCode get_up_adjacencies_2d(EntityHandle vid,
                                       std::vector<EntityHandle> &adjents);*/

    //! Given an edge, finds the faces incident on it.
    /** Given an edge, it first finds a matching half-edge corresponding to eid, and then
     * collects all the incident half-edges/faces via the sibhes map.
     * Parameters:
     * \param eid EntityHandle of the query edge
     * \param adjents Vector returning the incident faces
     * \param local_id By default false. If true, returns the local edge id's corresponding to the input edge
     * \param leids Vector returning local edge ids
    */

    ErrorCode get_up_adjacencies_2d(EntityHandle eid,
                                    std::vector<EntityHandle> &adjents,
                                    bool local_id = false,
                                    std::vector<int> * leids = NULL);

    //! Given a half-edge <fid, leid>, finds the faces incident on it.
    /**
     * Parameters:
     * \param fid EntityHandle of the containing face
     * \param leid local id of the edge w.r.t to the face
     * \param add_inent If true, adds the input fid into the returning vector of adjents.
     * \param adjents Vector returning the incident faces
     * \param local_id By default false. If true, returns the local edge id's as well.
     * \param leids Vector returning local edge ids
    */

    ErrorCode get_up_adjacencies_2d(EntityHandle fid,
                                    int leid,
                                    bool add_inent,
                                    std::vector<EntityHandle> &adjents,
                                    bool local_id = false,
                                    std::vector<int> * leids = NULL);

    //! Given an edge, finds edge-connected neighbor face
    /** Given an face, it gathers all the neighbor faces of each local edge of the face.
     * Parameters:
     * \param fid EntityHandle of the query face
     * \param adjents Vector returning neighbor faces
    */
       
    ErrorCode get_neighbor_adjacencies_2d(EntityHandle fid,
                                          std::vector<EntityHandle> &adjents);

    //! Given a range of faces, finds the total number of edges.
    
    int find_total_edges_2d(Range &faces);
    
    // 3D Maps and queries

    //! Given a range of cells, determines the sibling half-faces and stores them into SIBHFS_CID, SIBHFS_LFID tags.
    /** Compute all sibling half-faces for all half-faces in the given volume.
     * The sibling half-faces is defined in terms of the containing cell and the local id of the face w.r.t that cell.
     * That is, the map consists of two pieces of information: <EntityHandle cid, int lfid>
     * Parameters:
     * \param faces Range of cells
    */

    ErrorCode determine_sibling_halffaces( Range &cells);

    //! Given a range of cells, determines the incident half-faces and stores them into V2HF_CID, V2HF_LFID tags.
    /** Compute a map between a vertex and an incident half-face.
     * This map is not always required, but is essential for local neighborhood searching as it acts
     * like an anchor to start the search.
     * 	Parameters:
     * \param faces Range of cells
    */
    
    ErrorCode determine_incident_halffaces( Range &cells);
    
    //! Given a range of cells, tags all border vertices with a true value.
    /** Tag border vertices by using the sibhf_cid map. All vertices on half-faces with no sibling
     * half-faces are considered as border vertices.
     * Parameters:
     * \param cells Range of cells
     * \param isborder: A dense tag over all vertices of size 1. Value is true for a border vertex, otherwise is false.
    */

    ErrorCode determine_border_vertices( Range &cells,
                                         Tag isborder);

    //! Given an edge, finds the cells incident on it.
    /** Given an edge, it first finds a matching local edge in a cell corresponding to eid, and then
     * collects all the incident cells via the sibhfs map.
     * Parameters:
     * \param eid EntityHandle of the query edge
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local edge id's corresponding to the input edge
     * \param leids Vector returning local edge ids
    */

    ErrorCode get_up_adjacencies_edg_3d(EntityHandle eid,
                                        std::vector<EntityHandle> &adjents,
                                        bool local_id = false,
                                        std::vector<int> * leids = NULL);

    //! Given a local edge <cid, leid>, finds the cells incident on it.
    /** Given a local edge, it gathers all the incident cells via the sibhfs map.
     * Parameters:
     * \param cid EntityHandle of the cell containing the local edge
     * \param leid local edge id w.r.t the cell
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local edge id's corresponding to the input edge
     * \param leids Vector returning local edge ids
    */

    ErrorCode get_up_adjacencies_edg_3d(EntityHandle cid,
                                        int leid, std::vector<EntityHandle> &adjents,
                                        bool local_id = false,
                                        std::vector<int> * leids = NULL);

    //! Given an face, finds the cells incident on it.
    /** Given an face, it first finds a matching half-face in a cell corresponding to face, and then
     * collects all the incident cells via the sibhfs map.
     * Parameters:
     * \param fid EntityHandle of the query face
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local face id's corresponding to the input face
     * \param leids Vector returning local face ids
    */

    ErrorCode get_up_adjacencies_face_3d(EntityHandle fid,
                                         std::vector<EntityHandle> &adjents,
                                         bool local_id = false,
                                         std::vector<int>  * lfids = NULL);

    //! Given a local face <cid, lfid>, finds the cells incident on it.
    /** Given a local face, it gathers all the incident cells via the sibhfs map.
     * Parameters:
     * \param cid EntityHandle of the cell containing the local edge
     * \param lfid local face id w.r.t the cell
     * \param adjents Vector returning the incident cells
     * \param local_id By default false. If true, returns the local face id's corresponding to the input face
     * \param lfids Vector returning local face ids
    */

    ErrorCode get_up_adjacencies_face_3d(EntityHandle cid,
                                         int lfid,
                                         std::vector<EntityHandle> &adjents,
                                         bool local_id = false,
                                         std::vector<int> * lfids = NULL);
    
    //! Given a cell, finds face-connected neighbor cells
    /** Given a cell, it gathers all the neighbor cells of each local face of the cell.
     * Parameters:
     * \param cid EntityHandle of the query cell
     * \param adjents Vector returning neighbor cells
    */

    ErrorCode get_neighbor_adjacencies_3d(EntityHandle cid,
                                          std::vector<EntityHandle> &adjents);
    

  protected:

    Core * mb;

    HalfFacetRep();

    bool mInitAHFmaps;

    enum {
      MAX_VERTICES = 8,
      MAX_EDGES = 12,
      MAX_FACES = 6,
      MAX_VERTS_HF = 4,
      MAX_INCIDENT_HF = 4
    };

    Range _verts, _edges, _faces, _cells;
    Tag sibhvs_eid, sibhvs_lvid, v2hv_eid, v2hv_lvid;
    Tag sibhes_fid, sibhes_leid, v2he_fid, v2he_leid;
    Tag sibhfs_cid, sibhfs_lfid, v2hf_cid, v2hf_lfid;


    MESHTYPE get_mesh_type(int nverts, int nedges, int nfaces, int ncells);

    ErrorCode init_curve();
    ErrorCode init_surface();
    ErrorCode init_volume();

    ErrorCode deinit_curve();
    ErrorCode deinit_surface();
    ErrorCode deinit_volume();

    //! Contains the local information for 2D entities
    /** Given a face, find the face type specific information
     * Parameters:
     * \param face EntityHandle. Used to gather info about the type of face for which local info is required
     * \param nepf: Returns the number of vertices/edges for given face type.
    */

    int local_maps_2d(EntityHandle face);

    //! Contains the local information for 2D entities
    /** Given number of edges, returns local indices of next and previous local edges.
     * Parameters:
     * \param nepf: The number of vertices/edges for given face type.
     * \param next, prev: Local ids of next and previous edges w.r.t to the face
     *
     * Note: For 2D entities, the number of vertices and edges are same and each local edge is outgoing
     * corresponding to the local vertex, i.e,

            v2        v3 __e2__v2
            /\          |      |
        e2 /  \ e1    e3|      |e1
          /____\        |______|
        v0  e0  v1     v0  e0  v1
    */


    ErrorCode local_maps_2d(int nepf, int *next, int *prev);

    //! Given a half-edge as <he_fid,he_lid> , finds the half-edges incident on it and adds them
    //  to an input queue if it not already added.
    /** Given an half-edge, obtain all the incident half-edges via the sibhes map and add them to a given
     * queue of half-edges, if they do not already exist in the queue. This function is used to increment the
     * search space for finding a matching half-edge.
     * Parameters:
     * \param he_fid EntityHandle of query half-edge
     * \param he_lid Local id of query half-edge
     * \param queue_fid
     * \param queue_lid
     * \param qsize Array of faces and local edge ids. qsize is the current size of the queue_fid.
     * \param trackfaces Array containing faces. If fid of an incident half-edge doesn't belong to trackfaces,
     *  the half-edge is added to the queue.
     * \param tcount Current size of trackfaces
     */

    ErrorCode get_up_adjacencies_2d(EntityHandle he_fid,
                                    int he_lid,
                                    EntityHandle *queue_fid,
                                    int *queue_lid,
                                    int *qsize,
                                    EntityHandle *trackfaces,
                                    int *tcount);

    //! Given an edge, finds a matching half-edge in the surface.
    /** Given an edge eid, it first collects few half-edges belonging to one-ring neighborhood of
     * the starting vertex of the given edge, and then simultaneously searches and adds to the local list
     * of half-edges for searching, till it finds a matching half-edge.
     * Parameters:
     * \param eid EntityHandle of the query edge
     * \param hefid, helid: Returns the matching half-edge corresponding to the query edge.
    */

    bool find_matching_halfedge( EntityHandle eid,
                                 EntityHandle *hefid,
                                 int *helid);

    //! Gather half-edges to a queue of half-edges.
    /** Given a vertex vid, and a half-edge <he_fid,he_lid>, add another half-edge in the same face sharing the vertex
     * and add all incident half-edges to the queue via sibhes map.
     * Parameters:
     * \param vid EntityHandle of a vertex in a half-edge
     * \param <EntityHandle he_fid, int he_lid>:  Half-edge
     * \param queue_fid, queue_lid, count
     * \param trackfaces, tcount
    */

    ErrorCode gather_halfedges( EntityHandle vid,
                                EntityHandle he_fid,
                                int he_lid,
                                EntityHandle *queue_fid,
                                int *queue_lid,
                                int *qsize,
                                EntityHandle *trackfaces,
                                int *tcount);

    //! Obtains another half-edge belonging to the same face as the input half-edge
    /** It uses the local maps to find another half-edge that is either incident or outgoing depending
     * on vid and input half-edge
     * Parameters:
     * \param vid EntityHandle of a vertex in a half-edge
     * \param <EntityHandle he_fid, int he_lid>:  Half-edge
     * \param <EntityHandle he2_fid, int he2_lid>: Returns another half-edge in the same he_fid sharing vid.
    */

    ErrorCode another_halfedge( EntityHandle vid,
                                EntityHandle he_fid,
                                int he_lid,
                                EntityHandle *he2_fid,
                                int *he2_lid);

    //! Collect and compare to find a matching half-edge with the given edge connectivity.
    /** Given edge connectivity, compare to an input list of half-edges to find a matching half-edge
     * and add a list of half-edges belonging to the one-ring neighborhood to a queue till it finds a match.
     * Parameters:
     * \param edg_vert End vertices of an edge
     * \param queue_fid, queue_lid, count
     * \param trackfaces, tcount
     * \param <EntityHandle he_fid, int he_lid>: Returns matching half-edge
    */

    bool collect_and_compare(std::vector<EntityHandle> &edg_vert,
                             EntityHandle *queue_fid,
                             int *queue_lid,
                             int *qsize,
                             EntityHandle *trackfaces,
                             int *tcount,
                             EntityHandle *he_fid,
                             int *he_lid);


    //! The local maps for 3D entities.
    /** Types of 3D entities supported: tets, pyramid, prism, hex
        Determines the type from input "cell"

	_3d_numels:
	nvpc: Number of vertices per cell
	nepc: Number of edges per cell
	nfpc: Number of faces per cell

	_3d_numels:
	nvmax: The maximum number of vertices of all half-faces of the cell

	_3d_hf2v: Map half-face to vertices
	hf2v_num: Array storing number of vertices per half-face
	hf2v_map: Local ids of vertices of each half-face, stored in an array
	hf2v_idx: Starting index for each half-face to access vertices

	_3d_v2hf: Map vertex to half-face
	v2hf_num: Array storing number of incident half-faces per vertex
	v2hf_map: Local ids of incident half-faces, stored in an array
	v2hf_idx: Starting index for each vertex to access incident half-faces

	_3d_edges: Maps for edges
	e2v: Local edge to local vertices
	e2hf: Local edge to incident half-faces
	f2leid: Local edges for each half-faces

	_3d_lookup_leid:
	nvpc: Is an input #vertices per cell
	lookup_leid: Map between local vertex v0 to local vertex v1 storing the local edge id e = <v0,v1>
    */


    struct LocalMaps3D{
      short int num_verts_in_cell; // Number of vertices per cell
      short int num_edges_in_cell; // Number of edges per cell
      short int num_faces_in_cell; // Number of faces per cell

      int hf2v_num[MAX_FACES]; //
      int hf2v[MAX_FACES][MAX_VERTS_HF];

      int v2hf_num[MAX_VERTICES];
      int v2hf[MAX_VERTICES][MAX_INCIDENT_HF];

      int e2v[MAX_EDGES][2];
      int e2hf[MAX_EDGES][2];
      int f2leid[MAX_FACES][MAX_VERTS_HF];
      int lookup_leids[MAX_VERTICES][MAX_VERTICES];
    };

    static const LocalMaps3D lConnMap3D[4];

    int get_index_from_type(EntityHandle cid);

    //! Given an edge, finds a matching local edge in an incident cell.
    /** Find a local edge with the same connectivity as the input edge, belonging to an incident cell.
     * Parameters:
     * \param eid EntityHandle of the edge
     * \param cid Returns EntityHandle of the incident cell
     * \param leid Returns the local id of the edge corresponding to the input edge w.r.t the incident cell.
    */

    bool find_matching_implicit_edge_in_cell( EntityHandle eid,
                                              EntityHandle *cid,
                                              int *leid);

    //! Given a face, finds a matching local face in an incident cell.
    /** Find a local face with the same connectivity as the input face, belonging to an incident cell.
     * Parameters:
     * \param fid EntityHandle of the face
     * \param cid Returns EntityHandle of the incident cell
     * \param lfid Returns the local id of the face corresponding to the input face w.r.t the incident cell.
    */

    bool find_matching_halfface(EntityHandle fid,
                                EntityHandle *cid,
                                int *leid);

    bool find_match_in_array(EntityHandle ent,
                             EntityHandle *ent_list,
                             int count,
                             bool get_index = false,
                             int *index = NULL);

  };

} // namespace moab 

#endif

