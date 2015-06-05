/*! \file RefineSlabs.hpp
 * Refine a set of unstructured hexes using pillowing. One level, where the ideal is one hex into 8.
  For more finer refinements (more levels) a template based approach such as RefineSlabs.hpp would 
  be better.
 */


#ifndef REFINE_SLABS_HPP
#define REFINE_SLABS_HPP

#include <map>
// #include <unordered_map> // Scott Mitchell's mac compiler version doesn't seem to support this c++11 extention

#include "moab/Range.hpp"
#include "moab/HalfFacetRep.hpp"
#include "moab/CN.hpp"

namespace moab
{
  
  class Core;
  class HalfFacetRep;

  class RefineSlabs
  {
    
  public:

    RefineSlabs(Core *impl);
    
    ~RefineSlabs();

    // top level interface function, the only one    
    //! \brief Generate a refined mesh.
    /** Given a hex mesh and a subset of hexes, refine the subset conformally. 
        The ideal is to refine each hex into eight. 
        Inputs:
          Handles to the subset of hexes. (Hexes outside those are not changed.)
          Handles to the subset of quads of a surface mesh that should be considered internal to the set.
           (The surface mesh will be changed if any faces are specified.)
        Outputs:
          Handles to the new, refined set of hexes.
          Handles to the new, refined surface quads.
        Changed:
          The mesh database, referenced by impl in the constructor (mbImpl)
      */
    // I considered also moab Range; here a vector seems better
    typedef std::vector<EntityHandle> Entities;
    typedef std::vector<Entities> EntitiesVec;
    ErrorCode refine_mesh(Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities &fine_quads);

    // debug
    void register_entity_handles(EntityHandle *hexes, int num_hexes, EntityHandle *vertices, int num_vertices);
    bool is_registered_hex( EntityHandle hex );
    bool is_registered_vertex( EntityHandle vertex );

  protected:

    // ======= data members
    Core *mbImpl;
    HalfFacetRep *ahf, *refinement_ahf;
    ErrorCode new_refinement_ahf( size_t num_hexes_memory_estimate );

    // debugging
    // ensure that any coarse entity we request information from AHF/Core, was one of the coarse hexes passed in.
    // This allows range checking, and also (sometimes) if we passed in the wrong sort of index, since EntityType is a basic type.
    EntityHandle *registered_hexes, *registered_vertices;
    int registered_num_hexes, registered_num_vertices;

    // main routines
    // setup 
    ErrorCode initialize(); // companion to constructor
    ErrorCode initialize_refinement( Entities &coarse_hexes, Entities &coarse_quads );
    // mark the coarse_hexes as being the hexes we wish to refine
    // identify the boundary (surface) and interior nodes of the hexes
    ErrorCode mark_hex_nodes(Entities &coarse_hexes);
    // ditto for coarse quads. Recall the coarse quads are on the boundary of the coarse hexes and on the geometric surface.
    // Such quads are treated as if they were in the interior of the set, so they will not change, but new surface quads will be added around it.
    // If none are passed in, the surface mesh will not change.
    ErrorCode mark_surface_nodes(Entities &coarse_quads);
    // find all the slabs, subsets of coarse hexes, that will be pillowed one subset at a time.
    ErrorCode find_slabs( Entities &coarse_hexes, Entities &coarse_quads, EntitiesVec &slabs );
    // pillow all the slabs
    ErrorCode pillow_slabs( EntitiesVec &slabs );
    // pillow an individual slab
    ErrorCode pillow_slab( Entities &slab ); 
    // pillow the hexes of an individual slab
    void pillow_hexes( Entities &shrink_set, Entities &new_hexes );

    // AHF todo
    // remove the coarse hexes and quads from the AHF, and replace them with the refined fine_hexes and quad from the refinement_AHF
    ErrorCode replace_mesh( Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities &fine_quads ); // AHF todo
    // For a hex, replace one of its nodes with a new node. 
    // This is a hex-by-hex operation, so the connectivity between hexes might be ill-defined, and the AHF might be out of date,
    // while these sorts of operations are taking place. So, when we're done and the connectivity is well-defined again, call the update_AHF_connectivity function.
    void replace_node( EntityHandle chex, int node_lid, EntityHandle new_node); // AHF todo
    void udpate_AHF_connectivity(); // AHF todo


  protected:
    // ======= nested classes, little datastructures hanging off of entities, 
    // and interface functions to call AHF and Moab::Core routines to traverse the mesh

    // hexes and nodes are identified by an EntityHandle
    // edges and quads do not have EntityHandles
    // This is Scott Mitchells attempt to create a truly generic mesh entity for refining slabs
    // and a generic and clean way to access the adjacency calls needed
    class SlabEntity
    {
      // to tell if two entities are the same, call the is_equal function
    public:
      virtual ~SlabEntity() {}
    };
    // Since there is no such thing as a bad entity_handle flag, we implement a bool for indicating the entity is unasigned, etc
    // We initialize to -1 anyway, but don't count on this being a bad value
    static const EntityHandle bad_handle = -1;
    class SlabHex : public SlabEntity
    {
    public:
      bool operator==(const SlabHex& rhs) const;
    protected:
      EntityHandle entity_handle;

    };
    class SlabNode : public SlabEntity
    {
    public:
      bool operator==(const SlabNode& rhs) const;
    protected:
      EntityHandle entity_handle;

    };
    class SlabEdge : public SlabEntity
    {
    public:
      // == if nodes_match, regardless of directions_match. Is this the behavior we want?
      bool operator==(const SlabEdge& rhs) const;
      bool nodes_match( const SlabEdge &rhs ) const;
      bool directions_match( const SlabEdge &rhs ) const;
      // flip directions of the edge, by swapping head and tail node
      void flip(); 


      // store two vertices, in forward direction
      // store a hex and a local id and a direction
      int edge_lid;
      EntityHandle head_node, tail_node, hex;
      SlabEdge() : edge_lid(-1), head_node(bad_handle), tail_node(bad_handle), hex(bad_handle)  {}

      // the default copy constructors should be fine, but let's declare them explicitly 
      // just in case there is any compiler weirdness because SlabEntity is a base class with a virtual destructor

      // copy constructors; these are used by the vectors, for example
      // this could be sped up if desired, by using member initialization instead of assignment
      SlabEdge ( const SlabEdge & copy_me ) { assignment(copy_me); }
      SlabEdge ( SlabEdge & copy_me ) { assignment(copy_me); }
      // assignment operator =, should call the copy constructor
      void assignment( const SlabEdge &copy_me );
    };

    // unused
    // class SlabQuad : public SlabEntity
    // {
    // public:
    //   bool operator==(const SlabQuad& rhs) const
    //   {
    //     return (entity_handle == rhs.entity_handle);
    //   }
    //   // could be four vertices
    //   // could be a hex and a local id
    // };

    // global
    // true if m1 and m2 represent the same entity
    bool is_equal( const SlabEntity *m1, const SlabEntity *m2 ) const;

    class HexRefinement
    {
    public:
      Entities fine_hexes;
      int edge_refinement_level[12];

      HexRefinement() : 
      fine_hexes(0)
      {
        for (size_t i = 0; i < 12; ++i )
          edge_refinement_level[i] = 0;
      }
    };

    class HexCoarsening
    {
    public:
      EntityHandle coarse_hex;
      EntityHandle coarse_owner;
      HexCoarsening( EntityHandle chex ) : coarse_hex(chex), coarse_owner(chex) {}
      // copy constructor used
    };

    // data we hang off of both coarse and fine entities. 
    // IMO too complicated to deal with registering attributes, etc. For now, just use a map, and let the Moab experts migrate it later if efficiency is an issue
    class SlabData
    {
    public:
      // EntityHandle entity_handle;
      // no data needed, just its existence is enough, basically booleans we can assign to an entity
      bool is_coarse;
      enum Membership {EXTERNAL, BOUNDARY, INTERNAL};
      Membership membership; // membership wrt the refinement set 

      Membership shrink_membership; // membership wrt the current shrink set
      bool my_copy_good; // true if my_copy has been assigned
      EntityHandle my_copy; // points from coarse to fine mesh
      bool mini_me_good; // true if mini_me has been assigned
      EntityHandle mini_me; // points from a fine mesh node to its shrunk copy

      // todo memory optimization, coarse hexes have a refinement, fine hexes have a coarsening
      // make two separate derived classes
      HexRefinement *refinement;
      HexCoarsening *coarsening;

      void copy_data( SlabData *copy_me );

      SlabData() : is_coarse(true), membership(INTERNAL), shrink_membership(EXTERNAL), my_copy_good(false), my_copy(bad_handle), mini_me_good(false), mini_me(bad_handle), refinement(0), coarsening(0) {}
    };

    // return the HexRefinement of the coarse_hex
    HexRefinement *get_hex_refinement( EntityHandle coarse_hex );
    // create one if it doesn't exist
    HexRefinement *force_hex_refinement( EntityHandle coarse_hex );
    HexCoarsening *force_hex_coarsening( EntityHandle fine_hex );
    // add the fine_hex to the list of hexes that the coarse hex is divided into
    void add_refined_hex( EntityHandle coarse_hex, EntityHandle fine_hex );


    // get the nodes of the hex (quad), through moab core
    ErrorCode get_hex_nodes( EntityHandle hex, EntityHandle hex_nodes[8] );   
    void get_quad_nodes( EntityHandle hex, const EntityHandle hex_nodes[8], int face_lid, EntityHandle* quad_nodes );
    ErrorCode get_quad_nodes( EntityHandle quad, EntityHandle quad_nodes[4] );

    // get the hexes (quads) of the node, through ahf(true) or refinement_ahf(false)
    void get_all_hexes( EntityHandle node, Entities &hexes, bool is_coarse = true );
    void get_all_quads( EntityHandle node, Entities &quads, bool is_coarse = true );

    // convert the oriented edge defined by the hex, its local id in the hex, and the starting node, into a SlabEdge
    void get_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge );

    int get_edge_refinement_level(EntityHandle hex, int edge_lid) 
    {
      return force_hex_refinement(hex)->edge_refinement_level[edge_lid]; 
    }
    // given a node of a hex, return which node it is. e.g. 0, 1, ...8
    int get_hex_node_index( EntityHandle hex, EntityHandle node );

    // slab edge functions 

    // true if edge, and any of its possible representations are not present in edges
    bool unique( const std::vector< SlabEdge > &edges, const SlabEdge &edge );
    bool add_unique( std::vector< SlabEdge > &edges, const SlabEdge &edge );
    // given a slab edge, find an equivalent slab edge defined wrt the passed in hex
    bool get_matching_edge( EntityHandle hex, const SlabEdge &slab_edge, SlabEdge &match );
      // given a SlabEdge (edge and vertex) of a hex, find the other two edges that are opposite the edge in a quad of the hex
    void get_opp( const SlabEdge &edge, SlabEdge &opp1, SlabEdge &opp2 );
    // given a SlabEdge (edge and vertex) of a hex, find the other two edges sharing that vertex
    void get_adj( const SlabEdge &edge, SlabEdge &adj1, SlabEdge &adj2 );
    // get the three slab edges of the hex emanating from the given node
    void get_star_edges( EntityHandle hex, EntityHandle node, SlabEdge star[3] );

    // Find all the nearby slab edges, by traversing the mesh through ahf
    //   in contrast the above functions are for a single hex
    // ortho edges are the ones containing the vertex orthogonal to the edge in some hex
    // upper_slab_edges are the other ones containing the vertex
    // parallel_slab_edges are the face-opposite to the slab_edge
    void get_adjacent_slab_edges( const SlabEdge &slab_edge, std::vector< SlabEdge > &ortho_slab_edges, 
        std::vector< SlabEdge > &upper_slab_edges, std::vector< SlabEdge > &parallel_slab_edges );    

    typedef std::map< std::pair<EntityHandle, EntityHandle>, bool > EdgeDataMap;
    typedef EdgeDataMap::iterator EdgeDataIterator;
    EdgeDataMap edge_data_map;
    void reset_in_slab();
    bool any_in_slab( std::vector< SlabEdge > slab_edges );
    bool get_in_slab( const SlabEdge &slab_edge );
    void set_in_slab( const SlabEdge &slab_edge, bool new_value = true );

    // This is how we associate SlabData with specific EntityHandles
    // If efficiency becomes an issue, we could store the SlabData as attributes off of the actual entities
    typedef std::map<EntityHandle, SlabData*> SlabDataMap;
    //    typedef std::unordered_map<EntityHandle, SlabData*> SlabDataMap;
    typedef SlabDataMap::iterator SlabDataIterator;
    SlabDataMap slab_data_map;
    
    // get the slab_data associated with the entity_handle, 
    // which could be a coarse hex, coarse node, fine hex, or fine node
    SlabData *get_slab_data( EntityHandle entity_handle);
    // create one if it doesn't exist
    SlabData *force_slab_data( EntityHandle entity_handle );
    // This is a fine entity. Associate it with the passed-in coarse entity
    SlabData* set_coarse_entity( EntityHandle entity );
    // The first entity is the original. Associate it with the passed in copy of it.
    void set_copy( EntityHandle entity, EntityHandle copy );
    // Get the copy of the original entity
    bool get_copy( EntityHandle entity, EntityHandle &copy );

    // copy the coarse nodes into the fine ahf
    ErrorCode copy_nodes( Entities &coarse_hexes );
    // copy a coarse node into a fine one
    EntityHandle copy_node( EntityHandle coarse_node );

    ErrorCode create_node( EntityHandle node, EntityHandle &new_node );

    // copy the coarse hexes into the fine ahf
    ErrorCode copy_hexes( Entities &coarse_hexes );
    ErrorCode create_hex( EntityHandle fine_nodes[8], EntityHandle & new_hex  );

    // copy a fine node into another fine one, for pillowing a shrink set
    EntityHandle shrink_node( EntityHandle fine_node );

    // For a coarse node entity, associate it with the passed in fine node entity
    void set_fine_node( EntityHandle entity, EntityHandle fine );
    // For multiple coarse nodes
    void get_fine_nodes( EntityHandle *coarse_nodes, EntityHandle *fine_nodes, int num_nodes);
    // Get of above 
    bool get_fine_node( EntityHandle node, EntityHandle &fine_node );

    // code clarity (convenience) functions, making explicit the type of entity being dealt with
    void set_coarse_hex( EntityHandle hex )   { set_coarse_entity(hex); }
    void set_coarse_node( EntityHandle node ) { set_coarse_entity(node); }
    void set_coarse_quad( EntityHandle quad ) { set_coarse_entity(quad); }

    // true if the entity_handle is a coarse entity, false if a fine or other entity.
    bool is_coarse( EntityHandle entity_handle );
    // Get/set whether the entity is internal, external, or on the boundary of a set.
    // Here a set is a slab, or the coarse hexes being refined.
    SlabData::Membership get_membership( EntityHandle entity_handle );
    void set_membership( EntityHandle entity_handle, SlabData::Membership membership );
    // Here a set is a group (shrink set) of hexes about to be pillowed
    SlabData::Membership get_shrink_membership( EntityHandle entity_handle );
    void set_shrink_membership( EntityHandle entity_handle, SlabData::Membership membership );

    // mark the hexes as being in the pillow (shrink) set, determine whether nodes are interior or boundary
    void shrink_mark_slab( Entities &slab );
    void shrink_mark_coarse_slab( Entities &slab );
    void shrink_mark_fine_slab( Entities &slab, Entities &shrink_set);
    void remove_shrink_mark_slab( Entities &slab );
    // get the dimension of the geometric object that this mesh entity lies on. 
    // E.g. 3 if inside the volume, 2 if on its surface, 1 if on an edge of the surface, ...
    int get_geometry_dimension( EntityHandle entity_handle );

    // pick some unrefined edge of the hex to start growing a new slab
    bool find_seed_edge( EntityHandle hex, int &edge_lid, SlabEdge &slab_edge );
    void add_edge( const SlabEdge &slab_edge, Entities &slab );
    // traverse the sheet outward from the slab edge, adding hexes to the slab
    void extend_slab( const SlabEdge &slab_edge, Entities&slab );
    // True if the candidate slab_edge should be added to the slab? No if it is already refined, or already in the slab
    bool is_good_slab_edge( const SlabEdge &slab_edge );
    //   this version passes back the ortho and parallel edges, iff the return value is true
    bool is_good_slab_edge( const SlabEdge &slab_edge, std::vector<SlabEdge> &ortho, std::vector<SlabEdge> &upper, std::vector<SlabEdge> &parallel );
    //   this version makes the slab_edge out of the hex, its edge, and the orientation (node_01)
    bool is_good_slab_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge );
    // True if none of the slab edges have already been refined.
    bool none_refined( std::vector<SlabEdge> slab_edges );


  };
} //name space moab
#endif
