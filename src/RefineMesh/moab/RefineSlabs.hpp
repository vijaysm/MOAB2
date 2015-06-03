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

  protected:

    ErrorCode initialize();

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
      bool operator==(const SlabHex& rhs) const
      {
        return (entity_handle == rhs.entity_handle);
      }
    protected:
      EntityHandle entity_handle;

    };
    class SlabNode : public SlabEntity
    {
    public:
      bool operator==(const SlabNode& rhs) const
      {
        return (entity_handle == rhs.entity_handle);
      }
    protected:
      EntityHandle entity_handle;

    };
    class SlabEdge : public SlabEntity
    {
    public:
      bool operator==(const SlabEdge& rhs) const
      {
        return (head_node == rhs.head_node) &&
               (tail_node == rhs.tail_node);
              // hex and local id need not match
               // && (hex == rhs.hex);
      }
      // store two vertices, in forward direction
      // store a hex and a local id and a direction
      int edge_lid;
      EntityHandle head_node, tail_node, hex;
      SlabEdge() : edge_lid(-1), head_node(bad_handle), tail_node(bad_handle), hex(bad_handle)  {}

      bool nodes_match( const SlabEdge &rhs ) const
      {
        return ( head_node == rhs.head_node && tail_node == rhs.tail_node) 
            || ( head_node == rhs.tail_node && tail_node == rhs.head_node);
      }
      // -1 if opposite, +1 if same, 0 if nodes don't match
      bool directions_match( const SlabEdge &rhs ) const
      {
        return ( head_node == rhs.head_node && tail_node == rhs.tail_node);
      }
      void flip()
      {
        EntityHandle tmp = head_node;
        head_node = tail_node;
        tail_node = tmp;
      }

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

      void copy_data( SlabData *copy_me )
      {
        is_coarse = copy_me->is_coarse;
        membership = copy_me->membership;
        shrink_membership = copy_me->shrink_membership;
      }

      SlabData() : is_coarse(true), membership(INTERNAL), shrink_membership(EXTERNAL), my_copy_good(false), my_copy(bad_handle), mini_me_good(false), mini_me(bad_handle), refinement(0), coarsening(0) {}
    };

    HexRefinement *get_hex_refinement( EntityHandle coarse_hex )
    {
      SlabData *slab_data = get_slab_data( coarse_hex );
      if ( !slab_data )
        return 0;
      return slab_data->refinement;
    }
    HexRefinement *force_hex_refinement( EntityHandle coarse_hex )
    {
      SlabData *slab_data = get_slab_data( coarse_hex );
      assert( slab_data );
      if ( !slab_data->refinement )
        slab_data->refinement = new HexRefinement;
      return slab_data->refinement;
    }
    HexCoarsening *force_hex_coarsening( EntityHandle fine_hex )
    {
      SlabData *slab_data = get_slab_data( fine_hex );
      assert( slab_data );
      if ( !slab_data->coarsening )
        slab_data->coarsening = new HexCoarsening(fine_hex);
      return slab_data->coarsening;
    }

    void add_refined_hex( EntityHandle coarse_hex, EntityHandle fine_hex )
    {
      force_hex_refinement( coarse_hex )->fine_hexes.push_back( fine_hex );
      force_hex_coarsening( fine_hex );
    }

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
    bool unique( std::vector< SlabEdge > edges, SlabEdge edge );
    bool add_unique( std::vector< SlabEdge > edges, SlabEdge edge );
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
        std::vector< SlabEdge > &upper_slab_edges, std::vector< SlabEdge > parallel_slab_edges );    

    typedef std::map< std::pair<EntityHandle, EntityHandle>, bool > EdgeDataMap;
    typedef EdgeDataMap::iterator EdgeDataIterator;
    EdgeDataMap edge_data_map;
    void reset_in_slab();
    bool any_in_slab( std::vector< SlabEdge > slab_edges );
    bool get_in_slab( SlabEdge slab_edge );
    void set_in_slab( SlabEdge slab_edge, bool new_value = true );

    typedef std::map<EntityHandle, SlabData*> SlabDataMap;
//    typedef std::unordered_map<EntityHandle, SlabData*> SlabDataMap;
    typedef SlabDataMap::iterator SlabDataIterator;
    SlabDataMap slab_data_map;
    SlabData *get_slab_data( EntityHandle entity_handle) 
    { 
      SlabDataIterator it = slab_data_map.find(entity_handle);
      if ( it == slab_data_map.end() )
        return 0;
      return it->second;
    }
    SlabData *force_slab_data( EntityHandle entity_handle ) 
    { 
      return slab_data_map[entity_handle];
    }

    SlabData* set_coarse_entity( EntityHandle entity )
    {
      SlabData *slab_data = force_slab_data(entity);
      slab_data->is_coarse = true;
      return slab_data;
    }

    void set_copy( EntityHandle entity, EntityHandle copy )
    {      
      SlabData *slab_data = force_slab_data( entity );
      slab_data->my_copy_good = true;
      slab_data->my_copy = copy;

      SlabData *slab_data_copy = force_slab_data( copy );
      slab_data_copy->my_copy_good = true;
      slab_data_copy->my_copy = entity;
    }
    bool get_copy( EntityHandle entity, EntityHandle &copy )
    {
      SlabData * slab_data = get_slab_data( entity );
      if (slab_data && slab_data->my_copy_good)
      {
        copy = slab_data->my_copy;
        return true;
      }
      copy = bad_handle;
      return false;
    }
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

    void set_fine_node( EntityHandle entity, EntityHandle fine )
    {
      SlabData *slab_data = force_slab_data(entity);
      slab_data->mini_me_good = true;
      slab_data->mini_me = fine;

      SlabData *fine_slab = force_slab_data(fine);
      fine_slab->mini_me_good = true;
      fine_slab->mini_me = fine;
    }
    void get_fine_nodes( EntityHandle *coarse_nodes, EntityHandle *fine_nodes, int num_nodes);
    bool get_fine_node( EntityHandle node, EntityHandle &fine_node )
    {
      SlabData *slab_data = get_slab_data(node);
      if (slab_data && slab_data->mini_me_good)
      {
        fine_node = slab_data->mini_me;
        return true;
      }
      fine_node = bad_handle;
      return false;
    }

    void set_coarse_hex( EntityHandle hex )
    {
      set_coarse_entity(hex);
    }

    void set_coarse_node( EntityHandle node )
    {
      set_coarse_entity(node);
    }

    void set_coarse_quad( EntityHandle quad )
    {
      set_coarse_entity(quad);
    }
    bool get_coarse( EntityHandle entity_handle )
    {
      SlabData *slab_data = get_slab_data(entity_handle);
      if (!slab_data)
        return false;
      return slab_data->is_coarse;
    }
  
    SlabData::Membership get_membership( EntityHandle entity_handle )
    { 
      SlabData *slab_data = get_slab_data(entity_handle);
      if (!slab_data)
        return SlabData::EXTERNAL;
      return slab_data->membership;
    }
    void set_membership( EntityHandle entity_handle, SlabData::Membership membership )
    { 
      SlabData *slab_data = force_slab_data(entity_handle);
      slab_data->membership = membership;
      if (membership == SlabData::EXTERNAL)
        slab_data->is_coarse = false;
    }


    SlabData::Membership get_shrink_membership( EntityHandle entity_handle )
    { 
      SlabData *slab_data = get_slab_data(entity_handle);
      if (!slab_data)
        return SlabData::EXTERNAL;
      return slab_data->shrink_membership;
    }
    void set_shrink_membership( EntityHandle entity_handle, SlabData::Membership membership )
    { 
      SlabData *slab_data = get_slab_data(entity_handle);
      assert( slab_data );
      slab_data->shrink_membership = membership;
      if (membership == SlabData::EXTERNAL)
      {
        slab_data->mini_me_good = false;
        slab_data->mini_me = bad_handle;
      }
    }
    // mark the hexes as being in the pillow (shrink) set, determine whether nodes are interior or boundary
    void shrink_mark_slab( Entities &slab );
    void shrink_mark_coarse_slab( Entities &slab );
    void shrink_mark_fine_slab( Entities &slab, Entities &shrink_set);
    void remove_shrink_mark_slab( Entities &slab );
    // get the dimension of the geometric object that this mesh entity lies on
    int get_geometry_dimension( EntityHandle /*entity_handle*/ )
    {
      return 3; // zzyk AHF or MOAB or CGM todo
    }

    // pick some unrefined edge of the hex to start growing a new slab
    bool find_seed_edge( EntityHandle hex, int &edge_lid, SlabEdge &slab_edge );
    void add_edge( SlabEdge &slab_edge, Entities &slab );
    // traverse the sheet outward from the slab edge, adding hexes to the slab
    void extend_slab( SlabEdge slab_edge, Entities&slab );
    // True if the candidate slab_edge should be added to the slab? No if it is already refined, or already in the slab
    bool is_good_slab_edge( const SlabEdge &slab_edge );
    //   this version passes back the ortho and parallel edges, iff the return value is true
    bool is_good_slab_edge( const SlabEdge &slab_edge, std::vector<SlabEdge> &ortho, std::vector<SlabEdge> &upper, std::vector<SlabEdge> &parallel );
    //   this version makes the slab_edge out of the hex, its edge, and the orientation (node_01)
    bool is_good_slab_edge( EntityHandle hex, int edge_lid, int node_01, SlabEdge &slab_edge );
    // True if none of the slab edges have already been refined.
    bool none_refined( std::vector<SlabEdge> &slab_edges );

    // main routines
    ErrorCode initialize_refinement( Entities &coarse_hexes, Entities &coarse_quads );
    ErrorCode mark_hex_nodes(Entities &coarse_hexes);
    ErrorCode mark_surface_nodes(Entities &coarse_quads);
    ErrorCode find_slabs( Entities &coarse_hexes, Entities &coarse_quads, EntitiesVec &slabs );
    ErrorCode pillow_slabs( EntitiesVec &slabs );
    ErrorCode pillow_slab( Entities &slab ); // pillow an individual slab
    void pillow_hexes( Entities &shrink_set, Entities &new_hexes );

    // AHF todo
    ErrorCode replace_mesh( Entities &coarse_hexes, Entities &coarse_quads, Entities &fine_hexes, Entities &fine_quads ); // AHF todo
    void replace_node( EntityHandle chex, int node_lid, EntityHandle new_node); // AHF todo
    void udpate_AHF_connectivity(); // AHF todo



  protected:
    Core *mbImpl;
    HalfFacetRep *ahf, *refinement_ahf;

    ErrorCode new_refinement_ahf( size_t num_hexes_memory_estimate );

  };
} //name space moab
#endif
