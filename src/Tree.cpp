#include "moab/Tree.hpp"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"

#include <limits>

namespace moab 
{
    ErrorCode Tree::parse_common_options(FileOptions &options) 
    {
      double tmp_dbl;
      int tmp_int;
        // MAX_PER_LEAF: max entities per leaf; default = 6
      ErrorCode rval = options.get_int_option("MAX_PER_LEAF", tmp_int);
      if (MB_SUCCESS == rval) maxPerLeaf = std::max(tmp_int, 1);
      
        // MAX_DEPTH: max depth of the tree; default = 30
      rval = options.get_int_option("MAX_DEPTH", tmp_int);
      if (MB_SUCCESS == rval) maxDepth = tmp_int;
      if (maxDepth < 1) maxDepth = std::numeric_limits<unsigned>::max();

        // MIN_WIDTH: minimum width of box, used like a tolerance; default = 1.0e-10
      rval = options.get_real_option("MIN_WIDTH", tmp_dbl);
      if (MB_SUCCESS == rval) minWidth = tmp_dbl;

        // MESHSET_FLAGS: flags passed into meshset creation for tree nodes; should be a value from
        //          ENTITY_SET_PROPERTY (see Types.hpp); default = MESHSET_SET
      rval = options.get_int_option("MESHSET_FLAGS", tmp_int);
      if (MB_SUCCESS == rval && 0 <= tmp_int) meshsetFlags = (unsigned) tmp_int;
      else if (0 > tmp_int) return MB_FAILURE;

        // CLEAN_UP: if false, do not delete tree sets upon tree class destruction; default = true
      bool tmp_bool;
      rval = options.get_toggle_option("CLEAN_UP", true, tmp_bool);
      if (MB_SUCCESS == rval && !tmp_bool) cleanUp = false;

        // TAG_NAME: tag name to store tree information on tree nodes; default = "AKDTree"
      std::string tmp_str;
      rval = options.get_str_option("TAG_NAME", tmp_str);
      if (MB_SUCCESS == rval) boxTagName = tmp_str;

      return MB_SUCCESS;
    }

    static inline void box_accum( const CartVect& point,
                                  CartVect& bmin,
                                  CartVect& bmax )
    {
      for (unsigned j = 0; j < 3; ++j) {
        if (point[j] < bmin[j])
          bmin[j] = point[j];
        if (point[j] > bmax[j])
          bmax[j] = point[j];
      }
    }

    ErrorCode Tree::compute_bounding_box(Interface &iface, const Range& elems, CartVect &box_min, CartVect &box_max)
    {
      ErrorCode rval;
      box_min = CartVect(HUGE_VAL);
      box_max = CartVect(-HUGE_VAL);
      
      CartVect coords;
      EntityHandle const *conn, *conn2;
      int len, len2;
      Range::const_iterator i;
  
        // vertices
      const Range::const_iterator elem_begin = elems.lower_bound( MBEDGE );
      for (i = elems.begin(); i != elem_begin; ++i) {
        rval = iface.get_coords( &*i, 1, coords.array() );
        if (MB_SUCCESS != rval)
          return rval;
        box_accum( coords, box_min, box_max );
      }

        // elements with vertex-handle connectivity list
      const Range::const_iterator poly_begin = elems.lower_bound( MBPOLYHEDRON, elem_begin );
      std::vector<EntityHandle> dum_vector;
      for (i = elem_begin; i != poly_begin; ++i) {
        rval = iface.get_connectivity( *i, conn, len, true, &dum_vector);
        if (MB_SUCCESS != rval)
          return rval;

        for (int j = 0; j < len; ++j) {
          rval = iface.get_coords( conn+j, 1, coords.array() );
          if (MB_SUCCESS != rval)
            return rval;
          box_accum( coords, box_min, box_max );
        }
      }
  
        // polyhedra
      const Range::const_iterator set_begin  = elems.lower_bound( MBENTITYSET, poly_begin );
      for (i = poly_begin; i != set_begin; ++i) {
        rval = iface.get_connectivity( *i, conn, len, true );
        if (MB_SUCCESS != rval)
          return rval;

        for (int j = 0; j < len; ++j) {
          rval = iface.get_connectivity( conn[j], conn2, len2 );
          for (int k = 0; k < len2; ++k) {
            rval = iface.get_coords( conn2+k, 1, coords.array() );
            if (MB_SUCCESS != rval)
              return rval;
            box_accum( coords, box_min, box_max );
          }
        }
      }
  
        // sets
      CartVect tmin, tmax;
      for (i = set_begin; i != elems.end(); ++i) {
        Range tmp_elems;
        rval = iface.get_entities_by_handle(*i, tmp_elems);
        if (MB_SUCCESS != rval) return rval;
        rval = compute_bounding_box(iface, tmp_elems, tmin, tmax);
        if (MB_SUCCESS != rval) return rval;
      
        for (int j = 0; j < 3; ++j) {
          if (tmin[j] < box_min[j])
            box_min[j] = tmin[j];
          if (tmax[j] > box_max[j])
            box_max[j] = tmax[j];
        }
      }
  
      return MB_SUCCESS;
    }

    ErrorCode Tree::find_all_trees( Range& results )
    {
      return moab()->get_entities_by_type_and_tag( 0, MBENTITYSET, &boxTag, 0, 1, results );
    }

}
