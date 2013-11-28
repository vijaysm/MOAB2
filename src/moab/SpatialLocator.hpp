/**\file SpatialLocator.hpp
 * \class moab::SpatialLocator
 * \brief Tool to facilitate spatial location of a point in a mesh
 *
 * SpatialLocator facilitates searching for points in or performing ray traces on collections of mesh entities
 * in 2D or 3D.  This searching is facilitated by a tree-based decomposition of the mesh.  Various types
 * of trees are implemented in MOAB and can be used by this tool, but by default it uses AdaptiveKDTree
 * (see child classes of Tree for which others are available).  Parallel and serial searching are both 
 * supported.
 *
 * SpatialLocator can either cache the search results for points or pass back this information in arguments.  
 * Cached information is kept in locTable, indexed in the same order as search points passed in.  This information
 * consists of the entity handle containing the point and the parametric coordinates inside that element.
 * Information about the points searched, e.g. the entities from which those points are derived, can be stored
 * in the calling application if desired.
 *
 * In parallel, there is a separation between the proc deciding which points to search for (the "target" proc), 
 * and the proc locating the point in its local mesh (the "source" proc).  On the source proc, location 
 * information is cached in locTable, as in the serial case.  By default, this location information (handle and
 * parametric coords) is not returned to the target proc, since it would be of no use there.  Instead, the rank
 * of the source proc locating the point, and the index of that location info in the source proc's locTable, is
 * returned; this information is stored on the target proc in this class's parLocTable variable.  Again, 
 * information about the points searched should be stored in the calling application, if desired.
 *
 * This class uses the ElemEvaluator class for specification and evaluation of basis functions (used for computing
 * parametric coords within an entity).  See documentation and examples for that class for usage information.
 * 
 */

#ifndef MOAB_SPATIALLOCATOR_HPP
#define MOAB_SPATIALLOCATOR_HPP

#include "moab/Types.hpp"
#include "moab/Tree.hpp"
#include "moab/Range.hpp"
#include "moab/TupleList.hpp"

#include <string>
#include <vector>
#include <math.h>

namespace moab {

    class Interface;
    class ElemEvaluator;

    class SpatialLocator
    {
  public:
        /* constructor */
      SpatialLocator(Interface *impl, Range &elems, Tree *tree = NULL, ElemEvaluator *eval = NULL);

        /* destructor */
      virtual ~SpatialLocator();

        /* add elements to be searched */
      ErrorCode add_elems(Range &elems);
      
        /* get bounding box of this locator */
      ErrorCode get_bounding_box(BoundBox &box);
      
        /* locate a set of vertices, Range variant */
      ErrorCode locate_points(Range &vertices,
                              EntityHandle *ents, double *params, bool *is_inside = NULL,
                              double rel_tol = 0.0, double abs_tol = 0.0);
      
        /* locate a set of points */
      ErrorCode locate_points(const double *pos, int num_points,
                              EntityHandle *ents, double *params, bool *is_inside = NULL,
                              double rel_tol = 0.0, double abs_tol = 0.0);
      
        /* locate a set of vertices or entity centroids, storing results on TupleList in this class
         * Locate a set of vertices or entity centroids, storing the detailed results in member 
         * variable (TupleList) locTable (see comments on locTable for structure of that tuple).
         */
      ErrorCode locate_points(Range &ents,
                              double rel_tol = 0.0, double abs_tol = 0.0);
      
        /* locate a set of points, storing results on TupleList in this class
         * Locate a set of points, storing the detailed results in member variable (TupleList) locTable
         * (see comments on locTable for structure of that tuple).
         */
      ErrorCode locate_points(const double *pos, int num_points,
                              double rel_tol = 0.0, double abs_tol = 0.0);

#ifdef USE_MPI      
        /* locate a set of vertices or entity centroids, storing results on TupleList in this class
         * Locate a set of vertices or entity centroids, storing the detailed results in member 
         * variables (TupleList) locTable and parLocTable (see comments on locTable and parLocTable for 
         * structure of those tuples).
         */
      ErrorCode par_locate_points(Range &vertices,
                                  double rel_tol = 0.0, double abs_tol = 0.0);
      
        /* locate a set of points, storing results on TupleList in this class
         * Locate a set of points, storing the detailed results in member 
         * variables (TupleList) locTable and parLocTable (see comments on locTable and parLocTable for 
         * structure of those tuples).
         */
      ErrorCode par_locate_points(const double *pos, int num_points,
                                  double rel_tol = 0.0, double abs_tol = 0.0);
#endif

        /* return the tree */
      Tree *get_tree() {return myTree;}

        /* get the locTable
         */
      TupleList &loc_table() {return locTable;}
      
        /* get the locTable
         */
      const TupleList &loc_table() const {return locTable;}
      
        /* get the parLocTable
         */
      TupleList &par_loc_table() {return parLocTable;}
      
        /* get the parLocTable
         */
      const TupleList &par_loc_table() const {return parLocTable;}

        /* get elemEval */
      ElemEvaluator *elem_eval() {return elemEval;}
      
        /* get elemEval */
      const ElemEvaluator *elem_eval() const {return elemEval;}
      
        /* set elemEval */
      void elem_eval(ElemEvaluator *eval) {elemEval = eval;}
      
  private:

        /* locate a point */
      ErrorCode locate_point(const double *pos, 
                             EntityHandle &ent, double *params, bool *is_inside = NULL,
                             double rel_tol = 0.0, double abs_tol = 0.0);

        /* MOAB instance */
      Interface* mbImpl;

        /* elements being located */
      Range myElems;

        /* dimension of entities in locator */
      int myDim;
      
        /* tree used for location */
      Tree *myTree;
      
        /* element evaluator */
      ElemEvaluator *elemEval;

        /* whether I created the tree or not (determines whether to delete it or not on destruction) */
      bool iCreatedTree;

        /* \brief local locations table
         * This table stores detailed local location data results from locate_points, that is, location data
         * for points located on the local mesh.  Data is stored
         * in a TupleList, where each tuple consists of (p_i, hs_ul, r[3]_d), where
         *   p_i = (int) proc from which request for this point was made (0 if serial)
         *   hs_ul = (unsigned long) source entity containing the point
         *   r[3]_d = (double) parametric coordinates of the point in hs 
         */
      TupleList locTable;

        /* \brief parallel locations table
         * This table stores information about points located on a local or remote processor.  For 
         * points located on this processor's local mesh, detailed location data is stored in locTable.
         * For points located on remote processors, more communication is required to retrieve specific
         * location data (usually that information isn't useful on this processor).
         *
         * The tuple structure of this TupleList is (p_i, ri_i), where:
         *   p_i = (int) processor rank containing this point
         *   ri_i = (int) index into locTable on remote proc containing this point's location information
         * The indexing of parLocTable corresponds to that of the points/entities passed in.
         */
      TupleList parLocTable;
    };

    inline SpatialLocator::~SpatialLocator() 
    {
      if (iCreatedTree && myTree) delete myTree;
    }
    
    inline ErrorCode SpatialLocator::locate_point(const double *pos, 
                                                  EntityHandle &ent, double *params, bool *is_inside, 
                                                  double rel_tol, double abs_tol) 
    {
      return locate_points(pos, 1, &ent, params, is_inside, rel_tol, abs_tol);
    }

    inline ErrorCode SpatialLocator::get_bounding_box(BoundBox &box) 
    {
      return myTree->get_bounding_box(box);
    }

} // namespace moab 

#endif
