/**\file SpatialLocator.hpp
 * \class moab::SpatialLocator
 * \brief Tool to facilitate spatial location of a point in a mesh
 */

#ifndef MOAB_SPATIALLOCATOR_HPP
#define MOAB_SPATIALLOCATOR_HPP

#include "moab/Types.hpp"
#include "moab/Tree.hpp"
#include "moab/Range.hpp"

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
      
        /* locate a set of points */
      ErrorCode locate_points(const double *pos, int num_points,
                              EntityHandle *ents, double *params, 
                              double rel_tol = 0.0, double abs_tol = 0.0,
                              bool *is_inside = NULL);
      
        /* locate a point */
      ErrorCode locate_point(const double *pos, 
                             EntityHandle &ent, double *params, 
                             double rel_tol = 0.0, double abs_tol = 0.0,
                             bool *is_inside = NULL);

        /* return the tree */
      Tree *get_tree() {return myTree;}
      
  private:

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
    };

    inline SpatialLocator::~SpatialLocator() 
    {
      if (iCreatedTree && myTree) delete myTree;
    }
    
    inline ErrorCode SpatialLocator::locate_point(const double *pos, 
                                                  EntityHandle &ent, double *params, 
                                                  double rel_tol, double abs_tol,
                                                  bool *is_inside) 
    {
      return locate_points(pos, 1, &ent, params, rel_tol, abs_tol, is_inside);
    }

    inline ErrorCode SpatialLocator::get_bounding_box(BoundBox &box) 
    {
      return myTree->get_bounding_box(box);
    }

} // namespace moab 

#endif
