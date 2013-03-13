#ifndef ELEM_EVALUATOR_HPP
#define ELEM_EVALUATOR_HPP

#include <vector>

#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"
#include "moab/Matrix3.hpp"
#include "moab/CN.hpp"

namespace moab {
    typedef ErrorCode (*EvalFcn)(const double *params, const double *field, const int ndim, const int num_tuples, 
                          double *work, double *result);

    typedef ErrorCode (*ReverseEvalFcn)(const double *posn, const double *verts, const int nverts, const int ndim,
                                         const double tol, double *work, double *params, bool *is_inside);
        
    typedef ErrorCode (*JacobianFcn)(const double *params, const double *verts, const int nverts, const int ndim, 
                                     double *work, double *result);
        
    typedef ErrorCode (*IntegrateFcn)(const double *field, const double *verts, const int nverts, const int ndim,
                                      const int num_tuples, double *work, double *result);

    typedef ErrorCode (*InitFcn)(const double *verts, const int nverts, double *&work);

    class EvalSet
    {
  public:
        /** \brief Forward-evaluation of field at parametric coordinates */
      EvalFcn evalFcn;
        
        /** \brief Reverse-evaluation of parametric coordinates at physical space position */
      ReverseEvalFcn reverseEvalFcn;
        
        /** \brief Evaluate the jacobian at a specified parametric position */
      JacobianFcn jacobianFcn;
        
        /** \brief Forward-evaluation of field at parametric coordinates */
      IntegrateFcn integrateFcn;

        /** \brief Initialization function for an element */
      InitFcn initFcn;

        /** \brief Bare constructor */
      EvalSet() : evalFcn(NULL), reverseEvalFcn(NULL), jacobianFcn(NULL), integrateFcn(NULL), initFcn(NULL) {}

        /** \brief Constructor */
      EvalSet(EvalFcn eval, ReverseEvalFcn rev, JacobianFcn jacob, IntegrateFcn integ, InitFcn initf)
              : evalFcn(eval), reverseEvalFcn(rev), jacobianFcn(jacob), integrateFcn(integ), initFcn(initf)
          {}
        
        /** \brief Operator= */
      EvalSet &operator=(const EvalSet &eval) {
        evalFcn = eval.evalFcn;
        reverseEvalFcn = eval.reverseEvalFcn;
        jacobianFcn = eval.jacobianFcn;
        integrateFcn = eval.integrateFcn;
        initFcn = eval.initFcn;
        return *this;
      }

      static ErrorCode evaluate_reverse(EvalFcn eval, JacobianFcn jacob,
                                        const double *posn, const double *verts, const int nverts, 
                                        const int ndim, const double tol, double *work, double *params, 
                                        bool *inside);
      
    };

/**\brief Class facilitating local discretization-related functions
 * \class ElemEvaluator
 * This class implements discretization-related functionality operating
 * on data in MOAB.  A member of this class caches certain data about the element
 * it's currently operating on, but is not meant to be instantiated one-per-element,
 * but rather one-per-search (or other operation on a collection of elements).
 *
 * Actual discretization functionality is accessed through function pointers,
 * allowing applications to specialize the implementation of specific functions
 * while still using this class.
 *
 * This class depends on MOAB functionality for gathering entity-based data; the functions
 * it calls through function pointers depend only on POD (plain old data, or intrinsic data types).
 * This allows the use of other packages for serving these functions, without having to modify
 * them to get data through MOAB.  This should also promote efficiency, since in many cases they
 * will be able to read data from its native storage locations.
 */

    class ElemEvaluator {
  public:
        /** \brief Constructor */
      ElemEvaluator(Interface *impl, EntityHandle ent = 0, Tag tag = 0, int tag_dim = -1);

        /** \brief Evaluate cached tag at a given parametric location within the cached entity */
      ErrorCode eval(const double *params, double *result, int num_tuples = -1) const;
        
        /** \brief Evaluate cached tag at a given parametric location within the cached entity */
      ErrorCode eval(const CartVect &params, double *result, int num_tuples = -1) const;
        
        /** \brief Reverse-evaluate the cached entity at a given physical position */
      ErrorCode reverse_eval(const double *posn, double tol, double *params, bool *is_inside = NULL) const;
        
        /** \brief Reverse-evaluate the cached entity at a given physical position */
      ErrorCode reverse_eval(const CartVect &posn, double tol, CartVect &params, bool *is_inside = NULL) const;
        
        /** \brief Evaluate the jacobian of the cached entity at a given parametric location */
      ErrorCode jacobian(const double *params, double *result) const;
        
        /** \brief Evaluate the jacobian of the cached entity at a given parametric location */
      ErrorCode jacobian(const CartVect &params, Matrix3 &result) const;
        
        /** \brief Integrate the cached tag over the cached entity */
      ErrorCode integrate(double *result) const;
      
        /** \brief Return whether a physical position is inside the cached entity to within a tolerance */
      bool is_inside(const double *params, double tol) const;

        /** \brief Return whether a physical position is inside the cached entity to within a tolerance */
      bool is_inside(const CartVect &params, double tol) const;

        /** \brief Set the eval set for a given type entity */
      void set_eval_set(EntityType tp, const EvalSet &eval_set);
      
        /** \brief Get the eval set for a given type entity */
      inline EvalSet get_eval_set(EntityType tp) {return evalSets[tp];}
      
      inline EntityHandle get_ent_handle() const {return entHandle;};
      inline ErrorCode set_ent_handle(EntityHandle ent);

      inline const EntityHandle *get_vert_handles() const {return vertHandles;}
      inline int get_num_verts() const {return numVerts;}
      
      inline Tag get_tag_handle() const {return tagHandle;};
      inline ErrorCode set_tag_handle(Tag tag, int tag_dim = -1);
      inline ErrorCode set_tag(const char *tag_name, int tag_dim = -1);
      
      inline int get_tag_dim() const {return tagDim;};
      inline ErrorCode set_tag_dim(int dim);

      inline Interface *get_moab() {return mbImpl;}
      
  private:

        /** \brief Interface */
      Interface *mbImpl;
      
        /** \brief Entity handle being evaluated */
      EntityHandle entHandle;

        /** \brief Entity type */
      EntityType entType;

        /** \brief Entity dimension */
      int entDim;
            
        /** \brief Number of vertices cached here */
      int numVerts;

        /** \brief Cached copy of vertex handle ptr */
      const EntityHandle *vertHandles;
      
        /** \brief Cached copy of vertex positions */
      CartVect vertPos[CN::MAX_NODES_PER_ELEMENT];
      
        /** \brief Tag being evaluated */
      Tag tagHandle;

        /** \brief Whether tag is coordinates or something else */
      bool tagCoords;
      
        /** \brief Number of values in this tag */
      int numTuples;

        /** \brief Dimension of entities from which to grab tag */
      int tagDim;

        /** \brief Tag space */
      std::vector<unsigned char> tagSpace;
      
        /** \brief Evaluation methods for elements of various topologies */
      EvalSet evalSets[MBMAXTYPE];

        /** \brief Work space for element-specific data */
      double *workSpace;

    }; // class ElemEvaluator

    inline ElemEvaluator::ElemEvaluator(Interface *impl, EntityHandle ent, Tag tag, int tag_dim) 
            : mbImpl(impl), entHandle(0), entType(MBMAXTYPE), entDim(-1), numVerts(0), 
              vertHandles(NULL), tagHandle(0), tagCoords(false), numTuples(0), 
              tagDim(0), workSpace(NULL)
    {
      if (ent) set_ent_handle(ent);
      if (tag) set_tag_handle(tag, tag_dim);
    }
    
    inline ErrorCode ElemEvaluator::set_ent_handle(EntityHandle ent) 
    {
      entHandle = ent;
      if (workSpace) {
        delete [] workSpace;
        workSpace = NULL;
      }

      entType = mbImpl->type_from_handle(ent);
      entDim = mbImpl->dimension_from_handle(ent);

      ErrorCode rval = mbImpl->get_connectivity(ent, vertHandles, numVerts);
      if (MB_SUCCESS != rval) return rval;
      rval = mbImpl->get_coords(vertHandles, numVerts, vertPos[0].array());
      if (MB_SUCCESS != rval) return rval;
      if (tagHandle) {
        rval = set_tag_handle(tagHandle);
        if (MB_SUCCESS != rval) return rval;
      }

      if (evalSets[entType].initFcn) return (*evalSets[entType].initFcn)(vertPos[0].array(), numVerts, workSpace);
      return MB_SUCCESS;
    }
    
    inline ErrorCode ElemEvaluator::set_tag_handle(Tag tag, int tag_dim) 
    {
      ErrorCode rval = MB_SUCCESS;
      if (!tag && !tag_dim) {
        tagCoords = true;
        numTuples = 3;
        tagDim = 0;
        tagHandle = 0;
        return rval;
      }
      else if (tagHandle != tag) {
        tagHandle = tag;
        rval = mbImpl->tag_get_length(tagHandle, numTuples);
        if (MB_SUCCESS != rval) return rval;
        int sz;
        rval = mbImpl->tag_get_bytes(tag, sz);
        if (MB_SUCCESS != rval) return rval;
        tagSpace.reserve(CN::MAX_NODES_PER_ELEMENT*sz);
        tagCoords = false;
      }

      tagDim = (-1 == tag_dim ? entDim : tag_dim);
      
      if (entHandle) {
        if (0 == tagDim) {
          rval = mbImpl->tag_get_data(tagHandle, vertHandles, numVerts, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
        else if (tagDim == entDim) {
          rval = mbImpl->tag_get_data(tagHandle, &entHandle, 1, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
      }

      return rval;
    }

    inline ErrorCode ElemEvaluator::set_tag(const char *tag_name, int tag_dim) 
    {
      ErrorCode rval = MB_SUCCESS;
      Tag tag;
      if (!strcmp(tag_name, "COORDS")) {
        tagCoords = true;
        tagDim = 0;
        numTuples = 3;
        tagHandle = 0;
          // can return here, because vertex coords already cached when entity handle set
        return rval;
      }
      else {
        rval = mbImpl->tag_get_handle(tag_name, tag);
        if (MB_SUCCESS != rval) return rval;
      
        if (tagHandle != tag) {
          tagHandle = tag;
          rval = mbImpl->tag_get_length(tagHandle, numTuples);
          if (MB_SUCCESS != rval) return rval;
          int sz;
          rval = mbImpl->tag_get_bytes(tag, sz);
          if (MB_SUCCESS != rval) return rval;
          tagSpace.reserve(CN::MAX_NODES_PER_ELEMENT*sz);
          tagCoords = false;
        }

        tagDim = (-1 == tag_dim ? entDim : tag_dim);
      }
      
      if (entHandle) {
        if (0 == tagDim) {
          rval = mbImpl->tag_get_data(tagHandle, vertHandles, numVerts, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
        else if (tagDim == entDim) {
          rval = mbImpl->tag_get_data(tagHandle, &entHandle, 1, &tagSpace[0]);
          if (MB_SUCCESS != rval) return rval;
        }
      }

      return rval;
    }

    inline void ElemEvaluator::set_eval_set(EntityType tp, const EvalSet &eval_set) 
    {
      evalSets[tp] = eval_set;
      if (entHandle && evalSets[entType].initFcn) {
        ErrorCode rval = (*evalSets[entType].initFcn)(vertPos[0].array(), numVerts, workSpace);
        if (MB_SUCCESS != rval) throw rval;
      }
    }
    
    inline ErrorCode ElemEvaluator::eval(const double *params, double *result, int num_tuples) const 
    {
      assert(entHandle && MBMAXTYPE != entType);
      return (*evalSets[entType].evalFcn)(params, 
                                          (tagCoords ? (const double*) vertPos[0].array(): (const double*)&tagSpace[0]), 
                                          entDim, (-1 == num_tuples ? numTuples : num_tuples), 
                                          workSpace, result);
    }
        
    inline ErrorCode ElemEvaluator::reverse_eval(const double *posn, const double tol, double *params, bool *ins) const
    {
      assert(entHandle && MBMAXTYPE != entType);
      return EvalSet::evaluate_reverse(evalSets[entType].evalFcn, evalSets[entType].jacobianFcn, 
                                       posn, vertPos[0].array(), numVerts, entDim, tol, workSpace, 
                                       params, ins);
    }
        
      /** \brief Evaluate the jacobian of the cached entity at a given parametric location */
    inline ErrorCode ElemEvaluator::jacobian(const double *params, double *result) const
    {
      assert(entHandle && MBMAXTYPE != entType);
      return (*evalSets[entType].jacobianFcn)(params, vertPos[0].array(), numVerts, entDim, workSpace, result);
    }
        
      /** \brief Integrate the cached tag over the cached entity */
    inline ErrorCode ElemEvaluator::integrate(double *result) const
    {
      assert(entHandle && MBMAXTYPE != entType && (tagCoords || tagHandle));
      ErrorCode rval = MB_SUCCESS;
      if (!tagCoords) {
        if (0 == tagDim) rval = mbImpl->tag_get_data(tagHandle, vertHandles, numVerts, (void*)&tagSpace[0]);
        else rval = mbImpl->tag_get_data(tagHandle, &entHandle, 1, (void*)&tagSpace[0]);
        if (MB_SUCCESS != rval) return rval;
      }
      return (*evalSets[entType].integrateFcn)((tagCoords ? vertPos[0].array() : (const double *)&tagSpace[0]), 
                                               vertPos[0].array(), numVerts, entDim, numTuples, 
                                               workSpace, result);
    }

        /** \brief Evaluate cached tag at a given parametric location within the cached entity */
    inline ErrorCode ElemEvaluator::eval(const CartVect &params, double *result, int num_tuples) const 
    {
      return eval(params.array(), result, num_tuples);
    }
        
        /** \brief Reverse-evaluate the cached entity at a given physical position */
    inline ErrorCode ElemEvaluator::reverse_eval(const CartVect &posn, double tol, CartVect &params, bool *is_in) const
    {
      return reverse_eval(posn.array(), tol, params.array(), is_in);
    }
        
        /** \brief Evaluate the jacobian of the cached entity at a given parametric location */
    inline ErrorCode ElemEvaluator::jacobian(const CartVect &params, Matrix3 &result) const
    {
      return jacobian(params.array(), result.array());
    }
        
        /** \brief Return whether a physical position is inside the cached entity to within a tolerance */
    inline bool ElemEvaluator::is_inside(const CartVect &params, double tol) const
    {
      return (params[0]+tol >= -1.0 && params[0]-tol <= 1.0 &&
              params[1]+tol >= -1.0 && params[1]-tol <= 1.0 &&
              params[2]+tol >= -1.0 && params[2]-tol <= 1.0);
    }

      
} // namespace moab

#endif /*MOAB_ELEM_EVALUATOR_HPP*/
