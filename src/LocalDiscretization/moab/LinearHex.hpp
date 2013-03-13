#ifndef LINEAR_HEX_HPP
#define LINEAR_HEX_HPP
  /**\brief Shape function space for trilinear hexahedron, obtained by a pushforward of the canonical linear (affine) functions. */

#include "moab/ElemEvaluator.hpp"

namespace moab 
{
    
class LinearHex 
{
public:
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode evalFcn(const double *params, const double *field, const int ndim, const int num_tuples, 
                           double *work, double *result);
        
    /** \brief Reverse-evaluation of parametric coordinates at physical space position */
  static ErrorCode reverseEvalFcn(const double *posn, const double *verts, const int nverts, const int ndim,
                                  const double tol, double *work, double *params, bool *is_inside);
        
    /** \brief Evaluate the jacobian at a specified parametric position */
  static ErrorCode jacobianFcn(const double *params, const double *verts, const int nverts, const int ndim, 
                               double *work, double *result);
        
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode integrateFcn(const double *field, const double *verts, const int nverts, const int num_tuples, const int ndim, 
                                double *work, double *result);

  static EvalSet eval_set() 
      {
        return EvalSet(evalFcn, reverseEvalFcn, jacobianFcn, integrateFcn, (InitFcn)NULL);
      }
      
protected:
    /* Preimages of the vertices -- "canonical vertices" -- are known as "corners". */
  static const double corner[8][3];
  static const double gauss[1][2];
  static const unsigned int corner_count = 8;
  static const unsigned int gauss_count  = 1;
      
};// class LinearHex

} // namespace moab

#endif
