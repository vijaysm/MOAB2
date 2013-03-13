#ifndef SPECTRAL_QUAD_HPP
#define SPECTRAL_QUAD_HPP
  /*\brief Shape function space for spectral quad
   */

#include "EvalSet.hpp"

namespace moab 
{
    
class SpectralQuad
{
public:
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode evalFcn(const double *params, const double *field, const int ndim, const int num_tuples, 
                           double *work, double *result) const;
        
    /** \brief Reverse-evaluation of parametric coordinates at physical space position */
  static ErrorCode reverseEvalFcn(const double *posn, const double *verts, const int nverts, const int ndim,
                                  double *work, double *params) const;
        
    /** \brief Evaluate the jacobian at a specified parametric position */
  static ErrorCode jacobianFcn(const double *params, const double *verts, const int nverts, const int ndim, 
                               double *work, double *result) const;
        
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode integrateFcn(const double *field, const int num_tuples, double *work, double *result) const;

    /** \brief Initialize this EvalSet */
  static ErrorCode initFcn(const EntityHandle ent, double *&work) const;
      
  static EvalSet eval_set() 
      {
        return EvalSet(evalFcn, reverseEvalFcn, jacobianFcn, integrateFcn, initFcn);
      }
      
protected:
  static int _n;
  static real *_z[2];
  static lagrange_data _ld[2];
  static opt_data_2 _data; // we should use only 2nd component
  static real * _odwork;// work area

    // flag for initialization of data
  static bool _init;
  static real * _glpoints; // it is a space we can use to store gl positions for elements
    // on the fly; we do not have a tag yet for them, as in Nek5000 application
    // also, these positions might need to be moved on the sphere, for HOMME grids
    // do we project them or how do we move them on the sphere?
};// class SpectralQuad

} // namespace moab

#endif
