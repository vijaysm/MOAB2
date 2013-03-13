#ifndef SPECTRAL_HEX_HPP
#define SPECTRAL_HEX_HPP
  /**\brief Shape function space for spectral hexahedron
   */

#include "EvalSet.hpp"

namespace moab 
{
    
class SpectralHex 
{
public:
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode evalFcn(const double *params, const double *field, const int ndim, const int num_tuples, 
                           double *work, double *result);
        
    /** \brief Reverse-evaluation of parametric coordinates at physical space position */
  static ErrorCode reverseEvalFcn(const double *posn, const double *verts, const int nverts, const int ndim,
                                  double tol, double *work, double *params, bool *is_inside);
        
    /** \brief Evaluate the jacobian at a specified parametric position */
  static ErrorCode jacobianFcn(const double *params, const double *verts, const int nverts, const int ndim, 
                               double *work, double *result);
        
    /** \brief Forward-evaluation of field at parametric coordinates */
  static ErrorCode integrateFcn(const double *field, const double *verts, const int nverts, const int ndim,
                                const int num_tuples, double *work, double *result);

    /** \brief Initialize this EvalSet */
  static ErrorCode initFcn(const double *verts, const int nverts, double *&work);
      
  static EvalSet eval_set() 
      {
        return EvalSet(evalFcn, reverseEvalFcn, jacobianFcn, integrateFcn, initFcn);
      }
      
protected:
  static int _n;
  static real *_z[3];
  static lagrange_data _ld[3];
  static opt_data_3 _data;
  static real * _odwork;// work area
  static bool init_;
};// class SpectralHex

} // namespace moab

#endif
