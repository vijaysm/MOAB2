#include "moab/LinearQuad.hpp"
#include "moab/Matrix3.hpp"
#include "moab/Forward.hpp"

namespace moab 
{
    
    const double LinearQuad::corner[4][2] = {{ -1, -1},
                                             {  1, -1},
                                             {  1,  1},
                                             { -1,  1} };

      /* For each point, its weight and location are stored as an array.
         Hence, the inner dimension is 2, the outer dimension is gauss_count.
         We use a one-point Gaussian quadrature, since it integrates linear functions exactly.
      */
    const double LinearQuad::gauss[1][2] = { {  2.0,           0.0          } };

    ErrorCode LinearQuad::jacobianFcn(const double *params, const double *verts, const int /*nverts*/, const int /*ndim*/, 
                                      double *, double *result) 
    {
      Matrix3 *J = reinterpret_cast<Matrix3*>(result);
      *J = Matrix3(0.0);
      for (unsigned i = 0; i < 4; ++i) {
        const double   xi_p = 1 + params[0]*corner[i][0];
        const double  eta_p = 1 + params[1]*corner[i][1];
        const double dNi_dxi   = corner[i][0] * eta_p;
        const double dNi_deta  = corner[i][1] * xi_p;
        (*J)(0,0) += dNi_dxi   * verts[i*3+0];
        (*J)(1,0) += dNi_dxi   * verts[i*3+1];
        (*J)(0,1) += dNi_deta  * verts[i*3+0];
        (*J)(1,1) += dNi_deta  * verts[i*3+1];
      }
      (*J) *= 0.25;
      (*J)(2,2) = 1.0; /* to make sure the Jacobian determinant is non-zero */
      return MB_SUCCESS;
    }// LinearQuad::jacobian()

    ErrorCode LinearQuad::evalFcn(const double *params, const double *field, const int /*ndim*/, const int num_tuples, 
                                  double *, double *result) {
      for (int i = 0; i < num_tuples; i++) result[i] = 0.0;
      for (unsigned i = 0; i < 4; ++i) {
        const double N_i = (1 + params[0]*corner[i][0])
            * (1 + params[1]*corner[i][1]);
        for (int j = 0; j < num_tuples; j++) result[j] += N_i * field[i*num_tuples+j];
      }
      for (int i = 0; i < num_tuples; i++) result[i] *= 0.25;

      return MB_SUCCESS;
    }

    ErrorCode LinearQuad::integrateFcn(const double *field, const double *verts, const int nverts, const int ndim, 
                                       const int num_tuples, double *work, double *result) {
      double tmp_result[4];
      ErrorCode rval = MB_SUCCESS;
      for (int i = 0; i < num_tuples; i++) result[i] = 0.0;
      CartVect x;
      Matrix3 J;
      for(unsigned int j1 = 0; j1 < LinearQuad::gauss_count; ++j1) {
        x[0] = LinearQuad::gauss[j1][1];
        double w1 = LinearQuad::gauss[j1][0];
        for(unsigned int j2 = 0; j2 < LinearQuad::gauss_count; ++j2) {
          x[1] = LinearQuad::gauss[j2][1];
          double w2 = LinearQuad::gauss[j2][0];
          rval = evalFcn(x.array(), field, ndim, num_tuples, NULL, tmp_result);
          if (MB_SUCCESS != rval) return rval;
          rval = jacobianFcn(x.array(), verts, nverts, ndim, work, J[0]);
          if (MB_SUCCESS != rval) return rval;
          double tmp_det =  w1*w2*J.determinant();
          for (int i = 0; i < num_tuples; i++) result[i] += tmp_result[i]*tmp_det;
        }
      }
      return MB_SUCCESS;
    } // LinearHex::integrate_vector()

    ErrorCode LinearQuad::reverseEvalFcn(EvalFcn eval, JacobianFcn jacob, InsideFcn ins, 
                                         const double *posn, const double *verts, const int nverts, const int ndim,
                                         const double tol, double *work, double *params, bool *is_inside) 
    {
      return EvalSet::evaluate_reverse(eval, jacob, ins, posn, verts, nverts, ndim, tol, work, params, is_inside);
    } 

    bool LinearQuad::insideFcn(const double *params, const int ndim, const double tol) 
    {
      return EvalSet::inside_function(params, ndim, tol);
    }
    
} // namespace moab
