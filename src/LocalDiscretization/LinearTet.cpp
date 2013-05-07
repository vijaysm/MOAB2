#include "moab/LinearTet.hpp"
#include "moab/Forward.hpp"
#include <algorithm>

namespace moab 
{
    
    const double LinearTet::corner[4][3] = { {0,0,0},
                                             {1,0,0},
                                             {0,1,0},
                                             {0,0,1}};

    ErrorCode LinearTet::initFcn(const double *verts, const int /*nverts*/, double *&work) {
        // allocate work array as: 
        // work[0..8] = T
        // work[9..17] = Tinv
        // work[18] = detT
        // work[19] = detTinv
      assert(!work && verts);
      work = new double[20];
      Matrix3 *T = reinterpret_cast<Matrix3*>(work),
          *Tinv = reinterpret_cast<Matrix3*>(work+9);
      double *detT = work+18, *detTinv = work+19;
      
      *T = Matrix3(verts[1*3+0]-verts[0*3+0],verts[2*3+0]-verts[0*3+0],verts[3*3+0]-verts[0*3+0],
                   verts[1*3+1]-verts[0*3+1],verts[2*3+1]-verts[0*3+1],verts[3*3+1]-verts[0*3+1],
                   verts[1*3+2]-verts[0*3+2],verts[2*3+2]-verts[0*3+2],verts[3*3+2]-verts[0*3+2]);
      *Tinv = T->inverse();
      *detT = T->determinant();
      *detTinv = (0.0 == *detT ? HUGE : 1.0 / *detT);

      return MB_SUCCESS;
    }

    ErrorCode LinearTet::evalFcn(const double *params, const double *field, const int /*ndim*/, const int num_tuples, 
                                 double */*work*/, double *result) {
      assert(params && field && num_tuples > 0);
      std::vector<double> f0(num_tuples);
      std::copy(field, field+num_tuples, f0.begin());
      std::copy(field, field+num_tuples, result);
      
      
      for (unsigned i = 1; i < 4; ++i) {
        double p = 0.5*(params[i-1] + 1); // transform from -1 <= p <= 1 to 0 <= p <= 1
        for (int j = 0; j < num_tuples; j++)
          result[j] += (field[i*num_tuples+j]-f0[j])*p;
      }

      return MB_SUCCESS;
    }

    ErrorCode LinearTet::integrateFcn(const double *field, const double */*verts*/, const int nverts, const int /*ndim*/, const int num_tuples,
                                      double *work, double *result) 
    {
      assert(field && num_tuples > 0);
      std::fill(result, result+num_tuples, 0.0);
      for(int i = 0; i < nverts; ++i) {
        for (int j = 0; j < num_tuples; j++)
          result[j] += field[i*num_tuples+j];
      }
      double tmp = work[18]/24.0;
      for (int i = 0; i < num_tuples; i++) result[i] *= tmp;

      return MB_SUCCESS;
    }

    ErrorCode LinearTet::jacobianFcn(const double *, const double *, const int, const int , 
                                     double *work, double *result) 
    {
        // jacobian is cached in work array
      assert(work);
      std::copy(work, work+9, result);
      return MB_SUCCESS;
    }
    
    ErrorCode LinearTet::reverseEvalFcn(EvalFcn eval, JacobianFcn jacob, InsideFcn ins, 
                                        const double *posn, const double *verts, const int nverts, const int ndim,
                                        const double tol, double *work, double *params, bool *is_inside) 
    {
      assert(posn && verts);
      return evaluate_reverse(eval, jacob, ins, posn, verts, nverts, ndim, tol, work, params, is_inside);
    } 

    bool LinearTet::insideFcn(const double *params, const int , const double tol) 
    {
      return (params[0] >= -1.0-tol && params[1] >= -1.0-tol && params[2] >= -1.0-tol && 
              params[0] + params[1] + params[2] <= -1.0+tol);
      
    }
    
    ErrorCode LinearTet::evaluate_reverse(EvalFcn eval, JacobianFcn jacob, InsideFcn inside_f,
                                          const double *posn, const double *verts, const int nverts, 
                                          const int ndim, const double tol, double *work, 
                                          double *params, bool *inside) {
        // TODO: should differentiate between epsilons used for
        // Newton Raphson iteration, and epsilons used for curved boundary geometry errors
        // right now, fix the tolerance used for NR
      const double error_tol_sqr = tol*tol;
      CartVect *cvparams = reinterpret_cast<CartVect*>(params);
      const CartVect *cvposn = reinterpret_cast<const CartVect*>(posn);

        // find best initial guess to improve convergence
      CartVect tmp_params[] = {CartVect(-1,-1,-1), CartVect(1,-1,-1), CartVect(-1,1,-1), CartVect(-1,-1,1)};
      double resl = HUGE;
      CartVect new_pos, tmp_pos;
      ErrorCode rval;
      for (unsigned int i = 0; i < 4; i++) {
        rval = (*eval)(tmp_params[i].array(), verts, ndim, ndim, work, tmp_pos.array());
        if (MB_SUCCESS != rval) return rval;
        double tmp_resl = (tmp_pos-*cvposn).length_squared();
        if (tmp_resl < resl) {
          *cvparams = tmp_params[i];
          new_pos = tmp_pos;
          resl = tmp_resl;
        }        
      }

        // residual is diff between old and new pos; need to minimize that
      CartVect res = new_pos - *cvposn;
      Matrix3 J;
      rval = (*jacob)(cvparams->array(), verts, nverts, ndim, work, J[0]);
      double det = J.determinant();
      assert(det > std::numeric_limits<double>::epsilon());
      Matrix3 Ji = J.inverse(1.0/det);

      int iters=0;
        // while |res| larger than tol
      while (res % res > error_tol_sqr) {
        if(++iters>25)
          return MB_FAILURE;

          // new params tries to eliminate residual
        *cvparams -= Ji * res;

          // get the new forward-evaluated position, and its difference from the target pt
        rval = (*eval)(params, verts, ndim, ndim, work, new_pos.array());
        if (MB_SUCCESS != rval) return rval;
        res = new_pos - *cvposn;
      }

      if (inside)
        *inside = (*inside_f)(params, ndim, tol);

      return MB_SUCCESS;
    }// Map::evaluate_reverse()

} // namespace moab
