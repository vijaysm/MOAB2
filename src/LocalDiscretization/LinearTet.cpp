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
        for (int j = 0; j < num_tuples; j++)
          result[j] += (field[i*num_tuples+j]-f0[j])*params[i-1];
      }

      return MB_SUCCESS;
    }

    ErrorCode LinearTet::integrateFcn(const double *field, const double */*verts*/, const int nverts, const int num_tuples, const int /*ndim*/,
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
    
    ErrorCode LinearTet::reverseEvalFcn(const double *posn, const double *verts, const int nverts, const int ndim,
                                        const double tol, double *work, double *params, bool *is_inside) 
    {
      assert(posn && verts);
      return EvalSet::evaluate_reverse(evalFcn, jacobianFcn, posn, verts, nverts, ndim, tol, work, params, is_inside);
    } 

} // namespace moab
