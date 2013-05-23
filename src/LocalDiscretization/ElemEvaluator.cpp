#include <limits>

#include "moab/ElemEvaluator.hpp"
#include "moab/CartVect.hpp"
#include "moab/Matrix3.hpp"

// need to include eval set types here to support get_eval_set; alternative would be to have some
// type of registration, but we'd still need static registration for the built-in types
#include "moab/LinearQuad.hpp"
#include "moab/LinearTet.hpp"
#include "moab/LinearHex.hpp"
#include "moab/QuadraticHex.hpp"
//#include "moab/SpectralQuad.hpp"
//#include "moab/SpectralHex.hpp"

namespace moab { 
    ErrorCode EvalSet::evaluate_reverse(EvalFcn eval, JacobianFcn jacob, InsideFcn inside_f,
                                        const double *posn, const double *verts, const int nverts, 
                                        const int ndim, const double tol, double *work, 
                                        double *params, bool *inside) {
        // TODO: should differentiate between epsilons used for
        // Newton Raphson iteration, and epsilons used for curved boundary geometry errors
        // right now, fix the tolerance used for NR
      const double error_tol_sqr = tol*tol;
      CartVect *cvparams = reinterpret_cast<CartVect*>(params);
      const CartVect *cvposn = reinterpret_cast<const CartVect*>(posn);

        // initialize to center of element
      *cvparams = CartVect(-.4);
  
      CartVect new_pos;
        // evaluate that first guess to get a new position
      ErrorCode rval = (*eval)(cvparams->array(), verts, ndim, ndim, work, new_pos.array());
      if (MB_SUCCESS != rval) return rval;
      
        // residual is diff between old and new pos; need to minimize that
      CartVect res = new_pos - *cvposn;
      Matrix3 J;

      int iters=0;
        // while |res| larger than tol
      while (res % res > error_tol_sqr) {
        if(++iters>10)
          return MB_FAILURE;

          // get jacobian at current params
        rval = (*jacob)(cvparams->array(), verts, nverts, ndim, work, J[0]);
        double det = J.determinant();
        assert(det > std::numeric_limits<double>::epsilon());

          // new params tries to eliminate residual
        *cvparams -= J.inverse(1.0/det) * res;

          // get the new forward-evaluated position, and its difference from the target pt
        rval = (*eval)(params, verts, ndim, ndim, work, new_pos.array());
        if (MB_SUCCESS != rval) return rval;
        res = new_pos - *cvposn;
      }

      if (inside)
        *inside = (*inside_f)(params, ndim, tol);

      return MB_SUCCESS;
    }// Map::evaluate_reverse()

    bool EvalSet::inside_function(const double *params, const int ndims, const double tol) 
    {
      if (params[0] >= -1-tol && params[0] <= 1+tol &&
          (ndims < 2 || (params[1] >= -1-tol && params[1] <= 1+tol)) &&
          (ndims < 3 || (params[2] >= -1-tol && params[2] <= 1+tol))) 
        return true;
      else return false;
    }

        /** \brief Given type & #vertices, get an appropriate eval set */
    ErrorCode EvalSet::get_eval_set(EntityType tp, unsigned int num_vertices, EvalSet &eval_set) 
    {
      switch (tp) {
        case MBEDGE:
            break;
        case MBTRI:
            break;
        case MBQUAD:
            if (LinearQuad::compatible(tp, num_vertices, eval_set)) return MB_SUCCESS;
//            if (SpectralQuad::compatible(tp, num_vertices, eval_set)) return MB_SUCCESS;
            break;
        case MBTET:
            if (LinearTet::compatible(tp, num_vertices, eval_set)) return MB_SUCCESS;
            break;
        case MBHEX:
            if (LinearHex::compatible(tp, num_vertices, eval_set)) return MB_SUCCESS;
            if (QuadraticHex::compatible(tp, num_vertices, eval_set)) return MB_SUCCESS;
//            if (SpectralHex::compatible(tp, num_vertices, eval_set)) return MB_SUCCESS;
            break;
        default:
            break;
      }

      return MB_NOT_IMPLEMENTED;
    }
      
} // namespace moab
