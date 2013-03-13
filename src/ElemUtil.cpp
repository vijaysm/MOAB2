#include <iostream>
#include <limits>
#include <assert.h>

#include "ElemUtil.hpp"
#include "types.h"

namespace moab { 
namespace Element {

    bool Map::evaluate_reverse(const CartVect& x, CartVect &params, double tol, const CartVect& x0) const {
        // TODO: should differentiate between epsilons used for
        // Newton Raphson iteration, and epsilons used for curved boundary geometry errors
        // right now, fix the tolerance used for NR
      const double error_tol_sqr = tol*tol;
      double det;
      params = x0;
      CartVect delta = evaluate(params) - x;
      Matrix3 J;

      int iters=0;
      while (delta % delta > error_tol_sqr) {
        if(++iters>10)
          return false;

        J = jacobian(params);
        det = J.determinant();
        if (det < std::numeric_limits<double>::epsilon())
          return false;
        params -= J.inverse(1.0/det) * delta;
        delta = evaluate( params ) - x;
      }
      return true;
    }// Map::evaluate_reverse()

    const double LinearHex::corner[8][3] = { { -1, -1, -1 },
                                             {  1, -1, -1 },
                                             {  1,  1, -1 },
                                             { -1,  1, -1 },
                                             { -1, -1,  1 },
                                             {  1, -1,  1 },
                                             {  1,  1,  1 },
                                             { -1,  1,  1 } };

      /* For each point, its weight and location are stored as an array.
         Hence, the inner dimension is 2, the outer dimension is gauss_count.
         We use a one-point Gaussian quadrature, since it integrates linear functions exactly.
      */
    const double LinearHex::gauss[1][2] = { {  2.0,           0.0          } };

    CartVect LinearHex::evaluate( const CartVect& params ) const {
      CartVect x(0.0);
      for (unsigned i = 0; i < 8; ++i) {
        const double N_i = 
            (1 + params[0]*corner[i][0])
            * (1 + params[1]*corner[i][1])
            * (1 + params[2]*corner[i][2]);
        x += N_i * this->vertex[i];
      }
      x *= 0.125;
      return x;
    }// LinearHex::evaluate

    Matrix3 LinearHex::jacobian( const CartVect& params ) const {
      Matrix3 J(0.0);
      for (unsigned i = 0; i < 8; ++i) {
        const double   params_p = 1 + params[0]*corner[i][0];
        const double  eta_p = 1 + params[1]*corner[i][1];
        const double zeta_p = 1 + params[2]*corner[i][2];
        const double dNi_dparams   = corner[i][0] * eta_p * zeta_p;
        const double dNi_deta  = corner[i][1] *  params_p * zeta_p;
        const double dNi_dzeta = corner[i][2] *  params_p *  eta_p;
        J(0,0) += dNi_dparams   * vertex[i][0];
        J(1,0) += dNi_dparams   * vertex[i][1];
        J(2,0) += dNi_dparams   * vertex[i][2];
        J(0,1) += dNi_deta  * vertex[i][0];
        J(1,1) += dNi_deta  * vertex[i][1];
        J(2,1) += dNi_deta  * vertex[i][2];
        J(0,2) += dNi_dzeta * vertex[i][0];
        J(1,2) += dNi_dzeta * vertex[i][1];
        J(2,2) += dNi_dzeta * vertex[i][2];
      }
      return J *= 0.125;
    }// LinearHex::jacobian()

    void LinearHex::evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const {
      for (int i = 0; i < num_tuples; i++) eval[i] = 0.0;
      for (unsigned i = 0; i < 8; ++i) {
        const double N_i = (1 + params[0]*corner[i][0])
            * (1 + params[1]*corner[i][1])
            * (1 + params[2]*corner[i][2]);
        for (int j = 0; j < num_tuples; j++) eval[j] += N_i * field_values[i*num_tuples+j];
      }
      for (int i = 0; i < num_tuples; i++) eval[i] *= 0.125;
    }// LinearHex::evaluate_vector()

    void LinearHex::integrate_vector(const double *field_values, int num_tuples, double *integral) const {
      double tmp_integral[8];
      for (int i = 0; i < num_tuples; i++) integral[i] = 0.0;
      CartVect x;
      for(unsigned int j1 = 0; j1 < this->gauss_count; ++j1) {
        x[0] = this->gauss[j1][1];
        double w1 = this->gauss[j1][0];
        for(unsigned int j2 = 0; j2 < this->gauss_count; ++j2) {
          x[1] = this->gauss[j2][1];
          double w2 = this->gauss[j2][0];
          for(unsigned int j3 = 0; j3 < this->gauss_count; ++j3) {
            x[2] = this->gauss[j3][1];
            double w3 = this->gauss[j3][0];
            this->evaluate_vector(x,field_values, num_tuples, tmp_integral);
            double tmp_det =  w1*w2*w3*this->det_jacobian(x);
            for (int i = 0; i < num_tuples; i++) integral[i] += tmp_integral[i]*tmp_det;
          }
        }
      }
    }// LinearHex::integrate_vector()

    bool LinearHex::is_inside(const CartVect & params, double tol) const
    {
        // just look at the box+tol here
      return ( params[0]>=-1.-tol) && (params[0]<=1.+tol) &&
          ( params[1]>=-1.-tol) && (params[1]<=1.+tol) &&
          ( params[2]>=-1.-tol) && (params[2]<=1.+tol);
    }

      // those are not just the corners, but for simplicity, keep this name
      //
    const int QuadraticHex::corner[27][3] = {
        { -1, -1, -1 },
        {  1, -1, -1 },
        {  1,  1, -1 },  // corner nodes: 0-7
        { -1,  1, -1 },  // mid-edge nodes: 8-19
        { -1, -1,  1 },  // center-face nodes 20-25  center node  26
        {  1, -1,  1 },  //
        {  1,  1,  1 },
        { -1,  1,  1 }, //                    4   ----- 19   -----  7
        {  0, -1, -1 }, //                .   |                 .   |
        {  1,  0, -1 }, //            16         25         18      |
        {  0,  1, -1 }, //         .          |          .          |
        { -1,  0, -1 }, //      5   ----- 17   -----  6             |
        { -1, -1,  0 }, //      |            12       | 23         15
        {  1, -1,  0 }, //      |                     |             |
        {  1,  1,  0 }, //      |     20      |  26   |     22      |
        { -1,  1,  0 }, //      |                     |             |
        {  0, -1,  1 }, //     13         21  |      14             |
        {  1,  0,  1 }, //      |             0   ----- 11   -----  3
        {  0,  1,  1 }, //      |         .           |         .
        { -1,  0,  1 }, //      |      8         24   |     10
        {  0, -1,  0 }, //      |  .                  |  .
        {  1,  0,  0 }, //      1   -----  9   -----  2
        {  0,  1,  0 }, //
        { -1,  0,  0 },
        {  0,  0, -1 },
        {  0,  0,  1 },
        {  0,  0,  0 }
    };
      //QuadraticHex::QuadraticHex(const std::vector<CartVect>& vertices) : Map(vertices){};
    QuadraticHex::QuadraticHex():Map(0) {
    }

    double SH(const int i, const double params)
    {
      switch (i)
      {
        case -1: return (params*params-params)/2;
        case 0: return 1-params*params;
        case 1: return (params*params+params)/2;
        default: return 0.;
      }
    }
    double DSH(const int i, const double params)
    {
      switch (i)
      {
        case -1: return params-0.5;
        case 0: return -2*params;
        case 1: return params+0.5;
        default: return 0.;
      }
    }

    CartVect QuadraticHex::evaluate( const CartVect& params ) const
    {

      CartVect x(0.0);
      for (int i=0; i<27; i++)
      {
        const double sh= SH(corner[i][0], params[0])
            *SH(corner[i][1], params[1])
            *SH(corner[i][2], params[2]);
        x+=sh* vertex[i];
      }

      return x;
    }

    bool QuadraticHex::is_inside(const CartVect & params, double tol) const
    {// just look at the box+tol here
      return ( params[0]>=-1.-tol) && (params[0]<=1.+tol) &&
          ( params[1]>=-1.-tol) && (params[1]<=1.+tol) &&
          ( params[2]>=-1.-tol) && (params[2]<=1.+tol);
    }

    Matrix3  QuadraticHex::jacobian(const CartVect& params) const
    {
      Matrix3 J(0.0);
      for (int i=0; i<27; i++)
      {
        const double sh[3]={ SH(corner[i][0], params[0]),
                             SH(corner[i][1], params[1]),
                             SH(corner[i][2], params[2]) };
        const double dsh[3]={ DSH(corner[i][0], params[0]),
                              DSH(corner[i][1], params[1]),
                              DSH(corner[i][2], params[2]) };


        for (int j=0; j<3; j++)
        {
          J(j,0)+=dsh[0]*sh[1]*sh[2]*vertex[i][j]; // dxj/dr first column
          J(j,1)+=sh[0]*dsh[1]*sh[2]*vertex[i][j]; // dxj/ds
          J(j,2)+=sh[0]*sh[1]*dsh[2]*vertex[i][j]; // dxj/dt
        }
      }


      return J;
    }
    void QuadraticHex::evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const
    {
      for (int i = 0; i < num_tuples; i++) eval[i] = 0.0;
      for (int i=0; i<27; i++)
      {
        const double sh= SH(corner[i][0], params[0])
            *SH(corner[i][1], params[1])
            *SH(corner[i][2], params[2]);
        for (int j = 0; j < num_tuples; j++) 
          eval[j] += sh* field_values[i*num_tuples+j];
      }
    }

    void QuadraticHex::integrate_vector(const double *field_vertex_values, int num_tuples, double *integral) const
  {
  }

    const double LinearTet::corner[4][3] = { {0,0,0},
                                             {1,0,0},
                                             {0,1,0},
                                             {0,0,1}};

    LinearTet::LinearTet() : Map(0) {

    }// LinearTet::LinearTet()


    void LinearTet::set_vertices(const CartVect *v, int num_vs) {
      this->Map::set_vertices(v, num_vs);
      this->T = Matrix3(v[1][0]-v[0][0],v[2][0]-v[0][0],v[3][0]-v[0][0],
                        v[1][1]-v[0][1],v[2][1]-v[0][1],v[3][1]-v[0][1],
                        v[1][2]-v[0][2],v[2][2]-v[0][2],v[3][2]-v[0][2]);
      this->T_inverse = this->T.inverse();
      this->det_T = this->T.determinant();
      this->det_T_inverse = (0.0 == this->det_T ? HUGE : 1.0/this->det_T);
    }


    void LinearTet::evaluate_vector(const CartVect& params, const double *field_values, int num_tuples, double *eval) const {
      std::vector<double> f0(num_tuples);
      for (int j = 0; j < num_tuples; j++) {
        f0[j] = field_values[j];
        eval[j] = f0[j];
      }
      for (unsigned i = 1; i < 4; ++i) {
        for (int j = 0; j < num_tuples; j++)
          eval[j] += (field_values[i*num_tuples+j]-f0[j])*params[i-1];
      }
    }

    void LinearTet::integrate_vector(const double *field_values, int num_tuples, double *integral) const {
      for (int i = 0; i < num_tuples; i++) integral[i] = 0.0;
      for(unsigned int i = 0; i < 4; ++i) {
        for (int j = 0; j < num_tuples; j++)
          integral[j] += field_values[i*num_tuples+j];
      }
      double tmp = this->det_T/24.0;
      for (int i = 0; i < num_tuples; i++) integral[i] *= tmp;
    }

    bool LinearTet::is_inside(const CartVect & params, double tol) const
    {
        // linear tet space is a tetra with vertices (0,0,0), (1,0,0), (0,1,0), (0, 0, 1)
        // first check if outside bigger box, then below the plane x+y+z=1
      return ( params[0]>=-tol)  &&
          ( params[1]>=-tol)  &&
          ( params[2]>=-tol)  &&
          ( params[0]+params[1]+params[2] < 1.0+tol);
    }
      // SpectralHex

      // filescope for static member data that is cached
    int SpectralHex::_n;
    double *SpectralHex::_z[3];
    lagrange_data SpectralHex::_ld[3];
    opt_data_3 SpectralHex::_data;
    double * SpectralHex::_odwork;

    bool SpectralHex::_init = false;

    SpectralHex::SpectralHex() : Map(0)
    {
    }
      // the preferred constructor takes pointers to GL blocked positions
    SpectralHex::SpectralHex(int order, double * x, double *y, double *z) : Map(0)
    {
      Init(order);
      _xyz[0]=x; _xyz[1]=y; _xyz[2]=z;
    }
    SpectralHex::SpectralHex(int order) : Map(0)
    {
      Init(order);
      _xyz[0]=_xyz[1]=_xyz[2]=NULL;
    }
    SpectralHex::~SpectralHex()
    {
      if (_init)
        freedata();
      _init=false;
    }
    void SpectralHex::Init(int order)
    {
      if (_init && _n==order)
        return;
      if (_init && _n!=order)
      {
          // TODO: free data cached
        freedata();
      }
        // compute stuff that depends only on order
      _init = true;
      _n = order;
        //triplicates! n is the same in all directions !!!
      for(int d=0; d<3; d++){
        _z[d] = tmalloc(double, _n);
        lobatto_nodes(_z[d], _n);
        lagrange_setup(&_ld[d], _z[d], _n);
      }
      opt_alloc_3(&_data, _ld);

      unsigned int nf = _n*_n, ne = _n, nw = 2*_n*_n + 3*_n;
      _odwork = tmalloc(double, 6*nf + 9*ne + nw);
    }
    void SpectralHex::freedata()
    {
      for(int d=0; d<3; d++){
        free(_z[d]);
        lagrange_free(&_ld[d]);
      }
      opt_free_3(&_data);
      free(_odwork);
    }

    void SpectralHex::set_gl_points( double * x, double * y, double *z)
    {
      _xyz[0] = x;
      _xyz[1] = y;
      _xyz[2] = z;
    }
    CartVect SpectralHex::evaluate( const CartVect& params ) const
    {
        //piece that we shouldn't want to cache
      int d=0;
      for(d=0; d<3; d++){
        lagrange_0(&_ld[d], params[d]);
      }
      CartVect result;
      for (d=0; d<3; d++)
      {
        result[d] = tensor_i3(_ld[0].J,_ld[0].n,
                              _ld[1].J,_ld[1].n,
                              _ld[2].J,_ld[2].n,
                              _xyz[d],   // this is the "field"
                              _odwork);
      }
      return result;
    }
      // replicate the functionality of hex_findpt
    bool SpectralHex::evaluate_reverse(CartVect const & xyz, CartVect &params, double tol, const CartVect &init) const
    {
      params = init;
      
        //find nearest point
      double x_star[3];
      xyz.get(x_star);

      double r[3] = {0, 0, 0 }; // initial guess for parametric coords
      unsigned c = opt_no_constraints_3;
      double dist = opt_findpt_3(&_data, (const double **)_xyz, x_star, r, &c);
        // if it did not converge, get out with throw...
      if (dist > 0.9e+30)
        return false;
        //c tells us if we landed inside the element or exactly on a face, edge, or node
        // also, dist shows the distance to the computed point.
        //copy parametric coords back
      params = r;

      return is_inside(params, tol);
    }
    Matrix3  SpectralHex::jacobian(const CartVect& params) const
    {
      double x_i[3];
      params.get(x_i);
        // set the positions of GL nodes, before evaluations
      _data.elx[0]=_xyz[0];
      _data.elx[1]=_xyz[1];
      _data.elx[2]=_xyz[2];
      opt_vol_set_intp_3(&_data,x_i);
      Matrix3 J(0.);
        // it is organized differently
      J(0,0) = _data.jac[0]; // dx/dr
      J(0,1) = _data.jac[1]; // dx/ds
      J(0,2) = _data.jac[2]; // dx/dt
      J(1,0) = _data.jac[3]; // dy/dr
      J(1,1) = _data.jac[4]; // dy/ds
      J(1,2) = _data.jac[5]; // dy/dt
      J(2,0) = _data.jac[6]; // dz/dr
      J(2,1) = _data.jac[7]; // dz/ds
      J(2,2) = _data.jac[8]; // dz/dt
      return J;
    }
    void SpectralHex::evaluate_vector(const CartVect& params, const double *field, int num_tuples, double *eval) const
    {
        //piece that we shouldn't want to cache
      int d;
      for(d=0; d<3; d++){
        lagrange_0(&_ld[d], params[d]);
      }

      *eval = tensor_i3(_ld[0].J,_ld[0].n,
                        _ld[1].J,_ld[1].n,
                        _ld[2].J,_ld[2].n,
                        field,
                        _odwork);
    }
    void SpectralHex::integrate_vector(const double *field_values, int num_tuples, double *integral) const
    {
        // set the position of GL points
        // set the positions of GL nodes, before evaluations
      _data.elx[0]=_xyz[0];
      _data.elx[1]=_xyz[1];
      _data.elx[2]=_xyz[2];
      double params[3];
        //triple loop; the most inner loop is in r direction, then s, then t
      for (int l = 0; l < num_tuples; l++) integral[l] = 0.0;
        //double volume = 0;
      int index=0; // used fr the inner loop
      for (int k=0; k<_n; k++ )
      {
        params[2]=_ld[2].z[k];
          //double wk= _ld[2].w[k];
        for (int j=0; j<_n; j++)
        {
          params[1]=_ld[1].z[j];
            //double wj= _ld[1].w[j];
          for (int i=0; i<_n; i++)
          {
            params[0]=_ld[0].z[i];
              //double wi= _ld[0].w[i];
            opt_vol_set_intp_3(&_data,params);
            double wk= _ld[2].J[k];
            double wj= _ld[1].J[j];
            double wi= _ld[0].J[i];
            Matrix3 J(0.);
              // it is organized differently
            J(0,0) = _data.jac[0]; // dx/dr
            J(0,1) = _data.jac[1]; // dx/ds
            J(0,2) = _data.jac[2]; // dx/dt
            J(1,0) = _data.jac[3]; // dy/dr
            J(1,1) = _data.jac[4]; // dy/ds
            J(1,2) = _data.jac[5]; // dy/dt
            J(2,0) = _data.jac[6]; // dz/dr
            J(2,1) = _data.jac[7]; // dz/ds
            J(2,2) = _data.jac[8]; // dz/dt
            double bm = wk*wj*wi* J.determinant();
            for (int l = 0; l < num_tuples; l++)
              integral[l]+= bm*field_values[num_tuples*index+l];
              //volume +=bm;
          }
        }
      }
        //std::cout << "volume: " << volume << "\n";
    }
      // this is the same as a linear hex, although we should not need it
    bool SpectralHex::is_inside(const CartVect & params, double tol) const
    {
        // just look at the box+tol here
      return ( params[0]>=-1.-tol) && (params[0]<=1.+tol) &&
          ( params[1]>=-1.-tol) && (params[1]<=1.+tol) &&
          ( params[2]>=-1.-tol) && (params[2]<=1.+tol);
    }

      // SpectralHex

      // filescope for static member data that is cached
    int SpectralQuad::_n;
    double *SpectralQuad::_z[2];
    lagrange_data SpectralQuad::_ld[2];
    opt_data_2 SpectralQuad::_data;
    double * SpectralQuad::_odwork;
    double * SpectralQuad::_glpoints;
    bool SpectralQuad::_init = false;

    SpectralQuad::SpectralQuad() : Map(0)
    {
    }
      // the preferred constructor takes pointers to GL blocked positions
    SpectralQuad::SpectralQuad(int order, double * x, double *y, double *z) : Map(0)
    {
      Init(order);
      _xyz[0]=x; _xyz[1]=y; _xyz[2]=z;
    }
    SpectralQuad::SpectralQuad(int order) : Map(4)
    {
      Init(order);
      _xyz[0]=_xyz[1]=_xyz[2]=NULL;
    }
    SpectralQuad::~SpectralQuad()
    {
      if (_init)
        freedata();
      _init=false;
    }
    void SpectralQuad::Init(int order)
    {
      if (_init && _n==order)
        return;
      if (_init && _n!=order)
      {
          // TODO: free data cached
        freedata();
      }
        // compute stuff that depends only on order
      _init = true;
      _n = order;
        //duplicates! n is the same in all directions !!!
      for(int d=0; d<2; d++){
        _z[d] = tmalloc(double, _n);
        lobatto_nodes(_z[d], _n);
        lagrange_setup(&_ld[d], _z[d], _n);
      }
      opt_alloc_2(&_data, _ld);

      unsigned int nf = _n*_n, ne = _n, nw = 2*_n*_n + 3*_n;
      _odwork = tmalloc(double, 6*nf + 9*ne + nw);
      _glpoints = tmalloc (double, 3*nf);
    }

    void SpectralQuad::freedata()
    {
      for(int d=0; d<2; d++){
        free(_z[d]);
        lagrange_free(&_ld[d]);
      }
      opt_free_2(&_data);
      free(_odwork);
      free(_glpoints);
    }

    void SpectralQuad::set_gl_points( double * x, double * y, double *z)
    {
      _xyz[0] = x;
      _xyz[1] = y;
      _xyz[2] = z;
    }
    CartVect SpectralQuad::evaluate( const CartVect& params ) const
    {
        //piece that we shouldn't want to cache
      int d=0;
      for(d=0; d<2; d++){
        lagrange_0(&_ld[d], params[d]);
      }
      CartVect result;
      for (d=0; d<3; d++)
      {
        result[d] = tensor_i2(_ld[0].J,_ld[0].n,
                              _ld[1].J,_ld[1].n,
                              _xyz[d],
                              _odwork);
      }
      return result;
    }
      // replicate the functionality of hex_findpt
    bool SpectralQuad::evaluate_reverse(CartVect const & xyz, CartVect &params, double tol, const CartVect &init) const
    {
      params = init;

        //find nearest point
      double x_star[3];
      xyz.get(x_star);

      double r[2] = {0, 0 }; // initial guess for parametric coords
      unsigned c = opt_no_constraints_3;
      double dist = opt_findpt_2(&_data, (const double **)_xyz, x_star, r, &c);
        // if it did not converge, get out with throw...
      if (dist > 0.9e+30)
        throw Map::EvaluationError();
        //c tells us if we landed inside the element or exactly on a face, edge, or node
        // also, dist shows the distance to the computed point.
        //copy parametric coords back
      params = r;

      return is_inside(params, tol);
    }


    Matrix3  SpectralQuad::jacobian(const CartVect& /*params*/) const
    {
        // not implemented
      Matrix3 J(0.);
      return J;
    }


    void SpectralQuad::evaluate_vector(const CartVect& params, const double *field, int num_tuples, double *eval) const
    {
        //piece that we shouldn't want to cache
      int d;
      for(d=0; d<2; d++){
        lagrange_0(&_ld[d], params[d]);
      }

      *eval = tensor_i2(_ld[0].J,_ld[0].n,
                        _ld[1].J,_ld[1].n,
                        field,
                        _odwork);
    }
    void SpectralQuad:: integrate_vector(const double */*field_values*/, int /*num_tuples*/, double */*integral*/) const
    {
      // not implemented
    }
      // this is the same as a linear hex, although we should not need it
    bool SpectralQuad::is_inside(const CartVect & params, double tol) const
    {
        // just look at the box+tol here
      return ( params[0]>=-1.-tol) && (params[0]<=1.+tol) &&
          ( params[1]>=-1.-tol) && (params[1]<=1.+tol) ;
    }
      // something we don't do for spectral hex; we do it here because
      //       we do not store the position of gl points in a tag yet
    void SpectralQuad::compute_gl_positions()
    {
        // will need to use shape functions on a simple linear quad to compute gl points
        // so we know the position of gl points in parametric space, we will just compute those
        // from the 3d vertex position (corner nodes of the quad), using simple mapping
      assert (this->vertex.size()==4);
      static double corner_params[4][2]={ { -1., -1.},
                                          {  1., -1.},
                                          {  1.,  1.},
                                          { -1.,  1.} };
        // we will use the cached lobatto nodes in parametric space _z[d] (the same in both directions)
      int indexGL=0;
      int n2= _n*_n;
      for (int i=0; i<_n; i++)
      {
        double csi=_z[0][i];
        for (int j=0; j<_n; j++)
        {
          double eta = _z[1][j]; // we could really use the same _z[0] array of lobatto nodes
          CartVect pos(0.0);
          for (int k = 0; k < 4; k++) {
            const double N_k = (1 + csi*corner_params[k][0])
                * (1 + eta*corner_params[k][1]);
            pos += N_k * vertex[k];
          }
          pos *= 0.25;// these are x, y, z of gl points; reorder them
          _glpoints[indexGL] = pos[0]; // x
          _glpoints[indexGL+n2] = pos[1]; // y
          _glpoints[indexGL+2*n2] = pos[2]; // z
          indexGL++;
        }
      }
        // now, we can set the _xyz pointers to internal memory allocated to these!
      _xyz[0] =  &(_glpoints[0]);
      _xyz[1] =  &(_glpoints[n2]);
      _xyz[2] =  &(_glpoints[2*n2]);
    }
    void SpectralQuad::get_gl_points( double *& x, double *& y, double *& z, int & size)
    {
      x=  (double *)_xyz[0] ;
      y = (double *)_xyz[1] ;
      z = (double *)_xyz[2] ;
      size = _n*_n;
    }
    
}// namespace Element

} // namespace moab
