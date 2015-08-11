#include "moab/Solvers.hpp"
#include <iostream>
#include <assert.h>
#include <vector>
#include <limits>
#include <cmath>

namespace moab {

  /* This class implements the lowest level solvers required by polynomial fitting for high-order reconstruction.
   * An underlying assumption of the matrices passed are that they are given in a column-major form. So,
   * the first mrows values of V is the first column, and so on. This assumption is made because most of the
   *  operations are column-based in the current scenario.
   * */

  void Solvers::rescale_matrix(int mrows, int ncols, double *V, double *ts)
  {
    //This function rescales the input matrix using the norm of each column.
    double *v = new double[mrows];
    for (int i=0; i< ncols; i++)
      {
        for (int j=0; j<mrows; j++)
          v[j] = V[mrows*i+j];

        //Compute norm of the column vector
        double w = vec_2norm(v, mrows);

        if (abs(w)==0)
          ts[i] = 1;
        else
          {
            ts[i] = w;
            for (int j=0; j<mrows; j++)
              V[mrows*i+j] = V[mrows*i+j]/ts[i];
          }
      }
    delete [] v;
  }

  void Solvers::compute_qtransposeB(int mrows, int ncols, const double *Q, int bncols, double *bs)
  {
    for (int k=0; k<ncols; k++)
      {
        for (int j=0; j<bncols; j++)
          {
            double t2 = 0;
            for (int i=k; i<nrows; i++)
              t2 += Q[mrows*k+i]*bs[mrows*j+i];
            t2 = t2 + t2;

            for (int i=k; i<nrows; i++)
              bs[mrows*j+i] -= t2*Q[mrows*k+i];
          }
      }
  }

  void Solvers::qr_polyfit_safeguarded(double *V, int mrows, int ncols, double *D, int &rank)
  {
    double tol = 1e-8;
    rank = ncols;
    double *v = new double[mrows];

    for (int k=0; k<ncols; k++)
      {
        int nv = nrows-k+1;

        for (int j=0; j<nv; j++)
          v[j] = A[mrows*k + (j+k-1)];

        double t2=0;

        for (int j=0; j<nv; j++)
          t2 = t2 + v[j]*v[j];

        double t = sqrt(t2);
        double vnrm = 0;

        if (v[0] >=0)
          {
            vnrm = sqrt(2*(t2+v[0]*t));
            v[0] = v[0]+t;
          }
        else
          {
            vnrm = sqrt(2*(t2-v[0]*t));
            v[0] = v[0]-t;
          }

        if (vnrm>0)
          {
            for (int j=0; j<nv; j++)
              v[j] = v[j]/vnrm;
          }

        for(int j=k; j<ncols; j++)
          {
            t2 = 0;
            for (int i=0; i<nv; i++)
              t2 = t2 + v[i]*A[mrows*j+(i+k-1)];
            t2 = t2+t2;
            for (int i=0; i<nv; i++)
                A[mrows*j+(i+k-1)] = A[mrows*j+(i+k-1)] - t2*v[i];
          }

        D[k] = A[mrows*k+k];

        for (int i=0; i<nv; i++)
            A[mrows*k+(i+k-1)] = v[i];

        if ((abs(D[k])) < tol && (rank == ncols))
          {
            rank = k-1;
            break;
          }
      }

    delete [] v;
  }

  void Solvers::backsolve(int mrows, int ncols, double *R, int bncols, double *bs, double *ws)
  {
    for (int k=0; k< bncols; k++)
      {
        for (int j=ncols-1; j>=0; j--)
          {
            for (int i=j+1; j<ncols; j++)
              bs[mrows*k+j] = bs[mrows*k+j] - R[mrows*i+j]*bs[mrows*k+i];

            assert(R[mrows*j+j] != 0);

            bs[mrows*k+j] = bs[mrows*k+j]/R[mrows*j+j];

            for (int j=1; j<ncols; j++)
              bs[mrows*k+j] = bs[mrows*k+j]/ws[j];

          }
      }
  }

  void Solvers::backsolve_polyfit_safeguarded()
  {

  }

  void Solvers::vec_dotprod(const int len, const double* a, const double* b, double* c)
  {
    for(int i=0;i<len;++i){
        c[i] = a[i]*b[i];
      }
  }

  void Solvers::vec_scalarprod(const int len, const double* a, const double c, double *b)
  {
    for(int i=0;i<len;++i){
         b[i] = c*a[i];
      }
  }

  void Solvers::vec_crossprod(const double a[3], const double b[3], double (&c)[3])
  {
    c[0] =a[1]*b[2]-a[2]*b[1];
    c[1] =a[2]*b[0]-a[0]*b[2];
    c[2] =a[0]*b[1]-a[1]*b[0];
  }

  double Solvers::vec_innerprod(const int len, const double* a, const double* b)
  {
    double ans=0;
    for(int i=0;i<len;++i){
        ans += a[i]*b[i];
    }
    return ans;
  }

  double Solvers::vec_2norm(const int len, const double* a)
  {
    double w=0, s=0;
    for (int k=0; k<len; k++)
      w = std::max(w, abs(a[k]));

    if (w==0){
        return 0;
    }else{
        for (int k=0; k<len; k++){
          s += (a[k]/w)*(a[k]/w);
        }
        s=w*sqrt(s);
    }
    return s;
  }

  double Solvers::vec_normalize(const int len, const double* a, double* b)
  {
    double nrm=0,mx=0;
    for(int i=0;i<len;++i){
        mx = std::max(abs(a[i]),mx);
    }
    if(0==mx){
      for(int i=0;i<len){
        b[i] = 0;
      }
      return 0;
    }
    for(int i=0;i<len;++i){
        nrm += (a[i]/mx)*(a[i]/mx);
    }
    nrm = mx*sqrt(nrm);
    if(nrm==0){
        return nrm;
    }
    for(int i=0;i<len;++i){
        b[i] = a[i]/nrm;
    }
    return nrm;
  }

  void Solvers::vec_projoff(const int len, const double* a, const double* b, double* c)
  {
    //c = a-<a,b>b/<b,b>;
    double bnrm = vec_2norm(b,len);
    if (bnrm==0){
        for(int i=0;i<len;++i){
          c[i] = a[i];
        }
        return; 
    }
    double innerp = vec_innerprod(len,a,b)/bnrm;

    if(innerp==0){
      for(int i=0;i<len;++i){
        c[i] = a[i];
      }
      return;
    }

    for(int i=0;i<len;++i){
        c[i] = a[i]-innerp*b[i];
    }
  }

    void Solvers::vec_linear_operation(const int len, const double mu, const double* a, const double psi, const double* b, double* c)
    {
      for(int i=0;i<len;++i){
          c[i] = mu*a[i]+psi*b[i];
      }
    }

}