#ifndef SOLVERS_HPP
#define SOLVERS_HPP

namespace moab
{

  class Solvers
  {
    Solvers() {};
    ~Solvers() {};

    void rescale_matrix(int mrows, int ncols, double *V, double *ts);

    void gen_vander_bivar(const int npts,const double* us, const int degree, double* V);

    void compute_qtransposeB(int mrows, int ncols, const double *Q, int bncols, double *bs);

    void qr_polyfit_safeguarded(double *V, int mrows, int ncols, double *D, int &rank);

    void backsolve(int mrows, int ncols, double *R, int bncols, double *bs, double *ws);

    void backsolve_polyfit_safeguarded();

    void vec_dotprod(const int len, const double* a, const double* b, double* c);

    void vec_scalarprod(const int len, const double* a, const double c, double *b);

    void vec_crossprod(const double a[3], const double b[3], double (&c)[3]);

    double vec_innerprod(const int len, const double* a, const double* b);

    double vec_2norm(const int len, const double* a);

    double vec_normalize(const int len, const double* a, double* b);

    void vec_projoff(const int len, const double* a, const double* b, double* c);

    void vec_linear_operation(const int len, const double mu, const double* a, const double psi, const double* b, double* c);



  };

}
#endif
