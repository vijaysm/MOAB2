#ifndef SOLVERS_HPP
#define SOLVERS_HPP

namespace moab
{

  class Solvers
  {
    Solvers() {};
    ~Solvers() {};

    void rescale_matrix(int mrows, int ncols, double *V, double *ts);

    void compute_qtransposeB(int mrows, int ncols, const double *Q, int bncols, double *bs);

    void qr_polyfit_safeguarded(int mrows, int ncols, double *V, double *D, int &rank);

    void backsolve(int mrows, int ncols, double *R, int bncols, double *bs, double *ws);

    void backsolve_polyfit_safeguarded(int dim, int degree, bool interp, int mrows, int ncols, double *R, int bncols, double *bs, double *ws, double *degree_out);

    void vec_dotprod(const double* a, const double* b, const int len, double* c);

    void vec_scalarprod(const double* a, const int len, const double c, double *b);

    void vec_crossprod(const double a[3], const double b[3], double (&c)[3]);

    double vec_innerprod(const double* a, const double* b, const int len);

    double vec_2norm(const double* a, const int len);

    double vec_normalize(const double* a, const int len, double* b);

    void vec_projoff(const double* a, const double* b, const int len, double* c);

      void vec_linear_operation(const double mu, const double* a, const double psi, const double* b, const int len, double* c);



  };

}
#endif
