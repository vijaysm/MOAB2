#ifndef SOLVERS_HPP
#define SOLVERS_HPP
#include <vector>

namespace moab
{

  class Solvers
  {
    Solvers() {};
    ~Solvers() {};
public:
    //compute combinational number, n choose k, maximum output is std::numeric_limits<unsigned int>::max();
    //If overflows, return 0
    static unsigned int nchoosek(unsigned int n, unsigned int k);

    static unsigned int compute_numcols_vander_multivar(unsigned int kvars,unsigned int degree);

    static void gen_multivar_monomial_basis(const int kvars,const double* vars, const int degree, std::vector<double>& basis);

    static void gen_vander_multivar(const int npts,const int kvars, const double* us, const int degree, std::vector<double>& V);

    static void rescale_matrix(int mrows, int ncols, double *V, double *ts);

    static void compute_qtransposeB(int mrows, int ncols, const double *Q, int bncols, double *bs);

    static void qr_polyfit_safeguarded(int mrows, int ncols, double *V, double *D, int *rank);

    static void backsolve(int mrows, int ncols, double *R, int bncols, double *bs, double *ws);

    static void backsolve_polyfit_safeguarded();

    static void vec_dotprod(const int len, const double* a, const double* b, double* c);

    static void vec_scalarprod(const int len, const double* a, const double c, double* b);

    static void vec_crossprod(const double a[3], const double b[3], double (&c)[3]);

    static double vec_innerprod(const int len, const double* a, const double* b);

    static double vec_2norm(const int len, const double* a);

    static double vec_normalize(const int len, const double* a, double* b);

    static void vec_projoff(const int len, const double* a, const double* b, double* c);

    static void vec_linear_operation(const int len, const double mu, const double* a, const double psi, const double* b, double* c);

  };

}
#endif
