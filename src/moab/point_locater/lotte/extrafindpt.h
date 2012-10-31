#ifndef EXTRA_FINDPT_H
#define EXTRA_FINDPT_H
/* extra stuff that james didn't expose */

typedef struct {
  unsigned constraints;
  unsigned dn, d1, d2;
  real *x[3], *fdn[3];
} opt_face_data_3;

typedef struct {
  unsigned constraints;
  unsigned de, d1, d2;
  real *x[3], *fd1[3], *fd2[3];
} opt_edge_data_3;

typedef struct {
  unsigned constraints;
  real x[3], jac[9];
} opt_point_data_3;

typedef struct {
  lagrange_data *ld;
  unsigned size[4];
  const real *elx[3];
  opt_face_data_3 fd;
  opt_edge_data_3 ed;
  opt_point_data_3 pd;
  real *work;
  real x[3], jac[9];
} opt_data_3;


void opt_alloc_3(opt_data_3 *p, lagrange_data *ld);
void opt_free_3(opt_data_3 *p);
double opt_findpt_3(opt_data_3 *p, const real *const elx[3],
                           const real xstar[3], real r[3], unsigned *constr);
void opt_vol_set_intp_3(opt_data_3 *p, const real r[3]);

const unsigned opt_no_constraints_2 = 3+1;
const unsigned opt_no_constraints_3 = 9+3+1;

/* for 2d spectralQuad */
/*--------------------------------------------------------------------------

   2 - D

  --------------------------------------------------------------------------*/

typedef struct {
  unsigned constraints;
  unsigned de, d1;
  real *x[2], *fd1[2];
} opt_edge_data_2;

typedef struct {
  unsigned constraints;
  real x[2], jac[4];
} opt_point_data_2;

typedef struct {
  lagrange_data *ld;
  unsigned size[3];
  const real *elx[2];
  opt_edge_data_2 ed;
  opt_point_data_2 pd;
  real *work;
  real x[2], jac[4];
} opt_data_2;
void opt_alloc_2(opt_data_2 *p, lagrange_data *ld);
void opt_free_2(opt_data_2 *p);
double opt_findpt_2(opt_data_2 *p, const real *const elx[2],
                           const real xstar[2], real r[2], unsigned *constr);



#endif //EXTRA_FINDPT_H
