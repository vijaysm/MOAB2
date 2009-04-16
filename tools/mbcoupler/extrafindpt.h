
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


const unsigned opt_no_constraints_2 = 3+1;
const unsigned opt_no_constraints_3 = 9+3+1;
