#ifndef TYPES_H
#define TYPES_H
#include "float.h"

/* integer type to use for everything */
#if   defined(USE_LONG)
#  define INTEGER long
#elif defined(USE_LONG_LONG)
#  define INTEGER long long
#elif defined(USE_SHORT)
#  define INTEGER short
#else
#  define INTEGER int
#endif

/* when defined, use the given type for global indices instead of INTEGER */
#if   defined(USE_GLOBAL_LONG_LONG)
#  define GLOBAL_INT long long
#elif defined(USE_GLOBAL_LONG)
#  define GLOBAL_INT long
#else
#  define GLOBAL_INT long
#endif

/* floating point type to use for everything */
#if   defined(USE_FLOAT)
   typedef float real;
#  define floorr floorf
#  define ceilr  ceilf
#  define sqrtr  sqrtf
#  define fabsr  fabsf
#  define cosr   cosf
#  define sinr   sinf
#  define EPS   (128*FLT_EPSILON)
#  define PI 3.1415926535897932384626433832795028841971693993751058209749445923F
#elif defined(USE_LONG_DOUBLE)
   typedef long double real;
#  define floorr floorl
#  define ceilr  ceill
#  define sqrtr  sqrtl
#  define fabsr  fabsl
#  define cosr   cosl
#  define sinr   sinl
#  define EPS   (128*LDBL_EPSILON)
#  define PI 3.1415926535897932384626433832795028841971693993751058209749445923L
#else
   typedef double real;
#  define floorr floor
#  define ceilr  ceil
#  define sqrtr  sqrt
#  define fabsr  fabs
#  define cosr   cos
#  define sinr   sin
#  define EPS   (128*DBL_EPSILON)
#  define PI 3.1415926535897932384626433832795028841971693993751058209749445923
#endif

/* apparently uint and ulong can be defined already in standard headers */
#define uint uint_
#define ulong ulong_
#define sint sint_
#define slong slong_

typedef   signed INTEGER sint;
typedef unsigned INTEGER uint;
#undef INTEGER

#ifdef GLOBAL_INT
  typedef   signed GLOBAL_INT slong;
  typedef unsigned GLOBAL_INT ulong;
#else
  typedef sint slong;
  typedef uint ulong;
#endif

#endif

