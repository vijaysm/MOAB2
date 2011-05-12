#ifndef IREL_LASSO_HPP
#define IREL_LASSO_HPP

#include "iRel.h"

#define RETURN(CODE) ERROR((CODE), "")
#define RETURNR(CODE) ERRORR((CODE), "")

#define ERROR(CODE, MSG)                           \
  do {                                             \
    *ierr = LASSOI->set_last_error((CODE), (MSG)); \
    return;                                        \
  } while(false)

#define ERRORR(CODE, MSG) return LASSOI->set_last_error((CODE), (MSG))

#define CHK_ERROR(CODE)                            \
  do {                                             \
    *ierr = (CODE);                                \
    if ((CODE) != iBase_SUCCESS)                   \
      return;                                      \
  } while(false)

#define CHK_ERRORR(CODE)                           \
  do {                                             \
    if ((CODE) != iBase_SUCCESS)                   \
      return (CODE);                               \
  } while(false)

#endif
