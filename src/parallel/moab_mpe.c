#include "moab_mpe.h"

int MPE_Allocate_event()
{
  static int counter = MOAB_MPE_FIRST_EVENT;
  return counter++;
}
