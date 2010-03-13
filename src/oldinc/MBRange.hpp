#ifndef MBRange_HEADER
#define MBRange_HEADER

#include "MBTypes.h"

#include "moab/Range.hpp"

typedef moab::Range MBRange;
typedef moab::range_inserter mb_range_inserter;
using moab::intersect; //( const moab::Range&, const moab::Range& );
using moab::subtract;  //( const moab::Range&, const moab::Range& );
using moab::unite;     //( const moab::Range&, const moab::Range& );

#endif
