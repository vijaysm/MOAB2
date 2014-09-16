#ifdef WIN32
#pragma warning (disable : 4786)
#endif 

#include "RefTemplates.hpp"

namespace moab{

  const RefTemplates::refPatterns RefTemplates::refTemplates[2][MAX_DEGREE]=
  {
    //TRI
    {{1,0,0,3,4,{3,5},{{0.5,0},{0.5,0.5},{0,0.5}},{{0,3,5},{3,4,5},{3,1,4},{5,4,2}},{{0,0,2,2,0,0},{3,2,4,0,1,1},{0,0,0,0,2,0},{2,1,0,0,0,0}}},//deg2
     {},//deg3
     {} //deg5

    },

    //QUAD
    {}
  };


}//namesapce moab
  
