#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/VerdictWrapper.hpp"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "TestUtil.hpp"


std::string TestDir( STRINGIFY(MESHDIR) );

std::string filename = TestDir + "/mbtest1.g";

using namespace moab;

int main( int argc, char* argv[] )
{
  ErrorCode rval;
  Core moab_core;
  Interface* mb = &moab_core;
  if(argc > 1) if (argc > 1) filename = std::string(argv[1]);
  rval = mb->load_mesh( filename.c_str());
  if (MB_SUCCESS != rval) {
    std::cerr << "Error reading file: " << filename.c_str() << std::endl;
    return 1;
  }
  //std::cout << "loaded mesh file: " << filename.c_str() << std::endl;
  Range entities;
  rval = mb->get_entities_by_handle( 0, entities ); // all entities from the model
  if (MB_SUCCESS != rval) {
    std::cerr << "can't get entities from file " << filename.c_str() << std::endl;
    return 1;
  }
  VerdictWrapper vw(mb);
  // for size methods/quality, we need a size, to compute relative sizes and stuff
  rval = vw.set_size(1.0);
  for (Range::iterator eit=entities.begin(); eit!=entities.end(); eit++)
  {
    EntityHandle eh=*eit;
    EntityType etype=TYPE_FROM_HANDLE(eh);
    if (etype==MBVERTEX || etype>MBHEX)
      continue;
    for (int quality=0; quality<MB_QUALITY_COUNT; quality++)
    {
      QualityType q = (QualityType)quality;
      double qm;
      rval = vw.quality_measure(eh, q, qm);
      if (MB_NOT_IMPLEMENTED == rval)
        continue;
      if (MB_FAILURE == rval)
      {
        std::cerr << " failure for entity " << mb->list_entity(eh) << " quality " << vw.quality_name(q) << "\n";
        return 1;
      }
      if (MB_SUCCESS == rval)
      {
        std::cout << "Entity type " << (EntityType)mb->type_from_handle(eh) << " id:" << mb->id_from_handle(eh) << " quality:" <<
            vw.quality_name(q) << " : " << qm << "\n";
      }
    }
    // now compute all qualities for each entity handle
    std::map<QualityType, double> qualities;

    rval = vw.all_quality_measures(eh, qualities);
    if (MB_SUCCESS == rval)
    {
      mb->list_entity(eh);
      for (std::map<QualityType, double>::iterator mit=qualities.begin(); mit!=qualities.end(); mit++)
      {
        std::cout << "   " << vw.quality_name(mit->first) << " " << mit->second << " \n";
      }
    }

  }

  return 0;
}

