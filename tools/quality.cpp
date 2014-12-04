#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <stdio.h>
#include <iomanip>
#ifndef WIN32
#  include <sys/times.h>
#  include <limits.h>
#  include <unistd.h>
#endif
#include <time.h>
#ifdef USE_MPI
#  include "moab_mpi.h"
#endif
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#  include <termios.h>
#  include <sys/ioctl.h>
#endif
#include <math.h>
#include <assert.h>
#include <float.h>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/Interface.hpp"


#include "moab/VerdictWrapper.hpp"

using namespace moab;


static void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage: " << name << " <input_file> " << std::endl;
}

Core mb;


int main( int argc, char* argv[] )
{
  if (argc<2)
  {
    print_usage(argv[0], std::cerr);
    return 1;
  }
  if (MB_SUCCESS != mb.load_file(argv[1]) )
  {
    fprintf(stderr, "Error reading file: %s\n", argv[1] );
    return 1;
  }

  VerdictWrapper vw(&mb);
  vw.set_size(1.); // for relative size measures
  Range entities;
  ErrorCode rval = mb.get_entities_by_handle(0, entities);
  if (MB_SUCCESS!=rval )
  {
    fprintf(stderr, "Error getting entities from file %s\n", argv[1] );
    return 1;
  }
  for (EntityType et=MBEDGE; et<MBENTITYSET; et++)
  {
    Range typed=entities.subset_by_type(et);
    if (!typed.empty())
    {
      std::map<QualityType, double> qualities, minq, maxq;
      Range::iterator  it = typed.begin();
      rval = vw.all_quality_measures(*it, qualities);
      if (MB_SUCCESS!=rval )
      {
        fprintf(stderr, "Error getting quality for entity type %d with id %ld \n", et, mb.id_from_handle(*it) );
        return 1;
      }
      minq=qualities; maxq=qualities;
      it++;
      for ( ;it!=typed.end(); it++)
      {
        rval = vw.all_quality_measures(*it, qualities);
        if (MB_SUCCESS!=rval )
        {
          fprintf(stderr, "Error getting quality for entity type %d with id %ld \n", et, mb.id_from_handle(*it) );
          return 1;
        }
        std::map<QualityType, double>::iterator minit=minq.begin();
        std::map<QualityType, double>::iterator maxit=maxq.begin();
        for (std::map<QualityType, double>::iterator mit=qualities.begin();
            mit!=qualities.end();
            mit++, minit++, maxit++)
        {
          if (mit->second > maxit->second) maxit->second = mit->second;
          if (mit->second < minit->second) minit->second = mit->second;
        }
      }
      std::cout <<" \n\n   " <<  typed.size() << " entities of type " << et << "\n";
      std::cout <<std::setw(30) << "Quality Name" << std::setw(15) << "    MIN" << std::setw(15) << "  MAX" << "\n";
      std::map<QualityType, double>::iterator minit=minq.begin();
      for (std::map<QualityType, double>::iterator maxit=maxq.begin();
                  maxit!=maxq.end();
                  minit++, maxit++)
      {
        std::cout <<std::setw(30) << vw.quality_name(minit->first)
            << std::setw(15) << minit->second <<
               std::setw(15) << maxit->second << "\n";
      }


    }
  }
  return 0;
}


