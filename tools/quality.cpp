#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <stdio.h>
#include <iomanip>
#ifdef USE_MPI
  #include "moab_mpi.h"
  #include "moab/ParallelComm.hpp"
#endif
#if !defined(_MSC_VER) && !defined(__MINGW32__)
  #include <termios.h>
  #include <sys/ioctl.h>
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


int main( int argc, char* argv[] )
{
  int proc_id = 0, size = 1;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif
  if (argc<2 && 0==proc_id)
  {
    print_usage(argv[0], std::cerr);
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }
  Core mb;
  std::string read_options;
  if (size>1)
    read_options="PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS";

  if (MB_SUCCESS != mb.load_file(argv[1], 0, read_options.c_str() ) )
  {
    fprintf(stderr, "Error reading file: %s\n", argv[1] );
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  VerdictWrapper vw(&mb);
  vw.set_size(1.); // for relative size measures; maybe it should be user input
  Range entities;
  ErrorCode rval = mb.get_entities_by_handle(0, entities);
  if (MB_SUCCESS!=rval )
  {
    fprintf(stderr, "Error getting entities from file %s\n", argv[1] );
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }
  // get all edges and faces, force them to be created
  Range allfaces, alledges;
  Range cells = entities.subset_by_dimension(3);
  rval = mb.get_adjacencies(cells, 2, true, allfaces, Interface::UNION);
  if (MB_SUCCESS!=rval )
  {
    fprintf(stderr, "Error getting all faces" );
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }

  rval = mb.get_adjacencies(allfaces, 1, true, alledges, Interface::UNION);
  if (MB_SUCCESS!=rval )
  {
    fprintf(stderr, "Error getting all edges" );
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return 1;
  }
  entities.merge(allfaces);
  entities.merge(alledges);
  for (EntityType et=MBENTITYSET; et >= MBEDGE; et--)
  {
    int num_qualities = vw.num_qualities(et);
    if (!num_qualities)
      continue;
    Range owned=entities.subset_by_type(et);
    std::map<QualityType, double> qualities, minq, maxq;
    int ne_local = (int)owned.size();
    int ne_global = ne_local;

#ifdef USE_MPI
    int mpi_err;
    if (size>1)
    {
    // filter the entities not owned, so we do not process them more than once
      ParallelComm* pcomm = ParallelComm::get_pcomm(&mb, 0);
      Range current = owned;
      owned.clear();
      rval = pcomm->filter_pstatus(current, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned);
      if (rval != MB_SUCCESS)
      {
        MPI_Finalize();
        return 1;
      }
      ne_local = (int)owned.size();
      mpi_err = MPI_Reduce(&ne_local, &ne_global, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      if (mpi_err)
      {
        MPI_Finalize();
        return 1;
      }
    }
#endif
    if (ne_global>0)
    {
      if (ne_local>0)
      {
        Range::iterator  it = owned.begin();
        rval = vw.all_quality_measures(*it, qualities);
        if (MB_SUCCESS!=rval )
        {
          fprintf(stderr, "Error getting quality for entity type %d with id %ld \n", et, mb.id_from_handle(*it) );
  #ifdef USE_MPI
          MPI_Finalize();
  #endif
          return 1;
        }
        minq=qualities; maxq=qualities;
        it++;
        for ( ;it!=owned.end(); it++)
        {
          rval = vw.all_quality_measures(*it, qualities);
          if (MB_SUCCESS!=rval )
          {
            fprintf(stderr, "Error getting quality for entity type %d with id %ld \n", et, mb.id_from_handle(*it) );
  #ifdef USE_MPI
            MPI_Finalize();
  #endif
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
      }
      if (0==proc_id)
      {
        std::cout <<" \n\n   " <<  ne_global << " entities of type " <<  vw.entity_type_name(et) << "\n";
        std::cout <<std::setw(30) << "Quality Name" << std::setw(15) << "    MIN" << std::setw(15) << "  MAX" << "\n";
      }

      /*std::map<QualityType, double>::iterator minit;
      std::map<QualityType, double>::iterator maxit;
      if (ne_local > 0) {
        minit = minq.begin();
        maxit = maxq.begin();
      }*/
      QualityType quality_type=MB_EDGE_RATIO;
      for (int i=0; i<num_qualities; i++, quality_type=(QualityType)(quality_type+1))
      {
        while(! (vw.possible_quality(et, quality_type)) && quality_type<MB_QUALITY_COUNT)
          quality_type=(QualityType)(quality_type+1); // will get them in order
        const char * name_q = vw.quality_name(quality_type);
        double local_min, global_min;
        double local_max, global_max;
        if (ne_local>0)
        {
          local_min = global_min = minq[quality_type];
          local_max = global_max = maxq[quality_type];
        }
        else
        {
          local_min = global_min = 1.e38; // so this task has no entities of this type
          local_max = global_max = -1.e38;// it can get here only in parallel
        }
#ifdef USE_MPI
        mpi_err = MPI_Reduce(&local_min, &global_min, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
        if (mpi_err)
        {
          MPI_Finalize();
          return 1;
        }
        mpi_err = MPI_Reduce(&local_max, &global_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        if (mpi_err)
        {
          MPI_Finalize();
          return 1;
        }
#endif
        if (0==proc_id)
        {
          std::cout <<std::setw(30) << name_q
              << std::setw(15) << global_min <<
                 std::setw(15) << global_max << "\n";

        }

      }
    }
  }//etype
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}


