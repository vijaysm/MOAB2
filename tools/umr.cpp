#include <iostream>
#include <stdlib.h>
#include <vector>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include <float.h>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "../verdict/moab/VerdictWrapper.hpp"
#include "../RefineMesh/moab/NestedRefine.hpp"

#ifdef USE_MPI
#include "moab_mpi.h"
#include "moab/ParallelComm.hpp"
#include "MBParallelConventions.h"
#endif

/* Exit values */
#define SUCCESS 0
#define USAGE_ERROR 1
#define NOT_IMPLEMENTED 2

using namespace moab;

static void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage: " << name
         << " <options> <input_file> [<input_file2> ...]" << std::endl
         << "Options: " << std::endl
         << "\t-h             - Print this help text and exit." << std::endl
         << "\t-d             - Dimension of the mesh."<<std::endl
         << "\t-n             - Exact or a maximum number of levels for the hierarchy. Default 1." << std::endl
         << "\t-L             - Degree of refinement for each level. Pass an array or a number. It is mandatory to pass dimension and num_levels before to use this option. If this flag is not used, a default of deg 2 refinement is used. " << std::endl
         << "\t-V             - Pass a desired volume (absolute) . This will generate a hierarchy such that the maximum volume is reduced to the given amount approximately. The length of the hierarchy can be constrained further if a maximum number of levels is passed. It is mandatory to pass the dimension for this option. " << std::endl
         <<"\t-q             - Prints out the maximum volume of the mesh and exits the program. This option can be used as a guide to volume constrainted mesh hierarchy later. "<<std::endl
         << "\t-t             - Print out the time taken to generate hierarchy." << std::endl
         <<"\t-s             - Print out the mesh sizes of each level of the generated hierarchy." << std::endl
         << "\t-o             - Specify true for output files for the mesh levels of the hierarchy." << std::endl
         //<< "\t-O option      - Specify read option." << std::endl
#ifdef USE_MPI
         << "\t-p[0|1|2]      - Read in parallel[0], optionally also doing resolve_shared_ents (1) and exchange_ghosts (2)" << std::endl
#endif
    ;
  exit(USAGE_ERROR);
}

static void usage_error( const char* name )
{
  print_usage( name, std::cerr );
#ifdef USE_MPI
  MPI_Finalize();
#endif
  exit(USAGE_ERROR);
}

bool parse_id_list(const char* string, int dim, int nval, std::vector<int> &results );

bool make_opts_string( std::vector<std::string> options, std::string& opts );

ErrorCode get_degree_seq(Core &mb, EntityHandle fileset, int dim, double desired_vol, int &num_levels, std::vector<int> &level_degs);

ErrorCode get_max_volume(Core &mb, EntityHandle fileset, int dim, double &vmax);

int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  int num_levels = 0, dim=0;
  std::vector<int> level_degrees;
  bool do_flag = true;
  bool print_times = false, print_size=false, output = false;
  bool parallel = false, resolve_shared = false, exchange_ghosts = false;
  bool printusage = false, parselevels = true;
  bool qc_vol = false, only_quality = false; double cvol = 0 ;
  std::string infile;

  int i=1;
  while ( i<argc)
  {
    if (!argv[i][0] && 0==proc_id)
      {
        usage_error(argv[0]);
#ifdef USE_MPI
        MPI_Finalize();
#endif
      }

    if (do_flag && argv[i][0] == '-')
    {
      switch ( argv[i][1] )
      {
          // do flag arguments:
        case '-': do_flag = false;       break;
        case 't': print_times = true;    break;
        case 's': print_size = true; break;
        case 'h':
        case 'H': print_usage( argv[0], std::cerr ); printusage = true; break;
        case 'd': dim = atol(argv[i+1]); ++i; break;
        case 'n': num_levels = atol(argv[i+1]);   ++i; break;
        case 'L':  if (dim !=0 && num_levels !=0)
            { parselevels = parse_id_list(argv[i+1], dim, num_levels, level_degrees); ++i;}
          else
            {print_usage( argv[0], std::cerr ); printusage = true;}
          break;
        case 'V': qc_vol = true; cvol = strtod(argv[i+1], NULL); ++i; break;
        case 'q': only_quality = true; break;
        case 'o': output = true; break;
#ifdef USE_MPI
        case 'p':
            parallel = true;
            if (argv[i][2] == '1' || argv[i][2] == '2') resolve_shared = true;
            if (argv[i][2] == '2') exchange_ghosts = true;
            break;
#endif
        default:
          ++i;
          switch ( argv[i-1][1] )
            {
              //case 'O':  read_opts.push_back(argv[i]); break;
            default: std::cerr << "Invalid option: " << argv[i] << std::endl;
            }
        }
      ++i;
    }
    // do file names
    else {
      infile = argv[i] ;
      ++i;
    }
  }

  if (infile.empty() && !printusage)
    print_usage(argv[0], std::cerr);

  if (!parselevels)
    exit(USAGE_ERROR);

#ifdef USE_MPI
  parallel = true;
#endif

  ErrorCode error;
  Core moab;
  Interface *mb = &moab;
  EntityHandle fileset;

  //Create a fileset
  error = mb->create_meshset(MESHSET_SET, fileset);MB_CHK_ERR(error);

  //Set the read options for parallel file loading
  std::vector<std::string> read_opts, write_opts;
  std::string read_options, write_options;

  if (parallel && size > 1){
      read_opts.push_back("PARALLEL=READ_PART");
      read_opts.push_back("PARTITION=PARALLEL_PARTITION");
      if (resolve_shared) read_opts.push_back("PARALLEL_RESOLVE_SHARED_ENTS");
      if (exchange_ghosts) read_opts.push_back("PARALLEL_GHOSTS=3.0.1");

      if (output)
        write_opts.push_back("PARALLEL=WRITE_PART");
    }

  if (!make_opts_string(  read_opts,  read_options ) || !make_opts_string(  write_opts,  write_options ))
    {
#ifdef USE_MPI
      MPI_Finalize();
#endif
      return USAGE_ERROR;
    }

  //Load file
  error = mb->load_file(infile.c_str(), &fileset, read_options.c_str());MB_CHK_ERR(error);

  //Create the nestedrefine instance

#ifdef USE_MPI
  ParallelComm *pc = new ParallelComm(&moab, MPI_COMM_WORLD);
  NestedRefine uref(&moab,pc,fileset);
#else
  NestedRefine uref(&moab);
#endif

  std::vector<EntityHandle> lsets;

 // std::cout<<"Read input file"<<std::endl;

  if (only_quality)
    {
      double vmax;
      error = get_max_volume(moab, fileset, dim, vmax);MB_CHK_ERR(error);
#ifdef USE_MPI
      int rank = 0;
      MPI_Comm_rank(MPI_COMM_WORLD, &rank);
      if (rank==0)
        std::cout<<"Maximum global volume = "<<vmax<<std::endl;
#else
      std::cout<<"Maximum volume = "<<vmax<<std::endl;
#endif
      exit(SUCCESS);
    }

  //If a volume constraint is given, find an optimal degree sequence to reach the desired volume constraint.
  if (qc_vol){
      error = get_degree_seq(moab, fileset, dim, cvol, num_levels, level_degrees);MB_CHK_ERR(error);

      if (dim==0)
        print_usage(argv[0], std::cerr);
    }

  if (num_levels == 0)
    num_levels = 1;

  int *ldeg =  new int[num_levels];
  //std::cout<<"level_degrees.size = "<<level_degrees.size()<<std::endl;
  if (level_degrees.empty()){
      for (int l=0; l<num_levels; l++)
        ldeg[l] = 2;
    }
  else
    {
      for (int l=0; l<num_levels; l++)
        ldeg[l] = level_degrees[l];
    }

  std::cout<<"Starting hierarchy generation"<<std::endl;
  error = uref.generate_mesh_hierarchy( num_levels, ldeg, lsets);MB_CHK_ERR(error);

  if (print_times)
    {
      std::cout<<"Finished hierarchy generation in "<<uref.timeall.tm_total<<"  secs"<<std::endl;
      if (parallel)
        {
          std::cout<<"Time taken for refinement "<<uref.timeall.tm_refine<<"  secs"<<std::endl;
          std::cout<<"Time taken for resolving shared interface "<<uref.timeall.tm_presolve<<"  secs"<<std::endl;
        }
    }
  else
    std::cout<<"Finished hierarchy generation "<<std::endl;

  if (print_size)
    {
      Range all_ents, ents[4];
      error = mb->get_entities_by_handle(fileset, all_ents); MB_CHK_ERR(error);

      for (int k=0; k<4; k++)
        ents[k] = all_ents.subset_by_dimension(k);

      std::cout<<std::endl;
      if (qc_vol)
        {
          double volume;
          error = get_max_volume(moab, fileset, dim, volume);MB_CHK_ERR(error);
          std::cout<<"Mesh size for level 0"<<"  :: nverts = "<<ents[0].size()<<", nedges = "<<ents[1].size()<<", nfaces = "<<ents[2].size()<<", ncells = "<<ents[3].size()<<" :: Vmax = "<<volume<<std::endl;
        }
      else
        std::cout<<"Mesh size for level 0"<<"  :: nverts = "<<ents[0].size()<<", nedges = "<<ents[1].size()<<", nfaces = "<<ents[2].size()<<", ncells = "<<ents[3].size()<<std::endl;

      for (int l=0; l<num_levels; l++)
        {
          all_ents.clear();
          ents[0].clear(); ents[1].clear(); ents[2].clear(); ents[3].clear();
          error = mb->get_entities_by_handle(lsets[l+1], all_ents); MB_CHK_ERR(error);

          for (int k=0; k<4; k++)
            ents[k] = all_ents.subset_by_dimension(k);

         // std::cout<<std::endl;

          if (qc_vol)
            {
              double volume;
              error = get_max_volume(moab, lsets[l+1], dim, volume);MB_CHK_ERR(error);
              std::cout<<"Mesh size for level "<<l+1<<"  :: nverts = "<<ents[0].size()<<", nedges = "<<ents[1].size()<<", nfaces = "<<ents[2].size()<<", ncells = "<<ents[3].size()<<" :: Vmax = "<<volume<<std::endl;
            }
          else
          std::cout<<"Mesh size for level "<<l+1<<"  :: nverts = "<<ents[0].size()<<", nedges = "<<ents[1].size()<<", nfaces = "<<ents[2].size()<<", ncells = "<<ents[3].size()<<std::endl;
        }
    }


  if (output){
      for (int l=0; l<num_levels; l++){
          std::string::size_type idx1 = infile.find_last_of("\\/");
          std::string::size_type idx2 = infile.find_last_of(".");
          std::string file = infile.substr(idx1+1,idx2-idx1-1);
          std::stringstream out;
          if (parallel)
            out <<  "_ML_" <<l+1<<".h5m";
          else
            out <<  "_ML_" <<l+1<<".vtk";
          file = file + out.str();
          const char* output_file = file.c_str();
          mb->write_file(output_file, 0, write_options.c_str(), &lsets[l+1], 1);
        }
    }

  delete [] ldeg;

#ifdef USE_MPI
  MPI_Finalize();
#endif

  exit(SUCCESS);
}

ErrorCode get_degree_seq(Core &mb, EntityHandle fileset, int dim, double desired_vol, int &num_levels, std::vector<int> &level_degs)
{
  //Find max volume
  double vmax_global;
  ErrorCode error = get_max_volume(mb, fileset, dim, vmax_global);MB_CHK_ERR(error);

  int init_nl = num_levels;
  num_levels = 0; level_degs.clear();

  //Find a sequence that reduces the maximum volume to desired.
  //desired_vol should be a fraction or absolute value ?

  //double remV = vmax_global*desired_vol/dim;
  double Vdesired = desired_vol;
  double try_x;
  double remV = vmax_global;
  int degs[3][3] = {{5,3,2},{25,9,4},{0,27,8}};

  if (dim==1 || dim == 2)
    {
      while (remV -Vdesired>= 0){
          try_x = degs[dim-1][0];
          if ((remV/try_x - Vdesired)>=0)
            {
              level_degs.push_back(5);
              num_levels += 1;
              remV /= try_x;
            }
          else
            {
              try_x = degs[dim-1][1];
              if ((remV/try_x - Vdesired)>=0)
                {
                  level_degs.push_back(3);
                  num_levels += 1;
                  remV /= try_x;
                }
              else {
                  try_x = degs[dim-1][2];
                  if ((remV/try_x - Vdesired)>=0)
                    {
                      level_degs.push_back(2);
                      num_levels += 1;
                      remV /= try_x;
                    }
                  else
                    break;
                }
            }
        }
    }
  else
    {
      while (remV -Vdesired>= 0){
          try_x = degs[dim-1][1];
          if ((remV/try_x - Vdesired)>=0)
            {
              level_degs.push_back(3);
              num_levels += 1;
              remV /= try_x;
            }
          else
            {
              try_x = degs[dim-1][2];
              if ((remV/try_x - Vdesired)>=0)
                {
                  level_degs.push_back(2);
                  num_levels += 1;
                  remV /= try_x;
                }
              else
                break;
            }
        }
    }

  if (init_nl != 0 && init_nl < num_levels)
    {
      for (int i=level_degs.size(); i>= init_nl; i--)
        level_degs.pop_back();
      num_levels = init_nl;
    }

  return MB_SUCCESS;
}

ErrorCode get_max_volume(Core &mb,  EntityHandle fileset, int dim, double &vmax)
{
  ErrorCode error;
  VerdictWrapper vw(&mb);
  QualityType q;

  switch (dim) {
    case 1: q = MB_LENGTH; break;
    case 2: q = MB_AREA; break;
    case 3: q = MB_VOLUME; break;
    default: return MB_FAILURE; break;
    }

  //Get all entities of the highest dimension which is passed as a command line argument.
  Range allents, owned;
  error = mb.get_entities_by_handle(fileset, allents);MB_CHK_ERR(error);
  owned = allents.subset_by_dimension(dim);MB_CHK_ERR(error);

  //Get all owned entities
#ifdef USE_MPI
  int size = 1;
  MPI_Comm_size( MPI_COMM_WORLD, &size );
  int mpi_err;
  if (size>1)
    {
      // filter the entities not owned, so we do not process them more than once
      ParallelComm* pcomm = moab::ParallelComm::get_pcomm(&mb, 0);
      Range current = owned;
      owned.clear();
      error = pcomm->filter_pstatus(current, PSTATUS_NOT_OWNED, PSTATUS_NOT, -1, &owned);
      if (error != MB_SUCCESS)
        {
          MPI_Finalize();
          return MB_FAILURE;
        }
    }
#endif

  double vmax_local=0;
  //Find the maximum volume of an entity in the owned mesh
  for (Range::iterator it=owned.begin(); it != owned.end(); it++)
    {
      double volume;
      error = vw.quality_measure(*it, q, volume);MB_CHK_ERR(error);
      if (volume >vmax_local)
        vmax_local = volume;
    }

  //Get the global maximum
  double vmax_global = vmax_local;
#ifdef USE_MPI
  mpi_err = MPI_Reduce(&vmax_local, &vmax_global, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  if (mpi_err)
      {
        MPI_Finalize();
        return MB_FAILURE;
      }
#endif

  vmax = vmax_global;

  return MB_SUCCESS;
}


bool parse_id_list( const char* string, int dim, int nval, std::vector<int> &results )
{
  bool okay = true;
  char* mystr = strdup( string );
  for (const char* ptr = strtok(mystr, ","); ptr; ptr = strtok(0,","))
  {
    char* endptr;
    int val = strtol( ptr, &endptr, 0 );

    if (dim==1 || dim == 2){
        if(val != 2 && val != 3 && val != 5){
            std::cerr <<"Not a valid degree for the passed dimension" << std::endl;
            okay = false;
            break;
          }
      }
    else if (dim==3){
        if(val != 2 && val != 3){
            std::cerr <<"Not a valid degree for the passed dimension" << std::endl;
            okay = false;
            break;
          }
      }

    if (endptr == ptr || val <= 0) {
      std::cerr << "Not a valid id: " << ptr << std::endl;
      okay = false;
      break;
    }

    results.push_back(val);
  }

  if ((int)results.size() < nval)
    {
      for (int i=results.size(); i<=nval-1; i++)
        results.push_back(results[0]);
    }
  else if ((int)results.size() > nval)
    {
      for (int i=results.size(); i>nval; i--)
        results.pop_back();
    }

  free( mystr );
  return okay;
}

bool make_opts_string( std::vector<std::string> options, std::string& opts )
{
  opts.clear();
  if (options.empty())
    return true;

    // choose a separator character
  std::vector<std::string>::const_iterator i;
  char separator = '\0';
  const char* alt_separators = ";+,:\t\n";
  for (const char* sep_ptr = alt_separators; *sep_ptr; ++sep_ptr) {
    bool seen = false;
    for (i = options.begin(); i != options.end(); ++i)
      if (i->find( *sep_ptr, 0 ) != std::string::npos) {
        seen = true;
        break;
      }
    if (!seen) {
      separator = *sep_ptr;
      break;
    }
  }
  if (!separator) {
    std::cerr << "Error: cannot find separator character for options string" << std::endl;
    return false;
  }
  if (separator != ';') {
    opts = ";";
    opts += separator;
  }

    // concatenate options
  i = options.begin();
  opts += *i;
  for (++i; i != options.end(); ++i) {
    opts += separator;
    opts += *i;
  }

  return true;
}

