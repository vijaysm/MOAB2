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
#include "moab/VerdictWrapper.hpp"
#include "../RefineMesh/moab/NestedRefine.hpp"

#ifdef USE_MPI
  #include "moab_mpi.h"
  #include "moab/ParallelComm.hpp"
#endif

/* Exit values */
#define SUCCESS 0
#define USAGE_ERROR 1
#define NOT_IMPLEMENTED 2

using namespace moab;

double wtime() {
  double y = -1;
  struct timeval cur_time;
  gettimeofday(&cur_time, NULL);
  y = (double)(cur_time.tv_sec) + (double)(cur_time.tv_usec)*1.e-6;
  return (y);
}

static void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage: " << name
         << " <options> <input_file> [<input_file2> ...]" << std::endl
         << "Options: " << std::endl
         << "\t-h             - Print this help text and exit." << std::endl
         << "\t-d             - Dimension of the mesh."<<std::endl
         << "\t-n             - Exact or a maximum number of levels for the hierarchy. Default 1." << std::endl
         << "\t-L             - Degree of refinement for each level. Pass an array or a number. It is mandatory to pass dimension and num_levels before to use this option. If this flag is not used, a default of deg 2 refinement is used. " << std::endl
         << "\t-V             - Pass a quality constraint of maximum volume for the hierarchy. If a volume constraint is passed, then a corresponding maximum number of levels should be passed. A optimal degree sequence would be computed satisfying the constraints to generate the hierarchy. It is mandatory to pass dimension for this option. " << std::endl
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

ErrorCode get_degree_seq(Core &moab, EntityHandle fileset, int dim, double desired_vol, int &num_levels, std::vector<int> &level_degs);

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
  bool qc_vol = false; double cvol = 0 ;
  std::string infile;

 // for (int i = 1; i < argc; i++
  int i=1;
  while ( i<argc)
  {
    if (!argv[i][0])
      usage_error(argv[0]);

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


  ErrorCode error;
  Core moab;
  Interface *mb = &moab;
  EntityHandle fileset;

  //Create a fileset
  error = mb->create_meshset(MESHSET_SET, fileset);MB_CHK_ERR(error);

  //Set the read options for parallel file loading
  std::vector<std::string> read_opts, write_opts;
  std::string read_options, write_options;

  if (parallel){
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

  if (num_levels == 0)
    num_levels = 1;

  //If a volume constraint is given, find an optimal degree sequence to reach the desired volume constraint.
  if (qc_vol){
      exit(NOT_IMPLEMENTED);
      if (dim==0)
        print_usage(argv[0], std::cerr);
     // error = get_degree_seq(&moab, fileset, dim, cvol, &num_levels, level_degrees);MB_CHK_ERR(error);
    }

  int *ldeg =  new int[num_levels];

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
      std::cout<<"Mesh size for level 0"<<"  :: nverts = "<<ents[0].size()<<", nedges = "<<ents[1].size()<<", nfaces = "<<ents[2].size()<<", ncells = "<<ents[3].size()<<std::endl;

      for (int l=0; l<num_levels; l++)
        {
          all_ents.clear();
          ents[0].clear(); ents[1].clear(); ents[2].clear(); ents[3].clear();
          error = mb->get_entities_by_handle(lsets[l+1], all_ents); MB_CHK_ERR(error);

          for (int k=0; k<4; k++)
            ents[k] = all_ents.subset_by_dimension(k);

          std::cout<<std::endl;
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

ErrorCode get_degree_seq(Core &moab,EntityHandle fileset, int dim, double desired_vol, int &num_levels, std::vector<int> &level_degs)
{
  ErrorCode error;
  VerdictWrapper vw(&moab);

  Range entities;
  error = moab.get_entities_by_handle(fileset, entities);MB_CHK_ERR(error);

  Range edges, faces, cells;
  edges = entities.subset_by_dimension(1);
  faces = entities.subset_by_dimension(2);
  cells = entities.subset_by_dimension(3);

  if (!cells.empty())
    {

    }

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

