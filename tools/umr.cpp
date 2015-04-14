#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>

#include <math.h>
#include <assert.h>
#include <float.h>
#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/VerdictWrapper.hpp"
#include "moab/NestedRefine.hpp"
#ifdef USE_MPI
  #include "moab_mpi.h"
  #include "moab/ParallelComm.hpp"
#endif

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
         << "\t-n             - Exact or a maximum number of levels for the hierarchy." << std::endl
         << "\t-d             - Degree of refinement for each level. Pass an array to use different degrees for each level." << std::endl
         << "\t-qv           - Pass a quality constraint of maximum volume for the hierarchy. If a volume constraint is passed, then a corresponding maximum number of levels should be passed. A optimal degree sequence would be computed satisfying the constraints to generate the hierarchy." << std::endl
         << "\t-T             - Print out the time taken to generate hierarchy." << std::endl
         << "\t-o    - Specify true for output files for the mesh levels of the hierarchy." << std::endl
         //<< "\t-O option      - Specify read option." << std::endl
#ifdef USE_MPI
         << "\t-p[0|1|2]      - Read in parallel[0], optionally also doing resolve_shared_ents (1) and exchange_ghosts (2)" << std::endl
#endif
    ;
}

static void usage_error( const char* name )
{
  print_usage( name, std::cerr );
#ifdef USE_MPI
  MPI_Finalize();
#endif
  exit(USAGE_ERROR);
}


int main(int argc, char* argv[])
{
  int proc_id = 0, size = 1;
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &proc_id );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#endif

  int num_levels = 0;
  int *level_degrees;
  bool do_flag = true;
  bool print_times = false, output = false;
  bool parallel = false, resolve_shared = false, exhange_ghosts = false;
  bool print_usage = false;
  bool qc_vol = false; double cvol = 0 ;
  std::vector<std::string> file_list;

  for (i = 1; i < argc; i++)
  {
    if (!argv[i][0])
      usage_error(argv[0]);

    if (do_flag && argv[i][0] == '-')
    {
      switch ( argv[i][1] )
      {
          // do flag arguments:
        case '-': do_flag = false;       break;
        case 'T': print_times = true;    break;
        case 'h':
        case 'H': print_usage( argv[0], std::cerr ); printed_usage = true; break;
        case 'n': num_level = argv[i];   break;
        case 'd': parse_id_list(argv[i], num_level, level_degrees); break;
        case 'qv': qc_vol = true; cvol = argv[i]; break;
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
    }
    // do file names
    else {
      file_list.push_back( argv[i] );
    }
  }

  ErrorCode error;
  Core moab;
  Interface *mb = &moab;
  EntityHandle fileset;

  //Create a fileset
  error = mb->create_meshset(MESHSET_SET, fileset);MB_CHK_ERR(error);
  //Set the read options for parallel file loading
  std::string read_opts = NULL, write_opts = NULL;
  if (parallel){
      read_opts.push_back("PARALLEL=READ_PART");
      read_opts.push_back("PARTITION=PARALLEL_PARTITION");
      if (resolve_shared) read_opts.push_back("PARALLEL_RESOLVE_SHARED_ENTS");
      if (exchange_ghosts) read_opts.push_back("PARALLEL_GHOSTS=3.0.1");

      if (output)
        write_opts.push_back("PARALLEL=WRITE_PART");
    }
  //Load file
  error = mb->load_file(file_list[0], &fileset, read_opts.c_str());MB_CHK_ERR(error);

  //Create the nestedrefine instance
if (parallel){
    ParallelComm *pc = new ParallelComm(&moab, MPI_COMM_WORLD);
    NestedRefine uref(&moab,pc,fileset);
  }
else
  NestedRefine uref(&moab,);

//If a volume constraint is given, find an optimal degree sequence to reach the desired volume constraint.
Range entities;
error = mb->get_entities_by_handle(fileset, entities);MB_CHK_ERR(error);
if (qc_vol){
    int nl_new = 0;
    int *ldeg;
    get_degree_seq(&moab, fileset, &nl_new, ldeg);MB_CHK_ERR(error);


  }

std::vector<EntityHandle> lsets;
if (print_times)
  std::cout<<"Starting hierarchy generation"<<std::endl;

double time_start = wtime();
error = uref.generate_mesh_hierarchy( num_levels,level_degrees, set); MB_CHK_ERR(error);
double time_total = wtime() - time_start;

if (print_times)
  std::cout<<"Finished hierarchy generation in "<<time_total<<"  secs"<<std::endl;

if (output){
    for (int i=0; i<num_levels; i++){
        const char *output_file = tmpstr.c_str();
        mb->write_file(output_file, 0, write_opts, &lsets[l+1], 1);
      }
  }


}

void get_degree_seq(Core &moab,EntityHandle fileset, int &nl_new,int *ldeg)
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

}

void parse_id_list( const char* string, int nval, int *results )
{
  bool okay = true;
  char* mystr = strdup( string );
  for (const char* ptr = strtok(mystr, ","); ptr; ptr = strtok(0,","))
  {
    char* endptr;
    int val = strtol( ptr, &endptr, 0 );
    if (endptr == ptr || val <= 0) {
      std::cerr << "Not a valid id: " << ptr << std::endl;
      okay = false;
      break;
    }

    int val2 = val;
    if (*endptr == '-') {
      const char* sptr = endptr+1;
      val2 = strtol( sptr, &endptr, 0 );
      if (endptr == sptr || val2 <= 0) {
        std::cerr << "Not a valid id: " << sptr << std::endl;
        okay = false;
        break;
      }
      if (val2 < val) {
        std::cerr << "Invalid id range: " << ptr << std::endl;
        okay = false;
        break;
      }
    }

    if (*endptr) {
      std::cerr << "Unexpected character: " << *endptr << std::endl;
      okay = false;
      break;
    }

    for (; val <= val2; ++val)
      if (!results.insert( (int)val ).second)
        std::cerr << "Warning: duplicate Id: " << val << std::endl;

  }

  free( mystr );
  return okay;
}

