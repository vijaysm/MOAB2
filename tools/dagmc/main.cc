#include "GeometryQueryTool.hpp"
#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "MBCore.hpp"
#include "cgm2moab.hpp"
#include "cubfile.h"

#define GF_CUBIT_FILE_TYPE    "CUBIT"
#define GF_STEP_FILE_TYPE     "STEP"
#define GF_IGES_FILE_TYPE     "IGES"
#define GF_ACIS_TXT_FILE_TYPE "ACIS_SAT"
#define GF_ACIS_BIN_FILE_TYPE "ACIS_SAB"

/* Get the type of a file.
   Return value is one of the above constants
 */
const char* get_geom_file_type( const char* filename );
const char* get_geom_fptr_type( FILE* file );

int is_cubit_file( FILE* file );
int is_step_file( FILE* file );
int is_iges_file( FILE* file );
int is_acis_txt_file( FILE* file );
int is_acis_bin_file( FILE* file );

double parse_tol( char* argv[], int argc, int& i );
void usage(char *name);
void print_usage(char *name);
void help(char *name);

double DEFAULT_DISTANCE = 0.001;
double DEFAULT_LEN = 0.0;
int DEFAULT_NORM = 5;

int main( int argc, char* argv[] )
{
  char *name = argv[0];
  const char *file_type = NULL;
  
  const char* input_name = 0;
  const char* output_name = 0;

  double dist_tol = 0.001, len_tol = 0.0;
  int norm_tol = 5;
  bool actuate_attribs = true;
  
    // Process CL args
  bool process_options = true;
  for (int i = 1; i < argc; ++i) {
    if (!process_options || argv[i][0] != '-') {
      if (output_name) {
        std::cerr << "Unexpected argument: " << argv[i] << std::endl;
        usage(argv[0]);
      }
      else if (input_name)
        output_name = argv[i];
      else 
        input_name = argv[i];
      continue;
    }
    
    if (!argv[i][1] || argv[i][2]) {  // two chars long
      std::cerr << "Invalid option: " << argv[i] << std::endl;
      usage(argv[0]);
    }
    
    switch(argv[i][1]) {
      default:
        std::cerr << "Invalid option: " << argv[i] << std::endl;
        usage(argv[0]);
      case 'd':
        dist_tol = parse_tol( argv, argc, i );
        break;
      case 'D':
        dist_tol = 0.0;
        break;
      case 'n':
        norm_tol = (int)round(parse_tol( argv, argc, i ));
        break;
      case 'N':
        norm_tol = 0;
        break;
      case 'l':
        len_tol = parse_tol( argv, argc, i );
        break;
      case 'L':
        len_tol = 0;
        break;
      case 'a':
        actuate_attribs = true;
        break;
      case 'A':
        actuate_attribs = false;
        break;
      case 't':
        ++i;
        if (i == argc) {
          std::cerr << "Expected argument following '-t'" << std::endl;
          usage(argv[0]);
        }
        file_type = argv[i];
        break;
      case 'h':
        help(argv[0]);
      case '-':
        process_options = false;
        break;
    }
  }
  
  if (!output_name)
    usage(argv[0]);
        
  
    // Initialize CGM
  InitCGMA::initialize_cgma("ACIS");
  if (actuate_attribs) {
    CGMApp::instance()->attrib_manager()->set_all_auto_read_flags( actuate_attribs );
    CGMApp::instance()->attrib_manager()->set_all_auto_actuate_flags( actuate_attribs );
  }
  
  
    // Intitalize MOAB
  MBCore moab;
  MBInterface* iface = &moab;
  
    // Get CGM file type
  if (!file_type) {
    file_type = get_geom_file_type( input_name );
    if (!file_type) {
      std::cerr << input_name << " : unknown file type, try '-t'" << std::endl;
      exit(1);
    }
  }
  
    // If CUB file, extract ACIS geometry
  CubitStatus s;
  if (!strcmp( file_type, "CUBIT" )) {
    char* temp_name = tempnam( 0, "cgm" );
    if (!temp_name) {
      perror(name);
      exit(2);
    }
    
    FILE* cub_file = fopen( input_name, "r" );
    if (!cub_file) {
      free(temp_name);
      perror(input_name);
      exit(2);
    }
    
    FILE* tmp_file = fopen( temp_name, "w" );
    if (!tmp_file) {
      free(temp_name);
      perror(temp_name);
      exit(2);
    }
    
    int rval = cub_file_type( cub_file, tmp_file, CUB_FILE_ACIS );
    fclose( cub_file );
    fclose( tmp_file );
    if (rval) {
      remove( temp_name );
      free( temp_name );
      exit(2);
    }
    
    s = GeometryQueryTool::instance()->import_solid_model( temp_name, "ACIS_SAT" );
    remove( temp_name );
    free( temp_name );
  }
  else {
    s = GeometryQueryTool::instance()->import_solid_model( input_name, file_type );
  }
  if (CUBIT_SUCCESS != s) {
    std::cerr << "Failed to read '" << input_name << "' of type '" << file_type << "'" << std::endl;
    exit(2);
  }
  
    // copy geometry facets into mesh database
  if (!cgm2moab(iface, dist_tol, norm_tol, 
                len_tol, actuate_attribs)) {
    std::cerr << "Internal error copying geometry" << std::endl;
    exit(5);
  }
  
    // write mesh database
  MBErrorCode rval = iface->write_mesh( output_name );
  if (MB_SUCCESS != rval) {
    std::cerr << "Failed to write '" << output_name << "'" << std::endl;
    exit(2);
  }
  
  return 0;
}

// parse double CL arg
double parse_tol( char* argv[], int argc, int& i )
{
  ++i;
  if (i == argc) {
    std::cerr << "Expected value following option '" << argv[i-1] << "'" << std::endl;
    usage(argv[0]);
  }
  
  char* endptr = 0;
  double val = strtod( argv[i], &endptr );
  if (!endptr || *endptr || val < DBL_EPSILON) {
    std::cerr << "Invalid tolerance value for '" << argv[i-1] <<"': " << argv[i] << std::endl;
    usage(argv[0]);
  }
  
  return val;
}

void print_usage(char *name)
{
  std::cout << "Usage: " << std::endl;
  std::cout << name << 
      " [-d <tol>|-D] [-n <tol>|-N] [-l <tol>|-L] [-A|-a] [-t <type>] <input_file> <outupt_file>" <<
      std::endl;
  std::cout << name << " -h" << std::endl;
}

void usage(char *name)
{
  print_usage(name);
  exit(1);
}

void help(char *name)
{
  print_usage(name);
  std::cout << "\t-d  max. distance deviation between facet and geometry" << std::endl
            << "\t-D  no distance tolerance" << std::endl
            << "\t    (default:" << DEFAULT_DISTANCE << ")" <<std::endl
            << "\t-n  max. normal angle deviation (degrees) between facet and geometry" << std::endl
            << "\t    (default:" << DEFAULT_NORM << ")" <<std::endl
            << "\t-N  no normal tolerance" << std::endl
            << "\t-l  max. facet edge length" << std::endl
            << "\t-L  no facet edge length maximum" << std::endl
            << "\t    (default:" << DEFAULT_LEN << ")" << std::endl
            << "\t-a  force actuation of all CGM attributes" << std::endl
            << "\t-A  disable all CGM attributes" << std::endl
            << "\t-t  specify input file type (default is autodetect)" << std::endl
            << std::endl;
  exit(0);
}

const char* get_geom_file_type( const char* name )
{
  FILE* file;
  const char* result = 0;
  
  file = fopen( name, "r" );
  if (file) {
    result = get_geom_fptr_type( file );
    fclose( file );
  }
  
  return result;
}

const char* get_geom_fptr_type( FILE* file )
{
  static const char* CUBIT_NAME = GF_CUBIT_FILE_TYPE;
  static const char*  STEP_NAME = GF_STEP_FILE_TYPE;
  static const char*  IGES_NAME = GF_IGES_FILE_TYPE;
  static const char*   SAT_NAME = GF_ACIS_TXT_FILE_TYPE;
  static const char*   SAB_NAME = GF_ACIS_BIN_FILE_TYPE;
  
  if (is_cubit_file(file))
    return CUBIT_NAME;
  else if (is_step_file(file))
    return STEP_NAME;
  else if (is_iges_file(file))
    return IGES_NAME;
  else if (is_acis_bin_file(file))
    return SAB_NAME;
  else if (is_acis_txt_file(file))
    return SAT_NAME;
  else
    return 0;
}

int is_cubit_file( FILE* file )
{
  unsigned char buffer[4];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 4, 1, file) &&
         !memcmp(buffer, "CUBE", 4);
}

int is_step_file( FILE* file )
{
  unsigned char buffer[9];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 9, 1, file) &&
         !memcmp(buffer, "ISO-10303", 9);
}

int is_iges_file( FILE* file )
{
  unsigned char buffer[10];
  return !fseek(file, 72, SEEK_SET) &&
         fread(buffer, 10, 1, file) &&
         !memcmp(buffer, "S      1\r\n", 10);
}

int is_acis_bin_file( FILE* file )
{
  char buffer[15];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 15, 1, file) &&
         !memcmp(buffer, "ACIS BinaryFile", 9);
}

int is_acis_txt_file( FILE* file )
{
  char buffer[5];
  int version, length;
  
  if (fseek(file,0,SEEK_SET) || 
      2 != fscanf( file, "%d %*d %*d %*d %d ", &version, &length ))
    return 0;
    
  if (version < 1 || version >0xFFFF)
    return 0;
  
    // Skip appliation name
  if (fseek(file, length, SEEK_CUR))
    return 0;
    
    // Read length of version string followed by first 5 characters
  if (2 != fscanf(file, "%d %4s", &length, buffer))
    return 0;
    
  return !strcmp( buffer, "ACIS" );
}
