#include "GeometryQueryTool.hpp"
#include "InitCGMA.hpp"
#include "CGMApp.hpp"
#include "MBCore.hpp"
#include "cubfile.h"
#include "Tqdcfr.hpp"
#include "FileOptions.hpp"
#include "ReadNCDF.hpp"
#include "quads_to_tris.hpp"

#define GF_CUBIT_FILE_TYPE    "CUBIT"
#define GF_STEP_FILE_TYPE     "STEP"
#define GF_IGES_FILE_TYPE     "IGES"
#define GF_ACIS_TXT_FILE_TYPE "ACIS_SAT"
#define GF_ACIS_BIN_FILE_TYPE "ACIS_SAB"
#define GF_OCC_BREP_FILE_TYPE "OCC"

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
int is_occ_brep_file( FILE* file );

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
  const char* update_name = 0;
  const char* output_name = 0;
  const char* time_step = 0;
  double dist_tol = 0.001, len_tol = 0.0;
  int norm_tol = 5;
  bool actuate_attribs = true;
  
    // Process CL args
  bool process_options = true;
  if(argc < 5)
  {
    std::cerr << "Need a cub file, an update exodus file, an export h5m file and a time step." <<std::endl;
    exit(4);
  }

  for (int i = 1; i < argc; ++i) {
    if (!process_options || argv[i][0] != '-') {
      if(input_name && update_name && output_name)
        time_step = argv[i];
      else if (input_name && update_name)
        output_name = argv[i];
      else if(input_name) 
        update_name = argv[i];
      else
        input_name = argv[i];
      continue;
    }
  }   

    // Initialize CGM
  InitCGMA::initialize_cgma();
  if (actuate_attribs) {
    CGMApp::instance()->attrib_manager()->set_all_auto_read_flags( actuate_attribs );
    CGMApp::instance()->attrib_manager()->set_all_auto_actuate_flags( actuate_attribs );
  }
  
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
    
    MBCore *my_impl = new MBCore();
    Tqdcfr *my_tqd = new Tqdcfr(my_impl);
    ReadNCDF my_ex_reader(my_impl);
    MBEntityHandle file_set;
    char options[120] = "tdata=coord,";
    strcat(options, time_step);
    strcat(options,",set");
    FileOptions opts(options)  ;

    MBErrorCode result = my_tqd->load_file(input_name, file_set, opts, NULL, 0, 0);

    //opts = "tdata=coord, 100, sum, temp.exo";
    result =  my_ex_reader.load_file(update_name, file_set, opts, NULL, 0 , 0);

    // convert the quads to tris
    quads_to_tris( my_impl, file_set );

    result = my_impl->write_mesh( output_name );
    assert(!result);
  }
  return 0;
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
  static const char*  BREP_NAME = GF_OCC_BREP_FILE_TYPE;
  
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
  else if (is_occ_brep_file(file))
    return BREP_NAME;
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

int is_occ_brep_file( FILE* file )
{
  unsigned char buffer[6];
  return !fseek(file, 0, SEEK_SET) &&
         fread(buffer, 6, 1, file) &&
         !memcmp(buffer, "DBRep_", 6);
}
