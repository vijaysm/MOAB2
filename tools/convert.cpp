/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

// If Microsoft compiler, then WIN32
#if defined(_MSC_VER) || defined(__MINGW32__)
#  ifndef WIN32
#    define WIN32
#  endif
#endif

#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "MBReaderWriterSet.hpp"
#include <iostream>
#include <iomanip>
#include <set>
#include <cstdlib>
#include <algorithm>
#ifndef WIN32
#  include <sys/times.h>
#  include <limits.h>
#  include <unistd.h>
#endif
#include <time.h>
#ifdef USE_MPI
#  include <mpi.h>
#endif

/* Exit values */
#define USAGE_ERROR 1
#define READ_ERROR 2
#define WRITE_ERROR 3
#define OTHER_ERROR 4
#define ENT_NOT_FOUND 5

/* Tag to control ACIS dump from .cub file reader */
const char acis_dump_file_tag_name[] = "__ACISDumpFile";


void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage: " << name << 
    " [-a <sat_file>|-A] [-t] [subset options] [-f format] <input_file> <output_file>" << std::endl
    << "\t-f <format>    - Specify output file format" << std::endl
    << "\t-a <acis_file> - ACIS SAT file dumped by .cub reader (same as \"-o SAT_FILE=acis_file\"" << std::endl
    << "\t-A             - .cub file reader should not dump a SAT file (depricated default)" << std::endl
    << "\t-o option      - Specify write option." << std::endl
    << "\t-O option      - Specify read option." << std::endl
    << "\t-t             - Time read and write of files." << std::endl
    << "\t-g             - Enable verbose/debug output." << std::endl
    << "\t-h             - Print this help text and exit." << std::endl
    << "\t-l             - List available file formats and exit." << std::endl
    << "\t-I <dim>       - Generate internal entities of specified dimension." << std::endl
#ifdef USE_MPI
    << "\t-P             - Append processor ID to file name" << std::endl
#endif
    << "\t--             - treat all subsequent options as file names" << std::endl
    << "\t                 (allows file names beginning with '-')" << std::endl
    << "  subset options: " << std::endl
    << "\tEach of the following options should be followed by " << std::endl
    << "\ta list of ids.  IDs may be separated with commas.  " << std::endl
    << "\tRanges of IDs may be specified with a '-' between " << std::endl
    << "\ttwo values.  The list may not contain spaces." << std::endl
    << "\t-v  - volume" << std::endl
    << "\t-s  - surface" << std::endl
    << "\t-c  - curve" << std::endl
    << "\t-V  - vertex" << std::endl
    << "\t-m  - material set (block)" << std::endl
    << "\t-d  - Dirchlet set (nodeset)" << std::endl
    << "\t-n  - Neumann set (sideset)" << std::endl
    << "\tThe presence of one or more of the following flags limits " << std::endl
    << "\tthe exported mesh to only elements of the corresponding " << std::endl
    << "\tdimension.  Vertices are always exported." << std::endl
    << "\t-1  - edges " << std::endl
    << "\t-2  - tri, quad, polygon " << std::endl
    << "\t-3  - tet, hex, prism, etc. " << std::endl
    ;
}

void print_help( const char* name )
{
  std::cout << 
  " This program can be used to convert between mesh file\n"
  " formats, extract a subset of a mesh file to a separate\n"
  " file, or both.  The type of file to write is determined\n"
  " from the file extension (e.g. \".vtk\") protion of the\n"
  " output file name.\n"
  " \n"
  " While MOAB strives to export and import all data from\n"
  " each supported file format, most file formats do\n"
  " not support MOAB's entire data model.  Thus MOAB cannot\n"
  " guarantee lossless conversion for any file formats\n"
  " other than the native HDF5 representation.\n"
  "\n";
  
  print_usage( name, std::cout );
  exit(0);
}

void usage_error( const char* name )
{
  print_usage( name, std::cerr );
  exit(USAGE_ERROR);
} 


void list_formats( MBInterface* );
bool parse_id_list( const char* string, std::set<int>&  );
void print_id_list( const char*, std::ostream& stream, const std::set<int>& list );
void reset_times();
void write_times( std::ostream& stream );
void remove_entities_from_sets( MBInterface* gMB, MBRange& dead_entities, MBRange& empty_sets );
void remove_from_vector( std::vector<MBEntityHandle>& vect, const MBRange& ents_to_remove );
bool make_opts_string( std::vector<std::string> options, std::string& result );

int main(int argc, char* argv[])
{
  MBInterface* gMB;
  MBErrorCode result;
  MBRange range;

  int proc_id = 0;
  gMB = new MBCore();

  bool append_rank = false;
  int i, dim;
  bool dims[4] = {false, false, false, false};
  const char* format = NULL; // output file format
  const char* in = NULL;    // input file name
  const char* out = NULL;   // output file name
  bool verbose = false;
  std::set<int> geom[4], mesh[3];       // user-specified IDs 
  std::vector<MBEntityHandle> set_list; // list of user-specified sets to write
  std::vector<std::string> write_opts, read_opts;
  const char* const mesh_tag_names[] = { DIRICHLET_SET_TAG_NAME,
                                         NEUMANN_SET_TAG_NAME,
                                         MATERIAL_SET_TAG_NAME };
  const char* const geom_names[] = { "VERTEX",
                                     "CURVE",
                                     "SURFACE",
                                     "VOLUME" };
  
    // scan arguments
  bool do_flag = true;
  bool print_times = false;
  bool generate[] = { false, false, false };
  bool pval;
  for (i = 1; i < argc; i++)
  {
    if (!argv[i][0])
      usage_error(argv[0]);
      
    if (do_flag && argv[i][0] == '-')
    {
      if (!argv[i][1] || argv[i][2])
        usage_error(argv[0]);
      
      switch ( argv[i][1] )
      {
          // do flag arguments:
        case '-': do_flag = false;       break;
        case 'g': verbose = true;        break;
        case 't': print_times = true;    break;
        case 'A':                        break;
        case 'h': 
        case 'H': print_help( argv[0] ); break;
        case 'l': list_formats( gMB );   break;
#ifdef USE_MPI
        case 'P': append_rank = true;    break;
#endif
        case '1': case '2': case '3':
          dims[argv[i][1] - '0'] = true; break;
          // do options that require additional args:
        default: 
          ++i;
          if (i == argc || argv[i][0] == '-') {
            std::cerr << "Expected argument following " << argv[i-1] << std::endl;
            usage_error(argv[0]);
          }
          if (argv[i-1][1] == 'I') {
            int dim = atoi( argv[i] );
            if (dim < 1 || dim > 2) {
              std::cerr << "Invalid dimension value following -I" << std::endl;
              usage_error(argv[0]);
            }
            generate[dim] = true;
            continue;
          }
          pval = false;
          switch ( argv[i-1][1] )
          {
            case 'a': 
              read_opts.push_back( std::string("SAT_FILE=") + argv[i] );
              pval = true;
              break;
            case 'f': format = argv[i]; pval = true;              break;
            case 'o': write_opts.push_back(argv[i]); pval = true; break;
            case 'O':  read_opts.push_back(argv[i]); pval = true; break;
            case 'v': pval = parse_id_list( argv[i], geom[3] ); break;
            case 's': pval = parse_id_list( argv[i], geom[2] ); break;
            case 'c': pval = parse_id_list( argv[i], geom[1] ); break;
            case 'V': pval = parse_id_list( argv[i], geom[0] ); break;
            case 'm': pval = parse_id_list( argv[i], mesh[2] ); break;
            case 'n': pval = parse_id_list( argv[i], mesh[1] ); break;
            case 'd': pval = parse_id_list( argv[i], mesh[0] ); break;
            default: std::cerr << "Invalid option: " << argv[i] << std::endl;
          }
          
          if (!pval) {
            std::cerr << "Invalid flag or flag value: " << argv[i-1] << " " << argv[i] << std::endl;
            usage_error(argv[0]);
          }
      }
    }
      // do file names
    else if (!in)
      in = argv[i];
    else if (!out)
      out = argv[i];
    else  { // too many file names
      std::cerr << "Unexpexed argument: " << argv[i] << std::endl;
      usage_error(argv[0]);
    }
  }
  if (!in || !out) {
    std::cerr << "No output file name specified." << std::endl;
    usage_error(argv[0]);
  }
    
  std::string mod_out;
  if (append_rank) {
    char buffer[16];
    sprintf(buffer,".%d",proc_id);
    mod_out = out;
    mod_out += buffer;
    out = mod_out.c_str();
  }

    // construct options string from individual options
  std::string read_options, write_options;
  if (!make_opts_string(  read_opts,  read_options ) ||
      !make_opts_string( write_opts, write_options ))
    return USAGE_ERROR;
  
  
    // Read the input file.
  reset_times();
  MBEntityHandle read_meshset;
  result = gMB->load_file( in, read_meshset, read_options.c_str() );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to load \"" << in << "\"." << std::endl;
    std::cerr  << "Error code: " << gMB->get_error_string(result) << " (" << result << ")" << std::endl;
    std::string message;
    if (MB_SUCCESS == gMB->get_last_error(message) && !message.empty())
      std::cerr << "Error message: " << message << std::endl;
    return READ_ERROR;
  }
  std::cerr << "Read \"" << in << "\"" << std::endl;
  if (print_times) write_times( std::cerr );
  
    // Determine if the user has specified any geometry sets to write
  bool have_geom = false;
  for (dim = 0; dim <= 3; ++dim)
  {
    if (!geom[dim].empty())
      have_geom = true;
    if (verbose)
      print_id_list( geom_names[dim], std::cout, geom[dim] );
  }
  
    // True if the user has specified any sets to write
  bool have_sets = have_geom;
  
    // Get geometry tags
  MBTag dim_tag, id_tag;
  if (have_geom) 
  {
    int size;
    result = gMB->tag_get_handle( GLOBAL_ID_TAG_NAME, id_tag );
    if (MB_SUCCESS != result) 
    {
      std::cerr << "No ID tag defined."  << std::endl;
      have_geom = false;
    }
    else
    {
      result = gMB->tag_get_size( id_tag, size );
      if (MB_SUCCESS != result || size != sizeof(int))
      {
        std::cerr << "ID tag is invalid." << std::endl;
        have_geom = false;
      }
    }
    result = gMB->tag_get_handle( GEOM_DIMENSION_TAG_NAME, dim_tag );
    if (MB_SUCCESS != result) 
    {
      std::cerr << "No geometry tag defined."  << std::endl;
      have_geom = false;
    }
    else
    {
      result = gMB->tag_get_size( dim_tag, size );
      if (MB_SUCCESS != result || size != sizeof(int))
      {
        std::cerr << "Geometry tag is invalid." << std::endl;
        have_geom = false;
      }
    }
  }
  
    // Get geometry sets
  if ( have_geom ) 
  {
    int id_val;
    MBTag tags[] = { id_tag, dim_tag };
    const void* vals[] = { &id_val, &dim };
    for (dim = 0; dim <= 3; ++dim) 
    {
      int init_count = set_list.size();
      for (std::set<int>::iterator iter = geom[dim].begin(); iter != geom[dim].end(); ++iter) 
      {
        id_val = *iter;
        range.clear();
        result = gMB->get_entities_by_type_and_tag( 0, MBENTITYSET, tags, vals, 2, range );
        if (MB_SUCCESS != result || range.empty()) 
        {
          range.clear();
          std::cerr << geom_names[dim] << " " << id_val << " not found.\n";
        }
        std::copy( range.begin(), range.end(), std::back_inserter(set_list) );
      }
      
      if (verbose)
        std::cout << "Found " << (set_list.size()-init_count) << ' '
                  << geom_names[dim] << " sets" << std::endl;
    }
  }
  
    // Get mesh groupings
  for (i = 0; i < 3; ++i) 
  {
    if (verbose)
      print_id_list( mesh_tag_names[i], std::cout, mesh[i] );
    
    if (mesh[i].empty())
      continue;
    have_sets = true;
    
      // Get tag
    MBTag tag;
    result = gMB->tag_get_handle( mesh_tag_names[i], tag );
    if (MB_SUCCESS != result) 
    {
      std::cerr << "Tag not found: " << mesh_tag_names[i] << std::endl;
      continue;
    }
    int size;
    result = gMB->tag_get_size( tag, size );
    if (MB_SUCCESS != result || size != sizeof(int))
    {
      std::cerr << "Tag invalid: " << mesh_tag_names[i] << std::endl;
      continue;
    }

      // get entity sets
    int init_count = set_list.size();
    for (std::set<int>::iterator iter = mesh[i].begin(); iter != mesh[i].end(); ++iter) 
    {
      range.clear();
      const void* vals[] = { &*iter };
      result = gMB->get_entities_by_type_and_tag( 0, MBENTITYSET, &tag, vals, 1, range );
      if (MB_SUCCESS != result || range.empty()) 
      {
        range.clear();
        std::cerr << mesh_tag_names[i] << " " << *iter << " not found.\n";
      }
      std::copy( range.begin(), range.end(), std::back_inserter(set_list) );
    }
      
    if (verbose)
      std::cout << "Found " << (set_list.size()-init_count) << ' '
                << mesh_tag_names[i] << " sets" << std::endl;
  }
  
    // Check if output is limited to certain dimensions of elements
  bool bydim = false;
  for (dim = 1; dim < 4; ++dim)
    if (dims[dim])
      bydim = true;
  
    // Check conflicting input
  if (bydim) {
    if (generate[1] && !dims[1]) {
      std::cerr << "Warning: Request to generate 1D internal entities but not export them." << std::endl;
      generate[1] = false;
    } 
     if (generate[2] && !dims[2]) {
      std::cerr << "Warning: Request to generate 2D internal entities but not export them." << std::endl;
      generate[2] = false;
    } 
  }
 
    // Generate any internal entities
  if (generate[1] || generate[2]) {
    MBEntityHandle all_mesh = 0;
    const MBEntityHandle* sets = &all_mesh;
    int num_sets = 1;
    if (have_sets) {
      num_sets = set_list.size();
      sets = &set_list[0];
    }
    for (i = 0; i < num_sets; ++i) {
      MBRange dim3, dim2, adj;
      gMB->get_entities_by_dimension( sets[i], 3, dim3, true );
      if (generate[1]) {
        gMB->get_entities_by_dimension( sets[i], 2, dim2, true );
        gMB->get_adjacencies( dim3, 1, true, adj, MBInterface::UNION );
        gMB->get_adjacencies( dim2, 1, true, adj, MBInterface::UNION );
      }
      if (generate[2]) {
        gMB->get_adjacencies( dim3, 2, true, adj, MBInterface::UNION );
      }
      if (sets[i])
        gMB->add_entities( sets[i], adj );
    }
  }
      
    // Delete any entities not of the dimensions to be exported
  if (bydim) {
      // Get list of dead elements
    MBRange dead_entities , tmp_range;
    for (dim = 1; dim <= 3; ++dim) {
      if (dims[dim])
        continue;
      gMB->get_entities_by_dimension(0, dim, tmp_range );
      dead_entities.merge( tmp_range );
    }
      // Remove dead entities from all sets, and add all 
      // empty sets to list of dead entities.
    MBRange empty_sets;
    remove_entities_from_sets( gMB, dead_entities, empty_sets );
    while (!empty_sets.empty()) {
      if (!set_list.empty())
        remove_from_vector( set_list, empty_sets );
      dead_entities.merge( empty_sets );
      MBRange tmp_range;
      remove_entities_from_sets( gMB, empty_sets, tmp_range );
      empty_sets = subtract( tmp_range,  dead_entities );
    }
      // Destroy dead entities
    result = gMB->delete_entities( dead_entities );
  }
  
    // If user specified sets to write, but none were found, exit.
  if (have_sets && set_list.empty())
  {
    std::cerr << "Nothing to write." << std::endl;
    return ENT_NOT_FOUND;
  }
  
  if (verbose)
  {
    if (have_sets)
      std::cout << "Found " << set_list.size() 
              << " specified sets to write (total)." << std::endl;  
    else
      std::cout << "No sets specifed.  Writing entire mesh." << std::endl; 
  }  
  
    // Write the output file
  reset_times();
  if (have_sets) 
    result = gMB->write_file( out, format, write_options.c_str(), &set_list[0], set_list.size() );
  else
    result = gMB->write_file( out, format, write_options.c_str() );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to write \"" << out << "\"." << std::endl; 
    std::cerr  << "Error code: " << gMB->get_error_string(result) << " (" << result << ")" << std::endl;
    std::string message;
    if (MB_SUCCESS == gMB->get_last_error(message) && !message.empty())
      std::cerr << "Error message: " << message << std::endl;
    return WRITE_ERROR;
  }
  std::cerr << "Wrote \"" << out << "\"" << std::endl;
  if (print_times) write_times( std::cerr );

#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
}

bool parse_id_list( const char* string, std::set<int>& results )
{
  bool okay = true;
  char* mystr = strdup( string );
  for (const char* ptr = strtok(mystr, ","); ptr; ptr = strtok(0,","))
  {
    char* endptr;
    long val = strtol( ptr, &endptr, 0 );
    if (endptr == ptr || val <= 0) {
      std::cerr << "Not a valid id: " << ptr << std::endl;
      okay = false;
      break;
    }
    
    long val2 = val;
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

void print_id_list( const char* head, std::ostream& stream, const std::set<int>& list )
{
  stream << head << ": ";
  
  if (list.empty())
  {
    stream << "(none)" << std::endl;
    return;
  }
  
  int start, prev;
  std::set<int>::const_iterator iter = list.begin();
  start = prev = *(iter++);
  for (;;)
  {
    if (iter == list.end() || *iter != 1+prev) {
      stream << start;
      if (prev != start)
        stream << '-' << prev;
      if (iter == list.end())
        break;
      stream << ", ";
      start = *iter;
    }
    prev = *(iter++);
  }
  
  stream << std::endl;
}
    
    


void print_time( int clk_per_sec, const char* prefix, clock_t ticks, std::ostream& stream )
{
  ticks *= clk_per_sec/100;
  clock_t centi = ticks % 100;
  clock_t seconds = ticks / 100;
  if (seconds < 120)
  {
    stream << prefix << (ticks / 100) << "." << centi << "s" << std::endl;
  }
  else
  {
    clock_t minutes = (seconds / 60) % 60;
    clock_t hours = (seconds / 3600);
    seconds %= 60;
    if (hours)
      stream << hours << "h";
    if (minutes)
      stream << minutes << "m";
    if (seconds || centi)
      stream << seconds << "." << centi << "s";
    stream << " (" << (ticks/100) << "." << centi << "s)" << std::endl;
  }
}

clock_t usr_time, sys_time, abs_time;

#ifdef WIN32

void reset_times() 
{
  abs_time = clock();
}


void write_times( std::ostream& stream ) 
{
  clock_t abs_tm = clock();
  print_time( CLOCKS_PER_SEC, "  ", abs_tm - abs_time, stream );
  abs_time = abs_tm;
}

#else

void reset_times()
{
  tms timebuf;
  abs_time = times( &timebuf );
  usr_time = timebuf.tms_utime;
  sys_time = timebuf.tms_stime;
}

void write_times( std::ostream& stream )
{
  clock_t usr_tm, sys_tm, abs_tm;
  tms timebuf;
  abs_tm = times( &timebuf );
  usr_tm = timebuf.tms_utime;
  sys_tm = timebuf.tms_stime;
  print_time( sysconf(_SC_CLK_TCK), "  real:   ", abs_tm - abs_time, stream );
  print_time( sysconf(_SC_CLK_TCK), "  user:   ", usr_tm - usr_time, stream );
  print_time( sysconf(_SC_CLK_TCK), "  system: ", sys_tm - sys_time, stream );
  abs_time = abs_tm;
  usr_time = usr_tm;
  sys_time = sys_tm;
}

#endif

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


void list_formats( MBInterface* gMB )
{
  const char iface_name[] = "MBReaderWriterSet";
  MBErrorCode err;
  void* void_ptr = 0;
  MBReaderWriterSet* set;
  MBReaderWriterSet::iterator i;
  std::ostream& str = std::cout;
    
    // get MBReaderWriterSet
  err = gMB->query_interface( iface_name, &void_ptr );
  if (err != MB_SUCCESS || !void_ptr) {
    std::cerr << "Internal error:  Interface \"" << iface_name 
              << "\" not available.\n";
    exit(OTHER_ERROR);
  }
  set = (MBReaderWriterSet*)void_ptr;
  
    // get field with for format description
  size_t w = 0;
  for (i = set->begin(); i != set->end(); ++i)
    if (i->description().length() > w)
      w = i->description().length();
  
    // write table header
  str << "Format  " << std::setw(w) << std::left << "Description"
      << "  Read  Write  File Name Suffixes\n"
      << "------  " << std::setw(w) << std::setfill('-') << "" << std::setfill(' ')
      << "  ----  -----  ------------------\n";
      
    // write table data
  for (i = set->begin(); i != set->end(); ++i)
  {
    std::vector<std::string> ext;
    i->get_extensions( ext );
    str << std::setw(6) << i->name() << "  " 
        << std::setw(w) << std::left << i->description() << "  "
        << (i->have_reader() ?  " yes" :  "  no") << "  "
        << (i->have_writer() ? "  yes" : "   no") << " ";
    for (std::vector<std::string>::iterator j = ext.begin(); j != ext.end(); ++j)
      str << " " << *j;
    str << std::endl;
  }
  str << std::endl;
  
  gMB->release_interface( iface_name, void_ptr );
  exit(0);
}

void remove_entities_from_sets( MBInterface* gMB, MBRange& dead_entities, MBRange& empty_sets )
{
  empty_sets.clear();
  MBRange sets;
  gMB->get_entities_by_type( 0, MBENTITYSET, sets );
  for (MBRange::iterator i = sets.begin(); i != sets.end(); ++i) {
    MBRange set_contents;
    gMB->get_entities_by_handle( *i, set_contents, false );
    set_contents = intersect( set_contents, dead_entities );
    gMB->remove_entities( *i, set_contents );
    set_contents.clear();
    gMB->get_entities_by_handle( *i, set_contents, false );
    if (set_contents.empty())
      empty_sets.insert( *i );
  }
}

void remove_from_vector( std::vector<MBEntityHandle>& vect, const MBRange& ents_to_remove )
{
  MBRange::const_iterator i;
  std::vector<MBEntityHandle>::iterator j;
  for (i = ents_to_remove.begin(); i != ents_to_remove.end(); ++i) {
    j = std::find( vect.begin(), vect.end(), *i );
    if (j != vect.end())
      vect.erase( j );
  }
}
