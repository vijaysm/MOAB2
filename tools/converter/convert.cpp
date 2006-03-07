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

#include "MBCore.hpp"
#include "MBRange.hpp"
#include "MBTagConventions.hpp"
#include "MBReaderWriterSet.hpp"
#include <iostream>
#include <iomanip>
#include <set>
#ifndef WIN32
#  include <sys/times.h>
#  include <limits.h>
#include <unistd.h>
#endif
#include <time.h>

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
    " [-a <sat_file>|-A] [-t] [subset options] <input_file> <output_file>" << std::endl
    << "\t-a <acis_file> - ACIS SAT file dumped by .cub reader" << std::endl
    << "\t-A             - .cub file reader should not dump a SAT file" << std::endl
    << "\t-t             - Time read and write of files." << std::endl
    << "\t-g             - Enable verbose/debug output." << std::endl
    << "\t-h             - Print this help text and exit." << std::endl
    << "\t-l             - List available file formats and exit." << std::endl
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

int main(int argc, char* argv[])
{
  MBInterface* gMB;
  MBErrorCode result;
  MBRange range;
  
    // Get MB instance
  gMB = new MBCore();

  int i, dim;
  const char* in = NULL;    // input file name
  const char* out = NULL;   // output file name
  const char* acis = NULL;  // file to which to write geom data from .cub files.
  bool no_acis = false;
  bool verbose = false;
  std::set<int> geom[4], mesh[3];       // user-specified IDs 
  std::vector<MBEntityHandle> set_list; // list of user-specified sets to write
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
        case 'A': no_acis = true;        break;
        case 'h': 
        case 'H': print_help( argv[0] ); break;
        case 'l': list_formats( gMB );   break;
          // do options that require additional args:
        default: 
          ++i;
          if (i == argc || argv[i][0] == '-') {
            std::cerr << "Expected argument following " << argv[i-1] << std::endl;
            usage_error(argv[0]);
          }
          pval = false;
          switch ( argv[i-1][1] )
          {
            case 'a': acis = argv[i]; pval = true;              break;
            case 'v': pval = parse_id_list( argv[i], geom[3] ); break;
            case 's': pval = parse_id_list( argv[i], geom[2] ); break;
            case 'c': pval = parse_id_list( argv[i], geom[1] ); break;
            case 'V': pval = parse_id_list( argv[i], geom[0] ); break;
            case 'm': pval = parse_id_list( argv[i], mesh[2] ); break;
            case 'n': pval = parse_id_list( argv[i], mesh[1] ); break;
            case 'd': pval = parse_id_list( argv[i], mesh[0] ); break;
            default: std::cerr << "Invalid option: " << argv[i] << std::endl;
          }
          
          if (!pval)
            usage_error(argv[0]);
      }
    }
      // do file names
    else if (!in)
      in = argv[i];
    else if (!out)
      out = argv[i];
    else  // too many file names
      usage_error(argv[0]);
  }
  if (!in || !out)
    usage_error(argv[0]);
  
    // If requested, set mesh tag to indicate SAT file name to
    // dump from .cub file reader.
  if (no_acis || acis)
  {
    MBTag tag;
    result = gMB->tag_get_handle( acis_dump_file_tag_name, tag );
    if (result == MB_TAG_NOT_FOUND)
      result = gMB->tag_create( acis_dump_file_tag_name, sizeof(char*), MB_TAG_SPARSE, tag, NULL );
    else if(result != MB_SUCCESS)
    {
      std::cerr << "Failed to set ACIS dump file." << std::endl;
      return OTHER_ERROR;
    }

    result = gMB->tag_set_data( tag, 0, 0, &acis );
    if (result != MB_SUCCESS)
    {
      std::cerr << "Failed to set ACIS dump file." << std::endl;
      return OTHER_ERROR;
    }
    
    if (verbose)
    {
      if (no_acis)
        std::cout << "Disabling write of ACIS geometry data." << std::endl;
      else
        std::cout << "ACIS geometry data will be written to: " << acis << std::endl;
    }
  }
  
    // Read the input file.
  reset_times();
  result = gMB->load_mesh( in );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to load \"" << in << "\"." << std::endl; 
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
                << mesh_tag_names[dim] << " sets" << std::endl;
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
    result = gMB->write_mesh( out, &set_list[0], set_list.size() );
  else
    result = gMB->write_mesh( out );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to write \"" << out << "\"." << std::endl; 
    return WRITE_ERROR;
  }
  std::cerr << "Wrote \"" << out << "\"" << std::endl;
  if (print_times) write_times( std::cerr );

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
    
    


void print_time( const char* prefix, clock_t ticks, std::ostream& stream )
{
  ticks *= sysconf(_SC_CLK_TCK)/100;
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
  print_time( "  ", abs_tm - abs_time, stream );
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
  print_time( "  real:   ", abs_tm - abs_time, stream );
  print_time( "  user:   ", usr_tm - usr_time, stream );
  print_time( "  system: ", sys_tm - sys_time, stream );
  abs_time = abs_tm;
  usr_time = usr_tm;
  sys_time = sys_tm;
}

#endif

void list_formats( MBInterface* gMB )
{
  const char iface_name[] = "MBReaderWriterSet";
  MBErrorCode err;
  void* void_ptr = 0;
  MBReaderWriterSet* set;
  MBReaderWriterSet::iter_type i;
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
  str << std::setw(w) << std::left << "Format"
      << "  Read  Write  Extensions\n"
      << std::setw(w) << std::setfill('-') << "" << std::setfill(' ')
      << "  ----  -----  ----------\n";
      
    // write table data
  for (i = set->begin(); i != set->end(); ++i)
  {
    std::vector<std::string> ext;
    i->get_extensions( ext );
    str << std::setw(w) << std::left << i->description() << "  "
        << (i->have_reader() ?  " yes" :  "  no") << "  "
        << (i->have_writer() ? "  yes" : "   no") << " ";
    for (std::vector<std::string>::iterator j = ext.begin(); j != ext.end(); ++j)
      str << " " << *j;
    str << std::endl;
  }
  str << std::endl;
  
  exit(0);
}

