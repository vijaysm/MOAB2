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

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/ReaderWriterSet.hpp"
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
#  include "moab_mpi.h"
#endif
#include <stdio.h>

/* Exit values */
#define USAGE_ERROR 1
#define READ_ERROR 2
#define WRITE_ERROR 3
#define OTHER_ERROR 4
#define ENT_NOT_FOUND 5

using namespace moab;


void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage: " << name << 
    " [-t] [-h] [-g] [-l] [-O <read_option> [-O <reader_option> ...] ] <input_file> " << std::endl
    << "\t-t             - Time read of files." << std::endl
    << "\t-g             - Enable verbose/debug output." << std::endl
    << "\t-h             - Print this help text and exit." << std::endl
    << "\t-l             - List available file formats and exit." << std::endl
    << "\t-O option      - Specify read option." << std::endl
    << "\t--             - treat all subsequent options as file names" << std::endl
    << "\t                 (allows file names beginning with '-')" << std::endl
    ;
}

void print_help( const char* name )
{
  std::cout << 
    " This program reads a mesh file and reports times taken to read it.\n";
  print_usage( name, std::cout );
  exit(0);
}

void usage_error( const char* name )
{
  print_usage( name, std::cerr );
#ifdef USE_MPI
  MPI_Finalize();
#endif
  exit(USAGE_ERROR);
} 


void list_formats( Interface* );
void reset_times();
void write_times( std::ostream& stream );
bool make_opts_string( std::vector<std::string> options, std::string& result );

int main(int argc, char* argv[])
{
#ifdef USE_MPI
  MPI_Init(&argc,&argv);
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank); 
  MPI_Comm_size(MPI_COMM_WORLD, &size); 
#endif


  Interface* gMB;
  ErrorCode result;
  Range range;

  gMB = new Core();

  int i;
  const char* in = NULL;    // input file name
  bool verbose = false;
  std::vector<std::string> read_opts;
    
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
            case 'O':  read_opts.push_back(argv[i]); pval = true; break;
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
    else  { // too many file names
      std::cerr << "Unexpexed argument: " << argv[i] << std::endl;
      usage_error(argv[0]);
    }
  }
  if (!in) {
    std::cerr << "No input file name specified." << std::endl;
    usage_error(argv[0]);
  }
    

    // construct options string from individual options
  std::string read_options;
  if (!make_opts_string(  read_opts,  read_options ))
  {
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return USAGE_ERROR;
  }
  
    // Read the input file.
  reset_times();
  result = gMB->load_file( in, 0, read_options.c_str() );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to load \"" << in << "\"." << std::endl;
    std::cerr  << "Error code: " << gMB->get_error_string(result) << " (" << result << ")" << std::endl;
    std::string message;
    if (MB_SUCCESS == gMB->get_last_error(message) && !message.empty())
      std::cerr << "Error message: " << message << std::endl;
#ifdef USE_MPI
    MPI_Finalize();
#endif
    return READ_ERROR;
  }
#ifdef USE_MPI
  std::cout << "[proc " << rank << " of " << size << "]: ";
#endif
  std::cout << "Read \"" << in << "\"" << std::endl;
  if (print_times) write_times( std::cout );
  
#ifdef USE_MPI
  MPI_Finalize();
#endif
  return 0;
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


void list_formats( Interface* gMB )
{
  const char iface_name[] = "ReaderWriterSet";
  ErrorCode err;
  void* void_ptr = 0;
  ReaderWriterSet* set;
  ReaderWriterSet::iterator i;
  std::ostream& str = std::cout;
    
    // get ReaderWriterSet
  err = gMB->query_interface( iface_name, &void_ptr );
  if (err != MB_SUCCESS || !void_ptr) {
    std::cerr << "Internal error:  Interface \"" << iface_name 
              << "\" not available.\n";
    exit(OTHER_ERROR);
  }
  set = (ReaderWriterSet*)void_ptr;
  
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


