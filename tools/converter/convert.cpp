#include "MBCore.hpp"
#include <iostream>
#ifndef WIN32
#  include <sys/times.h>
#  include <limits.h>
#endif
#include <time.h>

/* Exit values */
#define USAGE_ERROR 1
#define READ_ERROR 2
#define WRITE_ERROR 3
#define OTHER_ERROR 4

/* Tag to control ACIS dump from .cub file reader */
const char acis_dump_file_tag_name[] = "__ACISDumpFile";


void usage_error( const char* name )
{
  std::cerr << "Usauge: " << name << 
    " [-d <sat_file>|-D] [-t] <input_file> <output_file>" << std::endl
    << "\t-d <acis_file> - ACIS SAT file dumped by .cub reader" << std::endl
    << "\t-D             - .cub file reader should not dump a SAT file" << std::endl
    << "\t-t             - Time read and write of files." << std::endl
    << "\t--             - treat all subsequent options as file names" << std::endl
    << "\t                 (allows file names beginning with '-')" << std::endl
    ;
  
  exit(USAGE_ERROR);
}


void reset_times();
void write_times( std::ostream& stream );

int main(int argc, char* argv[])
{
  MBInterface* gMB;
  MBErrorCode result;

    // scan arguments
  int i;
  const char* in = NULL;
  const char* out = NULL;
  const char* acis = NULL;
  bool do_acis = false;
  bool no_acis = false;
  bool do_flag = true;
  bool print_times = false;
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
        case '-': do_flag = false;    break;
        case 't': print_times = true; break;
        case 'd': do_acis = true;     break;
        case 'D': no_acis = true;     break;
        default: usage_error(argv[0]);
      }
    }
    else if (do_acis)
    {
      acis = argv[i];
      do_acis = false;
    }
    else if (!in)
      in = argv[i];
    else if (!out)
      out = argv[i];
    else
      usage_error(argv[0]);
  }
  if (do_acis || !in || !out)
    usage_error(argv[0]);
  
    // Get MB instance
  gMB = new MBCore();
  
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
  
    // Write the output file
  reset_times();
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
