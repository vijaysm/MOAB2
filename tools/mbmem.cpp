#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/CN.hpp"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <string>
#include <cstdio>
#include <cstdlib>

#ifndef _MSC_VER
#include <sys/times.h>
#include <sys/resource.h>
#endif

void usage( const char* argv0 )
{
  std::cerr << "Usage: " << argv0 << " [-h|-b|-k|-m] <filename> [<filename> ...]" << std::endl
            << "  -h : human readable units" << std::endl
            << "  -b : bytes" << std::endl
            << "  -k : kilobytes (1 kB == 1024 bytes)" << std::endl
            << "  -m : megabytes (1 MB == 1024 kB)" << std::endl
            << "  -g : gigabytes (1 GB == 1024 MB)" << std::endl
            << std::endl;
  std::exit(1);
}

enum Units { HUMAN, BYTES, KILOBYTES, MEGABYTES, GIGABYTES };
Units UNITS = HUMAN;

  // The core functionality of this example
void print_memory_stats( moab::Interface& mb );
  
  // main routine: read any specified files and call print_memory_stats
int main( int argc, char* argv[] )
{
  moab::Core mbcore;
  moab::Interface& mb = mbcore;
  moab::ErrorCode rval;
  bool no_more_flags = false;

    // load each file specified on command line
  for (int i = 1; i < argc; ++i) {
    if (!no_more_flags && argv[i][0] == '-') {
      if (!strcmp(argv[i],"-h"))
        UNITS = HUMAN;
      else if(!strcmp(argv[i],"-b"))
        UNITS = BYTES;
      else if(!strcmp(argv[i],"-k"))
        UNITS = KILOBYTES;
      else if(!strcmp(argv[i],"-m"))
        UNITS = MEGABYTES;
      else if(!strcmp(argv[i],"-g"))
        UNITS = GIGABYTES;
      else if(!strcmp(argv[i],"--"))
        no_more_flags = true;
      else
        usage(argv[0]);
    }
    else {
      rval = mb.load_file(argv[i]);

        // if file load failed, print some info and exit
      if (moab::MB_SUCCESS != rval) {
        std::string message;
        mb.get_last_error( message );
        std::cerr << mb.get_error_string(rval) << ": " << message << std::endl
                  << argv[i] << ": Failed to read file." << std::endl;
        return 1;
      }
    }
  }
  
    // print summary of MOAB's memory use
  print_memory_stats(mb);
  return 0;
}

  // struct to store memory stats
struct MemStats {
  unsigned long total_storage;
  unsigned long total_amortized;
  unsigned long entity_storage;
  unsigned long entity_amortized;
  unsigned long adjacency_storage;
  unsigned long adjacency_amortized;
  unsigned long tag_storage;
  unsigned long tag_amortized;
};

 // test if MemStats object indicates no memory
bool is_zero( const MemStats& stats );

  // populdate a MemStats structg by calling 
  // moab::Interface::estimated_memory_use
void get_mem_stats( moab::Interface& mb,
                    MemStats& data,
                    moab::EntityType type = moab::MBMAXTYPE );

  // Formatted string representation of memory size value
std::string memstr( unsigned long val );

  // Get string describing tag data type
std::string tag_type_string( moab::Interface& mb, moab::Tag tag );

  // Get string representation of tag storage type
std::string tag_storage_string( moab::Interface& mb, moab::Tag tag );

  // Center
std::string center( const char* str, size_t width );

void print_memory_stats( moab::Interface& mb )
{
  moab::ErrorCode rval;
  const char ANON_TAG_NAME[] = "(anonymous)";
  const int TYPE_WIDTH = 10;
  const int MEM_WIDTH = 7;
  const int MEM2_WIDTH = 2*MEM_WIDTH+1;
  const int MIN_TAG_NAME_WIDTH = strlen(ANON_TAG_NAME);
  const int DTYPE_WIDTH = 12;
  const int STORAGE_WIDTH = 8;

    // per-entity-type table header
  MemStats stats;
  std::cout.fill(' ');
  std::cout << std::left << std::setw(TYPE_WIDTH) << "Type" << ' '
            << center("Total",MEM2_WIDTH) << ' '
            << center("Entity",MEM2_WIDTH) << ' '
            << center("Adjacency",MEM2_WIDTH) << ' '
            << center("Tag",MEM2_WIDTH) << ' '
            << std::endl << std::setw(TYPE_WIDTH) << " ";
  for (int i = 0; i < 4; ++i)
    std::cout << ' ' << std::left << std::setw(MEM_WIDTH) << "Used"
              << ' ' << std::left << std::setw(MEM_WIDTH) << "Alloc";
  std::cout << std::endl;
  std::cout.fill('-');
  std::cout << std::setw(TYPE_WIDTH) << '-';
  for (int i = 0; i < 8; ++i)
    std::cout << ' ' << std::setw(MEM_WIDTH) << '-';
  std::cout.fill(' ');
  std::cout << std::endl;
  
    // per-entity-type memory use
  for (moab::EntityType t = moab::MBVERTEX; t != moab::MBMAXTYPE; ++t) {
    get_mem_stats( mb, stats, t );
    if (is_zero(stats)) continue;  // skip types with no allocated memory
    
    std::cout << std::left << std::setw(TYPE_WIDTH) << moab::CN::EntityTypeName(t) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.total_storage) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.total_amortized) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.entity_storage) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.entity_amortized) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.adjacency_storage) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.adjacency_amortized) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.tag_storage) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(stats.tag_amortized) << std::endl;
  }
  
    // get list of tags
  std::vector<moab::Tag> tags;
  std::vector<moab::Tag>::const_iterator ti;
  mb.tag_get_tags( tags );
  
    // figure out required field with to fit longest tag name 
  unsigned maxlen = MIN_TAG_NAME_WIDTH;
  for (ti = tags.begin(); ti != tags.end(); ++ti) {
    std::string name;
    rval = mb.tag_get_name( *ti, name );
    if (moab::MB_SUCCESS != rval)
      continue;
    if (name.size() > maxlen)
      maxlen = name.size();
  }
  
    // print header for per-tag data
  if (!tags.empty()) {
    std::cout.fill(' ');
    std::cout << std::endl
              << std::left << std::setw(maxlen) << "Tag Name" << ' '
              << std::left << std::setw(DTYPE_WIDTH) << "Type" << ' '
              << std::left << std::setw(STORAGE_WIDTH) << "Storage" << ' '
              << std::left << std::setw(MEM_WIDTH) << "Used" << ' '
              << std::left << std::setw(MEM_WIDTH) << "Alloc" << std::endl;
    std::cout.fill('-');
    std::cout << std::setw(maxlen) << '-' << ' '
              << std::setw(DTYPE_WIDTH) << '-' << ' '
              << std::setw(STORAGE_WIDTH) << '-' << ' '
              << std::setw(MEM_WIDTH) << '-' << ' '
              << std::setw(MEM_WIDTH) << '-' << std::endl;
    std::cout.fill(' ');
  }
              
    // print per-tag memory use
  for (ti = tags.begin(); ti != tags.end(); ++ti) {
    std::string name;
    rval = mb.tag_get_name( *ti, name );
    if (moab::MB_SUCCESS != rval || name.empty())
      name = ANON_TAG_NAME;
    
    unsigned long occupied, allocated;
    mb.estimated_memory_use( 0, 0, 0, 0, 0, 0, 0, 0, &*ti, 1, &occupied, &allocated );
    
    std::cout << std::left << std::setw(maxlen) << name << ' '
              << std::right << std::setw(DTYPE_WIDTH) << tag_type_string(mb,*ti) <<  ' '
              << std::right << std::setw(STORAGE_WIDTH) << tag_storage_string(mb,*ti) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(occupied) << ' '
              << std::right << std::setw(MEM_WIDTH) << memstr(allocated) << std::endl;
  }
  
    // print summary of overall memory use
  get_mem_stats( mb, stats );
  std::cout << std::endl
            << "TOTAL: (Used/Allocated)" << std::endl
            << "memory:    " << memstr(stats.total_storage) << "/" << memstr(stats.total_amortized) << std::endl
            << "entity:    " << memstr(stats.entity_storage) << "/" << memstr(stats.entity_amortized) << std::endl
            << "adjacency: " << memstr(stats.adjacency_storage) << "/" << memstr(stats.adjacency_amortized) << std::endl
            << "tag:       " << memstr(stats.tag_storage) << "/" << memstr(stats.tag_amortized) << std::endl
            << std::endl;

  std::FILE* filp = std::fopen("/proc/self/stat", "r");
  unsigned long vsize, rss;
  if (filp && 2 == std::fscanf(filp,
                "%*d " // pid
                "%*s " // comm
                "%*c " // state
                "%*d " // ppid
                "%*d " // pgrp
                "%*d " // session
                "%*d " // tty_nr
                "%*d " // tpgid
                "%*u " // flags
                "%*u " // minflt
                "%*u " // cminflt
                "%*u " // majflt
                "%*u " // cmajflt
                "%*u " // utime
                "%*u " // stime
                "%*d " // cutime
                "%*d " // cstime
                "%*d " // priority
                "%*d " // nice
                "%*d " // num_threads
                "%*d " // itrealvalue
                "%*u " // starttime
                "%lu " // vsize
                "%ld", // rss
                &vsize, &rss )) {
#ifndef _MSC_VER
    rss *= getpagesize();
#endif
    std::cout << std::endl << "SYSTEM:" 
              << std::endl << "Virtual memory:    " << memstr(vsize)
              << std::endl << "Resident set size: " << memstr(rss)
              << std::endl;
  }
  else {
#ifndef _MSC_VER
    struct rusage sysdata;
    if (getrusage( RUSAGE_SELF, &sysdata )) {
      std::cerr << "getrusage failed" << std::endl;
    }
    else {
      unsigned long rss = sysdata.ru_maxrss;
      rss *= getpagesize();
      std::cerr << std::endl << "SYSTEM:"
                << std::endl << "Resident set size: " << memstr(rss) 
                << std::endl;
    }
#endif
  }
}


bool is_zero( const MemStats& stats ) { return stats.total_amortized == 0; }

void get_mem_stats( moab::Interface& mb,
                    MemStats& data,
                    moab::EntityType type )
{
  if (type != moab::MBMAXTYPE) {
    moab::Range range;
    mb.get_entities_by_type( 0, type, range );
    mb.estimated_memory_use( range, 
                             &data.total_storage,
                             &data.total_amortized,
                             &data.entity_storage,
                             &data.entity_amortized,
                             &data.adjacency_storage,
                             &data.adjacency_amortized,
                             0, 0,
                             &data.tag_storage,
                             &data.tag_amortized );
  }
  else {
    mb.estimated_memory_use( 0, 0, 
                             &data.total_storage,
                             &data.total_amortized,
                             &data.entity_storage,
                             &data.entity_amortized,
                             &data.adjacency_storage,
                             &data.adjacency_amortized,
                             0, 0,
                             &data.tag_storage,
                             &data.tag_amortized );
  }
}

// rounded division
unsigned long rdiv( unsigned long num, unsigned long den )
{
  return (num + den/2) / den;
}

std::string memstr( unsigned long val )
{
  const unsigned long kb = 1024;
  const unsigned long mb = kb*kb;
  const unsigned long gb = kb*mb;
  const unsigned long tb = kb*gb;
  
  std::ostringstream s;
  if (UNITS == HUMAN) {
    if (val >= 10*tb)
      s << rdiv( val, tb ) << "TB";
    else if (val >= 10*gb)
      s << rdiv( val, gb ) << "GB";
    else if (val >= 10*mb)
      s << rdiv( val, mb ) << "MB";
    else if (val >= 10*kb)
      s << rdiv( val, kb ) << "kB";
    else if (val > 0)
      s << val << " B";
    else 
      s << "0  ";
  }
  else {
    unsigned long den;
    switch (UNITS) {
      case BYTES: den = 1; break;
      case KILOBYTES: den = kb; break;
      case MEGABYTES: den = mb; break;
      case GIGABYTES: den = gb; break;
      case HUMAN: break; // handled above, list here to suppress warning
    }
    
    s << rdiv( val, den );
  }
  return s.str();
}

std::string tag_type_string( moab::Interface& mb, moab::Tag tag )
{
  moab::ErrorCode rval;
  std::ostringstream s;
  
  moab::DataType type;
  rval = mb.tag_get_data_type( tag, type );
  if (moab::MB_SUCCESS != rval)
    return std::string();

  int typesize;
  std::string typestr;
  switch (type) {
    case moab::MB_TYPE_INTEGER:
      typestr = "int";
      typesize = sizeof(int);
      break;
    case moab::MB_TYPE_DOUBLE:
      typestr = "double";
      typesize = sizeof(double);
      break;
    case moab::MB_TYPE_HANDLE:
      typestr = "handle";
      typesize = sizeof(moab::EntityHandle);
      break;
    case moab::MB_TYPE_BIT:
      typesize = 1;
      typestr = "bits";
      break;
    case moab::MB_TYPE_OPAQUE:
      typesize = 1;
      typestr = "bytes";
      break;
    default:
      typesize = 1;
      typestr = "???";
      break;
  }

  int size;
  rval = mb.tag_get_size( tag, size );
  if (moab::MB_VARIABLE_DATA_LENGTH == rval) 
    s << "VAR " << typestr;
  else if (moab::MB_SUCCESS == rval) 
    s << size/typesize << " " << typestr;
  // else do nothing
  
  return s.str();
}

std::string tag_storage_string( moab::Interface& mb, moab::Tag tag )
{
  moab::ErrorCode rval;
  moab::TagType type;
  rval = mb.tag_get_type( tag, type );
  if (moab::MB_SUCCESS != rval)
    return std::string();

  switch (type) {
    case moab::MB_TAG_DENSE: return "dense"; 
    case moab::MB_TAG_SPARSE: return "sparse"; 
    case moab::MB_TAG_BIT: return "bit"; 
    default: return "(none)"; 
  }
}

 
std::string center( const char* str, size_t width )
{
  std::string text(str);
  if (text.size() >= width)
    return text;
 
  width -= text.size();
  if (1u == width) {
    text += " ";
    return text;
  }
 
  std::ostringstream s;
  s << std::setw(width/2) << ' ' << text << std::setw(width/2 + width%2) << ' ';
  return s.str();
}

