#include "RangeMap.hpp"
#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>

RangeMap<long,long,0> mapInstance;


bool insert( std::vector<long>& values );
bool verify( std::vector<long>& values );
bool print( std::vector<long>& values );
bool clear( std::vector<long>& values );
bool help();

bool do_command( char* cmd )
{
  std::vector<long> values;
  char *ptr, *end;
  const char space[] = " \t\r\n";
  char c = cmd[0];
  ++cmd;
  
  for (ptr = strtok( cmd, space ); ptr; ptr = strtok( 0, space )) {
    long val = strtol( ptr, &end, 0 );
    if (*end) {
      std::cerr << "Invalid integer value: " << ptr << std::endl;
      return false;
    }
    values.push_back( val );
  }
  
  switch (c) {
    case 'i': return insert( values );
    case 'v': return verify( values );
    case 'p': return print( values );
    case 'c': return clear( values );
    case 'h': return help( );
    default: 
      std::cerr << "Unknown command: " << c << std::endl;
      return false;
  }
}


bool insert( std::vector<long>& values )
{
  if (values.size() % 3) {
    std::cerr << "Invalid arguments for insertion.  "
                 "Insertion requires tuples of three values:  "
                 "{start_key,start_val,count}" << std::endl;
    return false;
  }
  for (unsigned i = 0; i < values.size(); i += 3) 
    if (mapInstance.end() == mapInstance.insert( values[i], values[i+1], values[i+2] ))
      std::cout << "Insertion of {" << values[i] << ", " 
                                    << values[i+1] << ", "
                                    << values[i+2] << "} FAILED" << std::endl;
  return true;
}

bool clear( std::vector<long>& values )
{
  if (values.size() % 2) {
    std::cerr << "(c)lear command expects ranges of values (or no values)" << std::endl;
    return false;
  }
  
  if(values.empty()) {
    mapInstance.clear();
    return true;
  }

  for (unsigned i = 0; i < values.size(); i += 2) 
    mapInstance.erase( values[i], values[i+1] );

  return true;
}

bool print( std::vector<long>& values )
{
  const int W = 20;
  if (values.empty()) {
    std::cout  << std::setw(W) << "Start Key" 
               << std::setw(W) << "Start Val" 
               << std::setw(W) << "Count" 
               << std::endl;
    RangeMap<long,long,0>::iterator i = mapInstance.begin();
    for (; i != mapInstance.end(); ++i) 
      std::cout << std::setw(W) << i->begin 
                << std::setw(W) << i->value 
                << std::setw(W) << i->count
                << std::endl;
  }
  else {
    for (unsigned i = 0; i < values.size(); ++i) 
      std::cout << std::setw(W) << values[i] << "->" 
                << std::setw(W) << mapInstance.find(values[i]) << std::endl;
  }
  return true;
}

bool verify( std::vector<long>& values )
{
  if (values.size() % 3) {
    std::cerr << "Invalid arguments for insertion.  "
                 "Insertion requires tuples of three values:  "
                 "{start_key,start_val,count}" << std::endl;
    return false;
  }
 
  bool result = true;
  RangeMap<long,long,0>::iterator i = mapInstance.begin();
  std::vector<long>::iterator j = values.begin();
  for (; i != mapInstance.end(); ++i) 
    if (j == values.end() ||
        i->begin != *j++ || 
        i->value != *j++ || 
        i->count != *j++)
      result = false;
      
  return result && j == values.end();
}

bool help()
{
  std::cout << "Commands: i p c v h" << std::endl
            << " i - insert list of values " << std::endl
            << " p - print entire map or lookup specified values" << std::endl
            << " c - clear map" << std::endl
            << " v - verify map contents" << std::endl
            << " h - this help" << std::endl;
  return true;
}

int main() {
  char buffer[1024];
  bool exit_val = 0;
  while (fgets(buffer, sizeof(buffer), stdin))
    exit_val += !do_command(buffer);
  return exit_val;
}

