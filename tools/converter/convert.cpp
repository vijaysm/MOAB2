#include "MBCore.hpp"
#include <iostream>


int main(int argc, char* argv[])
{
  MBInterface* gMB;
  MBErrorCode result;

  if (argc != 3)
  {
    std::cerr << "Usauge: " << argv[0] << " <input_file> <output_file>" << std::endl;
    return 1;
  }

  gMB = new MBCore();

  result = gMB->load_mesh( argv[1] );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to load \"" << argv[1] << "\"." << std::endl; 
    return 2;
  }
  
  result = gMB->write_mesh( argv[2] );
  if (MB_SUCCESS != result)
  { 
    std::cerr << "Failed to write \"" << argv[2] << "\"." << std::endl; 
    return 3;
  }
  
  return 0;
}

