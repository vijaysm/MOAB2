/*This unit test is for the uniform refinement capability based on AHF datastructures*/
#include <iostream>
#include "moab/Core.hpp"
#include "moab/NestedRefine.hpp"

using namespace moab;

#define STRINGIFY_(X) #X
#define STRINGIFY(X) STRINGIFY_(X)

int main(int argc, char *argv[])
{

  Core mb;
  Interface* mbImpl = &mb;
  ErrorCode error;

  if (argc==1)
    {
      std::cerr << "Usage: " << argv[0] << " [filename]" << std::endl;
      return 1;
    }

  const char *filename = argv[1];
  error = mbImpl->load_file(filename);
  if (error != MB_SUCCESS) return error;

  NestedRefine uref(&mb);

  int level_degrees[4] = {2,3,2,3};
  int num_levels = sizeof(level_degrees) / sizeof(int);
  EntityHandle *set = new EntityHandle[num_levels];

  std::cout<<"Starting hierarchy generation"<<std::endl;

  error = uref.generate_mesh_hierarchy(level_degrees, num_levels, set);
  if (error != MB_SUCCESS) return error;

  std::cout<<"Finished hierarchy generation"<<std::endl;

  std::stringstream file;
  file <<"mesh_hierarchy.vtk";
  std::string str = file.str();
  const char* output_file = str.c_str();
  error = mbImpl->write_file(output_file);
  if (error != MB_SUCCESS) return error;

  delete [] set;

  return 0;
}

