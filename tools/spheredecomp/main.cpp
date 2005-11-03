/*
 * Sphere decomp tool.  Meshes a group of spheres and the interstices between with
 * hex elements, by triangulating the vertices corresponding to sphere centers
 * and subdividing those tets.  For a description of the subdivision template used,
 * see comments in the subdivide_tet function below.
 */

#include "MBCore.hpp"
#include "SphereDecomp.hpp"
#include <iostream>

const char *SPHERE_RADII_TAG_NAME = "SPHERE_RADII";

#define RR if (MB_SUCCESS != result) return result

int main(int argc, char *argv[]) 
{
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <input_mesh> <output_mesh>" << std::endl;
    return 0;
  }
  
    // create MOAB
  MBInterface *mbImpl = new MBCore();
  
    // read in mesh
  MBErrorCode result = mbImpl->load_mesh(argv[1]); RR;

  MBTag sphere_radii_tag = 0;
  double dum_val = 0.1;
  result = mbImpl->tag_create(SPHERE_RADII_TAG_NAME, sizeof(double), 
                              MB_TAG_DENSE, MB_TYPE_DOUBLE, sphere_radii_tag, &dum_val); RR;

  SphereDecomp sd(mbImpl);
  
  result = sd.build_sphere_mesh(SPHERE_RADII_TAG_NAME); RR;
  
    // write mesh
  result = mbImpl->write_mesh(argv[2]); RR;
  
  return 0;
}

