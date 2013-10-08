/** @example LoadPartial.cpp
 * This example demonstrates how to load partially an h5m file
 * Sets in MOAB contain entities and have a tag on them to give meaning to the entities.
 * Tag name: PARALLEL_PARTITION_TAG_NAME "PARALLEL_PARTITION" is used for
 *   parallel partitions
 *
 * <b>Steps in this example </b>:
 *    -# Instantiate MOAB
 *    -# Get input mesh file name, already partitioned in 16 parts
 *    -# Load some parts; save
 *    -# Destroy the MOAB instance
 *
 *
 * <b> To compile: </b>
 *    make LoadPartial MOAB_DIR=<installdir> \n
 *
 * <b> To run: </b>
 *    -# LoadPartial <mesh-file> <int> \n
 *    -# LoadPartial (This uses the default <mesh-file>: <MOAB_SRC_DIR>/MeshFiles/unittest/mbtest2.g)
 */
#include <iostream>
#include <string>

// Include header for MOAB instance and range
#include "moab/Core.hpp"
#include "MBTagConventions.hpp"

int main(int argc, char **argv) {

  // instantiate & load a file
  moab::Interface *mb = new moab::Core();

  moab::ErrorCode rval;

  std::string file_name = "mpas_p16.h5m";
  int one = 1;

  if (argc>1)
    one = atoi(argv[1]);
  if (argc>2)
    file_name = argv[2];
  rval = mb->load_file(file_name.c_str(),
                       0, "DEBUG_IO=9;", PARALLEL_PARTITION_TAG_NAME, &one, 1);

  if (moab::MB_SUCCESS!=rval)
    std::cout << " failed to read\n";

  rval = mb->write_file("part.h5m");
  delete mb;
}
