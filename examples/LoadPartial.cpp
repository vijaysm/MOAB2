/** @example LoadPartial.cpp \n
 * \brief Load a part of a file  \n
 * <b>To run</b>: LoadPartial <file> <tag_name> <val1> <val2> ...\n
 *
 * In this example, it is shown how to load only a part of one file; the file must be organized in sets.
 * (cherry-picking only the sets we want)
 * The sets to load are identified by a tag name and the tag values for the sets of interest.
 * This procedure is used  when reading in parallel, as each processor will load only
 * its part of the file, identified either by partition or by material/block sets
 *  by default, this example will load parallel partition sets
 *  with values 1, 2, and 5 from ../MeshFiles/unittest/64bricks_1khex.h5m
 *  The example will always write the output to a file name part.h5m
 */

#include <iostream>
#include <vector>

// Include header for MOAB instance and tag conventions for
#include "moab/Core.hpp" 
#include "MBTagConventions.hpp"

int main(int argc, char **argv) {

    // instantiate & load a file 
    moab::Interface *mb = new moab::Core();

    moab::ErrorCode rval;
    if (argc <= 1) //
    {
      // the default file to load
      int set_tag_values[] = {1, 2, 5};
      int num_set_tag_values = 3;
      // this file is in the mesh files directory
      rval = mb->load_file("../MeshFiles/unittest/64bricks_1khex.h5m",
              0, 0, PARALLEL_PARTITION_TAG_NAME, set_tag_values, num_set_tag_values);
    }
    else
    {
      // first arg is input file, second is tag name, then are the tag values
      if (argc < 4)
      {
        std::cout<< " usage is " << argv[0] << " <file> <tag_name> <value1> <value2>  .. \n";
        return 0;
      }

      else
      {
        std::vector<int> vals(argc-3); // the first 3 args are exe, file, tagname; the rest are values
        for (int i=3; i<argc; i++)
          vals[i-3] = atoi(argv[i]);
        rval = mb->load_file(argv[1], 0, 0, argv[2], &vals[0], (int) vals.size() );
      }
    }
    if (moab::MB_SUCCESS!=rval)
      std::cout << " failed to read\n";
    rval = mb->write_file("part.h5m"); 
    if (moab::MB_SUCCESS!=rval)
      std::cout << " failed to write partial file.\n";
    else
      std::cout << " wrote successfully part.h5m.\n";
    delete mb;
    return 0;
} 
