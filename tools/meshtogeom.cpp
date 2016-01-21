#include <iostream>
#include <stdlib.h>
#include <vector>
#include <map>
#include <string>
#include <stdio.h>
#include <iomanip>
#include <fstream>
#include "moab/MOABConfig.h"

#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <termios.h>
#include <sys/ioctl.h>
#endif
#include <math.h>
#include <assert.h>
#include <float.h>

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "MBTagConventions.hpp"
#include "moab/Interface.hpp"


#include "moab/VerdictWrapper.hpp"
#include "iGeom.h"

using namespace moab;


static void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage: " << name << " <input_file> " << std::endl;
  stream << "This creates a textfile with material id followed by mesh/geom ratio for each material" << std::endl;
}


int main( int argc, char* argv[] )
{
  if (argc!=2)
    {
      print_usage(argv[0], std::cerr);
      return 1;
    }
  Core *mb = new Core();
  iGeom_Instance geom;
  char * out_file=NULL;
  std::ofstream ofile;

  if (argc==3)
    {
      std::cout<<" write verbose output to a CSV file " << argv[2] << "\n";
      out_file = argv[2];
      ofile.open (out_file);
    }
  if (MB_SUCCESS != mb->load_file(argv[1], 0 ) )
    {
      fprintf(stderr, "Error reading file: %s\n", argv[1] );
      return 1;
    }
  int result;
  std::string engine_opt = ";engine=ACIS";
  iGeom_newGeom( engine_opt.c_str(), &geom, &result, engine_opt.length() );
  if (iBase_SUCCESS != result) {
      printf("ERROR : can not instantiate engine\n");
      return 0;
    }
  iGeom_load(geom, argv[1], 0, &result, strlen(argv[1]), 0);
  if (iBase_SUCCESS != result) {
      printf("ERROR : can not load a geometry\n");
      return 0;
    }

  // now calculate geometric and mesh volume
  std::cout << "Computing MESHTOGEOM ratio" << std::endl;
  VerdictWrapper vw(mb);
  // Get all the materials in this mesh file
  moab::Tag mattag, gidtag;
  moab::Range matsets;
  mb->tag_get_handle("MATERIAL_SET", 1, MB_TYPE_INTEGER, mattag);
  mb->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER, gidtag);
  mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag, 0, 1, matsets );

  double mtot = 0.0, volume = 0.0;
  //loop thru all the materials
  moab::Range::iterator set_it;
  for (set_it = matsets.begin(); set_it != matsets.end(); set_it++)  {
      moab::EntityHandle this_set = *set_it;

      std::vector<moab::EntityHandle> geomsets_for_gid;
      mb->get_entities_by_type(this_set, MBENTITYSET, geomsets_for_gid);

      iBase_TagHandle igeom_gidtag;
      const char *tag1 = "GLOBAL_ID";
      int namelen = strlen(tag1);
      iGeom_getTagHandle(geom, tag1, &igeom_gidtag, &result, namelen);

      for(int volid = 0; volid < (int) geomsets_for_gid.size(); volid++){
          // get the gid of this volume
          int my_geom_id = 0;
          mb->tag_get_data(gidtag, &geomsets_for_gid[volid], 1, &my_geom_id);
          std::cout << " GID = " << my_geom_id << std::endl;
        }


      // get the id for this set
      int material_id;
      mb->tag_get_data(mattag, &this_set, 1, &material_id);

      std::vector<moab::EntityHandle> elems;
      mb->get_entities_by_dimension(this_set, 3, elems, true);
      // get all the elements in this material
      double mvol_temp = 0.0;
      for(unsigned int i=0; i<elems.size();i++){
          mvol_temp = 0.0;
          vw.quality_measure(elems[i], moab::MB_VOLUME, mvol_temp);
          mtot+=mvol_temp;
        }

      double meshtogeom = mtot/volume;
      std::cout << material_id << " has volume " << mtot << std::endl;

      moab::Tag mtog_tag;
      // now set the meshtogeom tag on this set
      mb->tag_get_handle( "MESHTOGEOM", 1, MB_TYPE_DOUBLE,
                          mtog_tag, MB_TAG_SPARSE|MB_TAG_CREAT );
      mb->tag_set_data(mtog_tag, &(*set_it), 1, &meshtogeom);

      elems.clear();
      mtot = 0.0;
    }
  delete mb;
  return 0;
}


