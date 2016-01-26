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
#include "iMesh.h"
#include "iRel.h"
#include "MBiMesh.hpp"


using namespace moab;


static void print_usage( const char* name, std::ostream& stream )
{
  stream << "Usage (Cubit file): " << name << " <input .cub file> <output mesh file with meshtogeom tag on material sets> <output text meshtogeom file>" << std::endl;
  stream << "\nUsage (OCC file): " << name << " <input OCC geometry file> <input mesh file> <output mesh file with meshtogeom tag on material sets> <output text meshtogeom file>" << std::endl;
}


int main( int argc, char* argv[] )
{
  int result;
  ErrorCode error = MB_SUCCESS;
  moab::Core *mb = new moab::Core();
  MBiMesh *mesh = new MBiMesh(mb);
  iMesh_Instance imesh = reinterpret_cast<iMesh_Instance>(mesh);
  iGeom_Instance geom;
  char * out_file=NULL;
  std::ofstream ofile;

  // Select engine based on CGM build with option
  std::string engine_opt;
#ifdef HAVE_OCC
  if (argc !=5){
    print_usage(argv[0], std::cerr);
    exit(0);
  }
  engine_opt = ";engine=OCC";
  out_file = argv[4];
  ofile.open (out_file);
  error = mb->load_file(argv[2], 0 ); MB_CHK_ERR(error);
#else
  if (argc !=4){
    print_usage(argv[0], std::cerr);
     exit(0);
  }
  engine_opt = ";engine=ACIS";
  out_file = argv[3];
  ofile.open (out_file);
  error = mb->load_file(argv[1], 0 ); MB_CHK_ERR(error);
#endif
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
  error = mb->tag_get_handle("MATERIAL_SET", 1, MB_TYPE_INTEGER, mattag); MB_CHK_ERR(error);
  error = mb->tag_get_handle("GLOBAL_ID", 1, MB_TYPE_INTEGER, gidtag); MB_CHK_ERR(error);
  error = mb->get_entities_by_type_and_tag( 0, MBENTITYSET, &mattag, 0, 1, matsets );MB_CHK_ERR(error);

  double mtot = 0.0, volume = 0.0;

  iRel_Instance assoc;
  iRel_PairHandle pair;

  /* initialize the Associate */
  iRel_create(0, &assoc, &result, 0);
  if(result!=0){std::cerr << "iRel failure " << result << std::endl; exit(1);}

  iRel_createPair(assoc, geom, iRel_ENTITY, iRel_IGEOM_IFACE, iRel_ACTIVE,
                  imesh, iRel_SET,  iRel_IMESH_IFACE, iRel_ACTIVE,
                  &pair, &result);
  if(result!=0){std::cerr << "iRel failure " << result << std::endl; exit(1);}

  iRel_inferAllRelations(assoc,pair, &result);
  if(result!=0){std::cerr << "iRel failure " << result << std::endl; exit(1);}

  //loop thru all the materials
  moab::Range::iterator set_it;
  for (set_it = matsets.begin(); set_it != matsets.end(); set_it++)  {
      moab::EntityHandle this_set = *set_it;

      std::vector<moab::EntityHandle> geomsets_for_gid;
      error = mb->get_entities_by_type(this_set, MBENTITYSET, geomsets_for_gid); MB_CHK_ERR(error);

      iBase_TagHandle igeom_gidtag;
      const char *tag1 = "GLOBAL_ID";
      int namelen = strlen(tag1);
      iGeom_getTagHandle(geom, tag1, &igeom_gidtag, &result, namelen);
      if(result!=0){std::cerr << "iGeom failure " << result << std::endl; exit(1);}

      for(unsigned int volid = 0; volid < geomsets_for_gid.size(); volid++){
          // get the gid of this volume
          int my_geom_id = 0;
          error = mb->tag_get_data(gidtag, &geomsets_for_gid[volid], 1, &my_geom_id); MB_CHK_ERR(error);

          iBase_EntityHandle ent2=NULL;
          iRel_getSetEntRelation(assoc, pair, (iBase_EntitySetHandle) geomsets_for_gid[volid], 1, &ent2, &result);
          if(result!=0){std::cerr << "iRel failure " << result << std::endl; exit(1);}

          double *myvol;
          int myvol_alloc, myvol_size;
          if(ent2 != NULL){
              iGeom_measure(geom, &ent2, 1, &myvol, &myvol_alloc, &myvol_size, &result);
              if(result!=0){std::cerr << "iGeom failure " << result << std::endl; exit(1);}
              volume+=*myvol;
            }
          else {
              std::cout << "Null entity resulting from iRel relation extraction." << std::endl;
            }
        }

      // get the id for this set
      int material_id;
      error = mb->tag_get_data(mattag, &this_set, 1, &material_id); MB_CHK_ERR(error);

      std::vector<moab::EntityHandle> elems;
      error = mb->get_entities_by_dimension(this_set, 3, elems, true); MB_CHK_ERR(error);
      // get all the elements in this material
      double mvol_temp = 0.0;
      for(unsigned int i=0; i<elems.size();i++){
          mvol_temp = 0.0;
          vw.quality_measure(elems[i], moab::MB_VOLUME, mvol_temp);
          mtot+=mvol_temp;
        }

      double meshtogeom = mtot/volume;
      std::cout << material_id << " has mesh volume " << mtot << " geometric volume " << volume << " ratio is " << meshtogeom << std::endl;
      ofile << "Material Id: "<<  material_id << " has mesh volume " << mtot << " geometric volume " << volume << " MESHTOGEOM ratio is " << meshtogeom << std::endl;

      moab::Tag mtog_tag;
      // now set the meshtogeom tag on this set
      error = mb->tag_get_handle( "MESHTOGEOM", 1, MB_TYPE_DOUBLE,
                                  mtog_tag, MB_TAG_SPARSE|MB_TAG_CREAT ); MB_CHK_ERR(error);
      error = mb->tag_set_data(mtog_tag, &(*set_it), 1, &meshtogeom); MB_CHK_ERR(error);

      elems.clear();
      mtot = 0.0;
      volume = 0.0;
    }
#ifdef HAVE_OCC
  error = mb->write_mesh(argv[3]); MB_CHK_ERR(error);
  std::cout << "wrote mesh file: " << argv[3] << " and meshtogeomfile: " << argv[4] << std::endl;
#else
  error = mb->write_mesh(argv[2]); MB_CHK_ERR(error);
  std::cout << "wrote mesh file: " << argv[2] << " and meshtogeomfile: " << argv[3] << std::endl;
#endif

  delete mb;
  ofile.close();
  return 0;
}


