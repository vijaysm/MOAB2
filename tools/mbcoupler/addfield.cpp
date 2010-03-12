#include <iostream>
#include "moab/Core.hpp"
#include <math.h>
#include <cstdlib>
using namespace std;


namespace moab {

//make a nice picture. tweak here
//217: x,y -> [-8, 8] /// z -> [-65, 65]
double physField(double x, double y, double z, double factor){
  double out;

    // 1/r^2 decay from {0,0,0}
    // tjt - changing to 1/r
  double scale = 1.0/factor;
  out = fabs(x*scale) + fabs(y*scale) + fabs(z*scale);
  out += 1e-1; // clamp
  out = 1/out;

  return out;
}


//get some sort of element center
void getHexPos(Interface *mbi, EntityHandle *hex, double &x, double &y, double &z){
  std::vector<EntityHandle> connect;
 
  mbi->get_connectivity(hex, 1, connect);
  double pos[3]={0,0,0};

  int numVerts = connect.size();
  for(int i=0; i<numVerts; i++){
    double tempPos[3];
    EntityHandle vert(connect[i]);

    mbi->get_coords(&vert, 1, tempPos);

    for(int j=0; j<3; j++){
      pos[j] += tempPos[j]/numVerts;
    }
  }
 
  x = pos[0];
  y = pos[1];
  z = pos[2];
}


void putElementField(Interface *mbi, char *tagname, double factor){
  Range hexes;

  mbi->get_entities_by_type(0, MBHEX, hexes);

  const double defVal = 0.;
  Tag fieldTag;
  mbi->tag_create(tagname, sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, fieldTag, &defVal);
 
  int numHexes = hexes.size();

  for(int i=0; i<numHexes; i++){
      //cout << hexes[i] << endl;
    EntityHandle hex = hexes[i];

    double x,y,z;
    getHexPos(mbi, &hex, x,y,z);

    double fieldValue =  physField(x,y,z, factor);

    mbi->tag_set_data(fieldTag, &hex, 1, &fieldValue);
  }

}


void putVertexField(Interface *mbi, char *tagname, double factor){
  Range verts;

  mbi->get_entities_by_type(0, MBVERTEX, verts);

  const double defVal = 0.;
  Tag fieldTag;
  mbi->tag_create(tagname, sizeof(double), MB_TAG_DENSE, MB_TYPE_DOUBLE, fieldTag, &defVal);
 
  int numVerts = verts.size();
  for(int i=0; i<numVerts; i++){
    EntityHandle vert = verts[i]; //?

    double vertPos[3];
    mbi->get_coords(&vert, 1, vertPos);

    double fieldValue =  physField(vertPos[0], vertPos[1], vertPos[2], factor);

    mbi->tag_set_data(fieldTag, &vert, 1, &fieldValue);
  }



}


//Using moab instead of imesh for dev speed
int main(int argc, char **argv){
  Interface *mbi = new Core();

  if (argc < 3){
    cout << "Usage: " << argv[0] << " <infile> <outfile> [factor]\n"
         << "Writes both vertex and element fields.\n";
    return 0;
  }

  mbi->load_mesh(argv[1]);

  double factor = 1.0;
  if (argc == 4) factor = atof(argv[3]);
  
  putVertexField(mbi, "vertex_field", factor);
  putElementField(mbi, "element_field", factor);

  ErrorCode result = mbi->write_mesh(argv[2]);
  if (MB_SUCCESS == result) cout << "wrote " << argv[2] << endl;
  else cout << "Failed to write " << argv[2] << endl;

    //  vector<double> coords;                                                                                                                               
    //  mbi->get_vertex_coordinates(coords);                                                                   
    //  double xavg = 0;                                                                                                                                           
    //  for (int i = 0; i < coords.size()/3; i++) xavg += coords[i];
    //  cout << xavg << endl;

  return 1;
}

} // namespace moab

