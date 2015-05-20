/*
 * create a helix mesh;
 *  input file : surface mesh in xz plane, with quads
 *  could be organized in edges/surfaces
 *  HelixSweep -i input.cub -p 10 -a 360
 */
#include "moab/Core.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/ReadUtilIface.hpp"
#include "moab/CartVect.hpp"


#include <time.h>
#include <iostream>
#include <vector>

using namespace moab;
using namespace std;

int main(int argc, char **argv)
{
  double pitch = 10.;
  double angle = 60.;
  int layers = 10;
  ProgOptions opts;

  string inputf =  "input.cub";
  opts.addOpt<string>("input,i", "Specify the input file name string (default input.cub)", &inputf);
  opts.addOpt<double>(string("pitch,p"), string("Total pitch (default=10.)"), &pitch);
  opts.addOpt<double>(string("angle,a"),  string("Total angle (default=60.)"), &angle);
  opts.addOpt<int>(string("layers,l"),  string("num layers (default=10)"), &layers);

  opts.parseCommandLine(argc, argv);

  Core mb;

  ErrorCode rval = mb.load_file(inputf.c_str()); MB_CHK_ERR(rval);

  ReadUtilIface* iface;
  rval = mb.query_interface(iface);MB_CHK_SET_ERR(rval, "Can't get reader interface");

  Range verts;
  rval = mb.get_entities_by_dimension(0, 0, verts); MB_CHK_ERR(rval);

  int num_newNodes = (layers+1) * (int)verts.size();

  vector<double*> arrays;
  EntityHandle startv;
  rval = iface->get_node_coords(3, num_newNodes, 0, startv, arrays);MB_CHK_SET_ERR(rval, "Can't get node coords");

  std::vector<double>  inico(3*verts.size());
  rval = mb.get_coords(verts, &inico[0]); MB_CHK_ERR(rval);

  // for each layer, compute an angle
  int nv = (int)verts.size();
  for (int i=0; i<=layers; i++)
  {
    double arad=angle/180*M_PI/layers*i;
    double deltaz = pitch/layers*angle/360.*i;
    for (int j=0; j<nv; j++)
    {
      double x=inico[3*j];

      double z=inico[3*j+2];
      arrays[0][i*nv+j] = x*cos(arad);
      arrays[1][i*nv+j] = x*sin(arad);
      arrays[2][i*nv+j] = z+deltaz;
    }
  }

  Range quads;
  rval = mb.get_entities_by_type(0, MBQUAD, quads); MB_CHK_ERR(rval);
  int num_hexas = layers * (int)quads.size();

  EntityHandle starte; // Connectivity
  EntityHandle* conn;
  //int num_v_per_elem = 8;

  rval = iface->get_element_connect(num_hexas, 8, MBHEX, 0, starte, conn);MB_CHK_SET_ERR(rval, "Can't get element connectivity");

  int ixh=0; // index in conn array
  for (Range::iterator qit=quads.begin(); qit!=quads.end(); qit++)
  {
    EntityHandle quad=*qit;
    const EntityHandle * conn4 =NULL;
    int nno;
    rval = mb.get_connectivity(quad, conn4, nno); MB_CHK_ERR(rval);
    if (nno!=4 )continue;
    CartVect pos[4];
    rval = mb.get_coords(conn4, 4, &(pos[0][0])); MB_CHK_ERR(rval);
    CartVect normal=(pos[1]-pos[0])*(pos[2]-pos[1]);
    bool positive = (normal % CartVect(0,1,0) > 0) ;

    int indices[4];
    for (int k=0; k<4; k++)
      indices[k] = verts.index(conn4[k]);

    if (positive)
    {
      for (int i=1; i<=layers; i++)
      {
        conn[ixh++]= startv + indices[0]+(i-1)*nv;
        conn[ixh++]= startv + indices[1]+(i-1)*nv;
        conn[ixh++]= startv + indices[2]+(i-1)*nv;
        conn[ixh++]= startv + indices[3]+(i-1)*nv;
        conn[ixh++]= startv + indices[0]+(i)*nv;
        conn[ixh++]= startv + indices[1]+(i)*nv;
        conn[ixh++]= startv + indices[2]+(i)*nv;
        conn[ixh++]= startv + indices[3]+(i)*nv;
      }
    }
    else
    {
      for (int i=1; i<=layers; i++)
      {
        conn[ixh++]= startv + indices[0]+(i-1)*nv;
        conn[ixh++]= startv + indices[3]+(i-1)*nv;
        conn[ixh++]= startv + indices[2]+(i-1)*nv;
        conn[ixh++]= startv + indices[1]+(i-1)*nv;
        conn[ixh++]= startv + indices[0]+(i)*nv;
        conn[ixh++]= startv + indices[3]+(i)*nv;
        conn[ixh++]= startv + indices[2]+(i)*nv;
        conn[ixh++]= startv + indices[1]+(i)*nv;
      }
    }

  }

  rval = mb.write_file("out.h5m");MB_CHK_ERR(rval);
  return 0;
}
