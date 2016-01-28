/** @example ExtrudePoly.cpp
 * Description: read a 2d mesh in plane, extrude to form  layers of prisms (polyhedra) \n
 *
 * To run: ./ExtrudePoly [meshfile] [outfile] [nlayers] [thickness_layer] \n
 */


#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include <iostream>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

// at most 20 edges per polygon
#define MAXEDGES  20
// Note: change the file name below to test a trivial "No such file or directory" error
string test_file_name = string(MESH_DIR) + string("/io/poly8-10.vtk");
string output = string("polyhedra.vtk");
int layers = 1;
double layer_thick = 1.0;

int main(int argc, char **argv)
{
  // Get MOAB instance
  Interface* mb = new (std::nothrow) Core;
  if (NULL == mb)
    return 1;

  // Need option handling here for input filename
  if (argc > 1) {
    // User has input a mesh file
    test_file_name = argv[1];
  }
  if (argc>2)
    output = argv[2];

  if (argc>3)
    layers = atoi(argv[3]);

  if (argc>4)
      layer_thick = atof(argv[4]);

  std::cout << "Run: " << argv[0] << " " << test_file_name << " " << output<<
      " " << layers << " " << layer_thick << "\n";
  // Load the mesh from vtk file
  ErrorCode rval = mb->load_mesh(test_file_name.c_str());MB_CHK_ERR(rval);

  // Get verts entities, by type
  Range verts;
  rval = mb->get_entities_by_type(0, MBVERTEX, verts);MB_CHK_ERR(rval);

  // Get faces, by dimension, so we stay generic to entity type
  Range faces;
  rval = mb->get_entities_by_dimension(0, 2, faces);MB_CHK_ERR(rval);
  cout << "Number of vertices is " << verts.size() << endl;
  cout << "Number of faces is " << faces.size() << endl;

  Range edges;
  // create all edges
  rval = mb->get_adjacencies(faces, 1, true, edges, Interface::UNION);MB_CHK_ERR(rval);

  cout << "Number of edges is " << edges.size() << endl;

  std::vector<double> coords;
  int nvPerLayer= (int)verts.size();
  coords.resize(3*nvPerLayer);

  rval = mb->get_coords(verts, &coords[0]); MB_CHK_ERR(rval);
  // create first vertices
  Range * newVerts = new Range [layers+1];
  newVerts[0] = verts; // just for convenience
  for (int ii=0; ii<layers; ii++)
  {
    for (int i=0;i<nvPerLayer; i++)
      coords[3*i+2] += layer_thick;

    rval = mb->create_vertices(&coords[0], nvPerLayer, newVerts[ii+1]); MB_CHK_ERR(rval);
  }
  // for each edge, we will create layers quads
  int nquads = edges.size()*layers;
  ReadUtilIface *read_iface;
  rval = mb->query_interface(read_iface);MB_CHK_SET_ERR(rval, "Error in query_interface");

  EntityHandle start_vert, start_elem, *connect;
       // Create quads
  rval = read_iface->get_element_connect(nquads, 4, MBQUAD, 0, start_elem, connect);MB_CHK_SET_ERR(rval, "Error in get_element_connect");
  int nedges = (int)edges.size();

  int indexConn=0;
  for (int j=0; j<nedges; j++)
  {
    EntityHandle edge=edges[j];

    const EntityHandle *conn2 = NULL;
    int num_nodes;
    rval = mb->get_connectivity( edge, conn2, num_nodes) ; MB_CHK_ERR(rval);
    if (2!=num_nodes)
      MB_CHK_ERR(MB_FAILURE);

    int i0= verts.index(conn2[0]);
    int i1= verts.index(conn2[1]);
    for  (int ii=0; ii<layers; ii++)
    {
      connect[indexConn++] = newVerts[ii][i0];
      connect[indexConn++] = newVerts[ii][i1];
      connect[indexConn++] = newVerts[ii+1][i1];
      connect[indexConn++] = newVerts[ii+1][i0];
    }
  }

  int nfaces= (int)faces.size();
  EntityHandle * allPolygons = new EntityHandle [nfaces*(layers+1)];
  for (int i=0; i<nfaces; i++)
    allPolygons[i] = faces[i];

  // vertices are parallel to the base vertices
  int indexVerts[MAXEDGES] = {0}; // polygons with at most MAXEDGES edges
  EntityHandle newConn[MAXEDGES] = {0};

  // edges will be used to determine the lateral faces of polyhedra (prisms)
  int indexEdges[MAXEDGES] = {0}; // index of edges in base polygon
  for (int j=0; j<nfaces; j++)
  {
    EntityHandle polyg=faces[j];

    const EntityHandle *connp = NULL;
    int num_nodes;
    rval = mb->get_connectivity( polyg, connp, num_nodes) ; MB_CHK_ERR(rval);

    for (int i=0; i<num_nodes; i++)
      indexVerts[i] = verts.index(connp[i]);

    for  (int ii=0; ii<layers; ii++)
    {
      // create a polygon on each layer
      for (int i=0; i<num_nodes; i++)
        newConn[i] = newVerts[ii+1] [indexVerts[i]]; // vertices in layer ii+1

      rval = mb->create_element(MBPOLYGON, newConn, num_nodes, allPolygons[nfaces*(ii+1)+j] ); MB_CHK_ERR(rval);

      // now create a polyhedra with top, bottom and lateral swept faces
    }
  }
  rval = mb->write_file(output.c_str());MB_CHK_ERR(rval);

  delete mb;

  return 0;
}
