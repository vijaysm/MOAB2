/*
 * mpas_mesh.cpp
 *
 *  Created on: Sep 29
 */

// simple mpas reader that should work in serial; read only the nodes and polygons
//
//
#include <iostream>
#include <sstream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>

#include "moab/Core.hpp"
#include "moab/Interface.hpp"
#include "moab/ReadUtilIface.hpp"
#include "MBTagConventions.hpp"
#include "CslamUtils.hpp"

#include "mpi.h"
#include "pnetcdf.h"

#define NCFUNC(func) ncmpi_ ## func

//! Collective I/O mode get
#define NCFUNCAG(func) ncmpi_get ## func ## _all

//! Independent I/O mode get
#define NCFUNCG(func) ncmpi_get ## func

//! Nonblocking get (request aggregation), used so far only for ucd mesh
#define NCFUNCREQG(func) ncmpi_iget ## func

#define NCDF_SIZE MPI_Offset
#define NCDF_DIFF MPI_Offset


#include <algorithm>
#include <string>
#include <assert.h>

#include <cmath>


using namespace moab;

#define INS_ID(stringvar, prefix, id) \
          sprintf(stringvar, prefix, id)

#define GET_DIM(ncdim, name, val)\
    {                            \
    int gdfail = NCFUNC(inq_dimid)(ncFile, name, &ncdim);          \
    if (NC_NOERR == gdfail) {                                             \
      MPI_Offset tmp_val;                                                   \
      gdfail = NCFUNC(inq_dimlen)(ncFile, ncdim, &tmp_val);                        \
      if (NC_NOERR != gdfail) {                                           \
        std::cout<<"couldn't get dimension length for" << name<< " \n"; \
        return 1;                                              \
      }                                                                 \
      else val = tmp_val;                                               \
    } else val = 0;}

#define GET_DIMB(ncdim, name, varname, id, val) \
          INS_ID(name, varname, id); \
          GET_DIM(ncdim, name, val);

#define GET_VAR(name, id, dims) \
    {                           \
    id = -1;\
    int gvfail = NCFUNC(inq_varid)( ncFile, name, &id);   \
    if (NC_NOERR == gvfail) {       \
    int ndims;\
    gvfail = NCFUNC(inq_varndims)( ncFile, id, &ndims);\
    if (NC_NOERR == gvfail) {\
    dims.resize(ndims);    \
    gvfail = NCFUNC(inq_vardimid)( ncFile, id, &dims[0]);}}}

#define GET_1D_INT_VAR(name, id, vals) \
    {GET_VAR(name, id, vals);  \
  if (-1 != id) {\
    NCDF_SIZE ntmp;\
    int ivfail = NCFUNC(inq_dimlen)(ncFile, vals[0], &ntmp);\
    vals.resize(ntmp);\
    NCDF_SIZE ntmp1 = 0;                                                           \
    ivfail = NCFUNC(get_vara_int_all)(ncFile, id, &ntmp1, &ntmp, &vals[0]);\
    if (NC_NOERR != ivfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable " <<name <<"\n";\
      return 1;}}}

#define GET_1D_DBL_VAR(name, id, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    std::cout << " var:" << name << " id: " << id << "\n"; \
    NCDF_SIZE ntmp;\
    int dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp);\
    vals.resize(ntmp);\
    NCDF_SIZE ntmp1 = 0;               \
    dvfail = NCFUNC(get_vara_double_all)( ncFile, id, &ntmp1, &ntmp, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

#define FILL_1D_DBL_VAR(name, id, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    MPI_Offset ntmp;\
    int dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp);\
    std::cout <<" var:" << name << " ntmp:" << ntmp << "\n"; \
    MPI_Offset ntmp1 = 0;                                                           \
    dvfail = NCFUNC(get_vara_double_all)( ncFile, id, &ntmp1, &ntmp, vals);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"FILL_1D_DBL_VAR:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

/*
int get_2d_flt_var(int ncFile, const char * name, int index, std::vector<float> & vals)
{
  int id;
  std::vector<int> dum_dims;
  GET_VAR(name, id, dum_dims);
  if (-1 != id) {
    size_t ntmp;
  int dvfail = nc_inq_dimlen(ncFile, dum_dims[1], &ntmp);
  vals.resize(ntmp);
  size_t ntmp1[2] = {index, 0};

  int ntimes;
  dvfail = nc_inq_dimlen(ncFile, dum_dims[0], &ntimes);
  if (index>=ntimes) {
    std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";
    return 1;
  }
  size_t count[2] ={1, ntmp};
  dvfail = nc_get_vara_float(ncFile, id, ntmp1, count, &vals[0]);
  if (NC_NOERR != dvfail) {
    std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";
    return 1;
  }
  return 0;
}
*/
/* get the variable along an index */
#define GET_2D_FLT_VAR(name, id, index, vals) \
    {std::vector<int> dum_dims;        \
  GET_VAR(name, id, dum_dims);\
  if (-1 != id) {\
    MPI_Offset ntmp, ntmp2;\
    int dvfail = NCFUNC(inq_dimlen)( ncFile, dum_dims[1], &ntmp);\
    dvfail = NCFUNC(inq_dimlen)(ncFile, dum_dims[0], &ntmp2);\
    vals.resize(ntmp);\
    MPI_Offset ntmp1[2] = {index, 0}; \
    if (index>=ntmp2) { \
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n"; \
      return 1; \
    } \
    MPI_Offset count[2] ={1, ntmp}; \
    dvfail = NCFUNC(get_vara_float)( ncFile, id, ntmp1, count, &vals[0]);\
    if (NC_NOERR != dvfail) {\
      std::cout<<"ReadNCDF:: Problem getting variable "<< name<<"\n";\
      return 1;}}}

int main(int argc, char ** argv)
{

  MPI_Init(&argc, &argv);
  int num_el = 50000;
  const char *data_file = "/home/iulian/source/myMoabRepo/MeshFiles/unittest/io/mpasx1.642.t.2.nc";
  if (argc>1)
  {
    data_file = argv[1];
  }
  std::cout << "num_el=" << num_el << "\n";

  Core moab;
  Interface & mb = moab;
  ErrorCode rval;

  int ncFile, temp_dim;        // netcdf/exodus file

  int fail = NCFUNC(open)(MPI_COMM_WORLD, data_file, 0, MPI_INFO_NULL, &ncFile);
  if (NC_NOWRITE != fail) {
    std::cout<<"ReadNCDF:: problem opening Netcdf/Exodus II "<<data_file <<"\n";
    return 1;
  }
  int maxEdges;
  GET_DIM(temp_dim, "maxEdges", maxEdges);
  std::cout << "maxEdges:" << maxEdges << "\n";

  int nVertices;
  GET_DIM(temp_dim, "nVertices", nVertices);
    std::cout << "nVertices:" << nVertices << "\n";

  int nCells;
  GET_DIM(temp_dim, "nCells", nCells);
  std::cout << "nCells:" << nCells << "\n";

  // see if the mesh from metadata makes sense
  // create polygons with the connectivity from conn array

  ReadUtilIface* readMeshIface;
  mb.query_interface( readMeshIface );
  if (NULL==readMeshIface)
    return 1;
  EntityHandle node_handle = 0;
  std::vector<double*> arrays;
  rval = readMeshIface->get_node_coords(3, nVertices,
      1, node_handle, arrays);
  if (MB_SUCCESS!= rval)
    return 1;

  std::vector<double> cor;

  FILL_1D_DBL_VAR("xVertex", temp_dim, &(arrays[0][0]));
  FILL_1D_DBL_VAR("yVertex", temp_dim, &(arrays[1][0]));
  FILL_1D_DBL_VAR("zVertex", temp_dim, &(arrays[2][0]));

  Range nodes(node_handle, node_handle+nVertices-1);
  std::vector<int> gids(nVertices);
  for (int j=0; j<nVertices; j++)
    gids[j]=j+1;
  Tag gid;
  int dum_val = 0;
  rval = mb.tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid, MB_TAG_DENSE | MB_TAG_CREAT, &dum_val);
  if (MB_SUCCESS != rval)
    return 1;

  rval = mb.tag_set_data(gid, nodes, &gids[0]);
  if (MB_SUCCESS!= rval)
       return 1;

  EntityHandle newSet;
  rval = mb.create_meshset(MESHSET_SET, newSet);
  if (MB_SUCCESS != rval)
    return 1;

  // so the nodes will be part
  mb.add_entities(newSet, nodes);


  std::vector<int>  nEdgesOnCell;
  GET_1D_INT_VAR("nEdgesOnCell", temp_dim, nEdgesOnCell);

  std::cout << " size of : nEdgesOnCell: " << nEdgesOnCell.size() << "\n";
  std::vector<int> numCellsWithNEdges(maxEdges, 0);
  int actualMaxEdges = 0;
  for (int j=0; j<nCells; j++)
  {
    int nEdgesThisCell = nEdgesOnCell[j];
    numCellsWithNEdges[nEdgesThisCell-1]++;
    if (actualMaxEdges<nEdgesThisCell)
      actualMaxEdges = nEdgesThisCell;
  }
  std::vector<EntityHandle> starthp(maxEdges, 0);
  std::vector<int> indexCellsWithNEdges(actualMaxEdges, 0);
  std::vector<EntityHandle*> conn_arr_local_cells_with_n_edges(actualMaxEdges);

  for (int j=2; j<maxEdges; j++)
  {
    std::cout << "cells with " << j+1 << " edges:" << numCellsWithNEdges[j] << "\n";
    if (numCellsWithNEdges[j]>0)
    {
      readMeshIface->get_element_connect(numCellsWithNEdges[j], j+1, MBPOLYGON, 0,
          starthp[j], conn_arr_local_cells_with_n_edges[j]);
    }
  }

  std::cout << " actualMaxEdges: " << actualMaxEdges << "\n";


  std::vector<int> vertices_on_cells(nCells * actualMaxEdges);
  NCDF_SIZE starts[2]={0, 0};
  NCDF_SIZE counts[2]={nCells, actualMaxEdges};
  std::vector<int> dum_dims;
  int verticesOnCellVarId;
  GET_VAR("verticesOnCell", verticesOnCellVarId, dum_dims);

  int nf = NCFUNCAG(_vara_int)(ncFile, verticesOnCellVarId, starts, counts, &vertices_on_cells[0]);
  if (nf!=0)
    return 1;

  for (int j = 1; j <= nCells; j++) {
    int nEdgesThisCell = nEdgesOnCell[j-1]; //
    int indxFill=indexCellsWithNEdges[nEdgesThisCell-1]*nEdgesThisCell;
    for (int k =0; k<nEdgesThisCell; k++)
    {
      conn_arr_local_cells_with_n_edges[nEdgesThisCell-1][indxFill+k]=
          (EntityHandle)vertices_on_cells[(j-1)*actualMaxEdges+k];
    }

    indexCellsWithNEdges[nEdgesThisCell-1]++;
    EntityHandle polygon=starthp[nEdgesThisCell-1]+indexCellsWithNEdges[nEdgesThisCell-1];
    mb.tag_set_data(gid, &polygon, 1, &j);
  }

  /*// create element tags with lat and lon, for centers
  Tag latTag=0;
  rval = mb.tag_get_handle("latCell", 1, MB_TYPE_DOUBLE, latTag, MB_TAG_DENSE | MB_TAG_CREAT);
  if (MB_SUCCESS != rval)
    return 1;
  Tag lonTag=0;
  rval = mb.tag_get_handle("lonCell", 1, MB_TYPE_DOUBLE, lonTag, MB_TAG_DENSE | MB_TAG_CREAT);
  if (MB_SUCCESS != rval)
    return 1;
  // do a deep copy or not?? for now it is OK :(
  std::vector<double> lonc;
  GET_1D_DBL_VAR("latCell", temp_dim, lonc);
  std::vector<double> latc;
  GET_1D_DBL_VAR("lonCell", temp_dim, latc);*/

  // now, we know that vertices are the same as handles, just set them dry :)


/*  // build a kd tree with the vertices
  EntityHandle tree_root;
  AdaptiveKDTree kd(&mb, true);
  rval = kd.build_tree(nodes, tree_root);
  if (MB_SUCCESS != rval)
    return 1;

  unsigned int  min_depth, max_depth;
  rval = kd.depth(tree_root, min_depth, max_depth);
  if (MB_SUCCESS != rval)
     return 1;
  std::cout << "min_depth, max_depth " << min_depth << " " << max_depth << "\n";
  // now, read the conn file created using spectral visu, and see how they fit
  // this can be directly imported to moab
  const char *myconn_file = "spec.R2.vtk";
  EntityHandle euler_set;
  rval = mb.create_meshset(MESHSET_SET, euler_set);
  if (MB_SUCCESS != rval)
    return 1;

  rval = mb.load_file(myconn_file, &euler_set);

  if (MB_SUCCESS != rval)
    return 1;

  mb.list_entities(&euler_set, 1);

  Range specQuads;
  rval = mb.get_entities_by_dimension(euler_set, 2, specQuads );
  if (MB_SUCCESS != rval)
     return 1;

  Range vertices;
  rval = mb.get_connectivity(specQuads, vertices);
  if (MB_SUCCESS != rval)
     return 1;

  // do a mapping, from position of vertices to the vertices in the kd tree.
  // find the closest vertex to each of this
  std::vector<EntityHandle> mappedTo(vertices.size());
  std::vector<double> mycoords(vertices.size()*3);
  rval = mb.get_coords(vertices, &mycoords[0]);
  double * ptr = &mycoords[0];
  size_t num_verts=vertices.size();
  for (size_t i=0; i<num_verts; i++, ptr+=3)
  {
    CartVect pos(ptr); // take 3 coordinates
    std::vector<EntityHandle> leaves;
    rval = kd.leaves_within_distance( tree_root,
                                          ptr,
                                          0.001,
                                        leaves);
    if (MB_SUCCESS != rval)
      return 1;
    Range closeVerts;
    for (std::vector<EntityHandle>::iterator vit = leaves.begin(); vit != leaves.end(); vit++)
    {
      rval= mb.get_entities_by_handle(*vit, closeVerts, Interface::UNION);
      if (moab::MB_SUCCESS != rval)
        return 1;
    }
    if (closeVerts.size()<1)
    {
      std::cout << "increase tolerance, no points close to " << pos << "\n";
      return 1;
    }
    std::vector<CartVect> coordsTarget(closeVerts.size());
    rval = mb.get_coords(closeVerts, &(coordsTarget[0][0]));
    if (MB_SUCCESS != rval)
      return 1;
    double minDist2=(pos-coordsTarget[0]).length_squared();
    EntityHandle closestVertex=closeVerts[0];
    for (unsigned int j=1; j<closeVerts.size(); j++)
    {
      double dist2=(pos-coordsTarget[j]).length_squared();
      if (minDist2>dist2)
      {
        closestVertex = closeVerts[j];
        minDist2=dist2;
      }
    }
    if (minDist2 > 0.00001)
    {
      std::cout << " problem with node " << vertices[i] << "  min dist2:" << minDist2 << "\n";
      return 1;
    }
    mappedTo [i] = closestVertex;
  }

  num_el = (int)specQuads.size();
  // tmp_ptr is of type int* and points at the same place as conn
  EntityHandle * conn = 0;

  EntityHandle elh;

  readMeshIface->get_element_connect(
          num_el,
          4,
          MBQUAD,
          1,
          elh,
          conn);

  EntityHandle firstVertMyMesh=vertices[0];
  for (int k=0; k<num_el; k++)
  {
    const EntityHandle * myconn=0;
    EntityHandle specElem = specQuads[k];
    int num_nodes=0;
    rval = mb.get_connectivity(specElem, myconn, num_nodes);
    if (MB_SUCCESS != rval || num_nodes !=4)
      return 1;

    int start_el = k*4;
    for (int j=0; j<4; j++)
       conn[start_el+j] = mappedTo[myconn[j]-firstVertMyMesh];
  }
  std::cout << " conn:" << conn[0] << " " << conn[1] << " " << conn[3]<< "\n";
  Range erange(elh, elh+num_el-1);

  mb.add_entities(newSet, erange);
  std::vector<int> gidels(num_el);
  Tag gid2;
  rval = mb.tag_get_handle("GLOBAL_ID_EL", 1, MB_TYPE_INTEGER,
        gid2, MB_TAG_SPARSE|MB_TAG_CREAT);

  if (MB_SUCCESS != rval)
      return 1;
  for (int k=0; k<num_el; k++)
    gidels[k]=k+1;
  mb.tag_set_data(gid2, erange, &gidels[0]);

  int times;
  GET_DIM(temp_dim, "time", times);
  Tag velotag;
  rval = mb.tag_get_handle("VELO", 3, MB_TYPE_DOUBLE,
            velotag, MB_TAG_DENSE|MB_TAG_CREAT);
  if (MB_SUCCESS!= rval)
    return 1;*/

 /* for (size_t tt=0; tt<56 (size_t)times; tt++)
  {
    // now, read velocities from file:
    // read the U850 and V850 variables
    std::vector<float> u850;
    GET_2D_FLT_VAR("U850", temp_dim, tt, u850);
    std::vector<float> v850;
    GET_2D_FLT_VAR("V850", temp_dim, tt, v850);

    std::cout << " U850:" << u850[0] << " " << u850[1] << " " << u850[5] << " "<< u850.size()<<"\n";
    std::cout << " V850:" << v850[0] << " " << v850[1] << " " << v850[5] << " "<< u850.size()<<"\n";
    // ok, use radius as 6371km; not needed

    std::vector<CartVect> velo850(ncol);

    std::stringstream fileName;
    fileName << "VELO0" <<  tt << ".h5m";
    std::cout << " read velocities at 850 for time:" << tt << "\n";


    for (int k=0; k<ncol; k++)
    {
      double latRad=lat[k]*conversion_factor;
      double lonRad=lon[k]*conversion_factor;
      CartVect U(-sin(lonRad), cos(lonRad), 0.);
      CartVect V(-sin(latRad)*cos(lonRad), -sin(latRad)*cos(lonRad), cos(latRad));
      velo850[k]=U*u850[k] +V*v850[k];
    }
    rval = mb.tag_set_data(velotag, nodes, &(velo850[0][0]));
    if (MB_SUCCESS!= rval)
      return 1;
    rval = mb.write_mesh(fileName.str().c_str(), &newSet, 1);
    if (MB_SUCCESS!= rval)
      return 1;
  }*/

  rval = mb.write_mesh("out.h5m");
  if (rval!=MB_SUCCESS)
    std::cout << " can't write file\n";
  std::cout << " wrote out.h5m\n";
  MPI_Finalize();

  return 0;
}
