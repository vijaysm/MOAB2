/*
 * PolyExoWriter.cpp
 *
 *  Created on: Aug 15, 2014
 */

#include <netcdfcpp.h>
#include "PolyExoWriter.hpp"


namespace moab {

PolyExoWriter::PolyExoWriter(Interface * moab): mb(moab){
  // TODO Auto-generated constructor stub

}

PolyExoWriter::~PolyExoWriter() {
  // TODO Auto-generated destructor stub
}
ErrorCode PolyExoWriter::write_poly_set(EntityHandle set, const char * filename)
{



  const int ParamFour = 4;
  const int ParamLenString = 33;

  Range polys;
  ErrorCode rval = mb->get_entities_by_dimension(set, 2, polys);
  // Determine the maximum number of nodes per element
  int nNodesPerElement = 0;
  for (Range::iterator pit=polys.begin(); pit!=polys.end(); pit++)
  {
    const EntityHandle * conn=0;
    int num_nodes;
    rval = mb->get_connectivity(*pit, conn, num_nodes);
    if (rval != MB_SUCCESS)
      return rval;
    if (nNodesPerElement < num_nodes)
      nNodesPerElement = num_nodes;
  }
  printf("Max nodes per element: %i", nNodesPerElement);

  // Output to a NetCDF Exodus file
  NcFile ncOut(filename, NcFile::Replace);

  // Random Exodus dimensions
  NcDim * dimLenString = ncOut.add_dim("len_string", ParamLenString);
  NcDim * dimLenLine = ncOut.add_dim("len_line", 81);
  NcDim * dimFour = ncOut.add_dim("four", ParamFour);
  NcDim * dimTime = ncOut.add_dim("time_step");
  NcDim * dimDimension = ncOut.add_dim("num_dim", 3);

  // Number of nodes
  Range verts;
  rval = mb->get_connectivity(polys, verts);
  if (rval != MB_SUCCESS)
    return rval;

  int nNodeCount = (int)verts.size();
  NcDim * dimNodes = ncOut.add_dim("num_nodes", nNodeCount);

  // Number of elements
  int nElementCount = (int)polys.size();
  NcDim * dimElements = ncOut.add_dim("num_elem", nElementCount);

  // Other dimensions
  NcDim * dimNumElementBlocks = ncOut.add_dim("num_el_blk", 1);
  NcDim * dimNumQARec = ncOut.add_dim("num_qa_rec", 1);
  NcDim * dimElementBlock1 = ncOut.add_dim("num_el_in_blk1", nElementCount);
  NcDim * dimNodesPerElement =
    ncOut.add_dim("num_nod_per_el1", nNodesPerElement);
  NcDim * dimAttBlock1 = ncOut.add_dim("num_att_in_blk1", 1);

  // Global attributes
  ncOut.add_att("api_version", 4.98f);
  ncOut.add_att("version", 4.98f);
  ncOut.add_att("floating_point_word_size", 8);
  ncOut.add_att("file_size", 0);

  char szTitle[128];
  sprintf(szTitle, "tempest(%s) 01/01/2013: 00:00:00", filename);
  ncOut.add_att("title", szTitle);

  // Time_whole (unused)
  ncOut.add_var("time_whole", ncDouble, dimTime);

  // QA records
  char szQARecord[ParamFour][ParamLenString] = {
    "Tempest", "13.0", "01/01/2013", "00:00:00"};

  NcVar * varQARecords =
    ncOut.add_var("qa_records", ncChar, dimNumQARec, dimFour, dimLenString);
  varQARecords->set_cur(0, 0, 0);
  varQARecords->put(&(szQARecord[0][0]), 1, 4, ParamLenString);

  // Coordinate names
  char szCoordNames[3][ParamLenString] = {"x", "y", "z"};

  NcVar * varCoordNames =
    ncOut.add_var("coor_names", ncChar, dimDimension, dimLenString);
  varCoordNames->set_cur(0, 0, 0);
  varCoordNames->put(&(szCoordNames[0][0]), 3, ParamLenString);

  // Element block names
  NcVar * varElementBlockNames =
    ncOut.add_var("eb_names", ncChar, dimNumElementBlocks, dimLenString);

  // Element map
  int * nElementMap = new int[nElementCount];
  for (int i = 0; i < nElementCount; i++) {
    nElementMap[i] = i+1;
  }

  NcVar * varElementMap =
    ncOut.add_var("elem_map", ncInt, dimElements);
  varElementMap->put(nElementMap, nElementCount);

  // red parent, blue parent; reuse the nElementMap array
  Tag tagr;
  rval = mb->tag_get_handle("RedParent", tagr);
  if (rval != MB_SUCCESS)
    return rval;
  rval = mb->tag_get_data(tagr, polys, nElementMap);
  NcVar * redParent =
    ncOut.add_var("red_parent", ncInt, dimElements);
  redParent->put(nElementMap, nElementCount);

  Tag tagb;
  rval = mb->tag_get_handle("BlueParent", tagb);
  if (rval != MB_SUCCESS)
    return rval;
  rval = mb->tag_get_data(tagb, polys, nElementMap);
  NcVar * blueParent =
    ncOut.add_var("blue_parent", ncInt, dimElements);
  blueParent->put(nElementMap, nElementCount);

  delete[] nElementMap;

  // Element block status
  int nOne = 1;

  NcVar * varElementBlockStatus =
    ncOut.add_var("eb_status", ncInt, dimNumElementBlocks);
  varElementBlockStatus->put(&nOne, 1);

  NcVar * varElementProperty =
    ncOut.add_var("eb_prop1", ncInt, dimNumElementBlocks);
  varElementProperty->put(&nOne, 1);
  varElementProperty->add_att("name", "ID");

  // Attributes
  double * dAttrib1 = new double[nElementCount];
  for (int i = 0; i < nElementCount; i++) {
    dAttrib1[i] = 1.0;
  }

  NcVar * varAttrib1 =
    ncOut.add_var("attrib1", ncDouble, dimElementBlock1, dimAttBlock1);
  varAttrib1->put(dAttrib1, nElementCount, 1);
  delete[] dAttrib1;

  // Face nodes (1-indexed)
  NcVar * varFaces =
    ncOut.add_var("connect1", ncInt, dimElementBlock1, dimNodesPerElement);

  varFaces->add_att("elem_type", "SHELL4");// this is not true, shell cannot have more than 4 nodes :)

  int * nConnect = new int[nNodesPerElement];
  for (int i = 0; i < nElementCount; i++) {
    EntityHandle polygon=polys[i];
    int num_nodes=0;
    const EntityHandle * conn = 0;
    rval = mb->get_connectivity(polygon, conn, num_nodes);
    if (rval != MB_SUCCESS)
      return rval;
    int k = 0;
    for (; k < num_nodes; k++) {
      nConnect[k] = verts.index(conn[k]) + 1;
    }
    for (; k < nNodesPerElement; k++) {
      nConnect[k] = nConnect[num_nodes-1];
    }

    varFaces->set_cur(i, 0);
    varFaces->put(nConnect, 1, nNodesPerElement);
  }

  delete[] nConnect;

  // Node list
  NcVar * varNodes =
    ncOut.add_var("coord", ncDouble, dimDimension, dimNodes);

  double * dCoord = new double[nNodeCount];

  double * coords = new double[3*nNodeCount];
  rval = mb->get_coords(verts, coords);
  if (rval != MB_SUCCESS)
    return rval;

  for (int i = 0; i < nNodeCount; i++) {
    dCoord[i] =coords[3*i];
  }
  varNodes->set_cur(0, 0);
  varNodes->put(dCoord, 1, nNodeCount);
  for (int i = 0; i < nNodeCount; i++) {
    dCoord[i] = coords[3*i+1];
  }
  varNodes->set_cur(1, 0);
  varNodes->put(dCoord, 1, nNodeCount);
  for (int i = 0; i < nNodeCount; i++) {
    dCoord[i] = coords[3*i+2];
  }
  varNodes->set_cur(2, 0);
  varNodes->put(dCoord, 1, nNodeCount);
  delete[] dCoord;
  delete[] coords;

  // Edge types
  NcVar * varEdgeTypes =
    ncOut.add_var("edge_type", ncInt,
      dimElementBlock1, dimNodesPerElement);

  int * nEdgeType = new int[nNodesPerElement];
  for (int i = 0; i < nElementCount; i++) {
    for (int k = 0; k < nNodesPerElement; k++) {
      nEdgeType[k] = static_cast<int>(0);// everything type 0 now; need to pass that with a tag
    }
    varEdgeTypes->set_cur(i, 0);
    varEdgeTypes->put(nEdgeType, 1, nNodesPerElement);
  }
  delete[] nEdgeType;

  return MB_SUCCESS;
}
} /* namespace moab */
