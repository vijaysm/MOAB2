/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/**
 * \file gttool_test.cpp
 *
 * \brief test geometrize method from GeomTopoTool
 *
 */
#include "moab/Core.hpp"
#include <iostream>

#include <assert.h>
#include <string.h>
#include <stdio.h>

#include "TestUtil.hpp"
#include "moab/GeomTopoTool.hpp"

using namespace moab;

std::string filename;
std::string ofile;
std::string ofile2;
std::string ofile3;
bool remove_output_file;
ErrorCode geometrize_test(Interface * mb, EntityHandle inputSet);

ErrorCode create_shell_test(Interface * mb);

ErrorCode duplicate_model_test(Interface * mb);

ErrorCode check_model_test(Interface * mb);

void handle_error_code(ErrorCode rv, int &number_failed, int &number_successful)
{
  if (rv == MB_SUCCESS) {
    std::cout << "Success";
    number_successful++;
  } else {
    std::cout << "Failure";
    number_failed++;
  }
}

int main(int argc, char *argv[])
{
  filename = TestDir + "/partBed.smf";
  ofile = "output.h5m";
  ofile2 = "shell.h5m";
  ofile3 = "shellCopy.h5m";

  remove_output_file = true;
  bool only_check = false;
  bool only_geometrize = false;


  if (argc == 1) {
    std::cout << "Using default input file and output files " << filename << " " << ofile <<
        " " << ofile2 <<  " " << ofile3 << std::endl;
  } else if (argc == 5) {
    filename = argv[1];
    ofile = argv[2];
    ofile2 = argv[3];
    ofile3 = argv[4];
    remove_output_file = false;
  } else if (argc==2) {
    ofile3 = argv[1];// check model only from file
    only_check = true;
    remove_output_file = false; // this is input now, do not delete it
  } else if (argc == 3) {
    filename = argv[1];
    ofile = argv[2];
    only_geometrize = true;
    remove_output_file = false;
  }
    else
  {
    std::cerr << "Usage: " << argv[0] << " [surface_mesh] [mbgeo_file] [shellfile] [copyshellfile] " << std::endl;
    return 1;
  }

  int number_tests_successful = 0;
  int number_tests_failed = 0;

  Core mbcore;
  Interface * mb = &mbcore;

  // use this as a tool to check a model
  if (only_check)
  {
    if (MB_SUCCESS==check_model_test( mb))
      std::cout << ofile3 << " passed gtt check\n";
    return 0;
  }

  mb->load_file(filename.c_str());

  //   FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed

  std::cout << "geometrize test: ";
  ErrorCode rval = geometrize_test( mb, 0); // just pass the root set
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  if (only_geometrize)
  {
    return number_tests_failed;
  }
  std::cout << "create shell test: ";
  rval = create_shell_test( mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "duplicate model test: ";
  rval = duplicate_model_test( mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "check_model test: ";
  rval = check_model_test( mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  return number_tests_failed;
}
ErrorCode geometrize_test(Interface * mb, EntityHandle inputSet)
{
  GeomTopoTool gtt(mb);
  EntityHandle outSet;
  ErrorCode rval=gtt.geometrize_surface_set(inputSet, outSet);MB_CHK_SET_ERR(rval, "Can't geometrize the set\n");

  std::cout<<"writing output file: " << ofile.c_str() << " ";
  rval=mb->write_file(ofile.c_str(), 0, 0, &outSet, 1);MB_CHK_SET_ERR(rval, "Can't write output file\n");
  if (remove_output_file)
  {
    remove(ofile.c_str());
  }
  return MB_SUCCESS;
}

ErrorCode create_shell_test(Interface * mb)
{
  // we should be able to delete mesh and create a model from scratch

  ErrorCode rval = mb->delete_mesh();MB_CHK_SET_ERR(rval, "Can't delete existing mesh\n");

  // create some vertices
  double coords [] = { 0, 0, 0,
                     1, 0, 0.1,
                     2, 0, 0,
                     3, 0, -0.1,
                     0, 1, 0,
                     1, 1, 0,
                     2, 1, 0,
                     3, 1, -0.1,
                     0, 2, 0,
                     1, 2, -0.1,
                     2, 2, -0.1,
                     3, 2, -0.2,
                     0, 0, 1,
                     1, 0, 0.9,
                     2, 0.1, 0.85,
                     3, 0.2, 0.8,
                     0, 0.1, 2,
                     1, 0.1, 2,
                     2.1, 0.2, 2.1,
                     3.1, 0.2, 2.1 };

  int nvert = 20;
  Range verts;
  rval = mb->create_vertices(coords, nvert, verts);MB_CHK_SET_ERR(rval, "Can't create vertices\n");

  EntityHandle connec [] = { 1, 2, 5,
                    5, 2, 6,
                    2, 3, 6,
                    6, 3, 7,
                    3, 4, 7,
                    7, 4, 8,
                    5, 6, 9,
                    9, 6, 10,
                    6, 7, 10,
                    10, 7, 11,
                    7, 8, 11,
                    11, 8, 12,  // first face, 1-12
                    13, 14, 1,
                    1, 14, 2,
                    14, 15, 2,
                    2, 15, 3,
                    15, 16, 3,
                    3, 16, 4,
                    17, 18, 13,
                    13, 18, 14,
                    18, 19, 14,
                    14, 19, 15,
                    19, 20, 15,
                    15, 20, 16 // second face, 13-24
                   };
  EntityHandle elem;
  int nbTri = sizeof(connec)/3/sizeof(EntityHandle);
  int i = 0;
  std::vector<EntityHandle> tris;
  for (i=0; i<nbTri; i++)
  {
    mb->create_element(MBTRI, &connec[3*i], 3, elem);
    tris.push_back(elem);
  }


  // create some edges too
  EntityHandle edges [] = { 1, 2,
                          2, 3,
                          3, 4,  // geo edge 1  1:3
                          4, 8,
                          8, 12,  // geo 2         4:5
                          12, 11,
                          11, 10,
                          10, 9, // geo 3 6:8
                          9, 5,
                          5, 1,  // geo 4  9:10
                          1, 13,
                          13, 17, // geo 5  11:12
                          17, 18,
                          18, 19,
                          19, 20, // geo 6  13:15
                          20, 16,
                          16, 4   // geo 7  16:17
                           };
  int nbEdges = sizeof(edges)/2/sizeof(EntityHandle);
  std::vector<EntityHandle> edgs;
  for (i=0; i<nbEdges; i++)
  {
    mb->create_element(MBEDGE, &edges[2*i], 2, elem);
    edgs.push_back(elem);
  }
  // create some sets, and create some ordered sets for edges
  EntityHandle face1, face2;
  rval = mb->create_meshset(MESHSET_SET, face1);MB_CHK_ERR(rval);
  rval = mb->add_entities(face1, &tris[0], 12);MB_CHK_ERR(rval);

  rval = mb->create_meshset(MESHSET_SET, face2);MB_CHK_ERR(rval);
  // next 12 triangles
  rval = mb->add_entities(face2, &tris[12], 12);MB_CHK_ERR(rval);

  // the orientation and senses need to be set for face edges
  moab::GeomTopoTool gTopoTool(mb, false);

  rval = gTopoTool.add_geo_set(face1, 2);MB_CHK_ERR(rval);

  rval = gTopoTool.add_geo_set(face2, 2);MB_CHK_ERR(rval);

  // create some edges
  EntityHandle edge[7]; //edge[0] has EH 1...
 ;
  for (i=0; i<7; i++)
  {
    rval = mb->create_meshset(MESHSET_ORDERED, edge[i]);MB_CHK_ERR(rval);
    rval = gTopoTool.add_geo_set(edge[i], 1);MB_CHK_ERR(rval);
  }

  // first 3 mesh edges...
  rval = mb->add_entities(edge[0], &edgs[0], 3);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[1], &edgs[3], 2);MB_CHK_ERR(rval); 
  rval = mb->add_entities(edge[2], &edgs[5], 3);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[3], &edgs[8], 2);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[4], &edgs[10], 2);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[5], &edgs[12], 3);MB_CHK_ERR(rval);
  rval = mb->add_entities(edge[6], &edgs[15], 2);MB_CHK_ERR(rval);

  // create some sets for vertices; also need to create some for parent/child relationships
  EntityHandle vertSets[6];// start from 0

  for (i=0; i<6; i++)
  {
    rval = mb->create_meshset(MESHSET_SET, vertSets[i]);MB_CHK_ERR(rval);
    rval = gTopoTool.add_geo_set(vertSets[i], 0);MB_CHK_ERR(rval);
  }

  EntityHandle v(1); // first vertex;
  rval = mb->add_entities(vertSets[0], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (4);
  rval = mb->add_entities(vertSets[1], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (9);
  rval = mb->add_entities(vertSets[2], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (12);
  rval = mb->add_entities(vertSets[3], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (17);
  rval = mb->add_entities(vertSets[4], &v, 1);MB_CHK_ERR(rval);
  v = EntityHandle (20);
  rval = mb->add_entities(vertSets[5], &v, 1);MB_CHK_ERR(rval);

  // need to add parent-child relations between sets
  // edge 1 : 1-2
  rval = mb ->add_parent_child( edge[0], vertSets[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[0], vertSets[1]);MB_CHK_ERR(rval);
  // edge 2 : 2-4
  rval = mb ->add_parent_child( edge[1], vertSets[1]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[1], vertSets[3]);MB_CHK_ERR(rval);

  // edge 3 : 4-3
  rval = mb ->add_parent_child( edge[2], vertSets[3]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[2], vertSets[2]);MB_CHK_ERR(rval);

  // edge 4 : 4-1
  rval = mb ->add_parent_child( edge[3], vertSets[2]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[3], vertSets[0]);MB_CHK_ERR(rval);

  // edge 5 : 1-5
  rval = mb ->add_parent_child( edge[4], vertSets[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[4], vertSets[4]);MB_CHK_ERR(rval);

  // edge 6 : 5-6
  rval = mb ->add_parent_child( edge[5], vertSets[4]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[5], vertSets[5]);MB_CHK_ERR(rval);

  // edge 7 : 6-2
  rval = mb ->add_parent_child( edge[6], vertSets[5]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( edge[6], vertSets[1]);MB_CHK_ERR(rval);

  // face 1: edges 1, 2, 3, 4
  rval = mb ->add_parent_child( face1, edge[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face1, edge[1]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face1, edge[2]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face1, edge[3]);MB_CHK_ERR(rval);

  // face 2: edges 1, 5, 6, 7
  rval = mb ->add_parent_child( face2, edge[0]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face2, edge[4]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face2, edge[5]);MB_CHK_ERR(rval);
  rval = mb ->add_parent_child( face2, edge[6]);MB_CHK_ERR(rval);

  // set senses !!
  std::vector<EntityHandle> faces;
  faces.push_back(face1); // the face1 has all edges oriented positively
  std::vector<int> senses;
  senses.push_back(moab::SENSE_FORWARD); //

  //faces.push_back(face1);
  //faces.push_back(face2);
  gTopoTool.set_senses(edge[1], faces, senses);
  gTopoTool.set_senses(edge[2], faces, senses);
  gTopoTool.set_senses(edge[3], faces, senses);

  faces[0]=face2; // set senses for edges for face2
  gTopoTool.set_senses(edge[4], faces, senses);
  gTopoTool.set_senses(edge[5], faces, senses);
  gTopoTool.set_senses(edge[6], faces, senses);

  // the only complication is edge1 (edge[0]), that has face1 forward and face 2 reverse
  faces[0]=face1;
  faces.push_back(face2);
  senses.push_back(moab::SENSE_REVERSE); // -1 is reverse; face 2 is reverse for edge1 (0)
  // forward == 0, reverse ==1
  gTopoTool.set_senses(edge[0], faces, senses);

  rval = mb->write_mesh(ofile2.c_str());MB_CHK_ERR(rval);

  rval = mb->delete_mesh();MB_CHK_ERR(rval);

  // now test loading it up
  rval = mb->load_file(ofile2.c_str());MB_CHK_ERR(rval);

  if (remove_output_file)
  {
    remove(ofile2.c_str());
  }
  // do some tests on geometry

  // it would be good to have a method on updating the geom topo tool
  // so we do not have to create another one
  moab::GeomTopoTool gTopoTool2(mb, true);// to find the geomsets
  Range ranges[5];
  rval =  gTopoTool2.find_geomsets(ranges);

  assert(MB_SUCCESS==rval);
  assert(ranges[0].size()==6);
  assert(ranges[1].size()==7);
  assert(ranges[2].size()==2);
  assert(ranges[3].size()==0);
  assert(ranges[4].size()==0);

  return MB_SUCCESS;
}

ErrorCode duplicate_model_test(Interface * mb)
{
  moab::GeomTopoTool gTopoTool2(mb, true);// to find the geomsets

  GeomTopoTool * newModel = NULL;
  ErrorCode rval = gTopoTool2.duplicate_model(newModel);
  if (NULL == newModel || rval!=MB_SUCCESS)
    return MB_FAILURE;

  Range ranges[5];
  rval = newModel->find_geomsets(ranges);MB_CHK_ERR(rval);

  assert(ranges[0].size()==6);
  assert(ranges[1].size()==7);
  assert(ranges[2].size()==2);
  assert(ranges[3].size()==0);

  // write the model to a test file
  EntityHandle rootModelSet = newModel->get_root_model_set();
  std::cout<<"writing duplicated model file: " << ofile3.c_str() << " ";
  rval=mb->write_file(ofile3.c_str(), 0, 0, &rootModelSet, 1);MB_CHK_SET_ERR(rval, "Can't write output files\n");

  delete newModel; // we are done with the new geom topo tool
  // do not delete yet the output file, delay after the next test
  /*if (remove_output_file)
  {
    remove(ofile3.c_str());
  }*/

  return MB_SUCCESS;
}

ErrorCode check_model_test(Interface * mb)
{
  ErrorCode rval = mb->delete_mesh();MB_CHK_SET_ERR(rval, "Can't delete existing mesh\n");
  
  rval = mb->load_file(ofile3.c_str());MB_CHK_ERR(rval);

  // do some tests on geometry
  // it would be good to have a method on updating the geom topo tool
  // so we do not have to create another one
  if (remove_output_file)
  {
    remove(ofile3.c_str());
  }
  moab::GeomTopoTool gTopoTool(mb, true);

  if (!gTopoTool.check_model())
    return MB_FAILURE;

  return MB_SUCCESS;
}

