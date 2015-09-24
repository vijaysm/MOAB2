/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/**
 * \file mbground_test.cpp
 *
 * \brief test mbground, a test for the FBEngine class, which is providing iGeom like methods to a MOAB db
 *
 */
#include "moab/Core.hpp"
#include <iostream>
#include <fstream>
#include <set>
#include <algorithm>
#include <vector>
#include <iterator>
#include <algorithm>
#include <iomanip>
#include <assert.h>
#include <string.h>
#include <math.h>
#include "moab/FBEngine.hpp"
#include "moab/GeomTopoTool.hpp"
#include "TestUtil.hpp"

std::string filename;
std::string filename_out;
std::string polyline_file_name;
double min_dot = 0.8;
bool keep_output;
int number_tests_successful = 0;
int number_tests_failed = 0;

using namespace moab;

ErrorCode split_test_across();
ErrorCode verify_split();

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
  filename = TestDir + "/PB.h5m";
  polyline_file_name = TestDir + "/polyline.txt";
  filename_out = "PB_ground.h5m";
  min_dot = 0.8;

  keep_output = false;

  if (argc == 1) {
    std::cout << "Using default input " << filename << " " << polyline_file_name     <<
        " " << min_dot<< " " << filename_out << std::endl;
    std::cout << "    default output file: " << filename_out << " will be deleted \n";
  } else if (argc == 5) {
    filename = argv[1];
    polyline_file_name = argv[2];
    min_dot = atof(argv[3]);
    filename_out = argv[4];
    keep_output = true;
  }
  else {
    std::cerr << "Usage: " << argv[0] << " [geom_filename] [polygon_file] [min_dot] [output_file]" << std::endl;
    return 1;
  }

  Core mbcore;
  Interface * mb = &mbcore;

  ErrorCode rval = mb->load_file(filename.c_str());
  MB_CHK_SET_ERR(rval,"failed to load input file");

  FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed

  if (!pFacet)
    return 1; // error

  // should the init be part of constructor or not?
  // this is where the obb tree is constructed, and smooth faceting initialized, too.
  rval = pFacet->Init();
  MB_CHK_SET_ERR(rval,"failed to initialize smoothing");

  delete pFacet;
  pFacet = NULL;

  // split_test_across
  std::cout << " split across test: ";
  rval = split_test_across();
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << " verify split ";
  rval = verify_split();
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";
  // when we are done, remove modified file if we want to
  if (!keep_output)
  {
    remove(filename_out.c_str());
  }
  return number_tests_failed;
}

  // this test will test a split like the one for grounding line
  // the first and last point of the polyline should be close to the
  // initial boundary of the face to be split
ErrorCode split_test_across()
{
  Core mbcore;
  Interface * mb = &mbcore;

  ErrorCode  rval = mb->load_file(filename.c_str());MB_CHK_SET_ERR(rval,"failed to load already modified file");

  FBEngine * pFacet = new FBEngine(mb, NULL, true);

  rval = pFacet->Init();MB_CHK_SET_ERR(rval,"failed to initialize smoothing");

  EntityHandle root_set;
  rval = pFacet->getRootSet(&root_set);MB_CHK_SET_ERR(rval, "ERROR : getRootSet failed!" );
  int top = 2; //  iBase_FACE;

  Range faces;
  rval = pFacet->getEntities(root_set, top, faces);MB_CHK_SET_ERR(rval,"Failed to get faces in split_test.");

  if (faces.size() !=1)
  {
    std::cout << "num faces model:" << faces.size() << "\n";
    return MB_FAILURE;//
  }
  // check only the second face

  EntityHandle second_face = faces[0];
  // use the polyPB.txt file to get the trimming polygon
   ;
  // read the file with the polygon user data

  std::ifstream datafile(polyline_file_name.c_str(), std::ifstream::in);
  if (!datafile) {
    std::cout << "can't read file\n";
    return MB_FAILURE;
  }
  //
  char temp[100];
  double direction[3];// normalized
  double gridSize;
  datafile.getline(temp, 100);// first line

  // get direction and mesh size along polygon segments, from file
  sscanf(temp, " %lf %lf %lf %lf ", direction, direction + 1, direction + 2,
      &gridSize);
  //NORMALIZE(direction);// just to be sure

  std::vector<double> xyz;
  while (!datafile.eof()) {
    datafile.getline(temp, 100);
    //int id = 0;
    double x, y, z;
    int nr = sscanf(temp, "%lf %lf %lf", &x, &y, &z);
    if (nr == 3) {
      xyz.push_back(x);
      xyz.push_back(y);
      xyz.push_back(z);
    }
  }
  int sizePolygon = (int)xyz.size()/3;
  if (sizePolygon < 2) {
    std::cerr << " Not enough points in the polygon" << std::endl;
    return MB_FAILURE;
  }

  EntityHandle newFace;// this test is with a "grounding" line
  // the second face should be the one that we want for test
  rval = pFacet->split_surface_with_direction(second_face, xyz, direction, /*closed*/0, /*min_dot */0.8, newFace);MB_CHK_ERR(rval);

  // save a new database, with 3 faces, eventually
  pFacet->delete_smooth_tags();
  delete pFacet;
  pFacet = NULL;// try not to write the obb tree

  rval = mb->write_file(filename_out.c_str());MB_CHK_SET_ERR(rval, "Writing mesh file failed\n");

  return rval;
}

ErrorCode verify_split()
{
  Interface * mb = new Core();

  ErrorCode  rval = mb->load_file(filename_out.c_str());MB_CHK_SET_ERR(rval, "Loading mesh file failed\n");

  moab::GeomTopoTool gTopoTool(mb, true);

  if (!gTopoTool.check_model())
    return MB_FAILURE;

  delete mb;
  return MB_SUCCESS;
}
