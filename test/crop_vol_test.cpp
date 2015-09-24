/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/*
 * crop_vol_test.cpp
 * \brief Driver to create a volume from 2 surfaces (top and bottom), a user defined polyline (contour)
 *   and a direction
 *
 */

#include "moab/Core.hpp"
#include <iostream>
#include <fstream>
#include <string.h>
#include "moab/FBEngine.hpp"
#include "moab/GeomTopoTool.hpp"
#include "TestUtil.hpp"

std::string filename_top;
std::string filename_bot;
std::string polygon_file_name;
double min_dot = 0.8;
std::string vol_file;

int number_tests_successful = 0;
int number_tests_failed = 0;

using namespace moab;

ErrorCode volume_test(FBEngine * pFacet);

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
  filename_bot = TestDir + "/BedCrop2.h5m";
  filename_top = TestDir + "/SECrop2.h5m";//"/polyPB.txt";
  polygon_file_name = TestDir + "/poly14.txt";
  vol_file = "volIce.h5m"; // output

  number_tests_successful = 0;
  bool remove_output = true;

  if (argc == 1) {
    std::cout << "Using default input files: " << filename_bot << " " << filename_top <<
        " " << polygon_file_name << " " << min_dot << " " << vol_file << std::endl;
    std::cout << "    default output file: " << vol_file << " will be deleted \n";
  } else if (argc == 6) {
    remove_output = false;
    filename_bot = argv[1];
    filename_top = argv[2];
    polygon_file_name = argv[3];
    min_dot = atof(argv[4]);
    vol_file = argv[5];
  }
  else
  {
    std::cout << " wrong number of arguments\n";
    return 1; // fail
  }
  Core mbcore;
  Interface * mb = &mbcore;

  ErrorCode rval = mb->load_file(filename_bot.c_str());
  MB_CHK_SET_ERR(rval,"failed to load bed file");

  rval = mb->load_file(filename_top.c_str());// load the second face (we know we have one face in each file)
  MB_CHK_SET_ERR(rval,"failed to load top file");

  FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed

  if (!pFacet)
    return 1; // error


  // should the init be part of constructor or not?
  // this is where the obb tree is constructed, and smooth faceting initialized, too.
  rval = pFacet->Init();
  MB_CHK_SET_ERR(rval,"failed to initialize smoothing");

  std::cout << "volume creation test: ";
  rval = volume_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  if (remove_output)
  {
    remove(vol_file.c_str());
  }
  return number_tests_failed;
}

ErrorCode volume_test (FBEngine * pFacet)
{
  // there should be 2 faces in the model
  std::ifstream datafile(polygon_file_name.c_str(), std::ifstream::in);
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
  if (sizePolygon < 3) {
    std::cerr << " Not enough points in the polygon" << std::endl;
    return MB_FAILURE;
  }
  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  MB_CHK_SET_ERR(rval, "ERROR : getRootSet failed!" );

  int top = 2; //  iBase_FACE;

  Range faces;
  rval = pFacet->getEntities(root_set, top, faces);
  MB_CHK_SET_ERR(rval,"Failed to get faces in volume_test.");


  if (2!=faces.size())
  {
    std::cerr << "expected 2 faces, top and bottom\n";
    return MB_FAILURE;
  }
  EntityHandle newFace1;// first test is with closed surface
  rval = pFacet->split_surface_with_direction(faces[0], xyz, direction,  /*closed*/1, min_dot, newFace1);
  MB_CHK_SET_ERR(rval,"Failed to crop first face.");

  EntityHandle newFace2;// first test is with closed surface
  rval = pFacet->split_surface_with_direction(faces[1], xyz, direction, /*closed*/1, min_dot, newFace2);
  MB_CHK_SET_ERR(rval,"Failed to crop second face.");

  // here we make the assumption that the edges, vertices are matching fine...
  EntityHandle volume;
  rval = pFacet->create_volume_with_direction(newFace1, newFace2, direction, volume);
  MB_CHK_SET_ERR(rval,"Failed to create volume.");

  Interface * mb = pFacet->moab_instance ();
  pFacet->delete_smooth_tags();
  // duplicate -- extract a new geom topo tool
  GeomTopoTool * gtt = pFacet-> get_gtt();
  GeomTopoTool * duplicate = NULL;
  std::vector<EntityHandle> gents;
  gents.push_back(volume);
  rval = gtt->duplicate_model(duplicate, &gents);
  MB_CHK_SET_ERR(rval,"Failed to extract volume.");
  EntityHandle newRootSet = duplicate->get_root_model_set() ;
  delete pFacet;
  pFacet = NULL;// try not to write the obb tree
  delete duplicate;
  rval = mb->write_file(vol_file.c_str(), NULL, NULL, &newRootSet, 1);

  return rval;

}
