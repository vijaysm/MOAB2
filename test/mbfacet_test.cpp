/**
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 */
/**
 * \file mbfacet_test.cpp
 *
 * \brief test mbfacet, a unit test for the FBEngine class, which is providing iGeom like methods to a MOAB db
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
std::string polygon_file_name;

std::string quads_file;

double min_dot = 0.8; // default

bool keep_output;
int number_tests_successful = 0;
int number_tests_failed = 0;

#define PROCESS_ERROR(A, B)  {if (A!=MB_SUCCESS) {  std::cout << B << std::endl; return 1; } }

#define CHECK( STR ) if (rval != MB_SUCCESS) return print_error( STR, rval, __FILE__, __LINE__ )

using namespace moab;

ErrorCode print_error(const char* desc, ErrorCode rval, const char* file,
    int line)
{
  std::cerr << "ERROR: " << desc << std::endl << "  Error code: " << rval
      << std::endl << "  At        : " << file << ':' << line << std::endl;
  return MB_FAILURE; // must always return false or CHECK macro will break
}

ErrorCode root_set_test(FBEngine * pFacet);
ErrorCode gentityset_test(FBEngine * pFacet);
ErrorCode geometry_evaluation_test(FBEngine * pFacet);
ErrorCode normals_test(FBEngine * pFacet);
ErrorCode ray_test(FBEngine * pFacet);
ErrorCode split_test(Interface * mb, FBEngine * pFacet);
ErrorCode check_split(Interface * mb);

ErrorCode split_quads_test();

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
  polygon_file_name = TestDir + "/polyPB.txt";
  filename_out = "PB_facet.h5m";
  quads_file = TestDir + "/quads.h5m";

  min_dot = 0.8;

  keep_output = false;
  bool only_cropping = false;
  if (argc == 1) {
    std::cout << "Using default input files: " << filename << " " << polygon_file_name <<
        " " << filename_out << " " << min_dot <<  " "<< quads_file << std::endl;
    std::cout << "    default output file: " << filename_out << " will be deleted \n";
  } else if (argc == 5) {
    only_cropping = true;
    filename = argv[1];
    polygon_file_name = argv[2];
    filename_out = argv[3];
    min_dot = atof(argv[4]);
  } else if (argc >=6) {
    filename = argv[1];
    polygon_file_name = argv[2];
    filename_out = argv[3];
    min_dot = atof(argv[4]);
    quads_file = argv[5];
    keep_output = true;
    std::cout << "Using input files: " << filename << " " << polygon_file_name <<
            " " << std::endl;
    std::cout << "    default output file: " << filename_out << " will be saved \n";
    std::cout << " min_dot: " << min_dot << " quads file:" << quads_file << "\n";

  } else {
    std::cerr << "Usage: " << argv[0] << " [geom_filename] [polygon_file] [output_file] [min_dot] [quads_file]" << std::endl;
    return 1;
  }

  Core mbcore;
  Interface * mb = &mbcore;

  ErrorCode rval = mb->load_file(filename.c_str());
  PROCESS_ERROR(rval,"failed to load input file");

  FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed

  if (!pFacet)
    return 1; // error

  // should the init be part of constructor or not?
  // this is where the obb tree is constructed, and smooth faceting initialized, too.
  rval = pFacet->Init();
  PROCESS_ERROR(rval,"failed to initialize smoothing");

  std::cout << "root set test: ";
  rval = root_set_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "gentity set test: ";
  rval = gentityset_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "geometry evals test: ";
  rval = geometry_evaluation_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "normal evals test: ";
  rval = normals_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "ray test: ";
  rval = ray_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "split with loop test ";
  rval = split_test(mb, pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  // pFacet has been deleted in split_test(), so we should mark it as NULL
  // Setting pFacet (value parameter) to NULL in split_test() does not work
  pFacet = NULL;

  std::cout << " check split: ";
  rval = check_split(mb);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  if (only_cropping)
    return number_tests_failed;// we do not remove the output file, either

  std::cout << " split quads test: ";
  rval = split_quads_test();
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  // when we are done, remove modified file if we want to
  if (!keep_output)
  {
    remove(filename_out.c_str());
  }
  return number_tests_failed;
}

ErrorCode root_set_test(FBEngine * pFacet)
{

  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );

  return MB_SUCCESS;
}

ErrorCode gentityset_test(FBEngine * pFacet)
{
  int num_type = 4;
  EntityHandle ges_array[4];
  int number_array[4];
  //int num_all_gentities_super = 0;
  int ent_type = 0; // iBase_VERTEX;

  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );

  // get the number of sets in the whole model
  int all_sets = 0;
  rval = pFacet->getNumEntSets(root_set, 0, &all_sets);
  CHECK( "Problem getting the number of all gentity sets in whole model." );// why do we count all sets


  // add gentities to entitysets by type
  for (; ent_type < num_type; ent_type++) {
    // initialize the entityset
    rval = pFacet->createEntSet(1, &ges_array[ent_type]);// 1 means isList, ordered
    CHECK( "Problem creating entityset." );

    // get entities by type in total "mesh"
    Range gentities;
    rval = pFacet->getEntities(root_set, ent_type, gentities);
    CHECK( "Failed to get gentities by type in gentityset_test." );

    // add gentities into gentity set
    rval = pFacet->addEntArrToSet(gentities, ges_array[ent_type]);
    CHECK( "Failed to add gentities in entityset_test." );

    // Check to make sure entity set really has correct number of entities in it
    rval = pFacet->getNumOfType(ges_array[ent_type], ent_type,
        &number_array[ent_type]);
    CHECK( "Failed to get number of gentities by type in entityset_test." );

    // compare the number of entities by type
    int num_type_gentity = gentities.size();

    if (number_array[ent_type] != num_type_gentity) {
      std::cerr << "Number of gentities by type is not correct" << std::endl;
      return MB_FAILURE;
    }

    // add to number of all entities in super set
    //num_all_gentities_super += num_type_gentity;
  }

  // make a super set having all entitysets
  EntityHandle super_set;
  rval = pFacet->createEntSet(1, &super_set);// 1 is list
  CHECK( "Failed to create a super set in gentityset_test." );

  for (int i = 0; i < num_type; i++) {
    rval = pFacet->addEntSet(ges_array[i], super_set);
    CHECK( "Failed to create a super set in gentityset_test." );
  }

  //----------TEST BOOLEAN OPERATIONS----------------//

  EntityHandle temp_ges1;
  rval = pFacet->createEntSet(1, &temp_ges1);
  CHECK( "Failed to create a super set in gentityset_test." );

  // Subtract
  // add all EDGEs and FACEs to temp_es1
  // get all EDGE entities
  Range gedges, gfaces, temp_gentities1;
  rval = pFacet->getEntities(ges_array[1], /* iBase_EDGE*/1, gedges);
  CHECK( "Failed to get gedge gentities in gentityset_test." );

  // add EDGEs to ges1
  rval = pFacet->addEntArrToSet(gedges, temp_ges1);
  CHECK( "Failed to add gedge gentities in gentityset_test." );

  // get all FACE gentities
  rval = pFacet->getEntities(ges_array[2], /*iBase_FACE*/2, gfaces);
  CHECK( "Failed to get gface gentities in gentityset_test." );

  // add FACEs to es1
  rval = pFacet->addEntArrToSet(gfaces, temp_ges1);
  CHECK( "Failed to add gface gentities in gentityset_test." );

  // subtract EDGEs
  rval = pFacet->gsubtract(temp_ges1, ges_array[1], temp_ges1);
  CHECK( "Failed to subtract gentitysets in gentityset_test." );

  rval = pFacet->getEntities(temp_ges1, 2, temp_gentities1);
  CHECK( "Failed to get gface gentities in gentityset_test." );

  if (gfaces.size() != temp_gentities1.size()) {
    std::cerr
        << "Number of entitysets after subtraction not correct \
             in gentityset_test."
        << std::endl;
    return MB_FAILURE;
  }

  // check there's nothing but gfaces in temp_ges1
  int num_gents;
  rval = pFacet->getNumOfType(temp_ges1, 1, &num_gents);
  CHECK( "Failed to get dimensions of gentities in gentityset_test." );
  if (0 != num_gents) {
    std::cerr << "Subtraction failed to remove all edges" << std::endl;
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

ErrorCode geometry_evaluation_test(FBEngine * pFacet)
{
  int i;
  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );

  int top = 0;// iBase_VERTEX;
  int num_test_top = 4; // iBase_ALL_TYPES;
  std::vector<std::vector<EntityHandle> > gentity_vectors(num_test_top);

  // fill the vectors of each topology entities
  //
  for (i = top; i < num_test_top; i++) {
    Range gentities;
    rval = pFacet->getEntities(root_set, i, gentities);
    CHECK("Failed to get gentities in eval tests.");

    gentity_vectors[i].resize(gentities.size());
    std::copy(gentities.begin(), gentities.end(), gentity_vectors[i].begin());
  }

  // check geometries in both directions
  double min[3], max[3], on[3];
  double near[3] = { .0, .0, .0 };
  std::vector<EntityHandle>::iterator vit;
  /*for (i = iBase_REGION; i >= iBase_VERTEX; i--) {
   if (i != iBase_EDGE) {*/
  for (i = 3; i >= 0; i--) {
    if (i != 1) {
      for (vit = gentity_vectors[i].begin(); vit != gentity_vectors[i].end(); ++vit) {
        EntityHandle this_gent = *vit;
        rval = pFacet->getEntBoundBox(this_gent, &min[0], &min[1], &min[2],
            &max[0], &max[1], &max[2]);
        CHECK("Failed to get bounding box of entity.");

        for (int j=0; j<3; j++)
          near[j] = (min[j]+max[j])/2.;
        rval = pFacet->getEntClosestPt(this_gent, near[0], near[1], near[2],
            &on[0], &on[1], &on[2]);
        CHECK("Failed to get closest point on entity.");
      }
    }
    else
    {
      // for edges, provide a little better help
      for (vit = gentity_vectors[i].begin(); vit != gentity_vectors[i].end(); ++vit) {
        EntityHandle this_gent = *vit;
        // we know that the edge is parametric, with par between 0 and 1

        // some random parameter

        rval = pFacet->getEntUtoXYZ(this_gent, 0.33, near[0], near[1], near[2]);
        CHECK("Failed to get a new point");

        std::cout << " entity of type " << i << " position:\n  "
                  << near[0] << " " << near[1] << " " << near[2] << "\n";
        rval = pFacet->getEntClosestPt(this_gent, near[0], near[1], near[2],
            &on[0], &on[1], &on[2]);
        CHECK("Failed to get closest point on entity.");
        std::cout << "   close by:  "
                          << on[0] << " " << on[1] << " " << on[2] << "\n";
      }
    }
  }


  return MB_SUCCESS;
}

//  test normals evaluations on the surface only
ErrorCode normals_test(FBEngine * pFacet)
{
  int i;
  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );

  int top = 0;//iBase_VERTEX;
  int num_test_top = 4; //iBase_ALL_TYPES;
  std::vector<std::vector<EntityHandle> > gentity_vectors(num_test_top);

  // fill the vectors of each topology entities
  //
  for (i = top; i < num_test_top; i++) {
    Range gentities;
    rval = pFacet->getEntities(root_set, i, gentities);
    CHECK("Failed to get gentities in eval tests.");

    gentity_vectors[i].resize(gentities.size());
    std::copy(gentities.begin(), gentities.end(), gentity_vectors[i].begin());
  }

  // check adjacencies in both directions
  double min[3], max[3];
  double normal[3] = { .0, .0, .0 };
  std::vector<EntityHandle>::iterator vit;
  for (i = 3/*iBase_REGION*/; i > 1 /*iBase_EDGE*/; i--) {
    for (vit = gentity_vectors[i].begin(); vit != gentity_vectors[i].end(); ++vit) {
      EntityHandle this_gent = *vit;
      rval = pFacet->getEntBoundBox(this_gent, &min[0], &min[1], &min[2],
          &max[0], &max[1], &max[2]);
      CHECK("Failed to get bounding box of entity.");

      rval = pFacet->getEntNrmlXYZ(this_gent, (max[0] + min[0]) / 2, (max[1]
          + min[1]) / 2, (max[2] + min[2]) / 2, &normal[0], &normal[1],
          &normal[2]);

      CHECK("Failed to get normal to the closest point.");
      std::cout << " entity of type " << i << " closest normal to center:\n  "
          << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
    }
  }

  return MB_SUCCESS;
}

//  ray test
ErrorCode ray_test(FBEngine * pFacet)
{

  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );

  int top = 2; //  iBase_FACE;

  Range faces;
  rval = pFacet->getEntities(root_set, top, faces);
  CHECK("Failed to get faces in ray_test.");

  // check only the first face

  // check adjacencies in both directions
  double min[3], max[3];

  EntityHandle first_face = faces[0];

  rval = pFacet->getEntBoundBox(first_face, &min[0], &min[1], &min[2], &max[0],
      &max[1], &max[2]);
  CHECK("Failed to get bounding box of entity.");

  // assume that the ray shot from the bottom of the box (middle) is a pretty good candidate
  // in z direction
  double x = (min[0] + max[0]) / 2, y = (min[1] + max[1]) / 2, z = min[2];
  std::vector<EntityHandle> intersect_entity_handles;
  std::vector<double> intersect_coords;
  std::vector<double> param_coords;
  rval = pFacet->getPntRayIntsct(x, y, z, // shot from
      0., 0., 1., // direction
      intersect_entity_handles,
      /*iBase_INTERLEAVED,*/
      intersect_coords, param_coords);

  CHECK("Failed to find ray intersections points ");
  for (unsigned int i = 0; i < intersect_entity_handles.size(); i++) {

    std::cout << " entity of type " << pFacet->moab_instance ()->type_from_handle(intersect_entity_handles[i]) <<
        " id from handle: " << pFacet->moab_instance ()->id_from_handle(intersect_entity_handles[i])  <<
        "\n" << intersect_coords[3 * i] << " " << intersect_coords[3 * i + 1] << " "
        << intersect_coords[3 * i + 2] << "\n" << " distance: " << param_coords[i] << "\n";
  }

  return MB_SUCCESS;
}

// this test is for creating 2 surfaces given a polyline and a direction for piercing.
// this test is for cropping a surface with a closed loop, no intersection
// with initial boundary
ErrorCode split_test(Interface * mb, FBEngine * pFacet)
{

  EntityHandle root_set;
  ErrorCode rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );

  int top = 2; //  iBase_FACE;

  Range faces;
  rval = pFacet->getEntities(root_set, top, faces);
  CHECK("Failed to get faces in split_test.");

  // check only the first face

  EntityHandle first_face = faces[0];
  // use the polyPB.txt file to get the trimming polygon
   ;
  // read the file with the polygon user data

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

  EntityHandle newFace;// first test is with closed surface
  rval = pFacet->split_surface_with_direction(first_face, xyz, direction,   /*closed*/1, min_dot, newFace);

  if (rval!=MB_SUCCESS)
    return rval;
  // save a new database
  pFacet->delete_smooth_tags();
  // get a new gtt tool

  // duplicate -- extract a new geom topo tool
  GeomTopoTool * gtt = pFacet-> get_gtt();
  GeomTopoTool * duplicate = NULL;
  std::vector<EntityHandle> gents;
  gents.push_back(newFace);
  rval = gtt->duplicate_model(duplicate, &gents);
  CHECK("Failed to extract surface.");
  EntityHandle newRootSet = duplicate->get_root_model_set() ;
  delete pFacet;
  pFacet = NULL;// try not to write the obb tree
  rval = mb->write_file(filename_out.c_str(), NULL, NULL, &newRootSet, 1);

  delete duplicate;

  return rval;
}

ErrorCode check_split(Interface * mb)
{
  // check loading the file in an empty db
  //delete pFacet;// should clean up the FBEngine
  ErrorCode rval = mb->delete_mesh();
  CHECK( "ERROR : delete mesh failed!" );

  rval = mb->load_file(filename_out.c_str());
  CHECK( "ERROR : can't load modified file!" );

  FBEngine * pFacet = new FBEngine(mb, NULL, true);// smooth facetting, no OBB tree passed

  // repeat tests on modified file

  // should the init be part of constructor or not?
  // this is where the obb tree is constructed, and smooth faceting initialized, too.
  rval = pFacet->Init();
  CHECK("failed to initialize smoothing");

  std::cout << "root set test: ";
  rval = root_set_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "gentity set test: ";
  rval = gentityset_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "geometry evals test: ";
  rval = geometry_evaluation_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "normal evals test: ";
  rval = normals_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "ray test: ";
  rval = ray_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  delete pFacet;
  pFacet = NULL;
  if (number_tests_failed>0)
    return MB_FAILURE;

  return MB_SUCCESS;
}

ErrorCode split_quads_test()
{
  Core mbcore;
  Interface * mb = &mbcore;

  ErrorCode  rval = mb->load_file(quads_file.c_str());
  CHECK("failed to load quads file");
  FBEngine * pFacet = new FBEngine(mb, NULL, true);

  rval = pFacet->Init();
  CHECK("failed to initialize smoothing");

  EntityHandle root_set;
  rval = pFacet->getRootSet(&root_set);
  CHECK( "ERROR : getRootSet failed!" );
  int top = 2; //  iBase_FACE;

  Range faces;
  rval = pFacet->getEntities(root_set, top, faces);
  CHECK("Failed to get faces in split_test.");

  if (faces.size() !=2)
  {
    std::cout << "num faces after loading quads model:" << faces.size() << "\n";
    return MB_FAILURE;//
  }


  // multiple tests
  std::cout << "gentity set test: ";
  rval = gentityset_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "geometry evals test: ";
  rval = geometry_evaluation_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "normal evals test: ";
  rval = normals_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  std::cout << "ray test: ";
  rval = ray_test(pFacet);
  handle_error_code(rval, number_tests_failed, number_tests_successful);
  std::cout << "\n";

  // export split triangles model, after deleting the smooth tags
  pFacet->delete_smooth_tags();

  delete pFacet;
  pFacet= NULL;
  std::string spl_file="q.split.h5m";
  rval = mb->write_file(spl_file.c_str(), 0, 0, &root_set, 1);
  CHECK("can't write result file");

  if (!keep_output)
  {
    remove(spl_file.c_str());
  }

  if (number_tests_failed>0)
    return MB_FAILURE;

  return MB_SUCCESS;
}

