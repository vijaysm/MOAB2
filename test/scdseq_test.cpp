/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Coroporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

//
// test the structured sequence stuff
//

#include "moab/ScdInterface.hpp"
#include "moab/Core.hpp"
#include "TestUtil.hpp"

#include <iostream>

using namespace moab;

void test_vertex_seq();
void test_element_seq();
void test_periodic_seq();
void test_parallel_partitions();

ErrorCode test_parallel_partition(int *gdims, int nprocs, int part_method);

ErrorCode check_vertex_sequence(const ScdBox *this_box, 
                                const int imin, const int jmin, const int kmin, 
                                const int imax, const int jmax, const int kmax, 
                                const EntityHandle this_start);
ErrorCode evaluate_vertex_box(ScdBox *this_box);
ErrorCode check_element_sequence(const ScdBox *this_box, 
                                 const HomCoord &min_params,
                                 const HomCoord &max_params,
                                 const EntityHandle this_start);
ErrorCode evaluate_element_sequence(ScdBox *this_box);
ErrorCode eseq_test1a(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max);
ErrorCode eseq_test1b(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max);
ErrorCode eseq_test1c(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max);
ErrorCode eseq_test2a(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max);
ErrorCode eseq_test2b(ScdInterface *scdi);
ErrorCode eseq_test2c(ScdInterface *scdi);
ErrorCode eseq_test2d(ScdInterface *scdi);

ErrorCode create_1d_3_sequences(ScdInterface *scdi,
                                HomCoord tmp_min, HomCoord tmp_max,
                                ScdBox **vseq, ScdBox **eseq);

ErrorCode create_2d_3_sequences(ScdInterface *scdi,
                                ScdBox **vseq, ScdBox **eseq);

ErrorCode create_2dtri_3_sequences(ScdInterface *scdi,
                                   const int int1, const int int2, const int int3,
                                   ScdBox **vseq, ScdBox **eseq);
ErrorCode create_3dtri_3_sequences(ScdInterface *scdi,
                                   const int int1, const int int2, const int int3, const int int4,
                                   ScdBox **vseq, ScdBox **eseq);
ErrorCode access_adjacencies(ScdBox *box);

void test_partition_methods();

void test_partition_method(ScdInterface::PartitionMethod pm);


// first comes general-capability code used by various tests; main and test functions
// come after these, starting with main
ErrorCode check_vertex_sequence(const ScdBox *this_box, 
                                const int imin, const int jmin, const int kmin, 
                                const int imax, const int jmax, const int kmax, 
                                const EntityHandle this_start) 
{
  ErrorCode result = MB_SUCCESS;
  
    // check data stored in sequence with that in arg list
  if (imin != this_box->box_min()[0] ||
      jmin != this_box->box_min()[1] ||
      kmin != this_box->box_min()[2] ||
      imax != this_box->box_max()[0] ||
      jmax != this_box->box_max()[1] ||
      kmax != this_box->box_max()[2]) {
    std::cout << "min/max params not correct for a sequence." << std::endl;
    result = MB_FAILURE;
  }

  HomCoord ijk1 = this_box->box_min(),
      ijk2 = this_box->box_max(),
      dijk = this_box->box_size();
  if (ijk2 - ijk1 + HomCoord(1,1,1,0) != dijk) {
    std::cout << "min_params/max_params/param_extents functions returned inconsistent results"
              << std::endl;
    result = MB_FAILURE;
  }
  
  if (this_start != this_box->start_vertex()) {
    std::cout << "Start handle for sequence wrong." << std::endl;
    result = MB_FAILURE;
  }
  
  return result;
}
  
ErrorCode check_element_sequence(const ScdBox *this_box, 
                                 const HomCoord &min_params,
                                 const HomCoord &max_params,
                                 const EntityHandle this_start) 
{
  ErrorCode result = MB_SUCCESS;
  
    // check data stored in sequence with that in arg list
  if (min_params.i() != this_box->box_min()[0] ||
      min_params.j() != this_box->box_min()[1] ||
      min_params.k() != this_box->box_min()[2] ||
      max_params.i() != this_box->box_max()[0] ||
      max_params.j() != this_box->box_max()[1] ||
      max_params.k() != this_box->box_max()[2]) {
    std::cout << "min/max params not correct for a sequence." << std::endl;
    result = MB_FAILURE;
  }

  if (this_box->box_max() - this_box->box_min() + HomCoord(1,1,1,0) != this_box->box_size()) {
    std::cout << "min_params/max_params/param_extents functions returned inconsistent results"
              << std::endl;
    result = MB_FAILURE;
  }
  
  if (this_start != this_box->start_element()) {
    std::cout << "Start handle for sequence wrong." << std::endl;
    result = MB_FAILURE;
  }

  if (this_box->boundary_complete() == false) {
    std::cout << "Element sequence didn't pass boundary_complete test." << std::endl;
    result = MB_FAILURE;
  }

  return result;
}
  
ErrorCode evaluate_vertex_sequence(ScdBox *this_box) 
{
  ErrorCode result = MB_SUCCESS;
  
    // first get the parametric extents
  HomCoord ijk1 = this_box->box_min(), ijk2 = this_box->box_max();

    // then the start vertex
  EntityHandle start_handle = this_box->start_vertex();
  
    // now evaluate all the vertices in forward and reverse
  EntityHandle tmp_handle, tmp_handle2;
  for (int i = ijk1[0]; i <= ijk2[0]; i++) {
    for (int j = ijk1[1]; j <= ijk2[1]; j++) {
      for (int k = ijk1[2]; k <= ijk2[2]; k++) {
          // compute what the vertex handle is supposed to be
        EntityHandle this_handle = start_handle + (i-ijk1[0]) + (j-ijk1[1])*(ijk2[0]-ijk1[0]+1) + 
            (k-ijk1[2])*(ijk2[1]-ijk1[1]+1)*(ijk2[0]-ijk1[0]+1);

          // get_vertex variants
        tmp_handle = this_box->get_vertex(i, j, k);

        if (this_box->get_vertex(i, j, k) != this_handle) {
          std::cout << "vertex seq: get_vertex(i, j, k) didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
        
        tmp_handle2 = this_box->get_vertex(HomCoord(i, j, k));
        if (tmp_handle2 != this_handle) {
          std::cout << "vertex seq: get_vertex(HomCoord) didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }

        int itmp, jtmp, ktmp;
        itmp = jtmp = ktmp = 0xdeadbeef;
        ErrorCode tmp_result = this_box->get_params(tmp_handle, itmp, jtmp, ktmp);
        if (MB_SUCCESS != tmp_result || i != itmp || j != jtmp || k != ktmp) {
          std::cout << "vertex seq: get_params didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "; itmp, jtmp, ktmp = " 
                    << itmp << ", " << jtmp << ", " << ktmp << std::endl;
          result = MB_FAILURE;
        }

        if (!this_box->contains(i, j, k) || !this_box->contains(HomCoord(i,j,k))) {
          std::cout << "vertex seq: contains didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
      }
    }
  }

  return result;
}

ErrorCode evaluate_element_sequence(ScdBox *this_box) 
{
  ErrorCode result = MB_SUCCESS;
  
    // first get the parametric extents
  HomCoord ijk1 = this_box->box_min(), ijk2 = this_box->box_max();

    // now evaluate all the vertices and elements in forward and reverse
  EntityHandle tmp_handle, tmp_handle2;
  int is_periodic_i = (this_box->is_periodic_i() ? 1 : 0),
      is_periodic_j = (this_box->is_periodic_j() ? 1 : 0);
  for (int i = ijk1[0]; i < ijk2[0]+is_periodic_i; i++) {
    for (int j = ijk1[1]; j < ijk2[1]+is_periodic_j; j++) {
      for (int k = ijk1[2]; k < ijk2[2]; k++) {

          // get_vertex variants
        tmp_handle = this_box->get_vertex(i, j, k);

        tmp_handle2 = this_box->get_vertex(HomCoord(i, j, k));
        if (tmp_handle2 != tmp_handle) {
          std::cout << "element seq: get_vertex(HomCoord) and get_vertex(i,j,k) didn't return" << std::endl
                    << "consistent results, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }

          // get element variants
        tmp_handle = this_box->get_element(i, j, k);
        
        tmp_handle2 = this_box->get_element(HomCoord(i,j,k));
        if (tmp_handle2 != tmp_handle) {
          std::cout << "element seq: get_element(HomCoord) and get_element(i,j,k) didn't return" << std::endl
                    << "consistent results, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
        
          // get_params
        int itmp, jtmp, ktmp;
        itmp = jtmp = ktmp = 0xdeadbeef;
        ErrorCode tmp_result = this_box->get_params(tmp_handle, itmp, jtmp, ktmp);
        if (MB_SUCCESS != tmp_result || i != itmp || j != jtmp || k != ktmp) {
          std::cout << "element seq: get_params didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "; itmp, jtmp, ktmp = " 
                    << itmp << ", " << jtmp << ", " << ktmp << std::endl;
          result = MB_FAILURE;
        }

        if (!this_box->contains(i, j, k) || !this_box->contains(HomCoord(i,j,k))) {
          std::cout << "element seq: contains didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
      }
    }
  }

  return result;
}

int main(int, char**) 
{
    // test partition methods
  RUN_TEST(test_parallel_partitions);
    
    // test creating and evaluating vertex sequences
  RUN_TEST(test_vertex_seq);

    // test creating and evaluating element sequences
  RUN_TEST(test_element_seq);

    // test periodic sequences
  RUN_TEST(test_periodic_seq);

}

void test_vertex_seq() 
{
  Core moab;
  ScdInterface *scdi;
  ErrorCode rval = moab.Interface::query_interface(scdi);
  CHECK_ERR(rval);
  
    // get the seq manager from gMB
  ScdBox *dum_box = NULL;
  
  
    // make a 1d sequence
  ErrorCode result = scdi->create_scd_sequence(HomCoord(-10, 0, 0), 
                                               HomCoord(10, 0, 0),
                                               MBVERTEX, 1, dum_box);
  CHECK_ERR(result);
  if (!dum_box || dum_box->start_vertex() == 0) CHECK_ERR(MB_FAILURE);

    // check sequence data
  result = check_vertex_sequence(dum_box, -10, 0, 0, 10, 0, 0, dum_box->start_vertex());
  CHECK_ERR(result);

    // evaluate that sequence for all possible values
  result = evaluate_vertex_sequence(dum_box);
  CHECK_ERR(result);

    // make a 2d sequence
  dum_box = NULL;
  result = scdi->create_scd_sequence(HomCoord(-10, -10, 0), 
                                     HomCoord(10, 10, 0),
                                     MBVERTEX, 1, dum_box);
  CHECK_ERR(result);

    // check sequence data
  result = check_vertex_sequence(dum_box, -10, -10, 0, 10, 10, 0, dum_box->start_vertex());
  CHECK_ERR(result);

    // evaluate that sequence for all possible values
  result = evaluate_vertex_sequence(dum_box);
  CHECK_ERR(result);

    // make a 3d sequence
  dum_box = NULL;
  result = scdi->create_scd_sequence(HomCoord(-10, -10, -10), 
                                     HomCoord(10, 10, 10),
                                     MBVERTEX, 1, dum_box);
  CHECK_ERR(result);
  if (!dum_box || dum_box->start_vertex() == 0) CHECK_ERR(MB_FAILURE);

    // check sequence data
  result = check_vertex_sequence(dum_box, -10, -10, -10, 10, 10, 10, dum_box->start_vertex());
  CHECK_ERR(result);

    // evaluate that sequence for all possible values
  result = evaluate_vertex_sequence(dum_box);
  CHECK_ERR(result);
}

void test_element_seq() 
{
  Core moab;
  ScdInterface *scdi;
  ErrorCode rval = moab.Interface::query_interface(scdi);
  CHECK_ERR(rval);

  HomCoord TEST_MIN_PARAMS(0,0,0);
  HomCoord TEST_BOX_MAX(11,5,2);

    // TEST 1: single vertex sequence blocks, unity mapping
  rval = eseq_test1a(scdi, TEST_MIN_PARAMS, TEST_BOX_MAX);
  CHECK_ERR(rval);
  rval = eseq_test1b(scdi, TEST_MIN_PARAMS, TEST_BOX_MAX);
  CHECK_ERR(rval);
  rval = eseq_test1c(scdi, TEST_MIN_PARAMS, TEST_BOX_MAX);
  CHECK_ERR(rval);

    // TEST 2: composites, 0d difference between element block "owning"
    // vertex block and element block sharing vertex block
  rval = eseq_test2a(scdi, TEST_MIN_PARAMS, TEST_BOX_MAX);
  CHECK_ERR(rval);
  rval = eseq_test2b(scdi);
  CHECK_ERR(rval);
  rval = eseq_test2c(scdi);
  CHECK_ERR(rval);
  rval = eseq_test2d(scdi);
  CHECK_ERR(rval);
}

ErrorCode eseq_test1a(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 1a: 1d single vertex seq block, min/max = (-10,0,0)/(10,0,0)
    // create vertex seq

    // first get 1d min/max by resetting j, k components
  tmp_min[1] = tmp_min[2] = tmp_max[1] = tmp_max[2] = 0;

  ScdBox *ebox, *vbox;
  ErrorCode result = scdi->create_scd_sequence(tmp_min, tmp_max,
                                               MBVERTEX, 1, vbox);
  CHECK_ERR(result);

    // now create the element sequence
  result = scdi->create_scd_sequence(tmp_min, tmp_max,
                                     MBEDGE, 1, ebox);
  CHECK_ERR(result);
  
    // add vertex seq to element seq
  result = ebox->add_vbox(vbox,tmp_min, tmp_min,
                          tmp_max, tmp_max, tmp_min, tmp_min);
  CHECK_ERR(result);
  
    // check/evaluate element sequence
  result = check_element_sequence(ebox, tmp_min, tmp_max, ebox->start_element());
  CHECK_ERR(result);
  
  result = evaluate_element_sequence(ebox);
  CHECK_ERR(result);
  
  return result;
}

ErrorCode eseq_test1b(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 1b: 2d single vertex seq block, min/max = (-10,-5,0)/(10,5,0)

    // first get 2d min/max by resetting k component
  tmp_min[2] = tmp_max[2] = 0;

    // get the seq manager from gMB
  ScdBox *ebox, *vbox;
  ErrorCode result = scdi->create_scd_sequence(tmp_min, tmp_max,
                                               MBVERTEX, 1, vbox);
  CHECK_ERR(result);

    // now create the element sequence
  result = scdi->create_scd_sequence(tmp_min, tmp_max,
                                     MBQUAD, 1, ebox);
  CHECK_ERR(result);

    // add vertex seq to element seq; first need to construct proper 3pt input (p1 is tmp_min)
  HomCoord p2(tmp_max.i(), tmp_min.j(), tmp_min.k());
  HomCoord p3(tmp_min.i(), tmp_max.j(), tmp_min.k());
  result = ebox->add_vbox(vbox, tmp_min, tmp_min,
                          p2, p2, p3, p3);
  CHECK_ERR(result);
  
  std::vector<EntityHandle> connect;
  EntityHandle dum_ent = ebox->start_element();
  result = scdi->impl()->get_connectivity(&dum_ent, 1, connect);
  CHECK_ERR(result);
  CHECK_EQUAL((unsigned int)connect.size(), (unsigned int)4);

    // check/evaluate element sequence
  result = check_element_sequence(ebox, tmp_min, tmp_max, ebox->start_element());
  CHECK_ERR(result);
  
  result = evaluate_element_sequence(ebox);
  CHECK_ERR(result);

  result = access_adjacencies(ebox);
  CHECK_ERR(result);
  
  return result;
}

ErrorCode eseq_test1c(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 1c: 3d single vertex seq block, min/max = (-10,-5,-1)/(10,5,1)

    // get the seq manager from gMB
  ScdBox *ebox, *vbox;
  ErrorCode result = scdi->create_scd_sequence(tmp_min, tmp_max,
                                               MBVERTEX, 1, vbox);
  CHECK_ERR(result);
  

    // now create the element sequence
  result = scdi->create_scd_sequence(tmp_min, tmp_max,
                                     MBHEX, 1, ebox);
  CHECK_ERR(result);
  
    // add vertex seq to element seq; first need to construct proper 3pt input (p1 is tmp_min)
  HomCoord p2(tmp_max.i(), tmp_min.j(), tmp_min.k());
  HomCoord p3(tmp_min.i(), tmp_max.j(), tmp_min.k());
  result = ebox->add_vbox(vbox, tmp_min, tmp_min, p2, p2, p3, p3);
  CHECK_ERR(result);
  
    // check/evaluate element sequence
  result = check_element_sequence(ebox, tmp_min, tmp_max, ebox->start_element());
  CHECK_ERR(result);
  
  result = evaluate_element_sequence(ebox);
  CHECK_ERR(result);

  result = access_adjacencies(ebox);
  CHECK_ERR(result);
  
  return result;
}

ErrorCode eseq_test2a(ScdInterface *scdi, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 2a: 1d composite block, 0d difference between owning/sharing blocks
    // create vertex seq
  
  ScdBox *vbox[3], *ebox[3];
  ErrorCode result = create_1d_3_sequences(scdi, tmp_min, tmp_max,
                                           vbox, ebox);
  CHECK_ERR(result);

    // whew; that's done; now check and evaluate

    // first check to make sure the parameter spaces tack onto one another
  if (ebox[0]->box_min() != tmp_min ||
      ebox[0]->box_max() != ebox[1]->box_min() ||
      ebox[1]->box_max() != ebox[2]->box_min())
    CHECK_ERR(MB_FAILURE);

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    result = check_element_sequence(ebox[i], ebox[i]->box_min(), ebox[i]->box_max(), 
                                    ebox[i]->start_element());
    CHECK_ERR(result);
  
    result = evaluate_element_sequence(ebox[i]);
    CHECK_ERR(result);
  }
  
  return result;
}

ErrorCode eseq_test2b(ScdInterface *scdi) 
{
    // TEST 2b: 2d composite block, 0d difference between owning/sharing blocks
    // create vertex seq
  
  ScdBox *ebox[3], *vbox[3];
  ErrorCode result = create_2d_3_sequences(scdi, vbox, ebox);
  CHECK_ERR(result);

    // whew; that's done; now check and evaluate

    // first check to make sure the parameter spaces tack onto one another
  if (ebox[0]->box_max() != HomCoord(ebox[1]->box_min().i(),ebox[1]->box_max().j(),0) ||
      ebox[1]->box_max() != HomCoord(ebox[2]->box_min().i(),ebox[2]->box_max().j(),0) ||
        // cheat on the i value of ebox0, since it's a periodic bdy (we don't check for that, so you
        // may get different values for that parameter)
      ebox[2]->box_max() != HomCoord(ebox[2]->box_max().i(),ebox[0]->box_max().j(),0))
    CHECK_ERR(MB_FAILURE);
    

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    result = check_element_sequence(ebox[i], ebox[i]->box_min(), ebox[i]->box_max(), 
                                    ebox[i]->start_element());
    CHECK_ERR(result);

    result = evaluate_element_sequence(ebox[i]);
    CHECK_ERR(result);
  }

  for (int i = 0; i < 3; i++) {
    result = access_adjacencies(ebox[i]);
    CHECK_ERR(result);
  }
  
  return result;
}

ErrorCode eseq_test2c(ScdInterface *scdi) 
{
    // TEST 2c: 2d composite block, 0d difference between owning/sharing blocks,
    // tri-valent shared vertex between the three blocks

    // interval settings: only 3 of them
  int int1 = 5, int2 = 15, int3 = 25;
  ScdBox *ebox[3], *vbox[3];
  ErrorCode result = create_2dtri_3_sequences(scdi, int1, int2, int3, vbox, ebox);
  CHECK_ERR(result);

    // whew; that's done; now check and evaluate

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    result = check_element_sequence(ebox[i], ebox[i]->box_min(), ebox[i]->box_max(), 
                                    ebox[i]->start_element());
    CHECK_ERR(result);
  
    result = evaluate_element_sequence(ebox[i]);
    CHECK_ERR(result);
  }
  
  for (int i = 0; i < 3; i++) {
    result = access_adjacencies(ebox[i]);
    CHECK_ERR(result);
  }

  return result;
}

ErrorCode eseq_test2d(ScdInterface *scdi) 
{
    // TEST 2d: 3d composite block, 0d difference between owning/sharing blocks,
    // tri-valent shared edge between the three blocks

    // interval settings: only 3 of them
  int int1 = 100, int2 = 100, int3 = 100, int4 = 100;
  ScdBox *ebox[3], *vbox[3];
  ErrorCode result = create_3dtri_3_sequences(scdi, int1, int2, int3, int4, vbox, ebox);
  CHECK_ERR(result);

    // whew; that's done; now check and evaluate

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    result = check_element_sequence(ebox[i], ebox[i]->box_min(), ebox[i]->box_max(), 
                                    ebox[i]->start_element());
    CHECK_ERR(result);
  
    result = evaluate_element_sequence(ebox[i]);
    CHECK_ERR(result);
  }
  
  for (int i = 0; i < 3; i++) {
    result = access_adjacencies(ebox[i]);
    CHECK_ERR(result);
  }

  return result;
}

void test_periodic_seq() 
{
  Core moab;
  ScdInterface *scdi;
  ErrorCode rval = moab.Interface::query_interface(scdi);
  CHECK_ERR(rval);
  HomCoord TEST_MIN_PARAMS(0,0,0);
  HomCoord TEST_BOX_MAX(11,5,2);

    // periodic in i
  ScdBox *new_box;
  rval = scdi->construct_box(TEST_MIN_PARAMS, TEST_BOX_MAX, NULL, 0, new_box, true, false);
  CHECK_ERR(rval);
  rval = evaluate_element_sequence(new_box);
  CHECK_ERR(rval);
  
    // periodic in j
  rval = scdi->construct_box(TEST_MIN_PARAMS, TEST_BOX_MAX, NULL, 0, new_box, false, true);
  CHECK_ERR(rval);
  rval = evaluate_element_sequence(new_box);
  CHECK_ERR(rval);
  
    // periodic in i and j
  rval = scdi->construct_box(TEST_MIN_PARAMS, TEST_BOX_MAX, NULL, 0, new_box, true, true);
  CHECK_ERR(rval);
  rval = evaluate_element_sequence(new_box);
  CHECK_ERR(rval);

    // 2d, periodic in i
  TEST_BOX_MAX[2] = 0;
  rval = scdi->construct_box(TEST_MIN_PARAMS, TEST_BOX_MAX, NULL, 0, new_box, true, false);
  CHECK_ERR(rval);
  rval = evaluate_element_sequence(new_box);
  CHECK_ERR(rval);
  
    // 2d, periodic in j
  rval = scdi->construct_box(TEST_MIN_PARAMS, TEST_BOX_MAX, NULL, 0, new_box, false, true);
  CHECK_ERR(rval);
  rval = evaluate_element_sequence(new_box);
  CHECK_ERR(rval);
  
    // 2d, periodic in i and j
  rval = scdi->construct_box(TEST_MIN_PARAMS, TEST_BOX_MAX, NULL, 0, new_box, true, true);
  CHECK_ERR(rval);
  rval = evaluate_element_sequence(new_box);
  CHECK_ERR(rval);
}

ErrorCode create_1d_3_sequences(ScdInterface *scdi,
                                HomCoord tmp_min, HomCoord tmp_max,
                                ScdBox **vbox, ScdBox **ebox) 
{
    // first get 1d min/max by resetting j, k components
  tmp_min[1] = tmp_min[2] = tmp_max[1] = tmp_max[2] = 0;

    // split the sequence in three in i parameter, at 1/2, 1/3 and 1/6, such that the 
    // total # vertices is the expected number from tmp_max - tmp_min)
  int idiff = (tmp_max[0] - tmp_min[0] + 1)/6;
  HomCoord vseq0_minmax[2] = {HomCoord(tmp_min), HomCoord(tmp_min[0]+3*idiff-1, tmp_max[1], tmp_max[2])};
  HomCoord vseq1_minmax[2] = {HomCoord(tmp_min), HomCoord(tmp_min[0]+2*idiff-1, tmp_max[1], tmp_max[2])};
  HomCoord vseq2_minmax[2] = {HomCoord(tmp_min), HomCoord(tmp_max[0]-5*idiff, tmp_max[1], tmp_max[2])};
  
    // create three vertex sequences
  vbox[0] = vbox[1] = vbox[2] = NULL;
  

    // first vertex sequence 
  ErrorCode result = scdi->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
                                               MBVERTEX, 1, vbox[0]);
  CHECK_ERR(result);

    // second vertex sequence 
  result = scdi->create_scd_sequence(vseq1_minmax[0], vseq1_minmax[1],
                                     MBVERTEX, 1, vbox[1]);
  CHECK_ERR(result);

    // third vertex sequence 
  result = scdi->create_scd_sequence(vseq2_minmax[0], vseq2_minmax[1],
                                     MBVERTEX, 1, vbox[2]);
  CHECK_ERR(result);

    // now create the three element sequences
  ebox[0] = ebox[1] = ebox[2] = NULL;
  
    // create the first element sequence
  result = scdi->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
                                     MBEDGE, 1, ebox[0]);
  CHECK_ERR(result);
  
    // add first vertex seq to first element seq, forward orientation, unity transform
  result = ebox[0]->add_vbox(vbox[0],
                             vseq0_minmax[0], vseq0_minmax[0],
                             vseq0_minmax[1], vseq0_minmax[1],
                             vseq0_minmax[0], vseq0_minmax[0]);
  CHECK_ERR(result);
  
    // create the second element sequence; make it use the second vseq in reverse and start
    // with the end vertex of the first sequence; parameterize it such that it tacks onto the 
    // end of the previous ebox
  result = scdi->create_scd_sequence(HomCoord(vseq0_minmax[1].i(), 0, 0),
                                     HomCoord(1+vseq0_minmax[1].i()+vseq1_minmax[1].i()-vseq1_minmax[0].i(), 0, 0),
                                     MBEDGE, 1, ebox[1]);
  CHECK_ERR(result);
  
    // add shared vertex from first vseq to this ebox; parameter space should be the same since
    // we're adding to that parameter space
  result = ebox[1]->add_vbox(vbox[0],
                             vseq0_minmax[0], vseq0_minmax[0],
                             vseq0_minmax[1], vseq0_minmax[1],
                             vseq0_minmax[0], vseq0_minmax[0],
                             true,
                             HomCoord(ebox[1]->box_min().i(),0,0),
                             HomCoord(ebox[1]->box_min().i(),0,0));
  CHECK_ERR(result);
  
    // add second vseq to this ebox, but reversed; parameter space should be such that the
    // last vertex in the second vseq occurs first, and so on
  result = ebox[1]->add_vbox(vbox[1],
                             vseq1_minmax[1], ebox[1]->box_min()+HomCoord::unitv[0],
                             vseq1_minmax[0], ebox[1]->box_max(),
                             vseq1_minmax[1], ebox[1]->box_min()+HomCoord::unitv[0]);
  CHECK_ERR(result);
  
    // create the third element sequence; make it use the third vseq (forward sense) and start
    // with the start vertex of the second vseq; parameterize it such that it tacks onto the 
    // end of the previous ebox
  result = scdi->create_scd_sequence(ebox[1]->box_max(),
                                     HomCoord(ebox[1]->box_max().i()+1+vseq2_minmax[1].i()-vseq2_minmax[0].i(),0,0),
                                     MBEDGE, 1, ebox[2]);
  CHECK_ERR(result);
  
    // add shared vertex from second vseq to this ebox; parameter space mapping such that we get
    // first vertex only of that vseq
  result = ebox[2]->add_vbox(vbox[1],
                             vseq0_minmax[0], ebox[2]->box_min(),
                             vseq0_minmax[0]+HomCoord::unitv[0], ebox[2]->box_min()-HomCoord::unitv[0],
                             vseq0_minmax[0], ebox[2]->box_min(),
                             true,
                             ebox[2]->box_min(),
                             ebox[2]->box_min());
  CHECK_ERR(result);
  
    // add third vseq to this ebox, forward orientation
  result = ebox[2]->add_vbox(vbox[2],
                             vseq2_minmax[0], ebox[2]->box_min()+HomCoord::unitv[0],
                             vseq2_minmax[1], ebox[2]->box_max(),
                             vseq1_minmax[0], ebox[2]->box_min()+HomCoord::unitv[0]);
  CHECK_ERR(result);

  return result;
}

ErrorCode create_2d_3_sequences(ScdInterface *scdi,
                                ScdBox **vbox, ScdBox **ebox) 
{
    // create 3 rectangular sequences attached end to end and back (periodic); sequences are 
    // assorted orientations, sequences have globally-consistent (periodic in i) parameter space

    // set vbox parametric spaces directly
  HomCoord vbox0_minmax[2] = {HomCoord(0,0,0), HomCoord(5,5,0)};
  HomCoord vbox1_minmax[2] = {HomCoord(-2,4,0), HomCoord(8,9,0)};
  HomCoord vbox2_minmax[2] = {HomCoord(0,0,0), HomCoord(8,5,0)};

    // create three vertex sequences
  vbox[0] = vbox[1] = vbox[2] = NULL;

    // first vertex sequence 
  ErrorCode result = scdi->create_scd_sequence(vbox0_minmax[0], vbox0_minmax[1],
                                               MBVERTEX, 1, vbox[0]);
  CHECK_ERR(result);

    // second vertex sequence 
  result = scdi->create_scd_sequence(vbox1_minmax[0], vbox1_minmax[1],
                                     MBVERTEX, 1, vbox[1]);
  CHECK_ERR(result);

    // third vertex sequence 
  result = scdi->create_scd_sequence(vbox2_minmax[0], vbox2_minmax[1],
                                     MBVERTEX, 1, vbox[2]);
  CHECK_ERR(result);

    // now create the three element sequences
  ebox[0] = ebox[1] = ebox[2] = NULL;
  
    // create the first element sequence
  result = scdi->create_scd_sequence(vbox0_minmax[0], vbox0_minmax[1],
                                     MBQUAD, 1, ebox[0]);
  CHECK_ERR(result);
  
    // add first vertex seq to first element seq, forward orientation, unity transform
  result = ebox[0]->add_vbox(vbox[0],
                               // p1: imin,jmin
                             vbox0_minmax[0], vbox0_minmax[0],
                               // p2: imax,jmin
                             HomCoord(vbox0_minmax[1].i(), vbox0_minmax[0].j(),0),
                             HomCoord(vbox0_minmax[1].i(), vbox0_minmax[0].j(),0),
                               // p3: imin,jmax
                             HomCoord(vbox0_minmax[0].i(), vbox0_minmax[1].j(),0),
                             HomCoord(vbox0_minmax[0].i(), vbox0_minmax[1].j(),0));
  
  CHECK_ERR(result);
  
    // create the second element sequence; make it use the right side of the first vbox; 
    // parameterize it such that it tacks onto imax of the previous ebox
  result = scdi->create_scd_sequence(HomCoord(ebox[0]->box_max().i(), ebox[0]->box_min().j(), 0),
                                     HomCoord(vbox0_minmax[1].i()+1+vbox1_minmax[1].i()-vbox1_minmax[0].i(), 
                                              ebox[0]->box_max().j(), 0),
                                     MBQUAD, 1, ebox[1]);
  CHECK_ERR(result);
  
    // add shared side from first vbox to this ebox; parameter space should be the same since
    // we're adding to that parameter space
  result = ebox[1]->add_vbox(vbox[0],
                               // p1: lower right of box 0, lower left of box 1
                             HomCoord(vbox0_minmax[1].i(), vbox0_minmax[0].j(), 0), ebox[1]->box_min(),
                               // p2: one up from p1
                             HomCoord(vbox0_minmax[1].i(), vbox0_minmax[0].j(), 0)+HomCoord::unitv[1], 
                             ebox[1]->box_min()+HomCoord::unitv[1],
                               // p3: one right of p1
                             HomCoord(vbox0_minmax[1].i(), vbox0_minmax[0].j(), 0)+HomCoord::unitv[0], 
                             ebox[1]->box_min()+HomCoord::unitv[0],
                               // set bb such that it's the right side of the vbox, left of local ebox
                             true,
                             ebox[1]->box_min(),
                             HomCoord(ebox[1]->box_min().i(),ebox[1]->box_max().j(),0));
  CHECK_ERR(result);
  
    // add second vbox to this ebox, with different orientation but all of it (no bb input)
  result = ebox[1]->add_vbox(vbox[1],
                               // p1: one right of top left of ebox1
                             vbox1_minmax[0], 
                             HomCoord(ebox[1]->box_min().i()+1, ebox[1]->box_max().j(), 0),
                               // p2: one right from p1
                             vbox1_minmax[0]+HomCoord::unitv[0], 
                             HomCoord(ebox[1]->box_min().i()+2, ebox[1]->box_max().j(), 0),
                               // p3: one down from p1
                             vbox1_minmax[0]+HomCoord::unitv[1], 
                             HomCoord(ebox[1]->box_min().i()+1, ebox[1]->box_max().j()-1, 0));
  CHECK_ERR(result);
  
    // create the third element sequence; make it use the third vbox (middle) as well as the side of the
    // second sequence (left) and a side of the first sequence (right); parameterize it such that it tacks onto the 
    // end of the previous ebox and the beginning of the 1st sequence (i.e. periodic in i)
  result = scdi->create_scd_sequence(HomCoord(ebox[1]->box_max().i(), ebox[1]->box_min().j(),0),
                                       // add one extra for each of left and right sides
                                     HomCoord(ebox[1]->box_max().i()+1+vbox2_minmax[1].i()-vbox2_minmax[0].i()+1,
                                              ebox[1]->box_max().j(),0),
                                     MBEDGE, 1, ebox[2]);
  CHECK_ERR(result);
  
    // add shared side from second vbox to this ebox; parameter space mapping such that we get
    // a side only of that vbox
  result = ebox[2]->add_vbox(vbox[1],
                               // p1: bottom left
                             vbox1_minmax[1], ebox[2]->box_min(),
                               // p2: one right from p1
                             vbox1_minmax[1]+HomCoord::unitv[0],
                             ebox[2]->box_min()+HomCoord::unitv[0],
                               // p3: one up
                             vbox1_minmax[1]-HomCoord::unitv[1],
                             ebox[2]->box_min()+HomCoord::unitv[1],
                               // bb input such that we only get left side of ebox parameter space
                             true,
                             ebox[2]->box_min(),
                             HomCoord(ebox[2]->box_min().i(), ebox[2]->box_max().j(), 0));
  CHECK_ERR(result);
  
    // add shared side from first vbox to this ebox; parameter space mapping such that we get
    // a side only of that vbox
  result = ebox[2]->add_vbox(vbox[0],
                               // p1: bottom right
                             vbox0_minmax[0], HomCoord(ebox[2]->box_max().i(), ebox[2]->box_min().j(),0),
                               // p2: one right from p1
                             vbox0_minmax[0]+HomCoord::unitv[0], 
                             HomCoord(ebox[2]->box_max().i()+1, ebox[2]->box_min().j(),0),
                               // p3: one up from p1
                             vbox0_minmax[0]+HomCoord::unitv[1],
                             HomCoord(ebox[2]->box_max().i(), ebox[2]->box_min().j()+1,0),
                               // bb input such that we only get left side of ebox parameter space
                             true,
                             HomCoord(ebox[2]->box_max().i(), ebox[2]->box_min().j(),0),
                             ebox[2]->box_max());
  CHECK_ERR(result);
  

    // add third vbox to this ebox
  result = ebox[2]->add_vbox(vbox[2],
                               // p1: top right and left one
                             vbox2_minmax[0], ebox[2]->box_max()-HomCoord::unitv[0],
                               // p2: one left of p1
                             vbox2_minmax[0]+HomCoord::unitv[0], ebox[2]->box_max()-HomCoord::unitv[0]*2,
                               // p3: one down from p1
                             vbox2_minmax[0]+HomCoord::unitv[1], ebox[2]->box_max()-HomCoord::unitv[0] -
                             HomCoord::unitv[1]);
  CHECK_ERR(result);

  return result;
}

ErrorCode create_2dtri_3_sequences(ScdInterface *scdi,
                                   const int int1, const int int2, const int int3,
                                   ScdBox **vbox, ScdBox **ebox) 
{
    // create 3 rectangular sequences arranged such that the all share a common (tri-valent) corner;
    // orient each region such that its origin is at the tri-valent corner and the k direction is
    // out of the page
    //
    // int1 and int2 controls the i and j intervals in region 0, int3 follows from that.

    // input is 3 interval settings controlling the 3 degrees of freedom on the interfacesp

    // set vbox parametric spaces directly from int1-3
    // use 0-based parameterization on vbox's just for fun, which means we'll have to transform into
    // ebox system
  HomCoord vbox0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,0)};
  HomCoord vbox1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3-1,int1,0)};
  HomCoord vbox2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2-1,int3-1,0)};

    // create three vertex sequences
  vbox[0] = vbox[1] = vbox[2] = NULL;

    // first vertex sequence 
  ErrorCode result = scdi->create_scd_sequence(vbox0_minmax[0], vbox0_minmax[1],
                                               MBVERTEX, 1, vbox[0]);
  CHECK_ERR(result);

    // second vertex sequence 
  result = scdi->create_scd_sequence(vbox1_minmax[0], vbox1_minmax[1],
                                     MBVERTEX, 1, vbox[1]);
  CHECK_ERR(result);

    // third vertex sequence 
  result = scdi->create_scd_sequence(vbox2_minmax[0], vbox2_minmax[1],
                                     MBVERTEX, 1, vbox[2]);
  CHECK_ERR(result);


    // now create the three element sequences

    // set ebox parametric spaces directly from int1-3
    // use 0-based parameterization on ebox's just for fun, which means we'll have to transform into
    // ebox system
  HomCoord ebox0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,0)};
  HomCoord ebox1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3,int1,0)};
  HomCoord ebox2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2,int3,0)};

  ebox[0] = ebox[1] = ebox[2] = NULL;
  
    // create the first element sequence
  result = scdi->create_scd_sequence(ebox0_minmax[0], ebox0_minmax[1], 
                                     MBQUAD, 1, ebox[0]);
  CHECK_ERR(result);
  
    // only need to add one vbox to this, unity transform
  result = ebox[0]->add_vbox(vbox[0],
                               // trick: if I know it's going to be unity, just input
                               // 3 sets of equivalent points
                             vbox0_minmax[0], vbox0_minmax[0],
                             vbox0_minmax[0], vbox0_minmax[0],
                             vbox0_minmax[0], vbox0_minmax[0]);
  
  CHECK_ERR(result);
  
    // create the second element sequence
  result = scdi->create_scd_sequence(ebox1_minmax[0], ebox1_minmax[1],
                                     MBQUAD, 1, ebox[1]);
  CHECK_ERR(result);
  
    // add shared side from first vbox to this ebox, with bb to get just the line
  result = ebox[1]->add_vbox(vbox[0],
                               // p1: origin in both systems
                             vbox0_minmax[0], ebox0_minmax[0],
                               // p2: one unit along the shared line (i in one, j in other)
                             vbox0_minmax[0]+HomCoord::unitv[0], ebox0_minmax[0]+HomCoord::unitv[1],
                               // p3: arbitrary
                             vbox0_minmax[0], ebox0_minmax[0],
                               // set bb such that it's the jmin side of vbox
                             true,
                             ebox[1]->box_min(),
                             HomCoord(ebox[1]->box_min().i(),ebox[1]->box_max().j(),0));
  CHECK_ERR(result);
  
    // add second vbox to this ebox, with different orientation but all of it (no bb input)
  result = ebox[1]->add_vbox(vbox[1],
                               // p1: origin/i+1 (vbox/ebox)
                             vbox1_minmax[0], ebox1_minmax[0]+HomCoord::unitv[0],
                               // p2: j+1 from p1
                             vbox1_minmax[0]+HomCoord::unitv[1], 
                             ebox1_minmax[0]+HomCoord::unitv[0]+HomCoord::unitv[1],
                               // p3: i+1 from p1
                             vbox1_minmax[0]+HomCoord::unitv[0], 
                             ebox[1]->box_min()+HomCoord::unitv[0]*2);
  CHECK_ERR(result);
  
    // create the third element sequence
  result = scdi->create_scd_sequence(ebox2_minmax[0], ebox2_minmax[1],
                                     MBQUAD, 1, ebox[2]);
  CHECK_ERR(result);
  
    // add shared side from second vbox to this ebox
  result = ebox[2]->add_vbox(vbox[1],
                               // p1: origin/j+1 (vbox/ebox)
                             vbox1_minmax[0], ebox[2]->box_min()+HomCoord::unitv[1],
                               // p2: i+1/j+2 (vbox/ebox)
                             vbox1_minmax[0]+HomCoord::unitv[0],
                             ebox[2]->box_min()+HomCoord::unitv[1]*2,
                               // p3: arbitrary
                             vbox1_minmax[0], ebox[2]->box_min()+HomCoord::unitv[1],
                               // bb input such that we only get one side of ebox parameter space
                             true,
                             ebox[2]->box_min()+HomCoord::unitv[1],
                             HomCoord(ebox[2]->box_min().i(), ebox[2]->box_max().j(), 0));
  CHECK_ERR(result);
  
    // add shared side from first vbox to this ebox
  result = ebox[2]->add_vbox(vbox[0],
                               // p1: origin/origin
                             vbox1_minmax[0], ebox2_minmax[0],
                               // p2: j+1/i+1
                             vbox1_minmax[0]+HomCoord::unitv[1], 
                             ebox2_minmax[0]+HomCoord::unitv[0], 
                               // p3: arbitrary
                             vbox1_minmax[0], ebox2_minmax[0],
                               // bb input such that we only get one side of ebox parameter space
                             true,
                             ebox2_minmax[0],
                             HomCoord(ebox2_minmax[1].i(), ebox2_minmax[0].j(),0));
  CHECK_ERR(result);

    // add third vbox to this ebox
  result = ebox[2]->add_vbox(vbox[2],
                               // p1: origin/i+1,j+1
                             vbox2_minmax[0], ebox[2]->box_min()+HomCoord::unitv[0]+HomCoord::unitv[1],
                               // p2: i+1 from p1
                             vbox2_minmax[0]+HomCoord::unitv[0], ebox[2]->box_min()+HomCoord::unitv[0]*2+HomCoord::unitv[1],
                               // p3: j+1 from p1
                             vbox2_minmax[0]+HomCoord::unitv[1], ebox[2]->box_min()+HomCoord::unitv[0]+HomCoord::unitv[1]*2);
  CHECK_ERR(result);

  return result;
}

ErrorCode create_3dtri_3_sequences(ScdInterface *scdi,
                                   const int int1, const int int2, const int int3, const int int4,
                                   ScdBox **vbox, ScdBox **ebox) 
{
    // create 3 brick sequences arranged such that the all share a common (tri-valent) edge;
    // orient each region similarly to the 2dtri_3_sequences test problem, swept into 3d in the 
    // positive k direction.  This direction is divided into int4 intervals
    //
    // int1 and int2 controls the i and j intervals in region 0, int3 follows from that; int4
    // divides the k axis

    // input is 4 interval settings controlling the 4 degrees of freedom on the interfacesp

    // set vbox parametric spaces directly from int1-4
    // use 0-based parameterization on vbox's just for fun, which means we'll have to transform into
    // ebox system
  HomCoord vbox0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,int4)};
  HomCoord vbox1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3-1,int1,int4)};
  HomCoord vbox2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2-1,int3-1,int4)};

    // create three vertex sequences
  vbox[0] = vbox[1] = vbox[2] = NULL;

    // first vertex sequence 
  ErrorCode result = scdi->create_scd_sequence(vbox0_minmax[0], vbox0_minmax[1],
                                               MBVERTEX, 1,
                                               vbox[0]);
  CHECK_ERR(result);

    // second vertex sequence 
  result = scdi->create_scd_sequence(vbox1_minmax[0], vbox1_minmax[1],
                                     MBVERTEX, 1,
                                     vbox[1]);
  CHECK_ERR(result);

    // third vertex sequence 
  result = scdi->create_scd_sequence(vbox2_minmax[0], vbox2_minmax[1],
                                     MBVERTEX, 1,
                                     vbox[2]);
  CHECK_ERR(result);


    // now create the three element sequences

    // set ebox parametric spaces directly from int1-4
    // use 0-based parameterization on ebox's just for fun, which means we'll have to transform into
    // ebox system
  HomCoord ebox0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,int4)};
  HomCoord ebox1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3,int1,int4)};
  HomCoord ebox2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2,int3,int4)};

  ebox[0] = ebox[1] = ebox[2] = NULL;
  
    // create the first element sequence
  result = scdi->create_scd_sequence(ebox0_minmax[0], ebox0_minmax[1], 
                                     MBHEX, 1,
                                     ebox[0]);
  CHECK_ERR(result);
  
    // only need to add one vbox to this, unity transform
  result = ebox[0]->add_vbox(vbox[0],
                               // trick: if I know it's going to be unity, just input
                               // 3 sets of equivalent points
                             vbox0_minmax[0], vbox0_minmax[0],
                             vbox0_minmax[0], vbox0_minmax[0],
                             vbox0_minmax[0], vbox0_minmax[0]);
  
  CHECK_ERR(result);
  
    // create the second element sequence
  result = scdi->create_scd_sequence(ebox1_minmax[0], ebox1_minmax[1],
                                     MBHEX, 1,
                                     ebox[1]);
  CHECK_ERR(result);
  
    // add shared side from first vbox to this ebox, with bb to get just the face
  result = ebox[1]->add_vbox(vbox[0],
                               // p1: origin in both systems
                             vbox0_minmax[0], ebox1_minmax[0],
                               // p2: one unit along the shared line (i in one, j in other)
                             vbox0_minmax[0]+HomCoord::unitv[0], ebox1_minmax[0]+HomCoord::unitv[1],
                               // p3: +k in both (not arbitrary, since interface is 2d)
                             vbox0_minmax[0]+HomCoord::unitv[2], ebox1_minmax[0]+HomCoord::unitv[2],
                               // set bb such that it's the jmin side of vbox
                             true,
                             ebox[1]->box_min(),
                             HomCoord(ebox[1]->box_min().i(),ebox[1]->box_max().j(),
                                      ebox[1]->box_max().k()));
  CHECK_ERR(result);
  
    // add second vbox to this ebox, with different orientation but all of it (no bb input)
  result = ebox[1]->add_vbox(vbox[1],
                               // p1: origin/i+1 (vbox/ebox)
                             vbox1_minmax[0], ebox1_minmax[0]+HomCoord::unitv[0],
                               // p2: j+1 from p1
                             vbox1_minmax[0]+HomCoord::unitv[1], 
                             ebox1_minmax[0]+HomCoord::unitv[0]+HomCoord::unitv[1],
                               // p3: i+1 from p1
                             vbox1_minmax[0]+HomCoord::unitv[0], 
                             ebox[1]->box_min()+HomCoord::unitv[0]*2);
  CHECK_ERR(result);
  
    // create the third element sequence
  result = scdi->create_scd_sequence(ebox2_minmax[0], ebox2_minmax[1],
                                     MBHEX, 1,
                                     ebox[2]);
  CHECK_ERR(result);
  
    // add shared side from second vbox to this ebox
  result = ebox[2]->add_vbox(vbox[1],
                               // p1: origin/j+1 (vbox/ebox)
                             vbox1_minmax[0], ebox[2]->box_min()+HomCoord::unitv[1],
                               // p2: i+1/j+2 (vbox/ebox)
                             vbox1_minmax[0]+HomCoord::unitv[0],
                             ebox[2]->box_min()+HomCoord::unitv[1]*2,
                               // p3: +k in both (not arbitrary, since interface is 2d)
                             vbox1_minmax[0]+HomCoord::unitv[2], ebox[2]->box_min()+HomCoord::unitv[1]+HomCoord::unitv[2],
                               // bb input such that we only get one side of ebox parameter space
                             true,
                             ebox[2]->box_min()+HomCoord::unitv[1],
                             HomCoord(ebox[2]->box_min().i(), ebox[2]->box_max().j(), 
                                      ebox[2]->box_max().k()));
  CHECK_ERR(result);
  
    // add shared side from first vbox to this ebox
  result = ebox[2]->add_vbox(vbox[0],
                               // p1: origin/origin
                             vbox0_minmax[0], ebox2_minmax[0],
                               // p2: j+1/i+1
                             vbox0_minmax[0]+HomCoord::unitv[1], 
                             ebox2_minmax[0]+HomCoord::unitv[0], 
                               // p3: +k in both (not arbitrary, since interface is 2d)
                             vbox0_minmax[0]+HomCoord::unitv[2], ebox[2]->box_min()+HomCoord::unitv[2],
                               // bb input such that we only get one side of ebox parameter space
                             true,
                             ebox2_minmax[0],
                             HomCoord(ebox2_minmax[1].i(), ebox2_minmax[0].j(),
                                      ebox2_minmax[1].k()));
  CHECK_ERR(result);

    // add third vbox to this ebox
  result = ebox[2]->add_vbox(vbox[2],
                               // p1: origin/i+1,j+1
                             vbox2_minmax[0], ebox[2]->box_min()+HomCoord::unitv[0]+HomCoord::unitv[1],
                               // p2: i+1 from p1
                             vbox2_minmax[0]+HomCoord::unitv[0], ebox[2]->box_min()+HomCoord::unitv[0]*2+HomCoord::unitv[1],
                               // p3: j+1 from p1
                             vbox2_minmax[0]+HomCoord::unitv[1], ebox[2]->box_min()+HomCoord::unitv[0]+HomCoord::unitv[1]*2);
  CHECK_ERR(result);

  return result;
}

ErrorCode access_adjacencies(ScdBox *box) 
{
    // access the adjacencies in this box in a few places
  HomCoord box_size = box->box_size(), box_min = box->box_min(), box_max = box->box_max();
  
  EntityHandle dum_entity;
  const EntityHandle *connect;
  int num_connect;
  ErrorCode rval;
  bool is_2d = (box_size.k() <= 1);

    // edges first; bottom:
  for (int dir = 0; dir < 3; dir++) {
      // don't do 3rd direction for 2d box
    if (2 == dir && is_2d) continue;
    
    rval = box->get_adj_edge_or_face(1, box_min.i(), box_min.j(), box_min.k(), dir, dum_entity, true);
    if (MB_SUCCESS != rval) return rval;
  
      // do a simple API call on that entity to make sure we can
    rval = box->sc_impl()->impl()->get_connectivity(dum_entity, connect, num_connect);
    if (MB_SUCCESS != rval) return rval;
  }
    // middle:
  for (int dir = 0; dir < 3; dir++) {
      // don't do 3rd direction for 2d box
    if (2 == dir && is_2d) continue;
    
    rval = box->get_adj_edge_or_face(1, box_min.i()+.5*box_size.i(), 
                                     box_min.j()+.5*box_size.j(), box_min.k()+.5*box_size.k(), 
                                     dir, dum_entity, true);
    if (MB_SUCCESS != rval) return rval;
  
      // do a simple API call on that entity to make sure we can
    rval = box->sc_impl()->impl()->get_connectivity(dum_entity, connect, num_connect);
    if (MB_SUCCESS != rval) return rval;
  }
  
    // top:
  for (int dir = 0; dir < 3; dir++) {
      // don't do 3rd direction for 2d box
    if (2 == dir && is_2d) continue;
    
    rval = box->get_adj_edge_or_face(1, 
                                     (box_max.i() == box_min.i() ? box_max.i() : box_max.i()-1),
                                     (box_max.j() == box_min.j() ? box_max.j() : box_max.j()-1),
                                     (box_max.k() == box_min.k() ? box_max.k() : box_max.k()-1),
                                     dir, dum_entity, true);
    if (MB_SUCCESS != rval) return rval;
  
      // do a simple API call on that entity to make sure we can
    rval = box->sc_impl()->impl()->get_connectivity(dum_entity, connect, num_connect);
    if (MB_SUCCESS != rval) return rval;
  }
  
  if (is_2d) return MB_SUCCESS;

    // now faces; bottom:
  for (int dir = 0; dir < 3; dir++) {
    rval = box->get_adj_edge_or_face(2, box_min.i(), box_min.j(), box_min.k(), dir, dum_entity, true);
    if (MB_SUCCESS != rval) return rval;
  
      // do a simple API call on that entity to make sure we can
    rval = box->sc_impl()->impl()->get_connectivity(dum_entity, connect, num_connect);
    if (MB_SUCCESS != rval) return rval;
  }
    // middle:
  for (int dir = 0; dir < 3; dir++) {
      // don't do 3rd direction for 2d box
    if (2 == dir && is_2d) continue;
    
    rval = box->get_adj_edge_or_face(2, box_min.i()+.5*box_size.i(), 
                                     box_min.j()+.5*box_size.j(), box_min.k()+.5*box_size.k(), 
                                     dir, dum_entity, true);
    if (MB_SUCCESS != rval) return rval;
  
      // do a simple API call on that entity to make sure we can
    rval = box->sc_impl()->impl()->get_connectivity(dum_entity, connect, num_connect);
    if (MB_SUCCESS != rval) return rval;
  }
  
    // top:
  for (int dir = 0; dir < 3; dir++) {
      // don't do 3rd direction for 2d box
    if (2 == dir && is_2d) continue;
    
    rval = box->get_adj_edge_or_face(2, 
                                     (box_max.i() == box_min.i() ? box_max.i() : box_max.i()-1),
                                     (box_max.j() == box_min.j() ? box_max.j() : box_max.j()-1),
                                     (box_max.k() == box_min.k() ? box_max.k() : box_max.k()-1),
                                     dir, dum_entity, true);
    if (MB_SUCCESS != rval) return rval;
  
      // do a simple API call on that entity to make sure we can
    rval = box->sc_impl()->impl()->get_connectivity(dum_entity, connect, num_connect);
    if (MB_SUCCESS != rval) return rval;
  }

  return MB_SUCCESS;
}

void test_parallel_partitions() 
{
  
  Core moab;
  ScdInterface *scdi;
  ErrorCode rval = moab.Interface::query_interface(scdi);
  CHECK_ERR(rval);
  int gdims[] = {0, 0, 0, 48, 40, 18};

    // test for various numbers of procs, powers of two
  int maxpow = 10;
#ifndef NDEBUG
  maxpow = 4;
#endif  
  for (int exp = 2; exp <= maxpow; exp += 2) {
    int nprocs = 0.1 + pow(2.0, (double)exp);
  
    // alljorkori
    ErrorCode rval = test_parallel_partition(gdims, nprocs, ScdInterface::ALLJORKORI);
    CHECK_ERR(rval);
    
    // alljkbal
    rval = test_parallel_partition(gdims, nprocs, ScdInterface::ALLJKBAL);
    CHECK_ERR(rval);
    
    // sqij
    rval = test_parallel_partition(gdims, nprocs, ScdInterface::SQIJ);
    CHECK_ERR(rval);
    
    // sqjk
    rval = test_parallel_partition(gdims, nprocs, ScdInterface::SQJK);
    CHECK_ERR(rval);
  }
}

ErrorCode test_parallel_partition(int *gdims, int nprocs, int part_method) 
{
  ErrorCode rval;
  int pto, bdy_ind_a[2], bdy_ind_b[2], rdims_a[6], rdims_b[6], facedims_a[6], facedims_b[6], ldims[6],
      pfrom;
  for (int p = 0; p < nprocs/2; p++) {
    rval = ScdInterface::compute_partition(part_method, nprocs, p,
                                           gdims, ldims);
    if (MB_SUCCESS != rval) continue;
    
    for (int k = -1; k <= 1; k++) {
      for (int j = -1; j <= 1; j++) {
        for (int i = -1; i <= 1; i++) {
          rval = ScdInterface::get_neighbor(p, nprocs, part_method,
                                            gdims, ldims,
                                            false, false,
                                            i, j, k, pto, bdy_ind_a, rdims_a, facedims_a);
          if (MB_SUCCESS != rval) return rval;
          if (-1 == pto) continue;

          if (facedims_a[0] < rdims_a[0] || facedims_a[0] > rdims_a[3] ||
              facedims_a[1] < rdims_a[1] || facedims_a[1] > rdims_a[4] ||
              facedims_a[2] < rdims_a[2] || facedims_a[2] > rdims_a[5]) CHECK_ERR(MB_FAILURE);
          
          
            // non-negative value of pto; check corresponding input from that proc to this
          rval = ScdInterface::get_neighbor(pto, nprocs, part_method,
                                            gdims, rdims_a, false, false, 
                                            -1*i, -1*j, -1*k, pfrom, bdy_ind_b, rdims_b, facedims_b);
          if (MB_SUCCESS != rval) return rval;
          for (int ind = 0; ind < 3; ind++)
            if (facedims_a[ind] < rdims_b[ind] || facedims_b[ind] > rdims_b[ind+3]) CHECK_ERR(MB_FAILURE);
          for (int ind = 0; ind < 6; ind++) {
            if (facedims_a[ind] != facedims_b[ind]) CHECK_ERR(MB_FAILURE);
            if (rdims_b[ind] != ldims[ind]) CHECK_ERR(MB_FAILURE);
          }
          
          if (bdy_ind_a[0] != bdy_ind_b[0] || bdy_ind_a[1] != bdy_ind_b[1]) CHECK_ERR(MB_FAILURE);
        } // i
      } // j
    } // k
  } // p

  return MB_SUCCESS;
}

