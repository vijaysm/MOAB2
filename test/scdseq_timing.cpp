/**
 * MOAB, a Mesh-Oriented datABase, is a software component for creating,
 * storing and accessing finite element mesh data.
 * 
 * Copyright 2004 Sandia Corporation.  Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 * 
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * 
 */

#include "ScdVertexSeq.hpp"
#include "ScdElementSeq.hpp"
#include "EntitySequenceManager.hpp"
#include "EntitySequence.hpp"
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

#include <iostream>
#include <time.h>

using namespace moab;

/*
// some timing/memory measurement includes; memory measurement
// doesn't work on linux, so they're commented out for now to avoid
// platform problems
#include <unistd.h>
#include <termios.h>
#include <sys/ioctl.h>
#include <ctype.h>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
*/

int create_3dtri_3_sequences(Core *gMB,
                             const int intervals,
                             EntityHandle *vstart, EntityHandle *estart);
int create_3dtri_ucd_sequences(Core *gMB,
                               const int intervals,
                               EntityHandle *vstart, EntityHandle *estart);
void print_time();



int main(int argc, char**argv) 
{
  int errors = 0;

    // first we need to make a new Core
  Core *gMB = new Core();

    // get the intervals
  if (argc < 2) {
    std::cout << "Usage: <scdseq_timing> <#intervals> " << std::endl
              << " where #intervals is the number of intervals on each side of each cube." << std::endl;
    return 0;
  }
  
  int intervals;
  sscanf(argv[1], "%d", &intervals);

  char do_option[8];
  
  bool do_scd = true, do_ucd = true;
  
  if (argc > 2) {
    sscanf(argv[2], "%s", do_option);
    if (do_option[0] == 'u') do_scd = false;
    else if (do_option[0] == 's') do_ucd = false;
    else {
      std::cout << "Didn't understand input; doing both scd and ucd." << std::endl;
    }
  }
  
  EntityHandle estart[3], vstart[3];
  int total_elements = intervals*intervals*intervals;
  std::vector<EntityHandle> connect;
  ErrorCode result;
  clock_t start, stop;
  float time;
  char inp[1];

    // wait for input to get memory reading
  std::cout << "Hit any key and return to continue...";
  std::cin >> inp;
  std::cout << std::endl;

  if (do_scd) {
    
      // create structured mesh
    errors = create_3dtri_3_sequences(gMB, intervals, vstart, estart);
    if (errors != 0) {
      std::cout << "Problem creating structured sequences." << std::endl;
      return errors;
    }

      // get connectivity

    start = clock();
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < total_elements; i++) {
        result = gMB->get_connectivity(estart[j]+i, connect);
        if (MB_SUCCESS != result) break;
        connect.clear();
      }
    }
    stop = clock();
    time = static_cast<float>(stop - start)/CLOCKS_PER_SEC;

    std::cout << "Time to get connectivity for scd mesh of " << 3*total_elements << " elements: " 
              << time << " seconds." << std::endl;

    print_time();
      // wait for input to get memory reading
    std::cout << "Hit any key and return to continue...";
    std::cin >> inp;
    std::cout << std::endl;
  
      // destroy this mesh
    delete gMB;
  }
  
  if (do_ucd) {

      // now do the same thing, only unstructured
    gMB = new Core();

      // create the elements
    errors = create_3dtri_ucd_sequences(gMB, intervals, vstart, estart);
    if (errors != 0) {
      std::cout << "Problem creating unstructured sequences." << std::endl;
      return errors;
    }

      // get connectivity
    std::vector<EntityHandle> connect;
    start = clock();
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < total_elements; i++) {
        result = gMB->get_connectivity(estart[j]+i, connect);
        if (MB_SUCCESS != result) break;
        connect.clear();
      }
    }
    stop = clock();
    time = static_cast<float>(stop - start)/CLOCKS_PER_SEC;

    std::cout << "Time to get connectivity for ucd mesh of " << 3*total_elements << " elements: " 
              << time << " seconds." << std::endl;

    print_time();
      // wait for input to get memory reading
    std::cout << "Hit any key and return to continue...";
    std::cin >> inp;
    std::cout << std::endl;

      // destroy this mesh
    delete gMB;
    
  }
}

int create_3dtri_3_sequences(Core *gMB,
                             const int intervals,
                             EntityHandle *vstart, EntityHandle *estart) 
{
    // create 3 brick esequences arranged such that the all share a common (tri-valent) edge;
    // orient each region similarly to the 2dtri_3_esequences test problem, swept into 3d in the 
    // positive k direction.  This direction is divided into intervals intervals
    //
    // intervals and intervals controls the i and j intervals in region 0, intervals follows from that; intervals
    // divides the k axis

    // input is 4 interval settings controlling the 4 degrees of freedom on the interfacesp
  int errors = 0;

    
  ScdVertexSeq *vseq[3];
  ScdElementSeq *eseq[3];
  
    // set vseq parametric spaces directly from intervals-4
    // use 0-based parameterization on vseq's just for fun, which means we'll have to transform into
    // eseq system
  HomCoord vseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(intervals,intervals,intervals)};
  HomCoord vseq1_minmax[2] = {HomCoord(0,0,0), HomCoord(intervals-1,intervals,intervals)};
  HomCoord vseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(intervals-1,intervals-1,intervals)};

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();

    // create three vertex sequences
  EntitySequence *dum_seq;
  vseq[0] = vseq[1] = vseq[2] = NULL;

    // first vertex sequence 
  ErrorCode result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
                                                     MBVERTEX, 1,
                                                     vstart[0], dum_seq);
  if (NULL != dum_seq) vseq[0] = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && vstart[0] != 0 && dum_seq != NULL && vseq[0] != NULL);

    // second vertex sequence 
  result = seq_mgr->create_scd_sequence(vseq1_minmax[0], vseq1_minmax[1],
                                        MBVERTEX, 1,
                                        vstart[1], dum_seq);
  if (NULL != dum_seq) vseq[1] = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && vstart[1] != 0 && dum_seq != NULL && vseq[1] != NULL);

    // third vertex sequence 
  result = seq_mgr->create_scd_sequence(vseq2_minmax[0], vseq2_minmax[1],
                                        MBVERTEX, 1,
                                        vstart[2], dum_seq);
  if (NULL != dum_seq) vseq[2] = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && vstart[2] != 0 && dum_seq != NULL && vseq[2] != NULL);


    // now create the three element sequences

    // set eseq parametric spaces directly from intervals-4
    // use 0-based parameterization on eseq's just for fun, which means we'll have to transform into
    // eseq system
  HomCoord eseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(intervals,intervals,intervals)};
  HomCoord eseq1_minmax[2] = {HomCoord(0,0,0), HomCoord(intervals,intervals,intervals)};
  HomCoord eseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(intervals,intervals,intervals)};

  eseq[0] = eseq[1] = eseq[2] = NULL;
  
    // create the first element sequence
  result = seq_mgr->create_scd_sequence(eseq0_minmax[0], eseq0_minmax[1], 
                                        MBHEX, 1,
                                        estart[0], dum_seq);
  if (NULL != dum_seq) eseq[0] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[0] != 0 && dum_seq != NULL && eseq[0] != NULL);
  
    // only need to add one vseq to this, unity transform
  result = eseq[0]->add_vsequence(vseq[0],
                                    // trick: if I know it's going to be unity, just input
                                    // 3 sets of equivalent points
                                  vseq0_minmax[0], vseq0_minmax[0],
                                  vseq0_minmax[0], vseq0_minmax[0],
                                  vseq0_minmax[0], vseq0_minmax[0]);
  
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add first vsequence to first element sequence in tri-composite 3d eseq." << std::endl;
    errors++;
  }
  
    // create the second element sequence
  result = seq_mgr->create_scd_sequence(eseq1_minmax[0], eseq1_minmax[1],
                                        MBHEX, 1,
                                        estart[1], dum_seq);
  if (NULL != dum_seq) eseq[1] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[1] != 0 && dum_seq != NULL && eseq[1] != NULL);
  
    // add shared side from first vseq to this eseq, with bb to get just the face
  result = eseq[1]->add_vsequence(vseq[0],
                                    // p1: origin in both systems
                                  vseq0_minmax[0], eseq0_minmax[0],
                                    // p2: one unit along the shared line (i in one, j in other)
                                  vseq0_minmax[0]+HomCoord::unitv[0], eseq0_minmax[0]+HomCoord::unitv[1],
                                    // p3: arbitrary
                                  vseq0_minmax[0], eseq0_minmax[0],
                                    // set bb such that it's the jmin side of vseq
                                  true,
                                  eseq[1]->min_params(),
                                  HomCoord(eseq[1]->min_params().i(),eseq[1]->max_params().j(),
                                           eseq[1]->max_params().k()));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to second element sequence in tri-composite 3d eseq." << std::endl;
    errors++;
  }
  
    // add second vseq to this eseq, with different orientation but all of it (no bb input)
  result = eseq[1]->add_vsequence(vseq[1],
                                    // p1: origin/i+1 (vseq/eseq)
                                  vseq1_minmax[0], eseq1_minmax[0]+HomCoord::unitv[0],
                                    // p2: j+1 from p1
                                  vseq1_minmax[0]+HomCoord::unitv[1], 
                                  eseq1_minmax[0]+HomCoord::unitv[0]+HomCoord::unitv[1],
                                    // p3: i+1 from p1
                                  vseq1_minmax[0]+HomCoord::unitv[0], 
                                  eseq[1]->min_params()+HomCoord::unitv[0]*2);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add second vseq to second element sequence in tri-composite 3d eseq." << std::endl;
    errors++;
  }
  
    // create the third element sequence
  result = seq_mgr->create_scd_sequence(eseq2_minmax[0], eseq2_minmax[1],
                                        MBHEX, 1,
                                        estart[2], dum_seq);
  if (NULL != dum_seq) eseq[2] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[2] != 0 && dum_seq != NULL && eseq[2] != NULL);
  
    // add shared side from second vseq to this eseq
  result = eseq[2]->add_vsequence(vseq[1],
                                    // p1: origin/j+1 (vseq/eseq)
                                  vseq1_minmax[0], eseq[2]->min_params()+HomCoord::unitv[1],
                                    // p2: i+1/j+2 (vseq/eseq)
                                  vseq1_minmax[0]+HomCoord::unitv[0],
                                  eseq[2]->min_params()+HomCoord::unitv[1]*2,
                                    // p3: arbitrary
                                  vseq1_minmax[0], eseq[2]->min_params()+HomCoord::unitv[1],
                                    // bb input such that we only get one side of eseq parameter space
                                  true,
                                  eseq[2]->min_params()+HomCoord::unitv[1],
                                  HomCoord(eseq[2]->min_params().i(), eseq[2]->max_params().j(), 
                                           eseq[2]->max_params().k()));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to third element sequence in tri-composite 3d eseq." << std::endl;
    errors++;
  }
  
    // add shared side from first vseq to this eseq
  result = eseq[2]->add_vsequence(vseq[0],
                                    // p1: origin/origin
                                  vseq1_minmax[0], eseq2_minmax[0],
                                    // p2: j+1/i+1
                                  vseq1_minmax[0]+HomCoord::unitv[1], 
                                  eseq2_minmax[0]+HomCoord::unitv[0], 
                                    // p3: arbitrary
                                  vseq1_minmax[0], eseq2_minmax[0],
                                    // bb input such that we only get one side of eseq parameter space
                                  true,
                                  eseq2_minmax[0],
                                  HomCoord(eseq2_minmax[1].i(), eseq2_minmax[0].j(),
                                           eseq2_minmax[1].k()));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add left shared vsequence to third element sequence in tri-composite 3d eseq." << std::endl;
    errors++;
  }

    // add third vseq to this eseq
  result = eseq[2]->add_vsequence(vseq[2],
                                    // p1: origin/i+1,j+1
                                  vseq2_minmax[0], eseq[2]->min_params()+HomCoord::unitv[0]+HomCoord::unitv[1],
                                    // p2: i+1 from p1
                                  vseq2_minmax[0]+HomCoord::unitv[0], eseq[2]->min_params()+HomCoord::unitv[0]*2+HomCoord::unitv[1],
                                    // p3: j+1 from p1
                                  vseq2_minmax[0]+HomCoord::unitv[1], eseq[2]->min_params()+HomCoord::unitv[0]+HomCoord::unitv[1]*2);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add third vseq to third element sequence in tri-composite 3d eseq." << std::endl;
    errors++;
  }

  return errors;
}

int create_3dtri_ucd_sequences(Core *gMB, const int intervals, 
                               EntityHandle *vstart, EntityHandle *estart) 
{
  
  ReadUtilIface* readMeshIface;
  gMB->query_interface(readMeshIface);
  
  int num_elements = intervals*intervals*intervals;
  int num_verts = (intervals+1)*(intervals+1)*(intervals+1);
  
  std::vector<double*> arrays;
  for (int i = 0; i < 3; i++) {
    readMeshIface->get_node_coords(3, num_verts,
                                   MB_START_ID, vstart[i], arrays);
    arrays.clear();
  }

  EntityHandle *conn[3];

    // allocate 3 arrays to initialize connectivity data
  for (int i = 0; i < 3; i++) {
    readMeshIface->get_element_connect(num_elements, 8, MBHEX, 1, estart[i], conn[i]);

    // now, initialize connectivity data to what it should be; just fudge for now
    for (int j = 0; j < num_elements*8; j++)
      conn[i][j] = vstart[i];
  }

  return 0;
}

void print_time() 
{
/*
  struct rusage r_usage;
  float utime, stime;
  getrusage(RUSAGE_SELF, &r_usage);

  if (r_usage.ru_maxrss == 0) {
      // this machine doesn't return rss - try going to /proc
      // print the file name to open
    char file_str[4096], dum_str[4096];
    int file_ptr = -1, file_len;
    file_ptr = open("/proc/self/stat", O_RDONLY);
    file_len = read(file_ptr, file_str, sizeof(file_str)-1);
    if (file_len == 0) return;
    close(file_ptr);
    file_str[file_len] = '\0';
      // read the preceeding fields and the ones we really want...
    int dum_int;
    unsigned int dum_uint, vm_size, rss;
    static int page_size = getpagesize();
    int num_fields = sscanf(file_str, 
                            "%d " // pid
                            "%s " // comm
                            "%c " // state
                            "%d %d %d %d %d " // ppid, pgrp, session, tty, tpgid
                            "%u %u %u %u %u " // flags, minflt, cminflt, majflt, cmajflt
                            "%d %d %d %d %d %d " // utime, stime, cutime, cstime, counter, priority
                            "%u %u " // timeout, itrealvalue
                            "%d " // starttime
                            "%u %u", // vsize, rss
                            &dum_int, 
                            dum_str, 
                            dum_str, 
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, &dum_uint, &dum_uint, &dum_uint,
                            &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, &dum_int, 
                            &dum_uint, &dum_uint, 
                            &dum_int,
                            &vm_size, &rss);
    if (num_fields == 24) {
      r_usage.ru_maxrss = rss/page_size;
      r_usage.ru_idrss = vm_size/page_size;
    }
  }
  utime = (float)r_usage.ru_utime.tv_sec +
    ((float)r_usage.ru_utime.tv_usec/1.e6);
  stime = (float)r_usage.ru_stime.tv_sec +
    ((float)r_usage.ru_stime.tv_usec/1.e6);
  static int pagesize = getpagesize();

  std::cout << "Execution time = " << utime+stime << ", max RSS = " 
            << r_usage.ru_maxrss*pagesize << " bytes." << std::endl;
  
*/

}
