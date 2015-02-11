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

// Cubit performance tests building mapped mesh with CubitNodes
// and CubitHexes.

// Different platforms follow different conventions for usage
#ifndef NT
#include <sys/resource.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include "CubitNode.hpp"
#include "NodeHex.hpp"

using namespace moab;

const double LENGTH = 1.0;

extern "C" 
{
  void __ctype_toupper() 
  {}
}

void print_time(const bool print_em, double &tot_time, double &utime, double &stime);

void build_coords(const int nelem, double *&coords) 
{
    // allocate the memory
  int numv = nelem+1;
  int numv_sq = numv*numv;
  int tot_numv = numv*numv*numv;
  coords = new double[3*tot_numv];

  double scale = LENGTH/nelem;
// use FORTRAN-like indexing
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))

  int idx;
  for (int i=0; i < numv; i++) {
    for (int j=0; j < numv; j++) {
      for (int k=0; k < numv; k++) {
        idx = VINDEX(i,j,k);
          // interleaved coordinate ordering
        coords[3*idx] = i*scale;
        coords[3*idx+1] = j*scale;
        coords[3*idx+2] = k*scale;
      }
    }
  }
}

void build_connect(const int nelem, int *&connect) 
{
    // allocate the memory
  int nume_tot = nelem*nelem*nelem;
  connect = new int[8*nume_tot];

  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int idx = 0;
  int vijk;
  for (int i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = VINDEX(i,j,k);
        connect[idx++] = vijk;
        connect[idx++] = vijk+1;
        connect[idx++] = vijk+1+numv;
        connect[idx++] = vijk+numv;
        connect[idx++] = vijk+numv*numv;
        connect[idx++] = vijk+1+numv*numv;
        connect[idx++] = vijk+1+numv+numv*numv;
        connect[idx++] = vijk+numv+numv*numv;
        assert(idx <= 8*nume_tot);
      }
    }
  }
}

int main(int argc, char* argv[])
{
  int nelem = 20;
  if (argc > 1)
  {
    sscanf(argv[1], "%d", &nelem);
  }
  std::cout << "number of elements: " << nelem << std::endl;

  // create a hex mesh in a 1x1x1 cube with nelem number of
  // elements along an edge
  int nnodes = nelem + 1;

  double *coords = NULL;
  int *connect = NULL;

    // build the coordinates
  build_coords(nelem, coords);
  
    // build the connectivity
  build_connect(nelem, connect);
  
  assert(NULL != connect && NULL != coords);

  double ttime0, ttime1, ttime2, ttime3, utime, stime;
  print_time(false, ttime0, utime, stime);
  int nodes_tot = nnodes*nnodes*nnodes;
  CubitNode** node_array = new CubitNode*[nodes_tot];
  int i;
  for (i=0; i < nodes_tot; i++)
    node_array[i] = new CubitNode(coords[3*i], coords[3*i+1], coords[3*i+2]);

  int nelem_tot = nelem*nelem*nelem;
  NodeHex** hex_array = new NodeHex*[nelem_tot];
  int j = 0;
  CubitNode*  conn[8];
  for (i = 0; i < nelem_tot; i++) {
    conn[0] = node_array[connect[j]];
    conn[1] = node_array[connect[j+1]];
    conn[2] = node_array[connect[j+2]];
    conn[3] = node_array[connect[j+3]];
    conn[4] = node_array[connect[j+4]];
    conn[5] = node_array[connect[j+5]];
    conn[6] = node_array[connect[j+6]];
    conn[7] = node_array[connect[j+7]];
    j += 8;
    
    hex_array[i] = new NodeHex(conn);
  }

  print_time(false, ttime1, utime, stime);

    // query element to vertex

  
  for (i = 0; i < nelem_tot; i++) {
    double centroid[3] = {0.0, 0.0, 0.0};
    hex_array[i]->hex_nodes(conn[0], conn[1], conn[2], conn[3], 
                            conn[4], conn[5], conn[6], conn[7]);
    for (j = 0; j < 8; j++) {
      centroid[0] += conn[j]->node_x();
      centroid[1] += conn[j]->node_y();
      centroid[2] += conn[j]->node_z();
    }
    
  }

  print_time(false, ttime2, utime, stime);

    // need to allocate & populate TDHexKnowing for these
  for (i = 0; i < nelem_tot; i++)
    hex_array[i]->add_hex_to_nodes();
    
  DLIList<CubitHex*> hexes;
  for (i = 0; i < nodes_tot; i++) {
    node_array[i]->all_hexes(hexes);
  }
    
  print_time(false, ttime3, utime, stime);

  std::cout << "CUBIT: nelem, construct, e_to_v query, v_to_e query = " 
            << nelem << ", " 
            << ttime1-ttime0 << ", " 
            << ttime2-ttime1 << ", " 
            << ttime3-ttime2 << " seconds" 
            << std::endl;
  return 0;
}

void print_time(const bool print_em, double &tot_time, double &utime, double &stime) 
{
  struct rusage r_usage;
  getrusage(RUSAGE_SELF, &r_usage);
  utime = (double)r_usage.ru_utime.tv_sec +
     ((double)r_usage.ru_utime.tv_usec/1.e6);
  stime = (double)r_usage.ru_stime.tv_sec +
     ((double)r_usage.ru_stime.tv_usec/1.e6);
  tot_time = utime + stime;
  if (print_em)
    std::cout << "User, system, total time = " << utime << ", " << stime 
              << ", " << tot_time << std::endl;
#ifndef LINUX
 std::cout << "Max resident set size = " << r_usage.ru_maxrss*4096 << " bytes" << std::endl;
 std::cout << "Int resident set size = " << r_usage.ru_idrss << std::endl;
#else
  system("ps o args,drs,rss | grep perf | grep -v grep");  // RedHat 9.0 doesnt fill in actual memory data 
#endif
}
