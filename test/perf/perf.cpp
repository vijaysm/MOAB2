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

// MOAB performance tests building mapped mesh with nodes and
// hexes created one at a time.  This also creates the node to hex adjacencies.
 
// Different platforms follow different conventions for usage
#if !defined(_MSC_VER) && !defined(__MINGW32__)
#include <sys/resource.h>
#else
#include <time.h>
#endif
#ifdef SOLARIS
extern "C" int getrusage(int, struct rusage *);
#ifndef RUSAGE_SELF
#include </usr/ucbinclude/sys/rusage.h>
#endif
#endif

#ifndef IS_BUILDING_MB
#define IS_BUILDING_MB
#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"
#include "VertexSequence.hpp"
#include "StructuredElementSeq.hpp"
#include "EntitySequence.hpp"
#include "SequenceManager.hpp"
#include "moab/HomXform.hpp"
#include "moab/SetIterator.hpp"

using namespace moab;

double LENGTH = 1.0;

void testA(const int nelem, const double *coords);
void testB(const int nelem, const double *coords, int *connect);
void testC(const int nelem, const double *coords);
void testD(const int nelem, const double *coords, int ver);
void testE(const int nelem, const double *coords, int *connect);
void print_time(const bool print_em, double &tot_time, double &utime, double &stime, long &imem, long &rmem);
void query_vert_to_elem();
void query_elem_to_vert();
void query_struct_elem_to_vert();
void query_elem_to_vert_direct();
void query_vert_to_elem_direct();
ErrorCode normalize_elems(double *coords);
void check_answers(const char *);

#define RC(msg) if (MB_SUCCESS != result) do {std::cout << "FAIL in " << msg << std::endl;return;} while (true)
#define RR(msg) if (MB_SUCCESS != result) do {std::cout << "FAIL in " << msg << std::endl;return result;} while (true)

void compute_edge(double *start, const int nelem,  const double xint,
                  const int stride) 
{
  for (int i = 1; i < nelem; i++) {
    start[i*stride] = start[0]+i*xint;
    start[nelem+1+i*stride] = start[nelem+1]+i*xint;
    start[2*(nelem+1)+i*stride] = start[2*(nelem+1)]+i*xint;
  }
}

void compute_face(double *a, const int nelem,  const double xint,
                  const int stride1, const int stride2) 
{
    // 2D TFI on a face starting at a, with strides stride1 in ada and stride2 in tse
  for (int j = 1; j < nelem; j++) {
    double tse = j * xint;
    for (int i = 1; i < nelem; i++) {
      double ada = i * xint;
      
      a[i*stride1+j*stride2] = (1.0 - ada)*a[i*stride1]
        + ada*a[i*stride1+nelem*stride2]
        + (1.0 - tse)*a[j*stride2]
        + tse*a[j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[0]
        - (1.0 - tse)*ada*a[nelem*stride1]
        - tse*(1.0 - ada)*a[nelem*stride2]
        - tse*ada*a[nelem*(stride1+stride2)];
      a[nelem+1+i*stride1+j*stride2] = (1.0 - ada)*a[nelem+1+i*stride1]
        + ada*a[nelem+1+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[nelem+1+j*stride2]
        + tse*a[nelem+1+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[nelem+1+0]
        - (1.0 - tse)*ada*a[nelem+1+nelem*stride1]
        - tse*(1.0 - ada)*a[nelem+1+nelem*stride2]
        - tse*ada*a[nelem+1+nelem*(stride1+stride2)];
      a[2*(nelem+1)+i*stride1+j*stride2] = (1.0 - ada)*a[2*(nelem+1)+i*stride1]
        + ada*a[2*(nelem+1)+i*stride1+nelem*stride2]
        + (1.0 - tse)*a[2*(nelem+1)+j*stride2]
        + tse*a[2*(nelem+1)+j*stride2+nelem*stride1]
        - (1.0 - tse)*(1.0 - ada)*a[2*(nelem+1)+0]
        - (1.0 - tse)*ada*a[2*(nelem+1)+nelem*stride1]
        - tse*(1.0 - ada)*a[2*(nelem+1)+nelem*stride2]
        - tse*ada*a[2*(nelem+1)+nelem*(stride1+stride2)];
    }
  }
}

void build_coords(const int nelem, double *&coords) 
{
  double ttime0, ttime1, utime1, stime1;
  long imem, rmem;
  print_time(false, ttime0, utime1, stime1, imem, rmem);
    // allocate the memory
  int numv = nelem+1;
  int numv_sq = numv*numv;
  int tot_numv = numv*numv*numv;
  coords = new double[3*tot_numv];

// use FORTRAN-like indexing
#define VINDEX(i,j,k) (i + (j*numv) + (k*numv_sq))
  int idx;
  double scale1, scale2, scale3;
    // use these to prevent optimization on 1-scale, etc (real map wouldn't have
    // all these equal)
  scale1 = LENGTH/nelem;
  scale2 = LENGTH/nelem;
  scale3 = LENGTH/nelem;

#ifdef REALTFI
    // use a real TFI xform to compute coordinates
    // compute edges
    // i (stride=1)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1);
  compute_edge(&coords[VINDEX(0,nelem,nelem)], nelem, scale1, 1);
    // j (stride=numv)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(0,0,nelem)], nelem, scale1, numv);
  compute_edge(&coords[VINDEX(nelem,0,nelem)], nelem, scale1, numv);
    // k (stride=numv^2)
  compute_edge(&coords[VINDEX(0,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(0,nelem,0)], nelem, scale1, numv_sq);
  compute_edge(&coords[VINDEX(nelem,nelem,0)], nelem, scale1, numv_sq);

    // compute faces
    // i=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, numv, numv_sq);
  compute_face(&coords[VINDEX(nelem,0,0)], nelem, scale1, numv, numv_sq);
    // j=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv_sq);
  compute_face(&coords[VINDEX(0,nelem,0)], nelem, scale1, 1, numv_sq);
    // k=0, nelem
  compute_face(&coords[VINDEX(0,0,0)], nelem, scale1, 1, numv);
  compute_face(&coords[VINDEX(0,0,nelem)], nelem, scale1, 1, numv);

    // initialize corner indices
  int i000 = VINDEX(0,0,0);
  int ia00 = VINDEX(nelem,0,0);
  int i0t0 = VINDEX(0,nelem,0);
  int iat0 = VINDEX(nelem,nelem,0);
  int i00g = VINDEX(0,0,nelem);
  int ia0g = VINDEX(nelem,0,nelem);
  int i0tg = VINDEX(0,nelem,nelem);
  int iatg = VINDEX(nelem,nelem,nelem);
  double cX, cY, cZ;
  int adaInts = nelem;
  int tseInts = nelem;
  int gammaInts = nelem;
  
  
  for (int i=1; i < nelem; i++) {
    for (int j=1; j < nelem; j++) {
      for (int k=1; k < nelem; k++) {
        idx = VINDEX(i,j,k);
        double tse = i*scale1;
        double ada = j*scale2;
        double gamma = k*scale3;
        double tm1 = 1.0 - tse;
        double am1 = 1.0 - ada;
        double gm1 = 1.0 - gamma;

        cX = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cY = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        cZ = gm1 *   (am1*(tm1*coords[i000] + tse*coords[i0t0])  +
                      ada*(tm1*coords[ia00] + tse*coords[iat0])) +
          gamma * (am1*(tm1*coords[i00g] + tse*coords[i0tg])  +
                   ada*(tm1*coords[ia0g] + tse*coords[iatg]));

        double *ai0k = &coords[VINDEX(k,0,i)];
        double *aiak = &coords[VINDEX(k,adaInts,i)];
        double *a0jk = &coords[VINDEX(k,j,0)];
        double *atjk = &coords[VINDEX(k,j,tseInts)];
        double *aij0 = &coords[VINDEX(0,j,i)];
        double *aijg = &coords[VINDEX(gammaInts,j,i)];
  
        coords[VINDEX(i,j,k)] = (   am1*ai0k[0] 
                                    + ada*aiak[0] 
                                    + tm1*a0jk[0] 
                                    + tse*atjk[0]
                                    + gm1*aij0[0] 
                                    + gamma*aijg[0] )/2.0 - cX/2.0;

        coords[nelem+1+VINDEX(i,j,k)] = (   am1*ai0k[nelem+1] 
                                            + ada*aiak[nelem+1] 
                                            + tm1*a0jk[nelem+1] 
                                            + tse*atjk[nelem+1]
                                            + gm1*aij0[nelem+1] 
                                            + gamma*aijg[nelem+1] )/2.0 - cY/2.0;

        coords[2*(nelem+1)+VINDEX(i,j,k)] = (   am1*ai0k[2*(nelem+1)] 
                                                + ada*aiak[2*(nelem+1)] 
                                                + tm1*a0jk[2*(nelem+1)] 
                                                + tse*atjk[2*(nelem+1)]
                                                + gm1*aij0[2*(nelem+1)] 
                                                + gamma*aijg[2*(nelem+1)] )/2.0 - cZ/2.0;
      }
    }
  }
  

#else
  for (int i=0; i < numv; i++) {
    for (int j=0; j < numv; j++) {
      for (int k=0; k < numv; k++) {
        idx = VINDEX(i,j,k);
          // blocked coordinate ordering
        coords[idx] = i*scale1;
        coords[tot_numv+idx] = j*scale2;
        coords[2*tot_numv+idx] = k*scale3;
      }
    }
  }
#endif
  print_time(false, ttime1, utime1, stime1, imem, rmem);
//  std::cout << "MOAB: TFI time = " << ttime1-ttime0 << " sec" 
//            << std::endl;
}

void build_connect(const int nelem, const EntityHandle /*vstart*/, int *&connect)
{
    // allocate the memory
  int nume_tot = nelem*nelem*nelem;
  connect = new int[8*nume_tot];

  EntityHandle vijk;
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int idx = 0;
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
        assert(i <= numv*numv*numv);
      }
    }
  }
}

Interface *gMB;
Tag pos_tag, pos2_tag;

void init() 
{
  gMB = new Core();
  double def_val[3] = {0.0, 0.0, 0.0};
  ErrorCode rval = gMB->tag_get_handle("position_tag", 3, MB_TYPE_DOUBLE, pos_tag,
                                       MB_TAG_DENSE | MB_TAG_CREAT, def_val);
  assert(MB_SUCCESS == rval);
  if (rval){} // empty line to remove compiler warning  
  rval = gMB->tag_get_handle("position2_tag", 3, MB_TYPE_DOUBLE, pos2_tag,
                             MB_TAG_DENSE | MB_TAG_CREAT, def_val);
  assert(MB_SUCCESS == rval);
}
  
int main(int argc, char* argv[])
{
  int nelem = 20;
  if (argc < 3) {
    std::cout << "Usage: " << argv[0] << " <ints_per_side> [A|B|C|D [1|2|3|4]|E]" << std::endl;
    return 1;
  }
  
  char which_test = '\0';
  int ver = 0;
  
  sscanf(argv[1], "%d", &nelem);
  if (argc >= 3) sscanf(argv[2], "%c", &which_test);
  if (argc >= 4) sscanf(argv[3], "%d", &ver);

  if (3 <= argc && which_test != 'A' && which_test != 'B' && which_test != 'C' && which_test != 'D'&& which_test != 'E') {
      std::cout << "Must indicate A or B, C, D or E for test." << std::endl;
      return 1;
  }

  if (4 <= argc && which_test == 'D' && (ver < 1 || ver > 4)) {
      std::cout << "Must indicate version 1, 2, 3, or 4 for test D." << std::endl;
      return 1;
  }
  
//  std::cout << "number of elements: " << nelem << "; test " 
//            << which_test << std::endl;

    // pre-build the coords array
  double *coords = NULL;
  build_coords(nelem, coords);
  assert(NULL != coords);

  int *connect = NULL;

    // test A: create structured mesh
  if ('\0' == which_test || 'A' == which_test) testA(nelem, coords);

  build_connect(nelem, 1, connect);

    // test B: create mesh using bulk interface
  if ('\0' == which_test || 'B' == which_test) testB(nelem, coords, connect);

    // test C: create mesh using individual interface
  if ('\0' == which_test || 'C' == which_test) testC(nelem, coords);
  
    // test D: query mesh using iterators
  if ('\0' == which_test || 'D' == which_test) testD(nelem, coords, ver);
  
    // test E: query mesh using direct access
  if ('\0' == which_test || 'E' == which_test) testE(nelem, coords, connect);
  
  return 0;
}

void query_elem_to_vert()
{
  Range all_hexes;
  ErrorCode result = gMB->get_entities_by_type(0, MBHEX, all_hexes); RC("query_elem_to_vert");
  const EntityHandle *connect;
  int num_connect;
  double dum_coords[24];
  for (Range::iterator eit = all_hexes.begin(); eit != all_hexes.end(); eit++) {
    result = gMB->get_connectivity(*eit, connect, num_connect); RC("query_elem_to_vert");
    result = gMB->get_coords(connect, num_connect, dum_coords); RC("query_elem_to_vert");

      // compute the centroid
    double centroid[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 8; j++) {
      centroid[0] += dum_coords[3*j+0];
      centroid[1] += dum_coords[3*j+1];
      centroid[2] += dum_coords[3*j+2];
    }
    centroid[0] *= 0.125; centroid[1] *= 0.125; centroid[2] *= 0.125;
    result = gMB->tag_set_data(pos_tag, &(*eit), 1, centroid); RC("query_elem_to_vert");
  }
}

void query_vert_to_elem()
{
  Range all_verts;
  std::vector<EntityHandle> neighbor_hexes;
  std::vector<double> neighbor_pos;
  double coords[3];
  neighbor_pos.resize(3*8); // average vertex will have 8 adjacent hexes
  ErrorCode result = gMB->get_entities_by_type(0, MBVERTEX, all_verts); RC("query_vert_to_elem");
  for (Range::iterator vit = all_verts.begin(); vit != all_verts.end(); vit++) {
    neighbor_hexes.clear();
    result = gMB->get_coords(&(*vit), 1, coords); RC("query_vert_to_elem");
    result = gMB->get_adjacencies(&(*vit), 1, 3, false, neighbor_hexes); RC("query_vert_to_elem");
    assert(neighbor_pos.size() >= 3*neighbor_hexes.size());
    result = gMB->tag_get_data(pos2_tag, &neighbor_hexes[0], neighbor_hexes.size(), &neighbor_pos[0]); RC("query_vert_to_elem");
    for (unsigned int i = 0; i < neighbor_hexes.size(); i++) {
      neighbor_pos[3*i] += coords[0];
      neighbor_pos[3*i+1] += coords[1];
      neighbor_pos[3*i+2] += coords[2];
    }

    result = gMB->tag_set_data(pos2_tag, &neighbor_hexes[0], neighbor_hexes.size(), &neighbor_pos[0]); RC("query_vert_to_elem");
  }

    // get all hexes and divide pos_tag by 8; reuse all_verts
  result = normalize_elems(coords); RC("query_vert_to_elem");
}

void query_elem_to_vert_iters(int chunk_size, bool check_valid, 
                              std::vector<EntityHandle> &connect, 
                              double *dum_coords, double *dum_pos)
{
  std::vector<EntityHandle> hexes;
  SetIterator *iter;
  ErrorCode result = gMB->create_set_iterator(0, MBHEX, -1, chunk_size, check_valid, iter); RC("query_elem_to_vert_iters");
  bool atend = false;
  while (!atend) {
    hexes.clear();
    result = iter->get_next_arr(hexes, atend); RC("query_elem_to_vert_iters");
    result = gMB->get_connectivity(&hexes[0], hexes.size(), connect); RC("query_elem_to_vert_iters");
    result = gMB->get_coords(&connect[0], connect.size(), dum_coords); RC("query_elem_to_vert_iters");
    result = gMB->tag_get_data(pos_tag, &hexes[0], hexes.size(), dum_pos); RC("query_elem_to_vert_iters");
    for (unsigned int i = 0; i < hexes.size(); i++) {
        // compute the centroid
      for (int j = 0; j < 8;j++) {
        dum_pos[3*i+0] += dum_coords[24*i+3*j];
        dum_pos[3*i+1] += dum_coords[24*i+3*j+1];
        dum_pos[3*i+2] += dum_coords[24*i+3*j+2];
      }
      dum_pos[3*i+0] *= 0.125; dum_pos[3*i+1] *= 0.125; dum_pos[3*i+2] *= 0.125;
    }
    result = gMB->tag_set_data(pos_tag, &hexes[0], hexes.size(), dum_pos); RC("query_elem_to_vert_iters");
  }
  
  delete iter;
}

void query_vert_to_elem_iters(int chunk_size, bool check_valid, std::vector<EntityHandle> &/*connect*/,
                              double *dum_coords, double *dum_pos)
{
  std::vector<EntityHandle> verts, neighbor_hexes;
  SetIterator *iter;
  ErrorCode result = gMB->create_set_iterator(0, MBVERTEX, -1, chunk_size, check_valid, iter); RC("query_vert_to_elem_iters");
  assert(MB_SUCCESS == result);
  bool atend = false;
  while (!atend) {
    verts.clear();
    result = iter->get_next_arr(verts, atend); RC("query_vert_to_elem_iters");
    result = gMB->get_coords(&verts[0], verts.size(), dum_coords); RC("query_vert_to_elem_iters");
    chunk_size = std::min((int)verts.size(), chunk_size);
    for (int i = 0; i < chunk_size; i++) {
      neighbor_hexes.clear();
      result = gMB->get_adjacencies(&verts[i], 1, 3, false, neighbor_hexes); 
      RC("query_vert_to_elem_iters");
      result = gMB->tag_get_data(pos2_tag, &neighbor_hexes[0], neighbor_hexes.size(), dum_pos); 
      RC("query_vert_to_elem_iters");
      for (unsigned int j = 0; j < neighbor_hexes.size(); j++) {
        dum_pos[3*j+0] += dum_coords[3*i+0];
        dum_pos[3*j+1] += dum_coords[3*i+1];
        dum_pos[3*j+2] += dum_coords[3*i+2];
      }
      result = gMB->tag_set_data(pos2_tag, &neighbor_hexes[0], neighbor_hexes.size(), dum_pos); 
      RC("query_vert_to_elem_iters");
    }
  }

  result = normalize_elems(dum_coords); RC("query_vert_to_elem_iters");
}

// normalize pos_tag on all elems by 1/8; coords assumed large enough to hold 3 doubles
ErrorCode normalize_elems(double *coords)
{
    // get all hexes and divide pos_tag by 8
  Range elems;
  ErrorCode result = gMB->get_entities_by_type(0, MBHEX, elems); RR("normalize");
  
  for (Range::iterator vit = elems.begin(); vit != elems.end(); vit++) {
    result = gMB->tag_get_data(pos2_tag, &(*vit), 1, coords); RR("normalize");
    coords[0] *= 0.125; coords[1] *= 0.125; coords[2] *= 0.125;
    result = gMB->tag_set_data(pos2_tag, &(*vit), 1, coords); RR("normalize");
  }

  return result;
}

void query_struct_elem_to_vert()
{
    // assumes brick mapped mesh with handles starting at zero
  Range all_hexes;
  ErrorCode result = gMB->get_entities_by_type(0, MBHEX, all_hexes); RC("query_struct_elem_to_vert");
  double dum_coords[24];
  std::vector<EntityHandle> connect;
  for (Range::iterator eit = all_hexes.begin(); eit != all_hexes.end(); eit++) {
    result = gMB->get_connectivity(&(*eit), 1, connect); RC("query_struct_elem_to_vert");
    result = gMB->get_coords(&connect[0], connect.size(), dum_coords); RC("query_struct_elem_to_vert");

    double centroid[3] = {0.0, 0.0, 0.0};
    for (int j = 0; j < 8; j++) {
      centroid[0] += dum_coords[3*j+0];
      centroid[1] += dum_coords[3*j+1];
      centroid[2] += dum_coords[3*j+2];
    }
    centroid[0] *= 0.125; centroid[1] *= 0.125; centroid[2] *= 0.125; 
    result = gMB->tag_set_data(pos_tag, &(*eit), 1, centroid); RC("query_struct_elem_to_vert");
  }
}

#if defined(_MSC_VER) || defined(__MINGW32__)
void print_time(const bool print_em, double &tot_time, double &utime, double &stime, long &imem, long &rmem) 
{
  utime = (double)clock() / CLOCKS_PER_SEC;
  if (print_em)
    std::cout << "Total wall time = " << utime << std::endl;
  tot_time = stime = 0;
  imem = rmem = 0;
}
#else
void print_time(const bool print_em, double &tot_time, double &utime, double &stime, long &imem, long &rmem) 
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
  if (print_em) {
    std::cout << "Max resident set size = " << r_usage.ru_maxrss << " kbytes" << std::endl;
    std::cout << "Int resident set size = " << r_usage.ru_idrss << std::endl;
  }
  imem = r_usage.ru_idrss;
  rmem = r_usage.ru_maxrss;
#else
  system("ps o args,drs,rss | grep perf | grep -v grep");  // RedHat 9.0 doesnt fill in actual memory data 
  imem = rmem = 0;
#endif
}
#endif

void testA(const int nelem, const double *coords) 
{
  double ttime0, ttime1, ttime2, ttime3, ttime4, utime, stime;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4;
  
  print_time(false, ttime0, utime, stime, imem0, rmem0);

    // make a 3d block of vertices
  EntitySequence *dum_seq = NULL;
  ScdVertexData *vseq = NULL;
  StructuredElementSeq *eseq = NULL;
  init();
  SequenceManager *seq_mgr = dynamic_cast<Core*>(gMB)->sequence_manager();
  HomCoord vseq_minmax[2] = {HomCoord(0,0,0), 
                             HomCoord(nelem, nelem, nelem)};
  EntityHandle vstart, estart;
  
  ErrorCode result = seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1],
                                                  MBVERTEX, 1, vstart, dum_seq); RC("testA");
  if (NULL != dum_seq) vseq = dynamic_cast<ScdVertexData*>(dum_seq->data());
  assert (MB_FAILURE != result && vstart != 0 && dum_seq != NULL && vseq != NULL);
    // now the element sequence
  result = seq_mgr->create_scd_sequence(vseq_minmax[0], vseq_minmax[1], 
                                        MBHEX, 1, estart, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<StructuredElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart != 0 && dum_seq != NULL && eseq != NULL);
  
    // only need to add one vseq to this, unity transform
    // trick: if I know it's going to be unity, just input 3 sets of equivalent points
  result = eseq->sdata()->add_vsequence(vseq, vseq_minmax[0], vseq_minmax[0], vseq_minmax[0], 
                               vseq_minmax[0], vseq_minmax[0], vseq_minmax[0]);
  assert(MB_SUCCESS == result);

    // set the coordinates of the vertices
  EntityHandle handle;
  int i;
  double dumv[3];
  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  for (i = 0, handle = vstart; i < num_verts; i++, handle++) {
    dumv[0] = coords[i]; dumv[1] = coords[num_verts+i]; dumv[2] = coords[2*num_verts+i]; 
    result = gMB->set_coords(&handle, 1, dumv);
    assert(MB_SUCCESS == result);
  }

  print_time(false, ttime1, utime, stime, imem1, rmem1);

    // query the mesh 2 ways
  query_struct_elem_to_vert();

  print_time(false, ttime2, utime, stime, imem2, rmem2);

  query_vert_to_elem();
  
  print_time(false, ttime3, utime, stime, imem3, rmem3);

#ifndef NDEBUG
  check_answers("A");
#endif  

  delete gMB;
  
  print_time(false, ttime4, utime, stime, imem4, rmem4);

  std::cout << "MOAB_scd:nelem,construct,e_to_v,v_to_e,after_dtor,total= " 
            << nelem << " "
            << ttime1-ttime0 << " " 
            << ttime2-ttime1 << " " 
            << ttime3-ttime2 << " " 
            << ttime4-ttime3 << " " 
            << ttime4-ttime0 << " seconds" 
            << std::endl;
  std::cout << "MOAB_scd_memory(rss):initial,after_construction,e-v,v-e,after_dtor= " 
            << rmem0 << " " << rmem1 << " " << rmem2 << " " << rmem3 << " " << rmem4 <<  " kb" << std::endl;
}

void testB(const int nelem, const double *coords, int *connect) 
{
  double ttime0, ttime1, ttime2, ttime3, ttime4, utime, stime;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4;
  
  print_time(false, ttime0, utime, stime, imem0, rmem0);

  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  int num_elems = nelem*nelem*nelem;
  EntityHandle vstart, estart;

  // get the read interface
  ReadUtilIface* readMeshIface;
  init();
  gMB->query_interface(readMeshIface);
  
  // create a sequence to hold the node coordinates
  // get the current number of entities and start at the next slot
  std::vector<double*> coord_arrays;
  ErrorCode result = readMeshIface->get_node_coords(3, num_verts, 1, vstart, coord_arrays); RC("testB");
  assert(MB_SUCCESS == result && 1 == vstart &&
         coord_arrays[0] && coord_arrays[1] && coord_arrays[2]);
    // memcpy the coordinate data into place
  memcpy(coord_arrays[0], coords, sizeof(double)*num_verts);
  memcpy(coord_arrays[1], &coords[num_verts], sizeof(double)*num_verts);
  memcpy(coord_arrays[2], &coords[2*num_verts], sizeof(double)*num_verts);
  
  EntityHandle *conn = 0;
  result = readMeshIface->get_element_connect(num_elems, 8, MBHEX, 1, estart, conn);
  assert(MB_SUCCESS == result);
  for (int i = 0; i < num_elems*8; i++)
    conn[i] = vstart + connect[i];

  delete [] connect;
  
  result = readMeshIface->update_adjacencies(estart, num_elems, 8, conn);
  assert(MB_SUCCESS == result);

  print_time(false, ttime1, utime, stime, imem1, rmem1);

    // query the mesh 2 ways
  query_elem_to_vert();

  print_time(false, ttime2, utime, stime, imem2, rmem2);

  query_vert_to_elem();
  
  print_time(false, ttime3, utime, stime, imem3, rmem3);

#ifndef NDEBUG
  check_answers("B");
#endif  

  delete gMB;

  print_time(false, ttime4, utime, stime, imem4, rmem4);

  std::cout << "MOAB_ucd_blocked:nelem,construct,e_to_v,v_to_e,after_dtor,total= " 
            << nelem << " "
            << ttime1-ttime0 << " " 
            << ttime2-ttime1 << " " 
            << ttime3-ttime2 << " " 
            << ttime4-ttime3 << " " 
            << ttime4-ttime0 << " seconds" 
            << std::endl;
  std::cout << "MOAB_ucdblocked_memory_(rss):initial,after_construction,e-v,v-e,after_dtor= " 
            << rmem0 << " " << rmem1 << " " << rmem2 << " " << rmem3 << " " << rmem4 <<  " kb" << std::endl;
}

void testC(const int nelem, const double *coords) 
{
  double ttime0, ttime1, ttime2, ttime3, ttime4, utime, stime;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4;
  
  print_time(false, ttime0, utime, stime, imem0, rmem0);

    // create the vertices; assume we don't need to keep a list of vertex handles, since they'll
    // be created in sequence
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int num_verts = numv*numv*numv;
  double dum_coords[3] = {coords[0], coords[num_verts], coords[2*num_verts]};
  EntityHandle vstart;

  init();
  ErrorCode result = gMB->create_vertex(dum_coords, vstart); RC("testC");
  assert(MB_SUCCESS == result && 1 == vstart);

  EntityHandle dum_vert, vijk;
  int i;
  for (i = 1; i < num_verts; i++) {
    dum_coords[0] = coords[i]; 
    dum_coords[1] = coords[num_verts+i]; 
    dum_coords[2] = coords[2*num_verts+i];
    result = gMB->create_vertex(dum_coords, dum_vert);
    assert(MB_SUCCESS == result);
  }

  EntityHandle dum_conn[8];
  for (i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = vstart+VINDEX(i,j,k);
        dum_conn[0] = vijk;
        dum_conn[1] = vijk+1;
        dum_conn[2] = vijk+1+numv;
        dum_conn[3] = vijk+numv;
        dum_conn[4] = vijk+numv*numv;
        dum_conn[5] = vijk+1+numv*numv;
        dum_conn[6] = vijk+1+numv+numv*numv;
        dum_conn[7] = vijk+numv+numv*numv;
        result = gMB->create_element(MBHEX, dum_conn, 8, dum_vert);
        assert(MB_SUCCESS == result);
      }
    }
  }

  print_time(false, ttime1, utime, stime, imem1, rmem1);

    // query the mesh 2 ways
  query_elem_to_vert();

  print_time(false, ttime2, utime, stime, imem2, rmem2);

  query_vert_to_elem();
  
  print_time(false, ttime3, utime, stime, imem3, rmem3);

#ifndef NDEBUG
  check_answers("C");
#endif  

  delete gMB;

  print_time(false, ttime4, utime, stime, imem4, rmem4);

  std::cout << "MOAB_ucd_indiv:nelem,construct,e_to_v,v_to_e,after_dtor,total= " 
            << nelem << " "
            << ttime1-ttime0 << " " 
            << ttime2-ttime1 << " " 
            << ttime3-ttime2 << " " 
            << ttime4-ttime3 << " " 
            << ttime4-ttime0 << " seconds" 
            << std::endl;
  std::cout << "MOAB_ucd_indiv_memory_(rss):initial,after_construction,e-v,v-e,after_dtor= " 
            << rmem0 << " " << rmem1 << " " << rmem2 << " " << rmem3 << " " << rmem4 <<  " kb" << std::endl;
}

void testD(const int nelem, const double *coords, int ver) 
{
  double ttime0, ttime1, ttime2, ttime3, ttime4, ttime5, ttime6, ttime7, ttime8, ttime9, ttime10, 
      utime, stime;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4, imem5, rmem5, imem6, rmem6,
      imem7, rmem7, imem8, rmem8, imem9, rmem9, imem10, rmem10;
  
  print_time(false, ttime0, utime, stime, imem0, rmem0);

    // create the vertices; assume we don't need to keep a list of vertex handles, since they'll
    // be created in sequence
  int numv = nelem + 1;
  int numv_sq = numv*numv;
  int num_verts = numv*numv*numv;
  std::vector<double> dum_coords(24), dum_pos(24);
  dum_coords[0] = coords[0]; dum_coords[1] = coords[num_verts]; dum_coords[2] = coords[2*num_verts];
  EntityHandle vstart;

  init();
  ErrorCode result = gMB->create_vertex(&dum_coords[0], vstart); RC("testD");
  assert(MB_SUCCESS == result && 1 == vstart);

  EntityHandle dum_vert, vijk;
  int i;
  for (i = 1; i < num_verts; i++) {
    dum_coords[0] = coords[i]; 
    dum_coords[1] = coords[num_verts+i]; 
    dum_coords[2] = coords[2*num_verts+i];
    result = gMB->create_vertex(&dum_coords[0], dum_vert);
    assert(MB_SUCCESS == result);
  }

  EntityHandle dum_conn[8];
  for (i=0; i < nelem; i++) {
    for (int j=0; j < nelem; j++) {
      for (int k=0; k < nelem; k++) {
        vijk = vstart+VINDEX(i,j,k);
        dum_conn[0] = vijk;
        dum_conn[1] = vijk+1;
        dum_conn[2] = vijk+1+numv;
        dum_conn[3] = vijk+numv;
        dum_conn[4] = vijk+numv*numv;
        dum_conn[5] = vijk+1+numv*numv;
        dum_conn[6] = vijk+1+numv+numv*numv;
        dum_conn[7] = vijk+numv+numv*numv;
        result = gMB->create_element(MBHEX, dum_conn, 8, dum_vert);
        assert(MB_SUCCESS == result);
      }
    }
  }

  print_time(false, ttime1, utime, stime, imem1, rmem1);

    // query the mesh 2 ways with !check_valid
  std::vector<EntityHandle> connect(8);
#ifndef NDEBUG
  // used only in debug mode
  double def_val[3] = {0.0, 0.0, 0.0};
#endif
  if (ver==0 || ver==1)
  {
    query_elem_to_vert_iters(1, false, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime2, utime, stime, imem2, rmem2);
    query_vert_to_elem_iters(1, false, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime3, utime, stime, imem3, rmem3);
#ifndef NDEBUG
    check_answers("D");
    result = gMB->tag_delete(pos_tag);
    result = gMB->tag_get_handle("position_tag", 3, MB_TYPE_DOUBLE, pos_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
    result = gMB->tag_delete(pos2_tag);
    result = gMB->tag_get_handle("position2_tag", 3, MB_TYPE_DOUBLE, pos2_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
#endif
  }

  if (ver==0 || ver==2)
  {
    if (ver!=0)
       print_time(false, ttime3, utime, stime, imem3, rmem3);
    query_elem_to_vert_iters(1, true, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime4, utime, stime, imem4, rmem4);
    query_vert_to_elem_iters(1, true, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime5, utime, stime, imem5, rmem5);
#ifndef NDEBUG
    check_answers("D");
    result = gMB->tag_delete(pos_tag);
    result = gMB->tag_get_handle("position_tag", 3, MB_TYPE_DOUBLE, pos_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
    result = gMB->tag_delete(pos2_tag);
    result = gMB->tag_get_handle("position2_tag", 3, MB_TYPE_DOUBLE, pos2_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
#endif
  }

  if (ver==0 || ver>=3)
  {
    dum_coords.resize(2400); dum_pos.resize(300);
  }
  if (ver==0 || ver==3)
  {
    if (ver!=0)
       print_time(false, ttime5, utime, stime, imem3, rmem3);
    query_elem_to_vert_iters(100, false, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime6, utime, stime, imem6, rmem6);
    query_vert_to_elem_iters(100, false, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime7, utime, stime, imem7, rmem7);
#ifndef NDEBUG
    check_answers("D");
    result = gMB->tag_delete(pos_tag);
    result = gMB->tag_get_handle("position_tag", 3, MB_TYPE_DOUBLE, pos_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
    result = gMB->tag_delete(pos2_tag);
    result = gMB->tag_get_handle("position2_tag", 3, MB_TYPE_DOUBLE, pos2_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
#endif
  }
  
  if (ver==0 || ver==4)
  {
    if (ver!=0)
       print_time(false, ttime7, utime, stime, imem3, rmem3);
    query_elem_to_vert_iters(100, true, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime8, utime, stime, imem8, rmem8);
    query_vert_to_elem_iters(100, true, connect, &dum_coords[0], &dum_pos[0]);
    print_time(false, ttime9, utime, stime, imem9, rmem9);
#ifndef NDEBUG
    check_answers("D");
    result = gMB->tag_delete(pos_tag);
    result = gMB->tag_get_handle("position_tag", 3, MB_TYPE_DOUBLE, pos_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
    result = gMB->tag_delete(pos2_tag);
    result = gMB->tag_get_handle("position2_tag", 3, MB_TYPE_DOUBLE, pos2_tag,
                                 MB_TAG_DENSE | MB_TAG_CREAT, def_val);
    assert(MB_SUCCESS == result);
#endif
  }

  if (ver>0 && ver<4)
    print_time(false, ttime9, utime, stime, imem9, rmem9);
  delete gMB;

  print_time(false, ttime10, utime, stime, imem10, rmem10);

  if (ver==0 || ver==1)
  {
    std::cout << "MOAB_ucd_iters_!check_valid_1:nelem,construct,e-v,v-e,after_dtor,total= "
              << nelem << " " << ttime1-ttime0 << " " << ttime2-ttime1 << " " << ttime3-ttime2 << " " << ttime10-ttime9 
              << " " << ttime3-ttime0 + ttime10-ttime9 << std::endl;
    std::cout << "MOAB_ucd_iters_memory_(rss)_!check_valid_1:initial,after_construction,e-v,v-e,after_dtor= " 
              << rmem0 << " " << rmem1 << " " << rmem2 << " " << rmem3 << " " << rmem10 << " kb" << std::endl;
  }
  if (ver==0 || ver==2)
  {
    std::cout << "MOAB_ucd_iters_check_valid_1:nelem,construct,e-v,v-e,after_dtor,total= "
              << nelem << " " << ttime1-ttime0 << " " << ttime4-ttime3 << " " << ttime5-ttime4 << " " << ttime10-ttime9 
              << " " << ttime1-ttime0 + ttime5-ttime3 + ttime10-ttime9 << std::endl;
    std::cout << "MOAB_ucd_iters_memory_(rss)_check_valid_1:initial,after_construction,e-v,v-e,after_dtor= " 
              << rmem0 << " " << rmem1 << " " << rmem2 << " " << rmem3 << " " << rmem10 << " kb" << std::endl;
  }
  if (ver==0 || ver==3)
  {
    std::cout << "MOAB_ucd_iters_!check_valid_100:nelem,construct,e-v,v-e,after_dtor,total= "
              << nelem << " " << ttime1-ttime0 << " " << ttime6-ttime5 << " " << ttime7-ttime6 << " " << ttime10-ttime9 
              << " " << ttime1-ttime0 + ttime7-ttime5 + ttime10-ttime9 << std::endl;
    std::cout << "MOAB_ucd_iters_memory_(rss)_!check_valid_100:initial,after_construction,e-v,v-e,after_dtor= " 
              << rmem0 << " " << rmem1 << " " << rmem6 << " " << rmem7 << " " << rmem10 << " kb" << std::endl;
  }
  if (ver==0 || ver==4)
  {
    std::cout << "MOAB_ucd_iters_check_valid_100:nelem,construct,e-v,v-e,after_dtor,total= "
              << nelem << " " << ttime1-ttime0 << " " << ttime8-ttime7 << " " << ttime9-ttime8 << " " << ttime10-ttime9 
              << " " << ttime1-ttime0 + ttime10-ttime7 << std::endl;
    std::cout << "MOAB_ucd_iters_memory_(rss)_check_valid_100:initial,after_construction,e-v,v-e,after_dtor= " 
              << rmem0 << " " << rmem1 << " " << rmem8 << " " << rmem9 << " " << rmem10 << " kb" << std::endl;
  }
}

void testE(const int nelem, const double *coords, int *connect) 
{
  double ttime0, ttime1, ttime2, ttime3, ttime4, ttime5, ttime6, utime, stime;
  long imem0, rmem0, imem1, rmem1, imem2, rmem2, imem3, rmem3, imem4, rmem4, imem5, rmem5, imem6, rmem6;
  
  print_time(false, ttime0, utime, stime, imem0, rmem0);

  int num_verts = (nelem + 1)*(nelem + 1)*(nelem + 1);
  int num_elems = nelem*nelem*nelem;
  EntityHandle vstart, estart;

  // get the read interface
  ReadUtilIface* readMeshIface;
  init();
  gMB->query_interface(readMeshIface);
  
  // create a sequence to hold the node coordinates
  // get the current number of entities and start at the next slot
  std::vector<double*> coord_arrays;
  ErrorCode result = readMeshIface->get_node_coords(3, num_verts, 1, vstart, coord_arrays); RC("testE");
  assert(MB_SUCCESS == result && 1 == vstart &&
         coord_arrays[0] && coord_arrays[1] && coord_arrays[2]);
    // memcpy the coordinate data into place
  memcpy(coord_arrays[0], coords, sizeof(double)*num_verts);
  memcpy(coord_arrays[1], &coords[num_verts], sizeof(double)*num_verts);
  memcpy(coord_arrays[2], &coords[2*num_verts], sizeof(double)*num_verts);
  
  EntityHandle *conn = 0;
  result = readMeshIface->get_element_connect(num_elems, 8, MBHEX, 1, estart, conn);
  assert(MB_SUCCESS == result);
  for (int i = 0; i < num_elems*8; i++)
    conn[i] = vstart + connect[i];

  delete [] connect;

  Range verts(vstart, vstart+num_verts-1), elems(estart, estart+num_elems-1);
  
  result = readMeshIface->update_adjacencies(estart, num_elems, 8, conn);
  assert(MB_SUCCESS == result);

  print_time(false, ttime1, utime, stime, imem1, rmem1);

    // query the mesh 2 ways
  query_elem_to_vert_direct();

  print_time(false, ttime2, utime, stime, imem2, rmem2);

  query_vert_to_elem_direct();
  
  print_time(false, ttime3, utime, stime, imem3, rmem3);

#ifndef NDEBUG
  check_answers("E");
#endif  

  query_elem_to_vert_direct();

  print_time(false, ttime4, utime, stime, imem4, rmem4);

  query_vert_to_elem_direct();
  
  print_time(false, ttime5, utime, stime, imem5, rmem5);

#ifndef NDEBUG
  check_answers("E");
#endif  

  delete gMB;

  print_time(false, ttime6, utime, stime, imem6, rmem6);

  std::cout << "MOAB_ucd_direct:nelem,construct,e_to_v,v_to_e,after_dtor,total= " 
            << nelem << " "
            << ttime1-ttime0 << " " 
            << ttime2-ttime1 << " " 
            << ttime3-ttime2 << " " 
            << ttime6-ttime5 << " " 
            << ttime3-ttime0 + ttime6-ttime5 << " seconds" 
            << std::endl;
  std::cout << "MOAB_ucd_direct_memory_(rss):initial,after_construction,e-v,v-e,after_dtor= " 
            << rmem0 << " " << rmem1 << " " << rmem2 << " " << rmem3 << " " << rmem6 <<  " kb" << std::endl;

  std::cout << "MOAB_ucd_direct2:nelem,construct,e_to_v,v_to_e,after_dtor,total= " 
            << nelem << " "
            << ttime1-ttime0 << " " 
            << ttime4-ttime3 << " " 
            << ttime5-ttime4 << " " 
            << ttime6-ttime5 << " " 
            << ttime1-ttime0 + ttime6-ttime3 << " seconds" 
            << std::endl;
  std::cout << "MOAB_ucd_direct2_memory_(rss):initial,after_construction,e-v,v-e,after_dtor= " 
            << rmem0 << " " << rmem1 << " " << rmem4 << " " << rmem5 << " " << rmem6 <<  " kb" << std::endl;
}

void query_elem_to_vert_direct()
{
  Range all_hexes, all_verts;
  ErrorCode result = gMB->get_entities_by_type(0, MBHEX, all_hexes);
  result = gMB->get_entities_by_type(0, MBVERTEX, all_verts);  RC("query_elem_to_vert_direct");
  EntityHandle *connect;
  int ecount, vcount, vpere;
  double *coords[3];

  result = gMB->connect_iterate(all_hexes.begin(), all_hexes.end(), connect, vpere, ecount);
  if (MB_SUCCESS != result || ecount != (int)all_hexes.size()) {
    std::cout << "FAILED in connect_iterate!" << std::endl;
    return;
  }
  result = gMB->coords_iterate(all_verts.begin(), all_verts.end(), coords[0], coords[1], coords[2], vcount);
  if (MB_SUCCESS != result || vcount != (int)all_verts.size()) {
    std::cout << "FAILED in coords_iterate!" << std::endl;
    return;
  }
  double *centroid;
  result = gMB->tag_iterate(pos_tag, all_hexes.begin(), all_hexes.end(), ecount, (void*&)centroid);
  if (MB_SUCCESS != result || ecount != (int)all_hexes.size()) {
    std::cout << "FAILED in connect_iterate!" << std::endl;
    return;
  }

  EntityHandle vstart = *all_verts.begin();
  for (int i = 0; i < ecount; i++) {
      // compute the centroid
    for (int j = 0; j < vpere; j++) {
      int vind = *connect-vstart;
      connect++;
      centroid[3*i+0] += coords[0][vind];
      centroid[3*i+1] += coords[1][vind];
      centroid[3*i+2] += coords[2][vind];
    }

      // now normalize
    centroid[3*i+0] *= 0.125; centroid[3*i+1] *= 0.125; centroid[3*i+2] *= 0.125;
  }
}

void query_vert_to_elem_direct()
{
  Range all_verts, tmp_ents;
  ErrorCode result = gMB->get_entities_by_type(0, MBVERTEX, all_verts);
  assert(MB_SUCCESS == result);

    // make sure vertex-element adjacencies are created
  result = gMB->get_adjacencies(&(*all_verts.begin()), 1, 3, false, tmp_ents);
  
  const std::vector<EntityHandle> **adjs;
  int count;
  result = gMB->adjacencies_iterate(all_verts.begin(), all_verts.end(), adjs, count);
  if (MB_SUCCESS != result && count != (int)all_verts.size()) {
    std::cout << "FAILED:adjacencies_iterate." << std::endl;
    return;
  }

  double *coords[3];
  result = gMB->coords_iterate(all_verts.begin(), all_verts.end(), coords[0], coords[1], coords[2], count);
  if (MB_SUCCESS != result || count != (int)all_verts.size()) {
    std::cout << "FAILED in coords_iterate!" << std::endl;
    return;
  }
    // get all hexes, then iterator over pos2_tag
  double *centroid;
  int ecount;
  tmp_ents.clear();
  result = gMB->get_entities_by_type(0, MBHEX, tmp_ents);
  assert(MB_SUCCESS == result);
  result = gMB->tag_iterate(pos2_tag, tmp_ents.begin(), tmp_ents.end(), ecount, (void*&)centroid);
  if (MB_SUCCESS != result || ecount != (int)tmp_ents.size()) {
    std::cout << "FAILED in tag_iterate!" << std::endl;
    return;
  }

  int i;
  Range::iterator vit;
  EntityHandle estart = *tmp_ents.begin();
  for (vit = all_verts.begin(), i = 0; vit != all_verts.end(); vit++, i++) {
    assert(adjs[i]);
    for (std::vector<EntityHandle>::const_iterator vit2 = adjs[i]->begin(); vit2 != adjs[i]->end(); vit2++)
      if (*vit >= estart) {
        int eind = *vit2 - estart;
        centroid[3*eind+0] += coords[0][i];
        centroid[3*eind+1] += coords[1][i];
        centroid[3*eind+2] += coords[2][i];
      }
  }

    // now normalize
  for (i = 0; i < (int)tmp_ents.size(); i++) {
    centroid[3*i+0] *= 0.125; centroid[3*i+1] *= 0.125; centroid[3*i+2] *= 0.125; 
  }
}

void check_answers(const char */*test_name*/)
{
  Range elems;
  ErrorCode result = gMB->get_entities_by_type(0, MBHEX, elems); RC("check_answers");
  
  double coords1[3], coords2[3], del[3];
  double diff = 0.0;
  for (Range::iterator vit = elems.begin(); vit != elems.end(); vit++) {
    result = gMB->tag_get_data(pos_tag, &(*vit), 1, coords1); RC("check_answers");
    result = gMB->tag_get_data(pos2_tag, &(*vit), 1, coords2); RC("check_answers");
    for (int i = 0; i < 3; i++) del[i] = fabs(coords1[i]-coords2[i]);
    if (del[0] || del[1] || del[2]) {
      double len2 = std::max(coords1[0]*coords1[0]+coords1[1]*coords1[1]+coords1[2]*coords1[2],
                             coords2[0]*coords2[0]+coords2[1]*coords2[1]+coords2[2]*coords2[2]),
          num = del[0]*del[0]+del[1]*del[1]+del[2]*del[2];
      if (len2 > 0.0) 
        diff = std::max(diff, num/sqrt(len2));
      else if (num > 0.0) 
        diff = sqrt(num);
    }
  }
  if (diff > 0.0) 
    std::cout << "Max relative difference = " << diff << std::endl;
}

