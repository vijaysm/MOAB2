//
// test the structured sequence stuff
//

#include "ScdVertexSeq.hpp"
#include "ScdElementSeq.hpp"
#include "EntitySequenceManager.hpp"
#include "EntitySequence.hpp"
#include "MBCore.hpp"

#include <iostream>

int test_vertex_seq(MBCore *gMB);
MBErrorCode check_vertex_sequence(const ScdVertexSeq *this_seq, 
                                 const int imin, const int jmin, const int kmin, 
                                 const int imax, const int jmax, const int kmax, 
                                 const MBEntityHandle this_start);
MBErrorCode evaluate_vertex_sequence(ScdVertexSeq *this_seq);

int test_element_seq(MBCore *gMB);
MBErrorCode check_element_sequence(const ScdElementSeq *this_seq, 
                                    const HomCoord &min_params,
                                    const HomCoord &max_params,
                                    const MBEntityHandle this_start);
MBErrorCode evaluate_element_sequence(ScdElementSeq *this_seq);
int eseq_test1a(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max);
int eseq_test1b(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max);
int eseq_test1c(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max);
int eseq_test2a(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max);
int eseq_test2b(MBCore *gMB);
int eseq_test2c(MBCore *gMB);
int eseq_test2d(MBCore *gMB);

int create_1d_3_sequences(MBCore *gMB,
                          HomCoord tmp_min, HomCoord tmp_max,
                          ScdVertexSeq **vseq, MBEntityHandle *vstart,
                          ScdElementSeq **eseq, MBEntityHandle *estart);

int create_2d_3_sequences(MBCore *gMB,
                          ScdVertexSeq **vseq, MBEntityHandle *vstart,
                          ScdElementSeq **eseq, MBEntityHandle *estart);

int create_2dtri_3_sequences(MBCore *gMB,
                             const int int1, const int int2, const int int3,
                             ScdVertexSeq **vseq, MBEntityHandle *vstart,
                             ScdElementSeq **eseq, MBEntityHandle *estart);
int create_3dtri_3_sequences(MBCore *gMB,
                             const int int1, const int int2, const int int3, const int int4,
                             ScdVertexSeq **vseq, MBEntityHandle *vstart,
                             ScdElementSeq **eseq, MBEntityHandle *estart);

// first comes general-capability code used by various tests; main and test functions
// come after these, starting with main
MBErrorCode check_vertex_sequence(const ScdVertexSeq *this_seq, 
                                   const int imin, const int jmin, const int kmin, 
                                   const int imax, const int jmax, const int kmax, 
                                   const MBEntityHandle this_start) 
{
  MBErrorCode result = MB_SUCCESS;
  
    // check data stored in sequence with that in arg list
  if (imin != this_seq->i_min() ||
      jmin != this_seq->j_min() ||
      kmin != this_seq->k_min() ||
      imax != this_seq->i_max() ||
      jmax != this_seq->j_max() ||
      kmax != this_seq->k_max()) {
    std::cout << "min/max params not correct for a sequence." << std::endl;
    result = MB_FAILURE;
  }

  int i1, j1, k1, i2, j2, k2, di, dj, dk;
  this_seq->min_params(i1, j1, k1);
  this_seq->max_params(i2, j2, k2);
  this_seq->param_extents(di, dj, dk);
  if (i2-i1+1 != di || j2-j1+1 != dj || k2-k1+1 != dk) {
    std::cout << "min_params/max_params/param_extents functions returned inconsistent results"
              << std::endl;
    result = MB_FAILURE;
  }
  
  if (this_start != this_seq->get_start_handle()) {
    std::cout << "Start handle for sequence wrong." << std::endl;
    result = MB_FAILURE;
  }
  
  return result;
}
  
MBErrorCode check_element_sequence(const ScdElementSeq *this_seq, 
                                    const HomCoord &min_params,
                                    const HomCoord &max_params,
                                    const MBEntityHandle this_start) 
{
  MBErrorCode result = MB_SUCCESS;
  
    // check data stored in sequence with that in arg list
  if (min_params.i() != this_seq->i_min() ||
      min_params.j() != this_seq->j_min() ||
      min_params.k() != this_seq->k_min() ||
      max_params.i() != this_seq->i_max() ||
      max_params.j() != this_seq->j_max() ||
      max_params.k() != this_seq->k_max()) {
    std::cout << "min/max params not correct for a sequence." << std::endl;
    result = MB_FAILURE;
  }

  int i1, j1, k1, i2, j2, k2, di, dj, dk;
  this_seq->min_params(i1, j1, k1);
  this_seq->max_params(i2, j2, k2);
  this_seq->param_extents(di, dj, dk);
  if (i2-i1+1 != di || j2-j1+1 != dj || k2-k1+1 != dk) {
    std::cout << "min_params/max_params/param_extents functions returned inconsistent results"
              << std::endl;
    result = MB_FAILURE;
  }
  
  if (this_start != this_seq->get_start_handle()) {
    std::cout << "Start handle for sequence wrong." << std::endl;
    result = MB_FAILURE;
  }

  if (this_seq->boundary_complete() == false) {
    std::cout << "Element sequence didn't pass boundary_complete test." << std::endl;
    result = MB_FAILURE;
  }

  return result;
}
  
MBErrorCode evaluate_vertex_sequence(ScdVertexSeq *this_seq) 
{
  MBErrorCode result = MB_SUCCESS;
  
    // first get the parametric extents
  int imin, jmin, kmin, imax, jmax, kmax, itmp, jtmp, ktmp;
  this_seq->min_params(imin, jmin, kmin);
  this_seq->max_params(imax, jmax, kmax);

    // then the start vertex
  MBEntityHandle start_handle = this_seq->get_start_handle();
  
    // now evaluate all the vertices in forward and reverse
  MBEntityHandle tmp_handle, tmp_handle2;
  for (int i = imin; i <= imax; i++) {
    for (int j = jmin; j <= jmax; j++) {
      for (int k = kmin; k <= kmax; k++) {
          // compute what the vertex handle is supposed to be
        MBEntityHandle this_handle = start_handle + (i-imin) + (j-jmin)*(imax-imin+1) + 
          (k-kmin)*(jmax-jmin+1)*(imax-imin+1);

          // get_vertex variants
        tmp_handle = this_seq->get_vertex(i, j, k);

        if (this_seq->get_vertex(i, j, k) != this_handle) {
          std::cout << "vertex seq: get_vertex(i, j, k) didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
        
        tmp_handle2 = this_seq->get_vertex(HomCoord(i, j, k));
        if (this_seq->get_vertex(HomCoord(i, j, k)) != this_handle) {
          std::cout << "vertex seq: get_vertex(HomCoord) didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }

        MBErrorCode tmp_result = this_seq->get_params(tmp_handle, itmp, jtmp, ktmp);
        if (MB_SUCCESS != tmp_result || i != itmp || j != jtmp || k != ktmp) {
          std::cout << "vertex seq: get_params didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "; itmp, jtmp, ktmp = " 
                    << itmp << ", " << jtmp << ", " << ktmp << std::endl;
          result = MB_FAILURE;
        }

        if (!this_seq->contains(i, j, k) || !this_seq->contains(HomCoord(i,j,k))) {
          std::cout << "vertex seq: contains didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
      }
    }
  }

  return result;
}

MBErrorCode evaluate_element_sequence(ScdElementSeq *this_seq) 
{
  MBErrorCode result = MB_SUCCESS;
  
    // first get the parametric extents
  int imin, jmin, kmin, imax, jmax, kmax, itmp, jtmp, ktmp;
  this_seq->min_params(imin, jmin, kmin);
  this_seq->max_params(imax, jmax, kmax);

    // now evaluate all the vertices and elements in forward and reverse
  MBEntityHandle tmp_handle, tmp_handle2;
  for (int i = imin; i < imax; i++) {
    for (int j = jmin; j < jmax; j++) {
      for (int k = kmin; k < kmax; k++) {

          // get_vertex variants
        tmp_handle = this_seq->get_vertex(i, j, k);

        tmp_handle2 = this_seq->get_vertex(HomCoord(i, j, k));
        if (tmp_handle2 != tmp_handle) {
          std::cout << "element seq: get_vertex(HomCoord) and get_vertex(i,j,k) didn't return" << std::endl
                    << "consistent results, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }

          // get element variants
        tmp_handle = this_seq->get_element(i, j, k);
        
        tmp_handle2 = this_seq->get_element(HomCoord(i,j,k));
        if (tmp_handle2 != tmp_handle) {
          std::cout << "element seq: get_element(HomCoord) and get_element(i,j,k) didn't return" << std::endl
                    << "consistent results, i, j, k = " << i << ", " 
                    << j << ", " << k << "." << std::endl;
          result = MB_FAILURE;
        }
        
          // get_params
        MBErrorCode tmp_result = this_seq->get_params(tmp_handle, itmp, jtmp, ktmp);
        if (MB_SUCCESS != tmp_result || i != itmp || j != jtmp || k != ktmp) {
          std::cout << "element seq: get_params didn't work, i, j, k = " << i << ", " 
                    << j << ", " << k << "; itmp, jtmp, ktmp = " 
                    << itmp << ", " << jtmp << ", " << ktmp << std::endl;
          result = MB_FAILURE;
        }

        if (!this_seq->contains(i, j, k) || !this_seq->contains(HomCoord(i,j,k))) {
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
  int errors = 0;

    // first we need to make a new MBCore
  MBCore *gMB = new MBCore();
  
    // test creating and evaluating vertex sequences
  errors += test_vertex_seq(gMB);

    // test creating and evaluating element sequences
  errors += test_element_seq(gMB);

  if (errors > 0)
    std::cout << errors << " errors found." << std::endl;
  else
    std::cout << "All tests passed." << std::endl;
    
  return errors;
}

int test_vertex_seq(MBCore *gMB) 
{
    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();
  
  int errors = 0;
  MBEntityHandle oned_start, twod_start, threed_start;
  MBEntitySequence *dum_seq = NULL;
  ScdVertexSeq *oned_seq = NULL, *twod_seq = NULL, *threed_seq = NULL;
  
    // make a 1d sequence
  MBErrorCode result = seq_mgr->create_scd_sequence(-10, 0, 0, 10, 0, 0,
                                                     MBVERTEX, 1,
                                                     oned_start, dum_seq);
  if (NULL != dum_seq) oned_seq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  if (MB_FAILURE == result || oned_start == 0 || dum_seq == NULL ||
      oned_seq == NULL) {
    std::cout << "Problems creating a 1d sequence." << std::endl;
    errors++;
    return errors;
  }

    // check sequence data
  result = check_vertex_sequence(oned_seq, -10, 0, 0, 10, 0, 0, oned_start);
  if (MB_SUCCESS != result) {
    std::cout << "Sequence data didn't check out for 1d sequence." << std::endl;
    errors++;
  }

    // evaluate that sequence for all possible values
  result = evaluate_vertex_sequence(oned_seq);
  if (MB_SUCCESS != result) {
    std::cout << "Trouble evaluating 1d sequence." << std::endl;
    errors++;    
  }

    // make a 2d sequence
  dum_seq = NULL;
  result = seq_mgr->create_scd_sequence(-10, -10, 0, 10, 10, 0,
                                        MBVERTEX, 1,
                                        twod_start, dum_seq);
  if (NULL != dum_seq) twod_seq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  if (MB_FAILURE == result || twod_start == 0 || dum_seq == NULL ||
      twod_seq == NULL) {
    std::cout << "Problems creating a 2d sequence." << std::endl;
    errors++;
    return errors;
  }

    // check sequence data
  result = check_vertex_sequence(twod_seq, -10, -10, 0, 10, 10, 0, twod_start);
  if (MB_SUCCESS != result) {
    std::cout << "Sequence data didn't check out for 2d sequence." << std::endl;
    errors++;
  }

    // evaluate that sequence for all possible values
  result = evaluate_vertex_sequence(twod_seq);
  if (MB_SUCCESS != result) {
    std::cout << "Trouble evaluating 2d sequence." << std::endl;
    errors++;    
  }

    // make a 3d sequence
  dum_seq = NULL;
  result = seq_mgr->create_scd_sequence(-10, -10, -10, 10, 10, 10,
                                        MBVERTEX, 1,
                                        threed_start, dum_seq);
  if (NULL != dum_seq) threed_seq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  if (MB_FAILURE == result || threed_start == 0 || dum_seq == NULL ||
      threed_seq == NULL) {
    std::cout << "Problems creating a 3d sequence." << std::endl;
    errors++;
    return errors;
  }

    // check sequence data
  result = check_vertex_sequence(threed_seq, -10, -10, -10, 10, 10, 10, threed_start);
  if (MB_SUCCESS != result) {
    std::cout << "Sequence data didn't check out for 3d sequence." << std::endl;
    errors++;
  }

    // evaluate that sequence for all possible values
  result = evaluate_vertex_sequence(threed_seq);
  if (MB_SUCCESS != result) {
    std::cout << "Trouble evaluating 3d sequence." << std::endl;
    errors++;    
  }

  return errors;
}

int test_element_seq(MBCore *gMB) 
{
  int errors = 0;
  HomCoord TEST_MIN_PARAMS(0,0,0);
  HomCoord TEST_MAX_PARAMS(23,11,5);

    // TEST 1: single vertex sequence blocks, unity mapping
  errors += eseq_test1a(gMB, TEST_MIN_PARAMS, TEST_MAX_PARAMS);
  errors += eseq_test1b(gMB, TEST_MIN_PARAMS, TEST_MAX_PARAMS);
  errors += eseq_test1c(gMB, TEST_MIN_PARAMS, TEST_MAX_PARAMS);

    // TEST 2: composites, 0d difference between element block "owning"
    // vertex block and element block sharing vertex block
  errors += eseq_test2a(gMB, TEST_MIN_PARAMS, TEST_MAX_PARAMS);
  errors += eseq_test2b(gMB);
  errors += eseq_test2c(gMB);
  errors += eseq_test2d(gMB);

  return errors;
}

int eseq_test1a(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 1a: 1d single vertex seq block, min/max = (-10,0,0)/(10,0,0)
    // create vertex seq

  int errors = 0;

    // first get 1d min/max by resetting j, k components
  tmp_min[1] = tmp_min[2] = tmp_max[1] = tmp_max[2] = 0;

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();
  
  MBEntityHandle oned_start;
  MBEntitySequence *dum_seq;
  ScdVertexSeq *oned_seq = NULL;
  MBErrorCode result = seq_mgr->create_scd_sequence(tmp_min, tmp_max,
                                                     MBVERTEX, 1,
                                                     oned_start, dum_seq);
  if (NULL != dum_seq) oned_seq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && oned_start != 0 && dum_seq != NULL && oned_seq != NULL);

    // now create the element sequence
  MBEntityHandle eseq_start;
  ScdElementSeq *eseq = NULL;
  result = seq_mgr->create_scd_sequence(tmp_min, tmp_max,
                                        MBEDGE, 1,
                                        eseq_start, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && eseq_start != 0 && dum_seq != NULL && eseq != NULL);
  
    // add vertex seq to element seq
  result = eseq->add_vsequence(oned_seq,
                               tmp_min, tmp_min,
                               tmp_max, tmp_max,
                               tmp_min, tmp_min);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add vsequence to 1d single-block element sequence." << std::endl;
    errors++;
  }
  
    // check/evaluate element sequence
  result = check_element_sequence(eseq, tmp_min, tmp_max, eseq_start);
  if (MB_SUCCESS != result) {
    std::cout << "1d single-block element sequence didn't pass check." << std::endl;
    errors++;
  }
  
  result = evaluate_element_sequence(eseq);
  if (MB_SUCCESS != result) {
    std::cout << "1d single-block element sequence didn't evaluate correctly." << std::endl;
    errors++;
  }

  return errors;
}

int eseq_test1b(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 1b: 2d single vertex seq block, min/max = (-10,-5,0)/(10,5,0)

  int errors = 0;

    // first get 2d min/max by resetting k component
  tmp_min[2] = tmp_max[2] = 0;

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();
  
  MBEntityHandle twod_start;
  MBEntitySequence *dum_seq;
  ScdVertexSeq *twod_seq = NULL;
  MBErrorCode result = seq_mgr->create_scd_sequence(tmp_min, tmp_max,
                                                     MBVERTEX, 1,
                                                     twod_start, dum_seq);
  if (NULL != dum_seq) twod_seq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && twod_start != 0 && dum_seq != NULL && twod_seq != NULL);

    // now create the element sequence
  MBEntityHandle eseq_start;
  ScdElementSeq *eseq = NULL;
  result = seq_mgr->create_scd_sequence(tmp_min, tmp_max,
                                        MBQUAD, 1,
                                        eseq_start, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && eseq_start != 0 && dum_seq != NULL && eseq != NULL);
  
    // add vertex seq to element seq; first need to construct proper 3pt input (p1 is tmp_min)
  HomCoord p2(tmp_max.i(), tmp_min.j(), tmp_min.k());
  HomCoord p3(tmp_min.i(), tmp_max.j(), tmp_min.k());
  result = eseq->add_vsequence(twod_seq,
                               tmp_min, tmp_min,
                               p2, p2, p3, p3);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add vsequence to 2d single-block element sequence." << std::endl;
    errors++;
  }
  
    // check/evaluate element sequence
  result = check_element_sequence(eseq, tmp_min, tmp_max, eseq_start);
  if (MB_SUCCESS != result) {
    std::cout << "2d single-block element sequence didn't pass check." << std::endl;
    errors++;
  }
  
  result = evaluate_element_sequence(eseq);
  if (MB_SUCCESS != result) {
    std::cout << "2d single-block element sequence didn't evaluate correctly." << std::endl;
    errors++;
  }

  return errors;
}

int eseq_test1c(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 1c: 3d single vertex seq block, min/max = (-10,-5,-1)/(10,5,1)

  int errors = 0;

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();
  
  MBEntityHandle threed_start;
  MBEntitySequence *dum_seq;
  ScdVertexSeq *threed_seq = NULL;
  MBErrorCode result = seq_mgr->create_scd_sequence(tmp_min, tmp_max,
                                                     MBVERTEX, 1,
                                                     threed_start, dum_seq);
  if (NULL != dum_seq) threed_seq = dynamic_cast<ScdVertexSeq*>(dum_seq);
  assert (MB_FAILURE != result && threed_start != 0 && dum_seq != NULL && threed_seq != NULL);

    // now create the element sequence
  MBEntityHandle eseq_start;
  ScdElementSeq *eseq = NULL;
  result = seq_mgr->create_scd_sequence(tmp_min, tmp_max,
                                        MBHEX, 1,
                                        eseq_start, dum_seq);
  if (NULL != dum_seq) eseq = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && eseq_start != 0 && dum_seq != NULL && eseq != NULL);
  
    // add vertex seq to element seq; first need to construct proper 3pt input (p1 is tmp_min)
  HomCoord p2(tmp_max.i(), tmp_min.j(), tmp_min.k());
  HomCoord p3(tmp_min.i(), tmp_max.j(), tmp_min.k());
  result = eseq->add_vsequence(threed_seq,
                               tmp_min, tmp_min, p2, p2, p3, p3);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add vsequence to 3d single-block element sequence." << std::endl;
    errors++;
  }
  
    // check/evaluate element sequence
  result = check_element_sequence(eseq, tmp_min, tmp_max, eseq_start);
  if (MB_SUCCESS != result) {
    std::cout << "3d single-block element sequence didn't pass check." << std::endl;
    errors++;
  }
  
  result = evaluate_element_sequence(eseq);
  if (MB_SUCCESS != result) {
    std::cout << "3d single-block element sequence didn't evaluate correctly." << std::endl;
    errors++;
  }

  return errors;
}

int eseq_test2a(MBCore *gMB, HomCoord tmp_min, HomCoord tmp_max) 
{
    // TEST 2a: 1d composite block, 0d difference between owning/sharing blocks
    // create vertex seq
  
  ScdVertexSeq *vseq[3];
  ScdElementSeq *eseq[3];
  MBEntityHandle vstart[3], estart[3];
  
  int errors = create_1d_3_sequences(gMB, tmp_min, tmp_max,
                                     vseq, vstart, eseq, estart);
  if (0 != errors) return errors;

    // whew; that's done; now check and evaluate

    // first check to make sure the parameter spaces tack onto one another
  if (eseq[0]->min_params() != tmp_min ||
      eseq[0]->max_params() != eseq[1]->min_params() ||
      eseq[1]->max_params() != eseq[2]->min_params()) {
    std::cout << "In composite 1d eseq, parameter spaces don't line up." << std::endl;
    errors++;
  }

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    MBErrorCode result = check_element_sequence(eseq[i], eseq[i]->min_params(), eseq[i]->max_params(), 
                                                 estart[i]);
    if (MB_SUCCESS != result) {
      std::cout << "1d composite element sequence " << i << " didn't pass check." << std::endl;
      errors++;
    }
  
    result = evaluate_element_sequence(eseq[i]);
    if (MB_SUCCESS != result) {
      std::cout << "1d composite element sequence " << i << " didn't evaluate correctly." << std::endl;
      errors++;
    }
  }
  
  return errors;
}

int eseq_test2b(MBCore *gMB) 
{
    // TEST 2b: 2d composite block, 0d difference between owning/sharing blocks
    // create vertex seq
  
  ScdVertexSeq *vseq[3];
  ScdElementSeq *eseq[3];
  MBEntityHandle vstart[3], estart[3];
  
  int errors = create_2d_3_sequences(gMB, vseq, vstart, eseq, estart);
  if (0 != errors) return errors;

    // whew; that's done; now check and evaluate

    // first check to make sure the parameter spaces tack onto one another
  if (eseq[0]->max_params() != HomCoord(eseq[1]->min_params().i(),eseq[1]->max_params().j(),0) ||
      eseq[1]->max_params() != HomCoord(eseq[2]->min_params().i(),eseq[2]->max_params().j(),0) ||
        // cheat on the i value of eseq0, since it's a periodic bdy (we don't check for that, so you
        // may get different values for that parameter)
      eseq[2]->max_params() != HomCoord(eseq[2]->max_params().i(),eseq[0]->max_params().j(),0)) {
    std::cout << "In composite 2d eseq, parameter spaces don't line up." << std::endl;
    errors++;
  }

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    MBErrorCode result = check_element_sequence(eseq[i], eseq[i]->min_params(), eseq[i]->max_params(), 
                                                 estart[i]);
    if (MB_SUCCESS != result) {
      std::cout << "2d composite element sequence " << i << " didn't pass check." << std::endl;
      errors++;
    }
  
    result = evaluate_element_sequence(eseq[i]);
    if (MB_SUCCESS != result) {
      std::cout << "2d composite element sequence " << i << " didn't evaluate correctly." << std::endl;
      errors++;
    }
  }
  
  return errors;
}

int eseq_test2c(MBCore *gMB) 
{
    // TEST 2c: 2d composite block, 0d difference between owning/sharing blocks,
    // tri-valent shared vertex between the three blocks

  ScdVertexSeq *vseq[3];
  ScdElementSeq *eseq[3];
  MBEntityHandle vstart[3], estart[3];

    // interval settings: only 3 of them
  int int1 = 5, int2 = 15, int3 = 25;
  int errors = create_2dtri_3_sequences(gMB, int1, int2, int3, vseq, vstart, eseq, estart);
  if (0 != errors) return errors;

    // whew; that's done; now check and evaluate

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    MBErrorCode result = check_element_sequence(eseq[i], eseq[i]->min_params(), eseq[i]->max_params(), 
                                                 estart[i]);
    if (MB_SUCCESS != result) {
      std::cout << "2d tri-composite element sequence " << i << " didn't pass check." << std::endl;
      errors++;
    }
  
    result = evaluate_element_sequence(eseq[i]);
    if (MB_SUCCESS != result) {
      std::cout << "2d tri-composite element sequence " << i << " didn't evaluate correctly." << std::endl;
      errors++;
    }
  }
  
  return errors;
}

int eseq_test2d(MBCore *gMB) 
{
    // TEST 2d: 3d composite block, 0d difference between owning/sharing blocks,
    // tri-valent shared edge between the three blocks

  ScdVertexSeq *vseq[3];
  ScdElementSeq *eseq[3];
  MBEntityHandle vstart[3], estart[3];

    // interval settings: only 3 of them
  int int1 = 100, int2 = 100, int3 = 100, int4 = 100;
  int errors = create_3dtri_3_sequences(gMB, int1, int2, int3, int4, vseq, vstart, eseq, estart);
  if (0 != errors) return errors;

    // whew; that's done; now check and evaluate

    // check/evaluate element sequences
  for (int i = 0; i < 3; i++) {
    MBErrorCode result = check_element_sequence(eseq[i], eseq[i]->min_params(), eseq[i]->max_params(), 
                                                 estart[i]);
    if (MB_SUCCESS != result) {
      std::cout << "3d tri-composite element sequence " << i << " didn't pass check." << std::endl;
      errors++;
    }
  
    result = evaluate_element_sequence(eseq[i]);
    if (MB_SUCCESS != result) {
      std::cout << "3d tri-composite element sequence " << i << " didn't evaluate correctly." << std::endl;
      errors++;
    }
  }
  
  return errors;
}

int create_1d_3_sequences(MBCore *gMB,
                          HomCoord tmp_min, HomCoord tmp_max,
                          ScdVertexSeq **vseq, MBEntityHandle *vstart,
                          ScdElementSeq **eseq, MBEntityHandle *estart) 
{
  int errors = 0;
  
    // first get 1d min/max by resetting j, k components
  tmp_min[1] = tmp_min[2] = tmp_max[1] = tmp_max[2] = 0;

    // split the sequence in three in i parameter, at 1/2, 1/3 and 1/6, such that the 
    // total # vertices is the expected number from tmp_max - tmp_min)
  int idiff = (tmp_max[0] - tmp_min[0] + 1)/6;
  HomCoord vseq0_minmax[2] = {HomCoord(tmp_min), HomCoord(tmp_min[0]+3*idiff-1, tmp_max[1], tmp_max[2])};
  HomCoord vseq1_minmax[2] = {HomCoord(tmp_min), HomCoord(tmp_min[0]+2*idiff-1, tmp_max[1], tmp_max[2])};
  HomCoord vseq2_minmax[2] = {HomCoord(tmp_min), HomCoord(tmp_max[0]-5*idiff, tmp_max[1], tmp_max[2])};
  
    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();

    // create three vertex sequences
  MBEntitySequence *dum_seq;
  vseq[0] = vseq[1] = vseq[2] = NULL;
  

    // first vertex sequence 
  MBErrorCode result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
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
  eseq[0] = eseq[1] = eseq[2] = NULL;
  
    // create the first element sequence
  result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
                                        MBEDGE, 1,
                                        estart[0], dum_seq);
  if (NULL != dum_seq) eseq[0] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[0] != 0 && dum_seq != NULL && eseq[0] != NULL);
  
    // add first vertex seq to first element seq, forward orientation, unity transform
  result = eseq[0]->add_vsequence(vseq[0],
                                  vseq0_minmax[0], vseq0_minmax[0],
                                  vseq0_minmax[1], vseq0_minmax[1],
                                  vseq0_minmax[0], vseq0_minmax[0]);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add first vsequence to first element sequence in composite 1d eseq." << std::endl;
    errors++;
  }
  
    // create the second element sequence; make it use the second vseq in reverse and start
    // with the end vertex of the first sequence; parameterize it such that it tacks onto the 
    // end of the previous eseq
  result = seq_mgr->create_scd_sequence(HomCoord(vseq0_minmax[1].i(), 0, 0),
                                        HomCoord(1+vseq0_minmax[1].i()+vseq1_minmax[1].i()-vseq1_minmax[0].i(), 0, 0),
                                        MBEDGE, 1,
                                        estart[1], dum_seq);
  if (NULL != dum_seq) eseq[1] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[1] != 0 && dum_seq != NULL && eseq[1] != NULL);
  
    // add shared vertex from first vseq to this eseq; parameter space should be the same since
    // we're adding to that parameter space
  result = eseq[1]->add_vsequence(vseq[0],
                                  vseq0_minmax[0], vseq0_minmax[0],
                                  vseq0_minmax[1], vseq0_minmax[1],
                                  vseq0_minmax[0], vseq0_minmax[0],
                                  true,
                                  HomCoord(eseq[1]->min_params().i(),0,0),
                                  HomCoord(eseq[1]->min_params().i(),0,0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to second element sequence in composite 1d eseq." << std::endl;
    errors++;
  }
  
    // add second vseq to this eseq, but reversed; parameter space should be such that the
    // last vertex in the second vseq occurs first, and so on
  result = eseq[1]->add_vsequence(vseq[1],
                                  vseq1_minmax[1], eseq[1]->min_params()+HomCoord::unitv[0],
                                  vseq1_minmax[0], eseq[1]->max_params(),
                                  vseq1_minmax[1], eseq[1]->min_params()+HomCoord::unitv[0]);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add second vseq to second element sequence in composite 1d eseq." << std::endl;
    errors++;
  }
  
    // create the third element sequence; make it use the third vseq (forward sense) and start
    // with the start vertex of the second vseq; parameterize it such that it tacks onto the 
    // end of the previous eseq
  result = seq_mgr->create_scd_sequence(eseq[1]->max_params(),
                                        HomCoord(eseq[1]->max_params().i()+1+vseq2_minmax[1].i()-vseq2_minmax[0].i(),0,0),
                                        MBEDGE, 1,
                                        estart[2], dum_seq);
  if (NULL != dum_seq) eseq[2] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[2] != 0 && dum_seq != NULL && eseq[2] != NULL);
  
    // add shared vertex from second vseq to this eseq; parameter space mapping such that we get
    // first vertex only of that vseq
  result = eseq[2]->add_vsequence(vseq[1],
                                  vseq0_minmax[0], eseq[2]->min_params(),
                                  vseq0_minmax[0]+HomCoord::unitv[0], eseq[2]->min_params()-HomCoord::unitv[0],
                                  vseq0_minmax[0], eseq[2]->min_params(),
                                  true,
                                  eseq[2]->min_params(),
                                  eseq[2]->min_params());
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to third element sequence in composite 1d eseq." << std::endl;
    errors++;
  }
  
    // add third vseq to this eseq, forward orientation
  result = eseq[2]->add_vsequence(vseq[2],
                                  vseq2_minmax[0], eseq[2]->min_params()+HomCoord::unitv[0],
                                  vseq2_minmax[1], eseq[2]->max_params(),
                                  vseq1_minmax[0], eseq[2]->min_params()+HomCoord::unitv[0]);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add third vseq to third element sequence in composite 1d eseq." << std::endl;
    errors++;
  }

  return errors;
}

int create_2d_3_sequences(MBCore *gMB,
                          ScdVertexSeq **vseq, MBEntityHandle *vstart,
                          ScdElementSeq **eseq, MBEntityHandle *estart) 
{
    // create 3 rectangular esequences attached end to end and back (periodic); vsequences are 
    // assorted orientations, esequences have globally-consistent (periodic in i) parameter space
  int errors = 0;

    // set vseq parametric spaces directly
  HomCoord vseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(5,5,0)};
  HomCoord vseq1_minmax[2] = {HomCoord(-2,4,0), HomCoord(8,9,0)};
  HomCoord vseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(8,5,0)};

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();

    // create three vertex sequences
  MBEntitySequence *dum_seq;
  vseq[0] = vseq[1] = vseq[2] = NULL;

    // first vertex sequence 
  MBErrorCode result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
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
  eseq[0] = eseq[1] = eseq[2] = NULL;
  
    // create the first element sequence
  result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
                                        MBQUAD, 1,
                                        estart[0], dum_seq);
  if (NULL != dum_seq) eseq[0] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[0] != 0 && dum_seq != NULL && eseq[0] != NULL);
  
    // add first vertex seq to first element seq, forward orientation, unity transform
  result = eseq[0]->add_vsequence(vseq[0],
                                    // p1: imin,jmin
                                  vseq0_minmax[0], vseq0_minmax[0],
                                    // p2: imax,jmin
                                  HomCoord(vseq0_minmax[1].i(), vseq0_minmax[0].j(),0),
                                  HomCoord(vseq0_minmax[1].i(), vseq0_minmax[0].j(),0),
                                    // p3: imin,jmax
                                  HomCoord(vseq0_minmax[0].i(), vseq0_minmax[1].j(),0),
                                  HomCoord(vseq0_minmax[0].i(), vseq0_minmax[1].j(),0));
  
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add first vsequence to first element sequence in composite 2d eseq." << std::endl;
    errors++;
  }
  
    // create the second element sequence; make it use the right side of the first vseq; 
    // parameterize it such that it tacks onto imax of the previous eseq
  result = seq_mgr->create_scd_sequence(HomCoord(eseq[0]->max_params().i(), eseq[0]->min_params().j(), 0),
                                        HomCoord(vseq0_minmax[1].i()+1+vseq1_minmax[1].i()-vseq1_minmax[0].i(), 
                                                 eseq[0]->max_params().j(), 0),
                                        MBQUAD, 1,
                                        estart[1], dum_seq);
  if (NULL != dum_seq) eseq[1] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[1] != 0 && dum_seq != NULL && eseq[1] != NULL);
  
    // add shared side from first vseq to this eseq; parameter space should be the same since
    // we're adding to that parameter space
  result = eseq[1]->add_vsequence(vseq[0],
                                  vseq0_minmax[0], vseq0_minmax[0],
                                  vseq0_minmax[1], vseq0_minmax[1],
                                  vseq0_minmax[0], vseq0_minmax[0],
                                    // set bb such that it's the right side of the vseq, left of local eseq
                                  true,
                                  eseq[1]->min_params(),
                                  HomCoord(eseq[1]->min_params().i(),eseq[1]->max_params().j(),0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to second element sequence in composite 2d eseq." << std::endl;
    errors++;
  }
  
    // add second vseq to this eseq, with different orientation but all of it (no bb input)
  result = eseq[1]->add_vsequence(vseq[1],
                                    // p1: one right of top left
                                  vseq1_minmax[0], 
                                  HomCoord(eseq[1]->min_params().i()+1, eseq[1]->max_params().j(), 0),
                                    // p2: one right from p1
                                  vseq1_minmax[0]+HomCoord::unitv[0], 
                                  HomCoord(eseq[1]->min_params().i()+2, eseq[1]->max_params().j(), 0),
                                    // p3: one down from p1
                                  vseq1_minmax[0]+HomCoord::unitv[1], 
                                  HomCoord(eseq[1]->min_params().i()+1, eseq[1]->max_params().j()-1, 0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add second vseq to second element sequence in composite 2d eseq." << std::endl;
    errors++;
  }
  
    // create the third element sequence; make it use the third vseq (middle) as well as the side of the
    // second sequence (left) and a side of the first sequence (right); parameterize it such that it tacks onto the 
    // end of the previous eseq and the beginning of the 1st sequence (i.e. periodic in i)
  result = seq_mgr->create_scd_sequence(HomCoord(eseq[1]->max_params().i(), eseq[1]->min_params().j(),0),
                                          // add one extra for each of left and right sides
                                        HomCoord(eseq[1]->max_params().i()+1+vseq2_minmax[1].i()-vseq2_minmax[0].i()+1,
                                                 eseq[1]->max_params().j(),0),
                                        MBEDGE, 1,
                                        estart[2], dum_seq);
  if (NULL != dum_seq) eseq[2] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[2] != 0 && dum_seq != NULL && eseq[2] != NULL);
  
    // add shared side from second vseq to this eseq; parameter space mapping such that we get
    // a side only of that vseq
  result = eseq[2]->add_vsequence(vseq[1],
                                    // p1: bottom left
                                  vseq1_minmax[1], eseq[2]->min_params(),
                                    // p2: one right from p1
                                  vseq1_minmax[1]+HomCoord::unitv[0],
                                  eseq[2]->min_params()+HomCoord::unitv[0],
                                    // p3: one up
                                  vseq1_minmax[1]-HomCoord::unitv[1],
                                  eseq[2]->min_params()+HomCoord::unitv[1],
                                    // bb input such that we only get left side of eseq parameter space
                                  true,
                                  eseq[2]->min_params(),
                                  HomCoord(eseq[2]->min_params().i(), eseq[2]->max_params().j(), 0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add left shared vsequence to third element sequence in composite 2d eseq." << std::endl;
    errors++;
  }
  
    // add shared side from first vseq to this eseq; parameter space mapping such that we get
    // a side only of that vseq
  result = eseq[2]->add_vsequence(vseq[0],
                                    // p1: bottom right
                                  vseq1_minmax[0], HomCoord(eseq[2]->max_params().i(), eseq[2]->min_params().j(),0),
                                    // p2: one right from p1
                                  vseq1_minmax[0]+HomCoord::unitv[0], 
                                  HomCoord(eseq[2]->max_params().i()+1, eseq[2]->min_params().j(),0),
                                    // p3: one up from p1
                                  vseq1_minmax[0]+HomCoord::unitv[1],
                                  HomCoord(eseq[2]->max_params().i(), eseq[2]->min_params().j()+1,0),
                                    // bb input such that we only get left side of eseq parameter space
                                  true,
                                  HomCoord(eseq[2]->max_params().i(), eseq[2]->min_params().j(),0),
                                  eseq[2]->max_params());
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add left shared vsequence to third element sequence in composite 2d eseq." << std::endl;
    errors++;
  }
  

    // add third vseq to this eseq
  result = eseq[2]->add_vsequence(vseq[2],
                                    // p1: top right and left one
                                  vseq2_minmax[0], eseq[2]->max_params()-HomCoord::unitv[0],
                                    // p2: one left of p1
                                  vseq2_minmax[0]+HomCoord::unitv[0], eseq[2]->max_params()-HomCoord::unitv[0]*2,
                                    // p3: one down from p1
                                  vseq2_minmax[0]+HomCoord::unitv[1], eseq[2]->max_params()-HomCoord::unitv[0] -
                                  HomCoord::unitv[1]);
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add third vseq to third element sequence in composite 2d eseq." << std::endl;
    errors++;
  }

  return errors;
}

int create_2dtri_3_sequences(MBCore *gMB,
                             const int int1, const int int2, const int int3,
                             ScdVertexSeq **vseq, MBEntityHandle *vstart,
                             ScdElementSeq **eseq, MBEntityHandle *estart) 
{
    // create 3 rectangular esequences arranged such that the all share a common (tri-valent) corner;
    // orient each region such that its origin is at the tri-valent corner and the k direction is
    // out of the page
    //
    // int1 and int2 controls the i and j intervals in region 0, int3 follows from that.

    // input is 3 interval settings controlling the 3 degrees of freedom on the interfacesp
  int errors = 0;

    // set vseq parametric spaces directly from int1-3
    // use 0-based parameterization on vseq's just for fun, which means we'll have to transform into
    // eseq system
  HomCoord vseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,0)};
  HomCoord vseq1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3-1,int1,0)};
  HomCoord vseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2-1,int3-1,0)};

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();

    // create three vertex sequences
  MBEntitySequence *dum_seq;
  vseq[0] = vseq[1] = vseq[2] = NULL;

    // first vertex sequence 
  MBErrorCode result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
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

    // set eseq parametric spaces directly from int1-3
    // use 0-based parameterization on eseq's just for fun, which means we'll have to transform into
    // eseq system
  HomCoord eseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,0)};
  HomCoord eseq1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3,int1,0)};
  HomCoord eseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2,int3,0)};

  eseq[0] = eseq[1] = eseq[2] = NULL;
  
    // create the first element sequence
  result = seq_mgr->create_scd_sequence(eseq0_minmax[0], eseq0_minmax[1], 
                                        MBQUAD, 1,
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
    std::cout << "Couldn't add first vsequence to first element sequence in tri-composite 2d eseq." << std::endl;
    errors++;
  }
  
    // create the second element sequence
  result = seq_mgr->create_scd_sequence(eseq1_minmax[0], eseq1_minmax[1],
                                        MBQUAD, 1,
                                        estart[1], dum_seq);
  if (NULL != dum_seq) eseq[1] = dynamic_cast<ScdElementSeq*>(dum_seq);
  assert (MB_FAILURE != result && estart[1] != 0 && dum_seq != NULL && eseq[1] != NULL);
  
    // add shared side from first vseq to this eseq, with bb to get just the line
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
                                  HomCoord(eseq[1]->min_params().i(),eseq[1]->max_params().j(),0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to second element sequence in tri-composite 2d eseq." << std::endl;
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
    std::cout << "Couldn't add second vseq to second element sequence in tri-composite 2d eseq." << std::endl;
    errors++;
  }
  
    // create the third element sequence
  result = seq_mgr->create_scd_sequence(eseq2_minmax[0], eseq2_minmax[1],
                                        MBQUAD, 1,
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
                                  HomCoord(eseq[2]->min_params().i(), eseq[2]->max_params().j(), 0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add shared vsequence to third element sequence in tri-composite 2d eseq." << std::endl;
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
                                  HomCoord(eseq2_minmax[1].i(), eseq2_minmax[0].j(),0));
  if (MB_SUCCESS != result) {
    std::cout << "Couldn't add left shared vsequence to third element sequence in tri-composite 2d eseq." << std::endl;
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
    std::cout << "Couldn't add third vseq to third element sequence in tri-composite 2d eseq." << std::endl;
    errors++;
  }

  return errors;
}

int create_3dtri_3_sequences(MBCore *gMB,
                             const int int1, const int int2, const int int3, const int int4,
                             ScdVertexSeq **vseq, MBEntityHandle *vstart,
                             ScdElementSeq **eseq, MBEntityHandle *estart) 
{
    // create 3 brick esequences arranged such that the all share a common (tri-valent) edge;
    // orient each region similarly to the 2dtri_3_esequences test problem, swept into 3d in the 
    // positive k direction.  This direction is divided into int4 intervals
    //
    // int1 and int2 controls the i and j intervals in region 0, int3 follows from that; int4
    // divides the k axis

    // input is 4 interval settings controlling the 4 degrees of freedom on the interfacesp
  int errors = 0;

    // set vseq parametric spaces directly from int1-4
    // use 0-based parameterization on vseq's just for fun, which means we'll have to transform into
    // eseq system
  HomCoord vseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,int4)};
  HomCoord vseq1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3-1,int1,int4)};
  HomCoord vseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2-1,int3-1,int4)};

    // get the seq manager from gMB
  EntitySequenceManager *seq_mgr = gMB->sequence_manager();

    // create three vertex sequences
  MBEntitySequence *dum_seq;
  vseq[0] = vseq[1] = vseq[2] = NULL;

    // first vertex sequence 
  MBErrorCode result = seq_mgr->create_scd_sequence(vseq0_minmax[0], vseq0_minmax[1],
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

    // set eseq parametric spaces directly from int1-4
    // use 0-based parameterization on eseq's just for fun, which means we'll have to transform into
    // eseq system
  HomCoord eseq0_minmax[2] = {HomCoord(0,0,0), HomCoord(int1,int2,int4)};
  HomCoord eseq1_minmax[2] = {HomCoord(0,0,0), HomCoord(int3,int1,int4)};
  HomCoord eseq2_minmax[2] = {HomCoord(0,0,0), HomCoord(int2,int3,int4)};

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

