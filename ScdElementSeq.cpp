#include "ScdElementSeq.hpp"
#include "ScdVertexSeq.hpp"
#include "MBInterface.hpp"
#include "MBReadUtilIface.hpp"
#include "MBCN.hpp"
#include "MBInternals.hpp"

ScdElementSeq::ScdElementSeq(EntitySequenceManager *seq_mgr,
                             MBEntityHandle start_handle,
                             const int imin, const int jmin, const int kmin,
                             const int imax, const int jmax, const int kmax) 
    : ElementEntitySequence(seq_mgr, start_handle,
                            (imax-imin)*(jmax-jmin)*(kmax-kmin),
                            MBCN::VerticesPerEntity(TYPE_FROM_HANDLE(start_handle)), 
                            true, false)
  
{
    // need to have meaningful parameters
  assert(imax >= imin && jmax >= jmin && kmax >= kmin);

    // correct num entities if necessary
  if (mNumEntities == 0) {
    int this_dim = MBCN::Dimension(TYPE_FROM_HANDLE(start_handle));
    if (this_dim == 1) mNumEntities = imax - imin;
    else if (this_dim == 2) mNumEntities = (imax - imin)*(jmax - jmin);

      // if neither of the previous two tests passed, it's an error
    else assert(false);

    assert(mNumEntities > 0);
  }
    
  elementParams[0] = HomCoord(imin, jmin, kmin);
  elementParams[1] = HomCoord(imax, jmax, kmax);
  elementParams[2] = HomCoord(1, 1, 1);
  
    // assign and compute parameter stuff
  dIJK[0] = elementParams[1][0] - elementParams[0][0] + 1;
  dIJK[1] = elementParams[1][1] - elementParams[0][1] + 1;
  dIJK[2] = elementParams[1][2] - elementParams[0][2] + 1;
  dIJKm1[0] = dIJK[0] - 1;
  dIJKm1[1] = dIJK[1] - 1;
  dIJKm1[2] = dIJK[2] - 1;
}

ScdElementSeq::~ScdElementSeq() 
{
}

bool ScdElementSeq::boundary_complete() const
{
    // test the bounding vertex sequences to see if they fully define the
    // vertex parameter space for this rectangular block of elements

  int p;
  std::vector<VertexSeqRef> minlist, maxlist;

    // pseudo code:
    // for each vertex sequence v:
  for (std::vector<VertexSeqRef>::const_iterator vseq = vertexSeqRefs.begin();
       vseq != vertexSeqRefs.end(); vseq++)
  {
    //   test min corner mincorner:
    bool mincorner = true;
    //   for each p = (i-1,j,k), (i,j-1,k), (i,j,k-1):
    for (p = 0; p < 3; p++) {

    //     for each vsequence v' != v:
      for (std::vector<VertexSeqRef>::const_iterator othervseq = vertexSeqRefs.begin();
           othervseq != vertexSeqRefs.end(); othervseq++) 
      {
        if (othervseq == vseq) continue;        
    //       if v.min-p contained in v'
        if ((*othervseq).contains((*vseq).minmax[0]-HomCoord::unitv[p])) {
    //         mincorner = false
          mincorner = false;
          break;
        }
      }
      if (!mincorner) break;
    }
  
    bool maxcorner = true;
    //   for each p = (i-1,j,k), (i,j-1,k), (i,j,k-1):
    for (p = 0; p < 3; p++) {

    //     for each vsequence v' != v:
      for (std::vector<VertexSeqRef>::const_iterator othervseq = vertexSeqRefs.begin();
           othervseq != vertexSeqRefs.end(); othervseq++) 
      {
        if (othervseq == vseq) continue;        
    //       if v.max+p contained in v'
        if ((*othervseq).contains((*vseq).minmax[1]+HomCoord::unitv[p])) {
    //         maxcorner = false
          maxcorner = false;
          break;
        }
      }
      if (!maxcorner) break;
    }

    //   if mincorner add to min corner list minlist
    if (mincorner) minlist.push_back(*vseq);
    //   if maxcorner add to max corner list maxlist
    if (maxcorner) maxlist.push_back(*vseq);
  }
  
    // 
    // if minlist.size = 1 & maxlist.size = 1 & minlist[0] = esequence.min &
    //         maxlist[0] = esequence.max+(1,1,1)
  if (minlist.size() == 1 && maxlist.size() == 1 &&
      minlist[0].minmax[0] == elementParams[0] && 
      maxlist[0].minmax[1] == elementParams[1])
      //   complete
    return true;
    // else

  return false;
}

