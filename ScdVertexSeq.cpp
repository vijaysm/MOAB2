#include "ScdVertexSeq.hpp"

    //! constructor
ScdVertexSeq::ScdVertexSeq(EntitySequenceManager* seq_mgr,
                           const MBEntityHandle start_vertex, 
                           const int imin, const int jmin, const int kmin,
                           const int imax, const int jmax, const int kmax) 
    : VertexEntitySequence(seq_mgr, start_vertex, 
                           (imax-imin+1)*(jmax-jmin+1)*(kmax-kmin+1), 
                           true)
{
    // need to have meaningful parameters
  assert(imax >= imin && jmax >= jmin && kmax >= kmin);
  
  vertexParams[0] = HomCoord(imin, jmin, kmin);
  vertexParams[1] = HomCoord(imax, jmax, kmax);
  vertexParams[2] = HomCoord(1,1,1);
  
  dIJK[0] = imax-imin+1; dIJK[1] = jmax-jmin+1; dIJK[2] = kmax-kmin+1;
  dIJKm1[0] = dIJK[0]-1;
  dIJKm1[1] = dIJK[1]-1;
  dIJKm1[2] = dIJK[2]-1;
}
  
MBEntityHandle ScdVertexSeq::get_unused_handle()
{
  assert(false);
  return 0;
}


