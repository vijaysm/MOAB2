#ifndef SCDELEMENTSEQ
#define SCDELEMENTSEQ

//
// Class: ScdElementSeq
//
// Purpose: represent a rectangular element of mesh
//
// A ScdElement represents a rectangular element of mesh, including both vertices and
// elements, and the parametric space used to address that element.  Vertex data,
// i.e. coordinates, may not be stored directly in the element, but the element returns
// information about the vertex handles of vertices in the element.  Vertex and element
// handles associated with the element are each contiguous.

#include "EntitySequence.hpp"
#include "HomXform.hpp"
#include "MBCN.hpp"
#include "ScdVertexSeq.hpp"
#include "MBInternals.hpp"

#include <vector>

  //! structure to hold references to bounding vertex blocks
class VertexSeqRef
{
private:
  HomCoord minmax[2];
  HomXform xform, invXform;
  ScdVertexSeq *srcSeq;
public:
  friend class ScdElementSeq;
  
  VertexSeqRef(const HomCoord &min, const HomCoord &max,
               const HomXform &tmp_xform, ScdVertexSeq *this_seq);
    
  bool contains(const HomCoord &coords) const;
};

class ScdElementSeq : public ElementEntitySequence
{

private:

    //! parameter min/max/stride, in homogeneous coords ijkh
  HomCoord elementParams[3];

    //! difference between max and min params plus one (i.e. # VERTICES in
    //! each parametric direction)
  int dIJK[3];
  
    //! difference between max and min params (i.e. # ELEMENTS in
    //! each parametric direction)
  int dIJKm1[3];

    //! bare constructor, so compiler doesn't create one for me
  ScdElementSeq();

    //! list of bounding vertex blocks
  std::vector<VertexSeqRef> vertexSeqRefs;

public:

    //! constructor
  ScdElementSeq(EntitySequenceManager *seq_mgr,
                MBEntityHandle start_handle,
                const int imin, const int jmin, const int kmin,
                const int imax, const int jmax, const int kmax);
  
  virtual ~ScdElementSeq();

    //! get handle of vertex at i, j, k
  MBEntityHandle get_vertex(const int i, const int j, const int k) const;
  
    //! get handle of vertex at homogeneous coords
  inline MBEntityHandle get_vertex(const HomCoord &coords) const;
  
    //! get handle of element at i, j, k
  MBEntityHandle get_element(const int i, const int j, const int k) const;
  
    //! get handle of element at homogeneous coords
  MBEntityHandle get_element(const HomCoord &coords) const;
  
    //! get min params for this element
  void min_params(HomCoord &coords) const;
  void min_params(int &i, int &j, int &k) const;
  const HomCoord &min_params() const;

    //! get max params for this element
  void max_params(HomCoord &coords) const;
  void max_params(int &i, int &j, int &k) const;
  const HomCoord &max_params() const;
  
    //! get the number of vertices in each direction, inclusive
  void param_extents(int &di, int &dj, int &dk) const;

    //! given a handle, get the corresponding parameters
  MBErrorCode get_params(const MBEntityHandle ehandle,
                          int &i, int &j, int &k) const;
  
    //! convenience functions for parameter extents
  int i_min() const {return (elementParams[0].hom_coord())[0];}
  int j_min() const {return (elementParams[0].hom_coord())[1];}
  int k_min() const {return (elementParams[0].hom_coord())[2];}
  int i_max() const {return (elementParams[1].hom_coord())[0];}
  int j_max() const {return (elementParams[1].hom_coord())[1];}
  int k_max() const {return (elementParams[1].hom_coord())[2];}

    //! test the bounding vertex sequences and determine whether they fully
    //! define the vertices covering this element block's parameter space
  bool boundary_complete() const;

    //! test whether this sequence contains these parameters
  bool contains(const int i, const int j, const int k) const;
  inline bool contains(const HomCoord &coords) const;
  
    // from parent class
  virtual MBEntityHandle get_unused_handle();

  virtual MBErrorCode get_connectivity(MBEntityHandle entity, 
                                       std::vector<MBEntityHandle>& connectivity,
                                       const bool topological_connectivity = false) const;

  virtual MBErrorCode get_connectivity(MBEntityHandle entity, 
                                       const MBEntityHandle*& connectivity,
                                       int &num_vertices,
                                       const bool topological_connectivity = false) const;

    //! get connectivity of an entity given entity's parameters
  inline MBErrorCode get_params_connectivity(const int i, const int j, const int k,
                                       std::vector<MBEntityHandle>& connectivity) const;
  
    //! add a vertex seq ref to this element sequence;
    //! if bb_input is true, bounding box (in eseq-local coords) of vseq being added 
    //! is input in bb_min and bb_max (allows partial sharing of vseq rather than the whole
    //! vseq); if it's false, the whole vseq is referenced and the eseq-local coordinates
    //! is computed from the transformed bounding box of the vseq
  MBErrorCode add_vsequence(ScdVertexSeq *vseq, 
                             const HomCoord &p1, const HomCoord &q1,
                             const HomCoord &p2, const HomCoord &q2,
                             const HomCoord &p3, const HomCoord &q3,
                             bool bb_input = false,
                             const HomCoord &bb_min = HomCoord::unitv[0],
                             const HomCoord &bb_max = HomCoord::unitv[0]);
  
  virtual MBErrorCode set_connectivity(MBEntityHandle entity, const MBEntityHandle *conn,
                                const int num_vertices);

  virtual void get_entities(MBRange& ) const;
  
  virtual MBErrorCode split(MBEntityHandle , 
                             MBEntitySequence*& );
  
  virtual MBErrorCode convert_realloc(bool& , bool& , bool& , 
                                      MBCore*, MBTag );
  
};

inline MBEntityHandle ScdElementSeq::get_element(const int i, const int j, const int k) const
{
  return mStartEntityHandle + (i-i_min()) + (j-j_min())*dIJKm1[0] + (k-k_min())*dIJKm1[0]*dIJKm1[1];
}
  
inline MBEntityHandle ScdElementSeq::get_element(const HomCoord &coord) const
{
  return get_element(coord.i(), coord.j(), coord.k());
}
  
  //! get min params for this element
inline void ScdElementSeq::min_params(HomCoord &coords) const
{
  coords = elementParams[0];
}

inline void ScdElementSeq::min_params(int &i, int &j, int &k) const
{
  i = elementParams[0].i();
  j = elementParams[0].j();
  k = elementParams[0].k();
}

inline const HomCoord &ScdElementSeq::min_params() const
{
  return elementParams[0];
}

//! get max params for this element
inline void ScdElementSeq::max_params(HomCoord &coords) const
{
  coords = elementParams[1];
}
  
inline void ScdElementSeq::max_params(int &i, int &j, int &k) const
{
  i = elementParams[1].i();
  j = elementParams[1].j();
  k = elementParams[1].k();
}

inline const HomCoord &ScdElementSeq::max_params() const
{
  return elementParams[1];
}

  //! get the number of vertices in each direction, inclusive
inline void ScdElementSeq::param_extents(int &di, int &dj, int &dk) const
{
  di = dIJK[0];
  dj = dIJK[1];
  dk = dIJK[2];
}

inline MBErrorCode ScdElementSeq::get_params(const MBEntityHandle ehandle,
                                              int &i, int &j, int &k) const
{
  if (TYPE_FROM_HANDLE(ehandle) != TYPE_FROM_HANDLE(mStartEntityHandle)) return MB_FAILURE;

  int hdiff = ehandle - mStartEntityHandle;

    // use double ?: test below because on some platforms, both sides of the : are
    // evaluated, and if dIJKm1[1] is zero, that'll generate a divide-by-zero
  k = (dIJKm1[1] > 0 ? hdiff / (dIJKm1[1] > 0 ? dIJKm1[0]*dIJKm1[1] : 1) : 0);
  j = (hdiff - (k*dIJKm1[0]*dIJKm1[1])) / dIJKm1[0];
  i = hdiff % dIJKm1[0];

  k += elementParams[0].k();
  j += elementParams[0].j();
  i += elementParams[0].i();

  return (ehandle >= mStartEntityHandle &&
          ehandle < mStartEntityHandle+number_entities() &&
          i >= i_min() && i <= i_max() &&
          j >= j_min() && j <= j_max() &&
          k >= k_min() && k <= k_max()) ? MB_SUCCESS : MB_FAILURE;
}

inline bool ScdElementSeq::contains(const int i, const int j, const int k) const 
{
  return contains(HomCoord(i, j, k));
}

inline bool ScdElementSeq::contains(const HomCoord &temp) const 
{
    // upper bound is < instead of <= because element params max is one less
    // than vertex params max
  return (temp >= elementParams[0] && temp < elementParams[1]);
}

inline MBEntityHandle ScdElementSeq::get_vertex(const int i, const int j, const int k) const
{
  return get_vertex(HomCoord(i,j,k));
}

inline MBEntityHandle ScdElementSeq::get_unused_handle() 
{
  return MB_FAILURE;
}

inline MBErrorCode ScdElementSeq::set_connectivity(MBEntityHandle , 
                                                    const MBEntityHandle *,
                                                    const int )
{
  return MB_FAILURE;
}

inline void ScdElementSeq::get_entities(MBRange& range) const
{
  range.insert(mStartEntityHandle, mStartEntityHandle+mNumEntities-1);
}
  
inline MBErrorCode ScdElementSeq::split(MBEntityHandle , 
                                        MBEntitySequence*& )
{
  return MB_FAILURE;
}

  // reallocated the sequence to hold extra/less nodes, pass in what you want, 
  // and will return whether it needed
  // reallocate space for those nodes
inline MBErrorCode ScdElementSeq::convert_realloc(bool& , bool& , bool& , 
                                                  MBCore*, MBTag )
{
  return MB_NOT_IMPLEMENTED;
}
  
inline bool VertexSeqRef::contains(const HomCoord &coords) const 
{
  return (minmax[0] <= coords && minmax[1] >= coords);
}

inline VertexSeqRef::VertexSeqRef(const HomCoord &this_min, const HomCoord &this_max,
                                  const HomXform &tmp_xform, ScdVertexSeq *this_seq)
    : xform(tmp_xform), invXform(tmp_xform.inverse()), srcSeq(this_seq)
{
  minmax[0] = HomCoord(this_min);
  minmax[1] = HomCoord(this_max); 
}

inline MBErrorCode ScdElementSeq::get_connectivity(MBEntityHandle entity, 
                                                   std::vector<MBEntityHandle>& connectivity,
                                                   const bool) const
{
    // get parameter values
  int i, j, k;
  MBErrorCode result = get_params(entity, i, j, k);
  if (MB_SUCCESS != result) return result;
  
  return get_params_connectivity(i, j, k, connectivity);
}

inline MBErrorCode ScdElementSeq::get_connectivity(MBEntityHandle , 
                                                   const MBEntityHandle*& ,
                                                   int &,
                                                   const bool) const 
{
    // this version of get_connectivity isn't supported yet!
  return MB_NOT_IMPLEMENTED;
}

inline MBEntityHandle ScdElementSeq::get_vertex(const HomCoord &coords) const
{
   for (std::vector<VertexSeqRef>::const_iterator it = vertexSeqRefs.begin();
        it != vertexSeqRefs.end(); it++) {
     if ((*it).minmax[0] <= coords && (*it).minmax[1] >= coords) {
         // first get the vertex block-local parameters
       HomCoord local_coords = coords / (*it).xform;
    
      // now get the vertex handle for those coords
       return (*it).srcSeq->get_vertex(local_coords);
     }
   }
   
     // got here, it's an error
   return 0;
}

inline MBErrorCode ScdElementSeq::add_vsequence(ScdVertexSeq *vseq, 
                                                 const HomCoord &p1, const HomCoord &q1,
                                                 const HomCoord &p2, const HomCoord &q2, 
                                                 const HomCoord &p3, const HomCoord &q3,
                                                 bool bb_input,
                                                 const HomCoord &bb_min,
                                                 const HomCoord &bb_max)
{
    // compute the transform given the vseq-local parameters and the mapping to
    // this element sequence's parameters passed in minmax
  HomXform M;
  M.three_pt_xform(p1, q1, p2, q2, p3, q3);
  
    // min and max in element seq's parameter system may not be same as those in 
    // vseq's system, so need to take min/max

  HomCoord minmax[2];
  if (bb_input) {
    minmax[0] = bb_min;
    minmax[1] = bb_max;
  }
  else {
    minmax[0] = vseq->min_params() * M;
    minmax[1] = vseq->max_params() * M;
  }
  
    // check against other vseq's to make sure they don't overlap
  for (std::vector<VertexSeqRef>::const_iterator vsit = vertexSeqRefs.begin();
       vsit != vertexSeqRefs.end(); vsit++) 
    if ((*vsit).contains(minmax[0]) || (*vsit).contains(minmax[1])) 
      return MB_FAILURE;
    
#define MIN(a,b) ((a <= b) ? a : b)
#define MAX(a,b) ((a >= b) ? a : b)
  HomCoord tmp_min(MIN(minmax[0].i(), minmax[1].i()), 
                   MIN(minmax[0].j(), minmax[1].j()), 
                   MIN(minmax[0].k(), minmax[1].k()));
  HomCoord tmp_max(MAX(minmax[0].i(), minmax[1].i()), 
                   MAX(minmax[0].j(), minmax[1].j()), 
                   MAX(minmax[0].k(), minmax[1].k()));

  
    // set up a new vertex sequence reference
  VertexSeqRef tmp_seq_ref(tmp_min, tmp_max, M, vseq);

    // add to the list
  vertexSeqRefs.push_back(tmp_seq_ref);
  
  return MB_SUCCESS;
}

inline MBErrorCode ScdElementSeq::get_params_connectivity(const int i, const int j, const int k,
                                                           std::vector<MBEntityHandle>& connectivity) const
{
  if (contains(i, j, k) == false) return MB_FAILURE;
  
  connectivity.push_back(get_vertex(i, j, k));
  connectivity.push_back(get_vertex(i+1, j, k));
  if (MBCN::Dimension(get_type()) < 2) return MB_SUCCESS;
  connectivity.push_back(get_vertex(i+1, j+1, k));
  connectivity.push_back(get_vertex(i, j+1, k));
  if (MBCN::Dimension(get_type()) < 3) return MB_SUCCESS;
  connectivity.push_back(get_vertex(i, j, k+1));
  connectivity.push_back(get_vertex(i+1, j, k+1));
  connectivity.push_back(get_vertex(i+1, j+1, k+1));
  connectivity.push_back(get_vertex(i, j+1, k+1));
  return MB_SUCCESS;
}

#endif
