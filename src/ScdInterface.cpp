#include "moab/ScdInterface.hpp"
#include "moab/Core.hpp"
#include "SequenceManager.hpp"
#include "EntitySequence.hpp"
#include "StructuredElementSeq.hpp"
#include "VertexSequence.hpp"
#include "ScdVertexData.hpp"
#ifdef USE_MPI
#  include "moab/ParallelComm.hpp"
#endif
#include "assert.h"
#include <iostream>
#include <functional>
#include "moab/TupleList.hpp"
#include "moab/gs.hpp"

#define ERRORR(rval, str) {if (MB_SUCCESS != rval)          \
      {std::cerr << str; return rval; }}

        
namespace moab 
{
    
ScdInterface::ScdInterface(Core *impl, bool boxes) 
        : mbImpl(impl), 
          searchedBoxes(false),
          boxPeriodicTag(0),
          boxDimsTag(0),
          globalBoxDimsTag(0),
          partMethodTag(0),
          boxSetTag(0)
{
  if (boxes) find_boxes(scdBoxes);
}

  // Destructor
ScdInterface::~ScdInterface() 
{
  std::vector<ScdBox*> tmp_boxes;
  tmp_boxes.swap(scdBoxes);

  for (std::vector<ScdBox*>::iterator rit = tmp_boxes.begin(); rit != tmp_boxes.end(); rit++)
    delete *rit;

  if (box_set_tag(false)) 
    mbImpl->tag_delete(box_set_tag());

}

Interface *ScdInterface::impl() const
{
  return mbImpl;
}

ErrorCode ScdInterface::find_boxes(std::vector<ScdBox*> &scd_boxes) 
{
  Range tmp_boxes;
  ErrorCode rval = find_boxes(tmp_boxes);
  if (MB_SUCCESS != rval) return rval;

  for (Range::iterator rit = tmp_boxes.begin(); rit != tmp_boxes.end(); rit++) {
    ScdBox *tmp_box = get_scd_box(*rit);
    if (tmp_box) scd_boxes.push_back(tmp_box);
    else rval = MB_FAILURE;
  }

  return rval;
}

ErrorCode ScdInterface::find_boxes(Range &scd_boxes) 
{
  ErrorCode rval = MB_SUCCESS;
  box_dims_tag();
  Range boxes;
  if (!searchedBoxes) {
    rval = mbImpl->get_entities_by_type_and_tag(0, MBENTITYSET, &boxDimsTag, NULL, 1, 
                                                boxes, Interface::UNION);
    searchedBoxes = true;
    if (!boxes.empty()) {
      scdBoxes.resize(boxes.size());
      rval = mbImpl->tag_get_data(boxSetTag, boxes, &scdBoxes[0]);
      ScdBox *dum = NULL;
      std::remove_if(scdBoxes.begin(), scdBoxes.end(), std::bind2nd(std::equal_to<ScdBox*>(), dum) ) ;
    }
  }

  for (std::vector<ScdBox*>::iterator vit = scdBoxes.begin(); vit != scdBoxes.end(); vit++)
    scd_boxes.insert((*vit)->box_set());

  return rval;
}

ScdBox *ScdInterface::get_scd_box(EntityHandle eh) 
{
  ScdBox *scd_box = NULL;
  if (!box_set_tag(false)) return scd_box;

  mbImpl->tag_get_data(box_set_tag(), &eh, 1, &scd_box);
  return scd_box;
}

ErrorCode ScdInterface::construct_box(HomCoord low, HomCoord high, double *coords, unsigned int num_coords,
                                      ScdBox *& new_box, bool is_periodic_i, bool is_periodic_j) 
{
    // create a rectangular structured mesh block
  ErrorCode rval;

  HomCoord tmp_size = high - low + HomCoord(1, 1, 1, 0);
  if ((tmp_size[1] && num_coords && (int)num_coords <= tmp_size[0]) ||
      (tmp_size[2] && num_coords && (int)num_coords <= tmp_size[0]*tmp_size[1]))
    return MB_FAILURE;

  rval = create_scd_sequence(low, high, MBVERTEX, 0, new_box);
  ERRORR(rval, "Trouble creating scd vertex sequence.");

  if (num_coords && coords) {
      // set the vertex coordinates
    double *xc, *yc, *zc;
    rval = new_box->get_coordinate_arrays(xc, yc, zc);
    ERRORR(rval, "Couldn't get vertex coordinate arrays.");

    unsigned int i = 0;
    for (int kl = low[2]; kl <= high[2]; kl++) {
      for (int jl = low[1]; jl <= high[1]; jl++) {
        for (int il = low[0]; il <= high[0]; il++) {
          xc[i] = coords[3*i];
          if (new_box->box_size()[1])
            yc[i] = coords[3*i+1];
          if (new_box->box_size()[2])
            zc[i] = coords[3*i+2];
          i++;
        }
      }
    }
  }

    // create element sequence
  SequenceManager *seq_mgr = mbImpl->sequence_manager();

  EntitySequence *tmp_seq;
  EntityHandle start_ent;

    // construct the sequence
  EntityType this_tp = MBHEX;
  if (1 >= tmp_size[2]) this_tp = MBQUAD;
  if (1 >= tmp_size[2] && 1 >= tmp_size[1]) this_tp = MBEDGE;
  rval = seq_mgr->create_scd_sequence(low, high, this_tp, 0, start_ent, tmp_seq, is_periodic_i, is_periodic_j);
  ERRORR(rval, "Trouble creating scd element sequence.");

  new_box->elem_seq(tmp_seq);
  new_box->start_element(start_ent);

    // add vertex seq to element seq, forward orientation, unity transform
  rval = new_box->add_vbox(new_box,
                             // p1: imin,jmin
                           low, low, 
                             // p2: imax,jmin
                           low + HomCoord(1, 0, 0),
                           low + HomCoord(1, 0, 0),
                             // p3: imin,jmax
                           low + HomCoord(0, 1, 0),
                           low + HomCoord(0, 1, 0));
  ERRORR(rval, "Error constructing structured element sequence.");

    // add the new hexes to the scd box set; vertices were added in call to create_scd_sequence
  Range tmp_range(new_box->start_element(), new_box->start_element() + new_box->num_elements() - 1);
  rval = mbImpl->add_entities(new_box->box_set(), tmp_range);
  ERRORR(rval, "Couldn't add new hexes to box set.");

  return MB_SUCCESS;
}


ErrorCode ScdInterface::create_scd_sequence(HomCoord low, HomCoord high, EntityType tp,
                                            int starting_id, ScdBox *&new_box,
                                            bool is_periodic_i, bool is_periodic_j)
{
  HomCoord tmp_size = high - low + HomCoord(1, 1, 1, 0);
  if ((tp == MBHEX && 1 >= tmp_size[2]) ||
      (tp == MBQUAD && 1 >= tmp_size[1]))
    return MB_TYPE_OUT_OF_RANGE;

  SequenceManager *seq_mgr = mbImpl->sequence_manager();

  EntitySequence *tmp_seq;
  EntityHandle start_ent, scd_set;

    // construct the sequence
  ErrorCode rval = seq_mgr->create_scd_sequence(low, high, tp, starting_id, start_ent, tmp_seq,
                                                is_periodic_i, is_periodic_j);
  if (MB_SUCCESS != rval) return rval;

    // create the set for this rectangle
  rval = create_box_set(low, high, scd_set);
  if (MB_SUCCESS != rval) return rval;

    // make the ScdBox
  new_box = new ScdBox(this, scd_set, tmp_seq);
  if (!new_box) return MB_FAILURE;

    // set the start vertex/element
  Range new_range;
  if (MBVERTEX == tp) {
    new_range.insert(start_ent, start_ent+new_box->num_vertices()-1);
  }
  else {
    new_range.insert(start_ent, start_ent+new_box->num_elements()-1);
  }

    // put the entities in the box set
  rval = mbImpl->add_entities(scd_set, new_range);
  if (MB_SUCCESS != rval) return rval;

    // tag the set with the box
  rval = mbImpl->tag_set_data(box_set_tag(), &scd_set, 1, &new_box);
  if (MB_SUCCESS != rval) return rval;

  return MB_SUCCESS;
}

ErrorCode ScdInterface::create_box_set(const HomCoord low, const HomCoord high,
                                       EntityHandle &scd_set, bool is_periodic_i, bool is_periodic_j) 
{
    // create the set and put the entities in it
  ErrorCode rval = mbImpl->create_meshset(MESHSET_SET, scd_set);
  if (MB_SUCCESS != rval) return rval;

    // tag the set with parameter extents
  int boxdims[6];
  for (int i = 0; i < 3; i++) boxdims[i] = low[i];
  for (int i = 0; i < 3; i++) boxdims[3+i] = high[i];
  rval = mbImpl->tag_set_data(box_dims_tag(), &scd_set, 1, boxdims);
  if (MB_SUCCESS != rval) return rval;

  int is_periodic[2] = {is_periodic_i, is_periodic_j};
  rval = mbImpl->tag_set_data(box_periodic_tag(), &scd_set, 1, is_periodic);
  if (MB_SUCCESS != rval) return rval;

  return rval;
}

Tag ScdInterface::box_periodic_tag(bool create_if_missing) 
{
  if (boxPeriodicTag || !create_if_missing) return boxPeriodicTag;

  ErrorCode rval = mbImpl->tag_get_handle("BOX_PERIODIC", 2, MB_TYPE_INTEGER, 
                                          boxPeriodicTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return boxPeriodicTag;
}

Tag ScdInterface::box_dims_tag(bool create_if_missing) 
{
  if (boxDimsTag || !create_if_missing) return boxDimsTag;

  ErrorCode rval = mbImpl->tag_get_handle("BOX_DIMS", 6, MB_TYPE_INTEGER, 
                                          boxDimsTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return boxDimsTag;
}

Tag ScdInterface::global_box_dims_tag(bool create_if_missing) 
{
  if (globalBoxDimsTag || !create_if_missing) return globalBoxDimsTag;

  ErrorCode rval = mbImpl->tag_get_handle("GLOBAL_BOX_DIMS", 6, MB_TYPE_INTEGER, 
                                          globalBoxDimsTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return globalBoxDimsTag;
}

Tag ScdInterface::part_method_tag(bool create_if_missing) 
{
  if (partMethodTag || !create_if_missing) return partMethodTag;

  ErrorCode rval = mbImpl->tag_get_handle("PARTITION_METHOD", 1, MB_TYPE_INTEGER, 
                                          partMethodTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return partMethodTag;
}

Tag ScdInterface::box_set_tag(bool create_if_missing) 
{
  if (boxSetTag || !create_if_missing) return boxSetTag;

  ErrorCode rval = mbImpl->tag_get_handle("__BOX_SET", sizeof(ScdBox*), MB_TYPE_OPAQUE,
                                          boxSetTag, MB_TAG_SPARSE|MB_TAG_CREAT);
  if (MB_SUCCESS != rval) return 0;
  return boxSetTag;
}

  //! Remove the box from the list on ScdInterface
ErrorCode ScdInterface::remove_box(ScdBox *box) 
{
  std::vector<ScdBox*>::iterator vit = std::find(scdBoxes.begin(), scdBoxes.end(), box);
  if (vit != scdBoxes.end()) {
    scdBoxes.erase(vit);
    return MB_SUCCESS;
  }
  else return MB_FAILURE;
}

  //! Add the box to the list on ScdInterface
ErrorCode ScdInterface::add_box(ScdBox *box) 
{
  scdBoxes.push_back(box);
  return MB_SUCCESS;
}

ErrorCode ScdInterface::get_boxes(std::vector<ScdBox*> &boxes) 
{
  std::copy(scdBoxes.begin(), scdBoxes.end(), std::back_inserter(boxes));
  return MB_SUCCESS;
}

ScdBox::ScdBox(ScdInterface *sc_impl, EntityHandle box_set,
               EntitySequence *seq1, EntitySequence *seq2) 
        : scImpl(sc_impl), boxSet(box_set), vertDat(NULL), elemSeq(NULL), startVertex(0), startElem(0),
          partMethod(-1)
{
  for (int i = 0; i < 6; i++) {
    boxDims[i] = 0;
    globalBoxDims[i] = 0;
  }  
  for (int i = 0; i < 2; i++) isPeriodic[i] = false;
  VertexSequence *vseq = dynamic_cast<VertexSequence *>(seq1);
  if (vseq) vertDat = dynamic_cast<ScdVertexData*>(vseq->data());
  if (vertDat) {
      // retrieve the parametric space
    for (int i = 0; i < 3; i++) {
      boxDims[i] = vertDat->min_params()[i];
      boxDims[3+i] = vertDat->max_params()[i];
    }
    startVertex = vertDat->start_handle();
  }
  else if (sc_impl->boxDimsTag) {
      // look for parametric space info on set
    ErrorCode rval = sc_impl->mbImpl->tag_get_data(sc_impl->boxDimsTag, &box_set, 1, boxDims);
    if (MB_SUCCESS == rval) {
      Range verts;
      sc_impl->mbImpl->get_entities_by_dimension(box_set, 0, verts);
      if (!verts.empty()) startVertex = *verts.begin();
    }
  }

  elemSeq = dynamic_cast<StructuredElementSeq *>(seq2);
  if (!elemSeq)
    elemSeq = dynamic_cast<StructuredElementSeq *>(seq1);

  if (elemSeq) {
    if (vertDat) {
        // check the parametric space to make sure it's consistent
      assert(elemSeq->sdata()->min_params() == HomCoord(boxDims, 3) && 
             (elemSeq->sdata()->max_params() + HomCoord(1, 1, 1)) == HomCoord(boxDims, 3));

    } 
    else {
        // get the parametric space from the element sequence
      for (int i = 0; i < 3; i++) {
        boxDims[i] = elemSeq->sdata()->min_params()[i];
        boxDims[3+i] = elemSeq->sdata()->max_params()[i];
      }
    }

    startElem = elemSeq->start_handle();

    elemSeq->is_periodic(isPeriodic);
  }
  else {
    Range elems;
    sc_impl->mbImpl->get_entities_by_dimension(box_set, (boxDims[2] == boxDims[5] ? 2 : 3), elems);
    if (!elems.empty()) startElem = *elems.begin();
    int dum[2];
    ErrorCode rval = sc_impl->mbImpl->tag_get_data(sc_impl->boxPeriodicTag, &box_set, 1, dum);
    if (MB_SUCCESS == rval) {
      isPeriodic[0] = (dum[0] ? true : false);
      isPeriodic[1] = (dum[1] ? true : false);
    }
  }

  assert(vertDat || elemSeq || 
         boxDims[0] != boxDims[3]|| boxDims[1] != boxDims[4]|| boxDims[2] != boxDims[5]);

  boxSize = HomCoord(boxDims+3, 3) - HomCoord(boxDims, 3) + HomCoord(1, 1, 1);
  boxSizeIJ = (boxSize[1] ? boxSize[1] : 1) * boxSize[0];
  boxSizeIM1 = boxSize[0]-(isPeriodic[0] ? 0 : 1);
  boxSizeIJM1 = (boxSize[1] ? (boxSize[1]-(isPeriodic[1] ? 0 : 1)) : 1) * boxSizeIM1;

  scImpl->add_box(this);
}

ScdBox::~ScdBox() 
{
    // reset the tag on the set
  ScdBox *tmp_ptr = NULL;
  if (boxSet) scImpl->mbImpl->tag_set_data(scImpl->box_set_tag(), &boxSet, 1, &tmp_ptr);
  scImpl->remove_box(this);
}

EntityHandle ScdBox::get_vertex_from_seq(int i, int j, int k) const
{
  assert(elemSeq);
  return elemSeq->get_vertex(i, j, k);
}

int ScdBox::box_dimension() const
{
  return (startElem ? scImpl->mbImpl->dimension_from_handle(startElem) : -1);
}

ErrorCode ScdBox::add_vbox(ScdBox *vbox,
                           HomCoord from1, HomCoord to1, 
                           HomCoord from2, HomCoord to2,
                           HomCoord from3, HomCoord to3,
                           bool bb_input,
                           const HomCoord &bb_min,
                           const HomCoord &bb_max)
{
  if (!vbox->vertDat) return MB_FAILURE;
  ScdVertexData *dum_data = dynamic_cast<ScdVertexData*>(vbox->vertDat);
  ErrorCode rval = elemSeq->sdata()->add_vsequence(dum_data, from1, to1, from2, to2, from3, to3,
                                                   bb_input, bb_min, bb_max);
  return rval;
}

bool ScdBox::boundary_complete() const
{
  return elemSeq->boundary_complete();
}

ErrorCode ScdBox::get_coordinate_arrays(double *&xc, double *&yc, double *&zc) 
{
  if (!vertDat) return MB_FAILURE;

  xc = reinterpret_cast<double*>(vertDat->get_sequence_data(0));
  yc = reinterpret_cast<double*>(vertDat->get_sequence_data(1));
  zc = reinterpret_cast<double*>(vertDat->get_sequence_data(2));
  return MB_SUCCESS;
}

ErrorCode ScdBox::get_coordinate_arrays(const double *&xc, const double *&yc, const double *&zc) const
{
  if (!vertDat) return MB_FAILURE;
  xc = reinterpret_cast<const double*>(vertDat->get_sequence_data(0));
  yc = reinterpret_cast<const double*>(vertDat->get_sequence_data(1));
  zc = reinterpret_cast<const double*>(vertDat->get_sequence_data(2));
  return MB_SUCCESS;
}

ErrorCode ScdBox::vert_dat(ScdVertexData *vert_dat)
{
  vertDat = vert_dat;
  return MB_SUCCESS;
}

ErrorCode ScdBox::elem_seq(EntitySequence *elem_seq)
{
  elemSeq = dynamic_cast<StructuredElementSeq*>(elem_seq);
  if (elemSeq) elemSeq->is_periodic(isPeriodic);

  if (isPeriodic[0])
    boxSizeIM1 = boxSize[0]-(isPeriodic[0] ? 0 : 1);
  if (isPeriodic[0] || isPeriodic[1])
    boxSizeIJM1 = (boxSize[1] ? (boxSize[1]-(isPeriodic[1] ? 0 : 1)) : 1) * boxSizeIM1;

  return (elemSeq ? MB_SUCCESS : MB_FAILURE);
}  

ErrorCode ScdBox::get_params(EntityHandle ent, HomCoord &ijkd) const 
{
    // check first whether this is an intermediate entity, so we know what to do
  int dimension = box_dimension();
  int this_dim = scImpl->impl()->dimension_from_handle(ent);

  if ((0 == this_dim && !vertDat) ||
      (this_dim && this_dim == dimension)) {
    assert(elemSeq);
    return elemSeq->get_params(ent, ijkd[0], ijkd[1], ijkd[2]);
  }

  else if (!this_dim && vertDat)
    return vertDat->get_params(ent, ijkd[0], ijkd[1], ijkd[2]);

  else return MB_NOT_IMPLEMENTED;
}

  //! Get the entity of specified dimension adjacent to parametric element
  /**
   * \param dim Dimension of adjacent entity being requested
   * \param i Parametric coordinates of cell being evaluated
   * \param j Parametric coordinates of cell being evaluated
   * \param k Parametric coordinates of cell being evaluated
   * \param dir Direction (0, 1, or 2), for getting adjacent edges (2d, 3d) or faces (3d) 
   * \param ent EntityHandle of adjacent entity
   * \param create_if_missing If true, creates the entity if it doesn't already exist
   */
ErrorCode ScdBox::get_adj_edge_or_face(int dim, int i, int j, int k, int dir, EntityHandle &ent,
                                       bool create_if_missing) const 
{
    // describe connectivity of sub-element in static array
    // subconnect[dim-1][dir][numv][ijk] where dimensions are:
    // [dim-1]: dim=1 or 2, so this is 0 or 1
    // [dir]: one of 0..2, for ijk directions in a hex
    // [numv]: number of vertices describing sub entity = 2*dim <= 4
    // [ijk]: 3 values for i, j, k
  int subconnect[2][3][4][3] = {
      {{{0, 0, 0}, {1, 0, 0}, {-1, -1, -1}, {-1, -1, -1}}, // i edge
       {{0, 0, 0}, {0, 1, 0}, {-1, -1, -1}, {-1, -1, -1}}, // j edge
       {{0, 0, 0}, {0, 0, 1}, {-1, -1, -1}, {-1, -1, -1}}}, // k edge

      {{{0, 0, 0}, {0, 1, 0}, {0, 1, 1}, {0, 0, 1}}, // i face
       {{0, 0, 0}, {1, 0, 0}, {1, 0, 1}, {0, 0, 1}}, // j face
       {{0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0}}}}; // k face

    // check proper input dimensions and lower bound
  if (dim < 1 || dim > 2 || i < boxDims[0] || j < boxDims[1] || k < boxDims[2]) 
    return MB_FAILURE;

    // now check upper bound; parameters must be < upper corner, since edges/faces
    // follow element parameterization, not vertex parameterization
  else if ((boxDims[3] != boxDims[0] && i > (is_periodic_i() ? boxDims[3] : boxDims[3]-1)) ||
           (boxDims[4] != boxDims[1] && j > (is_periodic_j() ? boxDims[4] : boxDims[4]-1)) ||
           (boxDims[5] != boxDims[2] && k >= boxDims[5])) return MB_FAILURE;

        // get the vertices making up this entity
  EntityHandle verts[4];
  for (int ind = 0; ind < 2*dim; ind++) {
    verts[ind] = get_vertex(i+subconnect[dim-1][dir][ind][0],
                            j+subconnect[dim-1][dir][ind][1],
                            k+subconnect[dim-1][dir][ind][2]);
    if (!verts[ind]) return MB_FAILURE;
  }
  
  Range ents;
  ErrorCode rval = scImpl->impl()->get_adjacencies(verts, 2*dim, dim, false, ents);
  if (MB_SUCCESS != rval) return rval;

  if (ents.size() > 1) return MB_FAILURE;
  
  else if (ents.size() == 1) {
    ent = *ents.begin();
  }
  else if (create_if_missing)
    rval = scImpl->impl()->create_element((1 == dim ? MBEDGE : MBQUAD), verts, 2*dim, ent);
    
  return rval;
}
    
ErrorCode ScdInterface::tag_shared_vertices(ParallelComm *pcomm, EntityHandle seth) 
{
#ifdef USE_MPI
    // first, look for box data on the set
  ScdBox *box = get_scd_box(seth);
  Range tmp_range;
  ErrorCode rval;
  if (!box) {
      // look for contained boxes
    rval = mbImpl->get_entities_by_type(seth, MBENTITYSET, tmp_range);
    if (MB_SUCCESS != rval) return rval;
    for (Range::iterator rit = tmp_range.begin(); rit != tmp_range.end(); rit++) {
      box = get_scd_box(*rit);
      if (box) break;
    }
  }
  
  if (!box) return MB_FAILURE;

    // check the # ents in the box against the num in the set, to make sure it's only 1 box;
    // reuse tmp_range
  tmp_range.clear();
  rval = mbImpl->get_entities_by_dimension(seth, box->box_dimension(), tmp_range);
  if (MB_SUCCESS != rval) return rval;
  if (box->num_elements() != (int)tmp_range.size()) return MB_FAILURE;
    
  const int *gdims = box->global_box_dims();
  if ((gdims[0] == gdims[3] && gdims[1] == gdims[4] && gdims[2] == gdims[5]) ||
      -1 == box->part_method()) return MB_FAILURE;

    // ok, we have a partitioned box; get the vertices shared with other processors
  std::vector<int> procs, offsets, shared_indices;
  rval = get_shared_vertices(pcomm, box, procs, offsets, shared_indices);
  if (MB_SUCCESS != rval) return rval;

    // post receives for start handles once we know how many to look for
  std::vector<MPI_Request> recv_reqs(procs.size(), MPI_REQUEST_NULL), 
      send_reqs(procs.size(), MPI_REQUEST_NULL);
  std::vector<EntityHandle> rhandles(4*procs.size()), shandles(4);
  for (unsigned int i = 0; i < procs.size(); i++) {
    int success = MPI_Irecv(&rhandles[4*i], 4*sizeof(EntityHandle),
                            MPI_UNSIGNED_CHAR, procs[i],
                            1, pcomm->proc_config().proc_comm(), 
                            &recv_reqs[i]);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }

    // send our own start handles
  shandles[0] = box->start_vertex();
  shandles[1] = 0;
  if (box->box_dimension() == 2) {
    shandles[2] = box->start_element();
    shandles[3] = 0;
  }
  else {
    shandles[2] = 0;
    shandles[3] = box->start_element();
  }
  for (unsigned int i = 0; i < procs.size(); i++) {
    int success = MPI_Isend(&shandles[0], 4*sizeof(EntityHandle), MPI_UNSIGNED_CHAR, procs[i], 
                            1, pcomm->proc_config().proc_comm(), &send_reqs[i]);
    if (success != MPI_SUCCESS) return MB_FAILURE;
  }
  
    // receive start handles and save info to a tuple list
  int incoming = procs.size();
  int p, j, k;
  MPI_Status status;
  TupleList shared_data;
  shared_data.initialize(1, 0, 2, 0, 
                         shared_indices.size()/2);
  shared_data.enableWriteAccess();

  j = 0; k = 0;
  while (incoming) {
    int success = MPI_Waitany(procs.size(), &recv_reqs[0], &p, &status);
    if (MPI_SUCCESS != success) return MB_FAILURE;
    unsigned int num_indices = (offsets[p+1]-offsets[p])/2;
    int *lh = &shared_indices[offsets[p]], *rh = lh + num_indices;
    for (unsigned int i = 0; i < num_indices; i++) {
      shared_data.vi_wr[j++] = procs[p];
      shared_data.vul_wr[k++] = shandles[0] + lh[i];
      shared_data.vul_wr[k++] = rhandles[4*p] + rh[i];
      shared_data.inc_n();
    }
    incoming--;
  }

    // sort by local handle
  TupleList::buffer sort_buffer;
  sort_buffer.buffer_init(shared_indices.size()/2);
  shared_data.sort(1, &sort_buffer);
  sort_buffer.reset();
  
    // process into sharing data
  std::map<std::vector<int>, std::vector<EntityHandle> > proc_nvecs;
  Range dum;
  rval = pcomm->tag_shared_verts(shared_data, proc_nvecs, dum, 0);
  if (MB_SUCCESS != rval) return rval;
  
    // create interface sets
  rval = pcomm->create_interface_sets(proc_nvecs, -1, -1);
  if (MB_SUCCESS != rval) return rval;

    // make sure buffers are allocated for communicating procs
  for (std::vector<int>::iterator pit = procs.begin(); pit != procs.end(); pit++)
    pcomm->get_buffers(*pit);


  shared_data.reset();  
#ifndef NDEBUG
  rval = pcomm->check_all_shared_handles();
  if (MB_SUCCESS != rval) return rval;
#endif
  
  return MB_SUCCESS;
  
#else
  return MB_FAILURE;
#endif
}

ErrorCode ScdInterface::get_neighbor_alljkbal(int nr, int np,
                                              const int *gdims, const int *ldims,
                                              bool periodic_i, bool periodic_j,
                                              int *dijk, 
                                              int &pto, int *bdy_ind, int *rdims, int *facedims) 
{
#ifdef USE_MPI
  if (dijk[0] != 0) {
    pto = -1;
    return MB_SUCCESS;
  }
  
  pto = -1;
  bdy_ind[0] = bdy_ind[1] = -1;
  int pj, pk; // pj, pk: # procs in j, k directions
  ErrorCode rval = compute_partition_alljkbal(np, nr, gdims, rdims, &pj);
  if (MB_SUCCESS != rval) return rval;
  pk = np / pj;
  assert(pj * pk == np);
  pto = -1;
  if ((1 == pk && dijk[2]) ||  // 1d in j means no neighbors with dk != 0
      (!(nr%pk) && -1 == dijk[2]) || // at -k bdy
      (nr%pk == pk-1 && 1 == dijk[2]) || // at +k bdy
      (nr < pk && -1 == dijk[1] && !periodic_j) ||  // down and not periodic
      (nr >= np-pk && 1 == dijk[1] && !periodic_j))  // up and not periodic
    return MB_SUCCESS;
    
  pto = nr;
  std::copy(ldims, ldims+6, facedims);
  
  if (0 != dijk[1]) {
    pto = (pto + dijk[1]*pk) % np;
    assert (pto >= 0 && pto < np);
    int dj = (gdims[4] - gdims[1]) / pj, extra = (gdims[4] - gdims[1]) % pj;
    if (-1 == dijk[1]) {
      facedims[4] = facedims[1];
      rdims[4] = rdims[1];
      rdims[1] -= dj;
      if (pto < extra) rdims[1]--;
    }
    else {
      facedims[1] = facedims[4];
      rdims[1] = rdims[4];
      rdims[4] += dj;
      if (pto < extra) rdims[4]++;
    }
  }
  if (0 != dijk[2]) {
    pto = (pto + dijk[2]) % np;
    assert (pto >= 0 && pto < np);
    facedims[2] = facedims[5] = (-1 == dijk[2] ? facedims[2] : facedims[5]);
    int dk = (gdims[5] - gdims[2]) / pk;
    if (-1 == dijk[2]) {
      facedims[5] = facedims[2];
      rdims[5] = rdims[2];
      rdims[2] -= dk;
    }
    else {
      facedims[2] = facedims[5];
      rdims[2] = rdims[5];
      rdims[5] += dk;
    }
  }

  if (dijk[1] == -1 && periodic_j && ldims[1] == gdims[1]) bdy_ind[1] = 1;

  assert(-1 == pto ||
         (rdims[0] >= gdims[0] && rdims[3] <= gdims[3] && 
          rdims[1] >= gdims[1] && rdims[4] <= gdims[4] && 
          rdims[2] >= gdims[2] && rdims[5] <= gdims[5] &&
          facedims[0] >= rdims[0] && facedims[3] <= rdims[3] &&
          facedims[1] >= rdims[1] && facedims[4] <= rdims[4] &&
          facedims[2] >= rdims[2] && facedims[5] <= rdims[5] &&
          facedims[0] >= ldims[0] && facedims[3] <= ldims[3] &&
          facedims[1] >= ldims[1] && facedims[4] <= ldims[4] &&
          facedims[2] >= ldims[2] && facedims[5] <= ldims[5]));
  
  return MB_SUCCESS;
#else
  return MB_FAILURE;
#endif  
}

ErrorCode ScdInterface::get_neighbor_sqij(int nr, int np,
                                          const int *gdims, const int *ldims,
                                          bool periodic_i, bool periodic_j,
                                          int *dijk, 
                                          int &pto, int *bdy_ind, int *rdims, int *facedims) 
{
#ifdef USE_MPI
  if (dijk[2] != 0) {
    pto = -1;
    return MB_SUCCESS;
  }
  
  pto = -1;
  bdy_ind[0] = bdy_ind[1] = -1;
  int pi, pj; // pi, pj: # procs in i, j directions
              // guess pi
  ErrorCode rval = compute_partition_sqij(np, nr, gdims, rdims, &pi);
  if (MB_SUCCESS != rval) return rval;
  pj = np / pi;
  assert(pi * pj == np);
  pto = -1;
  if ((!(nr%pi) && -1 == dijk[0] && !periodic_i) ||  // left and not periodic
      ((nr%pi) == pi-1 && 1 == dijk[0] && !periodic_i) ||  // right and not periodic
      (!(nr/pi) && -1 == dijk[1] && !periodic_j) || // bottom and not periodic
      (nr/pi == pj-1 && 1 == dijk[1] && !periodic_j))  // top and not periodic
    return MB_SUCCESS;
  
  std::copy(ldims, ldims+6, facedims);
  pto = nr;
  int dj = (gdims[4] - gdims[1]) / pj, jextra = (gdims[4] - gdims[1]) % dj,
      di = (gdims[3] - gdims[0]) / pi, iextra = (gdims[3] - gdims[0]) % di;
  
  if (0 != dijk[0]) {
    pto = (pto + dijk[0]) % np;
    assert (pto >= 0 && pto < np);
    if (-1 == dijk[0]) {
      facedims[3] = facedims[0];
      rdims[3] = rdims[0];
      rdims[0] -= di;
      if (pto%pi < iextra) rdims[0]--;
    }
    else {
      facedims[0] = facedims[3];
      rdims[0] = rdims[3];
      rdims[3] += di;
      if (pto%pi < iextra) rdims[3]++;
    }
  }
  if (0 != dijk[1]) {
    pto = (pto + dijk[1]*pi) % np;
    assert (pto >= 0 && pto < np);
    if (-1 == dijk[1]) {
      facedims[4] = facedims[1];
      rdims[4] = rdims[1];
      rdims[1] -= dj;
      if (pto/pi < jextra) rdims[1]--;
    }
    else {
      facedims[1] = facedims[4];
      rdims[1] = rdims[4];
      rdims[4] += dj;
      if (pto/pi < jextra) rdims[4]++;
    }
  }

  if (dijk[0] == -1 && periodic_i && ldims[0] == gdims[0]) bdy_ind[0] = 1;
  if (dijk[1] == -1 && periodic_j && ldims[1] == gdims[1]) bdy_ind[1] = 1;

  assert(-1 == pto ||
         (rdims[0] >= gdims[0] && rdims[3] <= gdims[3] && 
          rdims[1] >= gdims[1] && rdims[4] <= gdims[4] && 
          rdims[2] >= gdims[2] && rdims[5] <= gdims[5] &&
          facedims[0] >= rdims[0] && facedims[3] <= rdims[3] &&
          facedims[1] >= rdims[1] && facedims[4] <= rdims[4] &&
          facedims[2] >= rdims[2] && facedims[5] <= rdims[5] &&
          facedims[0] >= ldims[0] && facedims[3] <= ldims[3] &&
          facedims[1] >= ldims[1] && facedims[4] <= ldims[4] &&
          facedims[2] >= ldims[2] && facedims[5] <= ldims[5]));

  return MB_SUCCESS;
#else
  return MB_FAILURE;
#endif  
}

ErrorCode ScdInterface::get_neighbor_sqjk(int nr, int np,
                                          const int *gdims, const int *ldims,
                                          bool periodic_i, bool periodic_j,
                                          int *dijk, 
                                          int &pto, int *bdy_ind, int *rdims, int *facedims) 
{
#ifdef USE_MPI
  if (dijk[0] != 0) {
    pto = -1;
    return MB_SUCCESS;
  }
  
  pto = -1;
  bdy_ind[0] = bdy_ind[1] = -1;
  int pj, pk; // pj, pk: # procs in j, k directions
  ErrorCode rval = compute_partition_sqjk(np, nr, gdims, rdims, &pj);
  if (MB_SUCCESS != rval) return rval;
  pk = np / pj;
  assert(pj * pk == np);
  pto = -1;
  if ((!(nr%pj) && -1 == dijk[1] && !periodic_j) ||  // down and not periodic
      ((nr%pj) ==  pj-1 && 1 == dijk[1] && !periodic_j) ||  // up and not periodic
      (!(nr/pj) && -1 == dijk[2]) || // k- bdy 
      ((nr/pj) == pk-1 && 1 == dijk[2])) // k+ bdy
    return MB_SUCCESS;
    
  std::copy(ldims, ldims+6, facedims);
  pto = nr;
  int dj = (gdims[4] - gdims[1]) / pj, jextra = (gdims[4] - gdims[1]) % dj,
      dk = (gdims[5] - gdims[2]) / pk, kextra = (gdims[5] - gdims[2]) % dk;
  
  if (0 != dijk[1]) {
    pto = (pto + dijk[1]) % np;
    assert (pto >= 0 && pto < np);
    if (-1 == dijk[1]) {
      facedims[4] = facedims[1];
      rdims[4] = rdims[1];
      rdims[1] -= dj;
      if (pto%pj < jextra) rdims[1]--;
    }
    else {
      facedims[1] = facedims[4];
      rdims[1] = rdims[4];
      rdims[4] += dj;
      if (pto%pj < jextra) rdims[4]++;
    }
  }
  if (0 != dijk[2]) {
    pto = (pto + dijk[2]*pj) % np;
    assert (pto >= 0 && pto < np);
    if (-1 == dijk[2]) {
      facedims[5] = facedims[2];
      rdims[5] = rdims[2];
      rdims[2] -= dk;
      if (pto/pj < kextra) rdims[2]--;
    }
    else {
      facedims[2] = facedims[5];
      rdims[2] = rdims[5];
      rdims[5] += dk;
      if (pto/pj < kextra) rdims[5]++;
    }
  }

  if (dijk[1] == -1 && periodic_j && ldims[1] == gdims[1]) bdy_ind[1] = 1;

  assert(-1 == pto ||
         (rdims[0] >= gdims[0] && rdims[3] <= gdims[3] && 
          rdims[1] >= gdims[1] && rdims[4] <= gdims[4] && 
          rdims[2] >= gdims[2] && rdims[5] <= gdims[5] &&
          facedims[0] >= rdims[0] && facedims[3] <= rdims[3] &&
          facedims[1] >= rdims[1] && facedims[4] <= rdims[4] &&
          facedims[2] >= rdims[2] && facedims[5] <= rdims[5] &&
          facedims[0] >= ldims[0] && facedims[3] <= ldims[3] &&
          facedims[1] >= ldims[1] && facedims[4] <= ldims[4] &&
          facedims[2] >= ldims[2] && facedims[5] <= ldims[5]));

  return MB_SUCCESS;
#else
  return MB_FAILURE;
#endif  
}

ErrorCode ScdInterface::get_neighbor_alljorkori(int nr, int np,
                                                const int *gdims, const int *ldims,
                                                bool periodic_i, bool periodic_j,
                                                int *dijk, 
                                                int &pto, int *bdy_ind, int *rdims, int *facedims) 
{
#ifdef USE_MPI
  int ind = -1;
  for (int i = 0; i < 3; i++)
    if (gdims[i] != ldims[i] || gdims[i+3] != ldims[i+3]) ind = i;
  if (!dijk[ind]) {
    pto = -1;
    return MB_SUCCESS;
  }

  pto = -1;
  bdy_ind[0] = bdy_ind[1] = -1;
  if ((dijk[0] && dijk[1] && dijk[2]) ||
      (fabs(dijk[0]) + fabs(dijk[1]) + fabs(dijk[2]) != 1)) return MB_SUCCESS;
  
  ErrorCode rval = MB_SUCCESS;
  std::copy(ldims, ldims+6, facedims);
  
  if (ldims[ind] == gdims[ind] &&
      ind < 2 && dijk[ind] == -1 && 
      ((periodic_i && 0 == ind) || (periodic_j && 1 == ind))) 
    bdy_ind[ind] = 1;

  if (-1 == dijk[ind] && nr) {
      // actual left neighbor
    pto = nr-1;
    facedims[ind+3] = facedims[ind];
    rval = compute_partition(ALLJORKORI, np, pto, gdims, rdims);
  }
  else if (1 == dijk[ind] && nr < np-1) {
      // actual right neighbor
    pto = nr+1;
    facedims[ind] = facedims[ind+3];
    rval = compute_partition(ALLJORKORI, np, pto, gdims, rdims);
  }
  else if (-1 == dijk[ind] && !nr && 
           ((0 == ind && periodic_i) || (1 == ind && periodic_j))) {
      // left across periodic bdy
    pto = np - 1;
    bdy_ind[ind] = 1;
    facedims[ind+3] = facedims[ind];
    rval = compute_partition(ALLJORKORI, np, pto, gdims, rdims);
  }
  else if (1 == dijk[ind] && nr == np-1 && 
           ((0 == ind && periodic_i) || (1 == ind && periodic_j))) {
      // right across periodic bdy
    pto = 0;
    bdy_ind[ind] = 1;
    facedims[ind] = facedims[ind+3];
    rval = compute_partition(ALLJORKORI, np, pto, gdims, rdims);
  }

  assert(-1 == pto ||
         (rdims[0] >= gdims[0] && rdims[3] <= gdims[3] && 
          rdims[1] >= gdims[1] && rdims[4] <= gdims[4] && 
          rdims[2] >= gdims[2] && rdims[5] <= gdims[5] &&
          facedims[0] >= rdims[0] && facedims[3] <= rdims[3] &&
          facedims[1] >= rdims[1] && facedims[4] <= rdims[4] &&
          facedims[2] >= rdims[2] && facedims[5] <= rdims[5] &&
          facedims[0] >= ldims[0] && facedims[3] <= ldims[3] &&
          facedims[1] >= ldims[1] && facedims[4] <= ldims[4] &&
          facedims[2] >= ldims[2] && facedims[5] <= ldims[5]));

  return rval;
#else
  return MB_FAILURE;
#endif  
}
  
  //! get shared vertices for alljorkori partition scheme
ErrorCode ScdInterface::get_shared_vertices(ParallelComm *pcomm, ScdBox *box, 
                                            std::vector<int> &procs,
                                            std::vector<int> &offsets, std::vector<int> &shared_indices) 
{
#ifdef USE_MPI
    // get index of partitioned dimension
  const int *ldims = box->box_dims();
  ErrorCode rval;
  int ijkrem[6], ijkface[6];

  for (int k = -1; k <= 1; k ++) {
    for (int j = -1; j <= 1; j ++) {
      for (int i = -1; i <= 1; i ++) {
        if (!i && !j && !k) continue;
        int pto, bdy_ind[2];
        rval = get_neighbor(pcomm->proc_config().proc_rank(), pcomm->proc_config().proc_size(),
                            box->part_method(),
                            box->global_box_dims(), box->box_dims(), 
                            box->is_periodic_i(), box->is_periodic_j(),
                            i, j, k, pto, bdy_ind, ijkrem, ijkface);
        if (MB_SUCCESS != rval) return rval;
        if (-1 != pto) {
          assert(std::find(procs.begin(), procs.end(), pto) == procs.end());
          procs.push_back(pto);
          offsets.push_back(shared_indices.size());
          rval = get_indices(bdy_ind, ldims, ijkrem, ijkface, shared_indices);
          if (MB_SUCCESS != rval) return rval;
        }
      }
    }
  }

  offsets.push_back(shared_indices.size());

  return MB_SUCCESS;
#else
  return MB_FAILURE;
#endif
}

} // namespace moab
