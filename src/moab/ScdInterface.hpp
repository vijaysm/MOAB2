/** \file ScdInterface.hpp
 */
#ifndef SCD_INTERFACE
#define SCD_INTERFACE

#include "moab/Interface.hpp"
#include "moab/HomXform.hpp"

#include <iostream>
#include <vector>    
#include "assert.h"

namespace moab {

class StructuredElementSeq;
class EntitySequence;
class ScdVertexData;
class EntitySequence;
class ScdBox;
class Core;
class ParallelComm;

/** \class ScdInterface ScdInterface.hpp "moab/ScdInterface.hpp"
 * \brief A structured mesh interface for MOAB-based data
 *
 * Structured mesh in MOAB is created and accessed through the ScdInterface and ScdBox classes.
 * 
 * \section Construction Construction and Representation
 * Structured mesh can be constructed in one of two ways.  First, a rectangular block of mesh, 
 * both vertices and edges/quads/hexes, can be created in one shot, using the construct_box method.
 * In this case, there are single sequences of vertices/entities.  The second method for creating
 * structured mesh is to create the structured blocks of vertices and elements separately.  In
 * this case, different blocks of elements can share blocks of vertices, and each block of
 * elements has its own independent parametric space.  The algorithms behind this representation
 * are described in T. Tautges, "MOAB-SD: Integrated structured and unstructured mesh representation",
 * Eng. w Comp, vol 20 no. 3.
 * 
 * Structured mesh is represented in MOAB down at the element sequence level, which is something
 * applications don't see.  In addition, when structured mesh is created, entity sets are also
 * created and tagged with information about the parametric space.  In particular, the BOX_DIMS
 * tag is used to indicate the lower and upper corners in parametric space (this
 * tag is integer size 6).  Structured mesh blocks are also available through ScdBox class objects
 * returned by ScdInterface.  These class objects should be treated only as references to the 
 * structured mesh blocks; that is, the structured mesh referenced by these objects is not deleted
 * when the ScdBox instance is destroyed.  Functions for destroying the actual mesh are are available 
 * on this class, though.
 *
 * Structured mesh blocks are returned in the form of ScdBox class objects.  Each ScdBox instance
 * represents a rectangular block of vertices and possibly elements (edges, quads, or hexes).  The
 * edge/quad/hex entity handles for a ScdBox are guaranteed to be contiguous, starting at a starting value
 * which is also available through the ScdBox class.  However, vertex handles may or may not be
 * contiguous, depending on the construction method.  The start vertex handle is also available from
 * the ScdBox class.
 *
 * \section Parameters Parametric Space
 *
 * Each structured box has a parametric (ijk) space, which can be queried through the ScdBox interface.
 * For non-periodic boxes, the edge/quad/hex parameter bounds are one less in each dimension than that
 * of the vertices.  Entity handles are allocated in column-major order, that is, with the i parameter 
 * varying fastest, then j, then k.
 *
 * Boxes can be periodic in i, or j, or both i and j.  If only i or j is periodic, the corresponding mesh
 * is a ring or annular cylinder; if both i and j are periodic, the corresponding mesh is an annular
 * torus.  A box cannot be periodic in all three parameters.  If i and/or j is periodic, the parameter extent
 * in the/each periodic direction is equal to that of the vertices in that direction.
 *
 * \section Adjs Adjacent Entities
 * This interface supports parametric access to intermediate-dimension entities, e.g. adjacent faces
 * and edges in a 3d mesh.  In this case, a direction parameter is added, to identify the parametric
 * direction of the entities being requested.  For example, to retrieve the faces adjacent to a hex
 * with parameters ijk, in the i parametric direction, you would use the parameters ijk0.  These 
 * intermediate entities are not stored in a structured representation, but their parametric positions
 * can be evaluated based on their adjacencies to higher-dimensional entities.  Thanks to Milad Fatenejad
 * for the thinking behind this.
 *
 * \section Evaluation Evaluation
 * The ScdBox class provides functions for evaluating the mesh based on the ijk parameter space.
 * These functions are inlined where possible, for efficiency.
*/

class ScdInterface 
{
public:
  friend class ScdBox;
  
    //! Constructor
    /** Constructor; if find_boxes is true, this will search for entity sets marked as
     * structured blocks, based on the BOX_DIMS tag.  Structured mesh blocks will be stored
     * in this interface class for future retrieval.  Structured mesh blocks created through
     * this interface will also be stored here.
     * \param impl MOAB instance
     * \param find_boxes If true, search all the entity sets, caching the structured mesh blocks
     */
  ScdInterface(Core *impl, bool find_boxes = false);
  
    // Destructor
  ~ScdInterface();

    //! Return the MOAB Interface instance *
  Interface *impl() const;
  
    //! Construct new structured mesh box, including both vertices and elements
    /** Parameter range
     * for vertex box is [low-high], for elements is [low-high).  Construct quads by passing
     * in low[2] == high[2], and edges by passing in low[1] == high[1] and low[2] == high[2].
     * The result is passed back in a ScdBox*, which is a *reference* to the box of structured mesh.
     * That is, the actual mesh is retained in MOAB when the ScdBox is destroyed.  To actually destroy
     * the mesh, call the destroy_mesh function on the ScdBox object first, before destroying it.
     * \param low Lower corner in parameter space
     * \param high Higher corner in parameter space
     * \param coords Coordinates of vertices, in column-major order; if NULL, no coords are set
     * \param num_coords Number of coordinate values; if zero, no coords are set
     * \param new_box Reference to box of structured mesh
     * \param is_periodic_i True if box is periodic in i direction
     * \param is_periodic_j True if box is periodic in j direction
     */
  ErrorCode construct_box(HomCoord low, HomCoord high, double *coords, unsigned int num_coords,
                          ScdBox *& new_box, bool is_periodic_i = false, bool is_periodic_j = false);

    //! Create a structured sequence of vertices, quads, or hexes
    /** Starting handle for the sequence is available from the returned ScdBox.  
     * If creating a structured quad or hex box, subsequent calls must be made to ScdBox::add_vbox, 
     * until all the vertices have been filled in for the box.
     * \param low Lower corner of structured box
     * \param high Higher corner of structured box
     * \param type EntityType, one of MBVERTEX, MBEDGE, MBQUAD, MBHEX
     * \param starting_id Requested start id of entities
     * \param new_box Reference to the newly created box of entities
     * \param is_periodic_i True if box is periodic in i direction
     * \param is_periodic_j True if box is periodic in j direction
     */
  ErrorCode create_scd_sequence(HomCoord low, HomCoord high, EntityType type,
                                int starting_id, ScdBox *&new_box, 
                                bool is_periodic_i = false, bool is_periodic_j = false);

    //! Return all the structured mesh blocks in this MOAB instance, as ScdBox objects
    /** Return the structured blocks in this MOAB instance.  If these were not searched for
     * at instantiation time, then the search is done now.
     * \param boxes Vector of ScdBox objects representing structured mesh blocks
     */
  ErrorCode find_boxes(std::vector<ScdBox*> &boxes);

    //! Return all the structured mesh blocks in this MOAB instance, as entity set handles
    /** Return the structured blocks in this MOAB instance.  If these were not searched for
     * at instantiation time, then the search is done now.
     * \param boxes Range of entity set objects representing structured mesh blocks
     */
  ErrorCode find_boxes(Range &boxes);

    //! Return all the structured mesh blocks known by ScdInterface (does not search)
    /** Return the structured blocks in this ScdInterface instance.  Does not search for new boxes,
     * just returns the contents of the list.
     * \param boxes Structured boxes
     */
  ErrorCode get_boxes(std::vector<ScdBox*> &boxes);

    //! Return the tag marking the lower and upper corners of boxes
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_dims_tag(bool create_if_missing = true);

    //! Return the tag marking the global lower and upper corners of boxes
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag global_box_dims_tag(bool create_if_missing = true);

    //! Return the tag marking the partitioning method used to partition the box in parallel
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag part_method_tag(bool create_if_missing = true);

    //! Return the tag marking whether box is periodic in i and j
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_periodic_tag(bool create_if_missing = true);

    //! Return the tag marking the ScdBox for a set
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_set_tag(bool create_if_missing = true);

    //! Return the ScdBox corresponding to the entity set passed in
    /** If the entity isn't a structured box set, NULL is returned.
     * \param eh Entity whose box is being queried
     */
  ScdBox *get_scd_box(EntityHandle eh);

    //! Partition method enumeration
  enum PartitionMethod {ALLJORKORI=0, ALLJKBAL, SQIJ, SQJK};
  
    //! Compute a partition of structured parameter space
    /** Compute a partition of structured parameter space
     * \param part_method Partition method; should be from PartitionMethod enum, or -1
     * \param np #procs
     * \param nr Rank of this proc
     * \param iMin/jMax Global minimum/maximum i parameter
     * \param jMin/jMax Global minimum/maximum j parameter
     * \param kMin/kMax Global minimum/maximum k parameter
     * \param iMin/jMax Local minimum/maximum i parameter
     * \param jMin/jMax Local minimum/maximum j parameter
     * \param kMin/kMax Local minimum/maximum k parameter
     */
  static ErrorCode compute_partition(int part_method,
                                     int np, int nr,
                                     const int *gijk, int *lijk);
  
  static ErrorCode get_neighbor(int pfrom, int np, int part_method,
                                const int *gdims, const int *ldims,
                                bool periodic_i, bool periodic_j,
                                int di, int dj, int dk,
                                int &pto, int *bdy_ind, int *rdims, int *facedims);
  
    //! Tag vertices with sharing data for parallel representations
    /** Given the ParallelComm object to use, tag the vertices shared with other processors
     */
  ErrorCode tag_shared_vertices(ParallelComm *pcomm, EntityHandle seth);
  
protected:
    //! Remove the box from the list on ScdInterface
  ErrorCode remove_box(ScdBox *box);
  
    //! Add the box to the list on ScdInterface
  ErrorCode add_box(ScdBox *box);
  
private:
    //! Create an entity set for a box, and tag with the parameters
    /** \param low Lower corner parameters for this box
     * \param high Upper corner parameters for this box
     * \param scd_set Entity set created
     * \param is_periodic_i True if box is periodic in i direction
     * \param is_periodic_j True if box is periodic in j direction
     */
  ErrorCode create_box_set(const HomCoord low, const HomCoord high,
                           EntityHandle &scd_set,
                           bool is_periodic_i = false, bool is_periodic_j = false);
  
    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by partitioning j, k, or i only.
     * If j is greater than #procs, partition that, else k, else i.
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_alljorkori(int np, int nr,
                                                       const int *gijk, int *lijk);
  
    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by partitioning j, and possibly k,
     * seeking square regions of jk space
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_alljkbal(int np, int nr, const int *gijk, int *lijk,
                                                     int *pjp = NULL);

    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by seeking square ij partitions
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_sqij(int np, int nr, const int *gijk, int *lijk,
                                                 int *pip = NULL);
  
    //! Compute a partition of structured parameter space
    /** Partitions the structured parametric space by seeking square jk partitions
     * For description of arguments, see ScdInterface::compute_partition.
     */
  inline static ErrorCode compute_partition_sqjk(int np, int nr, const int *gijk, int *lijk,
                                                 int *pjp = NULL);

    //! Get vertices shared with other processors
    /** Shared vertices returned as indices into each proc's handle space
     * \param box Box used to get parametric space info
     * \param procs Procs this proc shares vertices with
     * \param offsets Offsets into indices list for each proc
     * \param shared_indices local/remote indices of shared vertices
     */
  static ErrorCode get_shared_vertices(ParallelComm *pcomm, ScdBox *box, std::vector<int> &procs,
                                       std::vector<int> &offsets, std::vector<int> &shared_indices);

  static ErrorCode get_indices(const int *bdy_ind, const int *ldims, int *rdims, int *face_dims, 
                               std::vector<int> &shared_indices);
  
  static ErrorCode get_neighbor_alljorkori(int pfrom, int np,
                                           const int *gdims, const int *ldims,
                                           bool periodic_i, bool periodic_j,
                                           int *dijk, 
                                           int &pto, int *bdy_ind, int *rdims, int *facedims);
  
  static ErrorCode get_neighbor_alljkbal(int pfrom, int np,
                                         const int *gdims, const int *ldims,
                                         bool periodic_i, bool periodic_j,
                                         int *dijk, 
                                         int &pto, int *bdy_ind, int *rdims, int *facedims);
  
  static ErrorCode get_neighbor_sqij(int pfrom, int np,
                                     const int *gdims, const int *ldims,
                                     bool periodic_i, bool periodic_j,
                                     int *dijk, 
                                     int &pto, int *bdy_ind, int *rdims, int *facedims);
  
  static ErrorCode get_neighbor_sqjk(int pfrom, int np,
                                     const int *gdims, const int *ldims,
                                     bool periodic_i, bool periodic_j,
                                     int *dijk, 
                                     int &pto, int *bdy_ind, int *rdims, int *facedims);
  
  static int gtol(const int *gijk, int i, int j, int k);

  //! interface instance
  Core *mbImpl;

    //! whether we've searched the database for boxes yet
  bool searchedBoxes;
  
    //! structured mesh blocks; stored as ScdBox objects, can get sets from those
  std::vector<ScdBox*> scdBoxes;

    //! tag representing whether box is periodic in i and j
  Tag boxPeriodicTag;

    //! tag representing box lower and upper corners
  Tag boxDimsTag;

    //! tag representing global lower and upper corners
  Tag globalBoxDimsTag;

    //! tag representing partition method
  Tag partMethodTag;

    //! tag pointing from set to ScdBox
  Tag boxSetTag;
  
};

class ScdBox 
{
  friend class ScdInterface;
  
public:

    //! Destructor
  ~ScdBox();

    //! Return the ScdInterface responsible for this box
  inline ScdInterface *sc_impl() const;

    //! Add a vertex box to this box
    /* Add a vertex box to the element sequence referenced by this box.  The passed in vbox must
     * be a vertex box, with parametric extents no larger than that of this box.  This vbox is
     * oriented to this box by matching parameters from1-from3 in vbox to to1-to3 in this box.
     * If bb_input is true, only the part of the vertex sequence between bb_min and bb_max is referenced
     * \param vbox The vertex box being added to this box
     * \param from1 1st reference point on vbox
     * \param to1 1st reference point on this box
     * \param from2 2nd reference point on vbox
     * \param to2 2nd reference point on this box
     * \param from3 3rd reference point on vbox
     * \param to3 3rd reference point on this box
     * \param bb_input If true, subsequent parameters list extents of vbox referenced
     * \param bb_min Lower corner of rectangle referenced
     * \param bb_max Upper corner of rectangle referenced
     */
  ErrorCode add_vbox(ScdBox *vbox,
                     HomCoord from1, HomCoord to1, 
                     HomCoord from2, HomCoord to2,
                     HomCoord from3, HomCoord to3,
                     bool bb_input = false,
                     const HomCoord &bb_min = HomCoord::unitv[0],
                     const HomCoord &bb_max = HomCoord::unitv[0]);

    //! Return whether this box has all its vertices defined
    /** Tests whether vertex boxs added with add_vbox have completely defined the vertex parametric 
     * space for this box.
     *
     */
  bool boundary_complete() const;

    //! Return highest topological dimension of box
  inline int box_dimension() const;
  
    //! Starting vertex handle for this box
  inline EntityHandle start_vertex() const;
  
    //! Starting entity handle for this box
    /** If this is a vertex box, the start vertex handle is returned.
     */
  inline EntityHandle start_element() const;

    //! Return the number of elements in the box
    /* Number of elements is (boxSize[0]-1)(boxSize[1]-1)(boxSize[2]-1)
     */
  inline int num_elements() const;
  
    //! Return the number of vertices in the box
    /* Number of vertices is boxSize[0] * boxSize[1] * boxSize[2]
     */
  inline int num_vertices() const;
  
    //! Return the parametric coordinates for this box
    /**
     * \return IJK parameters of lower and upper corners
     */
  inline const int *box_dims() const;
  
    //! Return the parametric coordinates for the entire mesh
    /**
     * \return IJK parameters of lower and upper corners of global parametric space
     */
  inline const int *global_box_dims() const;
  
    //! Set the parametric coordinates for the entire mesh
    /**
     * \param box_dims Lower and upper corners of global parametric space
     */
  inline void set_global_box_dims(int *box_dims);
  
    //! Get the partMethod member
  inline int part_method() const;
  
    //! Set the partMethod member
  inline void part_method(int method);
  
    //! Return the lower corner parametric coordinates for this box
  inline HomCoord box_min() const;
  
    //! Return the upper corner parametric coordinates for this box
  inline HomCoord box_max() const;
  
    //! Return the parameter extents for this box
  inline HomCoord box_size() const;
  
    //! Return the parametric extents for this box
    /**
     * \param ijk IJK extents of this box
     */
  inline void box_size(int *ijk) const;
  
    //! Return the parametric extents for this box
    /**
     * \param i I extent of this box
     * \param j J extent of this box
     * \param k K extent of this box
     */
  void box_size(int &i, int &j, int &k) const;
  
    //! Get the element at the specified coordinates
    /**
     * \param ijk Parametric coordinates being evaluated
     */
  EntityHandle get_element(HomCoord ijk) const;
  
    //! Get the element at the specified coordinates
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  EntityHandle get_element(int i, int j = 0, int k = 0) const;
  
    //! Get the vertex at the specified coordinates
    /**
     * \param ijk Parametric coordinates being evaluated
     */
  EntityHandle get_vertex(HomCoord ijk) const;
  
    //! Get the vertex at the specified coordinates
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  EntityHandle get_vertex(int i, int j = 0, int k = 0) const;
  
    //! Get parametric coordinates of the specified entity
    /** This function returns MB_ENTITY_NOT_FOUND if the entity is not
     * in this ScdBox.
     * \param ent Entity being queried
     * \param i Parametric coordinates returned
     * \param j Parametric coordinates returned
     * \param k Parametric coordinates returned
     * \param dir Parametric coordinate direction returned (in case of getting adjacent
     *            edges (2d, 3d) or faces (3d); not modified otherwise
     */
  ErrorCode get_params(EntityHandle ent, int &i, int &j, int &k, int &dir) const;
  
    //! Get parametric coordinates of the specified entity, intermediate entities not allowed (no dir parameter)
    /** This function returns MB_ENTITY_NOT_FOUND if the entity is not
     * in this ScdBox, or MB_FAILURE if the entity is an intermediate-dimension entity.
     * \param ent Entity being queried
     * \param i Parametric coordinates returned
     * \param j Parametric coordinates returned
     * \param k Parametric coordinates returned
     */
  ErrorCode get_params(EntityHandle ent, int &i, int &j, int &k) const;
  
    //! Get parametric coordinates of the specified entity
    /** This function returns MB_ENTITY_NOT_FOUND if the entity is not
     * in this ScdBox.
     * \param ent Entity being queried
     * \param ijkd Parametric coordinates returned (including direction, in case of 
     *            getting adjacent edges (2d, 3d) or faces (3d))
     */
  ErrorCode get_params(EntityHandle ent, HomCoord &ijkd) const;
  
    /** \brief Get the adjacent edge or face at a parametric location
     * This function gets the left (i=0), front (j=0), or bottom (k=0) edge or face for a parametric element.
     * Left, front, or bottom is indicated by dir = 0, 1, or 2, resp.  All edges and faces in a structured
     * mesh block can be accessed using these parameters.
     * \param dim Dimension of adjacent entity being requested
     * \param i Parametric coordinates of cell being evaluated
     * \param j Parametric coordinates of cell being evaluated
     * \param k Parametric coordinates of cell being evaluated
     * \param dir Direction (0, 1, or 2), for getting adjacent edges (2d, 3d) or faces (3d) 
     * \param ent Entity returned from this function
     * \param create_if_missing If true, creates the entity if it doesn't already exist
     */
  ErrorCode get_adj_edge_or_face(int dim, int i, int j, int k, int dir, EntityHandle &ent,
                                 bool create_if_missing = true) const;

    //! Return whether the box contains the parameters passed in
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  bool contains(int i, int j, int k) const;

    //! Return whether the box contains the parameters passed in
    /**
     * \param i Parametric coordinates being evaluated
     * \param j Parametric coordinates being evaluated
     * \param k Parametric coordinates being evaluated
     */
  bool contains(const HomCoord ijk) const;

    //! Set/Get the entity set representing the box
  void box_set(EntityHandle this_set);
  EntityHandle box_set();

    //! Get coordinate arrays for vertex coordinates for a structured block
    /** Returns error if there isn't a single vertex sequence associated with this structured block
     * \param xc X coordinate array pointer returned
     * \param yc Y coordinate array pointer returned
     * \param zc Z coordinate array pointer returned
     */
  ErrorCode get_coordinate_arrays(double *&xc, double *&yc, double *&zc);
  
    //! Get read-only coordinate arrays for vertex coordinates for a structured block
    /** Returns error if there isn't a single vertex sequence associated with this structured block
     * \param xc X coordinate array pointer returned
     * \param yc Y coordinate array pointer returned
     * \param zc Z coordinate array pointer returned
     */
  ErrorCode get_coordinate_arrays(const double *&xc, const double *&yc, const double *&zc) const;

    //! Return whether box is periodic in i
    /** Return whether box is periodic in i
     * \return True if box is periodic in i direction
     */
  bool is_periodic_i() const;
  
    //! Return whether box is periodic in j
    /** Return whether box is periodic in j
     * \return True if box is periodic in j direction
     */
  bool is_periodic_j() const;
  
    //! Return whether box is periodic in i and j
    /** Return whether box is periodic in i and j
     * \param is_periodic_ij Non-zero if periodic in i [0] or j [1]
     */
  void is_periodic(bool is_periodic_ij[2]) const;
  
private:
    //! Constructor
    /** Create a structured box instance; this constructor is private because it should only be called
     * from ScdInterface, a friend class.  This constructor takes two sequences, one of which can be
     * NULL.  If both sequences come in non-NULL, the first should be a VertexSequence* corresponding to
     * a structured vertex sequence and the second should be a StructuredElementSeq*.  If the 2nd is NULL,
     * the first can be either of those types.  The other members of this class are taken from the sequences
     * (e.g. parametric space) or the box set argument.  Tags on the box set should be set from the caller.
     * \param sc_impl A ScdInterface instance
     * \param box_set Entity set representing this rectangle
     * \param seq1 An EntitySequence (see ScdBox description)
     * \param seq2 An EntitySequence (see ScdBox description), or NULL
     */
  ScdBox(ScdInterface *sc_impl, EntityHandle box_set,
         EntitySequence *seq1, EntitySequence *seq2 = NULL);

    //! function to get vertex handle directly from sequence
    /** \param i Parameter being queried
     * \param j Parameter being queried
     * \param k Parameter being queried
     */
  EntityHandle get_vertex_from_seq(int i, int j, int k) const;

    //! set the vertex sequence
  ErrorCode vert_dat(ScdVertexData *vert_dat);
  
    //! get the vertex sequence
  ScdVertexData *vert_dat() const;
  
    //! set the element sequence
  ErrorCode elem_seq(EntitySequence *elem_seq);
  
    //! get the element sequence
  StructuredElementSeq *elem_seq() const;
  
    //! Set the starting vertex handle for this box
  void start_vertex(EntityHandle startv);
  
    //! Set the starting entity handle for this box
  void start_element(EntityHandle starte);

    //! interface instance
  ScdInterface *scImpl;

    //! entity set representing this box
  EntityHandle boxSet;

    //! vertex sequence this box represents, if there's only one, otherwise they're
    //! retrieved from the element sequence
  ScdVertexData *vertDat;

    //! element sequence this box represents
  StructuredElementSeq *elemSeq;
  
    //! starting vertex handle for this box
  EntityHandle startVertex;
  
    //! starting element handle for this box
  EntityHandle startElem;

    //! lower and upper corners
  int boxDims[6];

    //! lower and upper corners of global box
  int globalBoxDims[6];

    //! partition method used to partition global parametric space
  int partMethod;
  
    //! is periodic in i or j
  bool isPeriodic[2];
  
    //! parameter extents
  HomCoord boxSize;

    //! convenience parameters, (boxSize[1]-1)*(boxSize[0]-1) and boxSize[0]-1
  int boxSizeIJ;
  int boxSizeIJM1;
  int boxSizeIM1;
  
};

inline ErrorCode ScdInterface::compute_partition(int part_method, int np, int nr,
                                                 const int *gijk, int *lijk)
{
  ErrorCode rval = MB_SUCCESS;
  switch (part_method) {
    case ALLJORKORI:
    case -1:
        rval = compute_partition_alljorkori(np, nr, gijk, lijk);
        break;
    case ALLJKBAL:
        rval = compute_partition_alljkbal(np, nr, gijk, lijk);
        break;
    case SQIJ:
        rval = compute_partition_sqij(np, nr, gijk, lijk);
        break;
    case SQJK:
        rval = compute_partition_sqjk(np, nr, gijk, lijk);
        break;
  }

  return rval;
}

inline ErrorCode ScdInterface::compute_partition_alljorkori(int np, int nr, const int *gijk, int *lijk) 
{
    // partition *the elements* over the parametric space; 1d partition for now, in the j, k, or i
    // parameters
  if (gijk[4] - gijk[1] > np) {
    int dj = (gijk[4] - gijk[1]) / np;
    int extra = (gijk[4] - gijk[1]) % np;
    lijk[1] = gijk[1] + nr*dj + 
        std::min(nr, extra);
    lijk[4] = lijk[1] + dj + (nr < extra ? 1 : 0);

    lijk[2] = gijk[2]; lijk[5] = gijk[5];
    lijk[0] = gijk[0]; lijk[3] = gijk[3];
  }
  else if (gijk[5] - gijk[2] > np) {
    int dk = (gijk[5] - gijk[2]) / np;
    int extra = (gijk[5] - gijk[2]) % np;
    lijk[2] = gijk[2] + nr*dk + 
        std::min(nr, extra);
    lijk[5] = lijk[2] + dk + (nr < extra ? 1 : 0);

    lijk[1] = gijk[1]; lijk[4] = gijk[4];
    lijk[0] = gijk[0]; lijk[3] = gijk[3];
  }
  else if (gijk[3] - gijk[0] > np) {
    int di = (gijk[3] - gijk[0]) / np;
    int extra = (gijk[3] - gijk[0]) % np;
    lijk[0] = gijk[0] + nr*di + 
        std::min(nr, extra);
    lijk[3] = lijk[0] + di + (nr < extra ? 1 : 0);

    lijk[2] = gijk[2]; lijk[5] = gijk[5];
    lijk[1] = gijk[1]; lijk[4] = gijk[4];
  }
  else {
      // Couldn't find a suitable partition...
    return MB_FAILURE;
  }

  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_alljkbal(int np, int nr, const int *gijk, int *lijk, int *njp) 
{
    // improved, possibly 2-d partition
  std::vector<double> kfactors;
  kfactors.push_back(1);
  int K = gijk[5] - gijk[2];
  for (int i = 2; i < K; i++) 
    if (!(K%i) && !(np%i)) kfactors.push_back(i);
  kfactors.push_back(K);
  
    // compute the ideal nj and nk
  int J = gijk[4] - gijk[1];
  double njideal = sqrt(((double)(np*J))/((double)K));
  double nkideal = (njideal*K)/J;
  
  int nk, nj;
  if (nkideal < 1.0) {
    nk = 1;
    nj = np;
  }
  else {
    std::vector<double>::iterator vit = std::lower_bound(kfactors.begin(), kfactors.end(), nkideal);
    if (vit == kfactors.begin()) nk = 1;
    else nk = (int)*(--vit);
    nj = np / nk;
  }

  int dk = K / nk;
  int dj = J / nj;
  
  lijk[2] = (nr % nk) * dk;
  lijk[5] = lijk[2] + dk;
  
  int extra = J % nj;
  
  lijk[1] = gijk[1] + (nr / nk) * dj + std::min(nr / nk, extra);
  lijk[4] = lijk[1] + dj + (nr / nk < extra ? 1 : 0);

  lijk[0] = gijk[0];
  lijk[3] = gijk[3];

  if (njp) *njp = nj;
  
  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_sqij(int np, int nr, const int *gijk, int *lijk,
                                                      int *pip) 
{
    // square IxJ partition

  if (np == 1) {
      if (lijk) {
	lijk[0] = gijk[0];
	lijk[3] = gijk[3];
	lijk[1] = gijk[1];
	lijk[4] = gijk[4];
	lijk[2] = gijk[2];
	lijk[5] = gijk[5];
      }
    }
  else {
    std::vector<double> pfactors, ppfactors;
    for (int i = 2; i <= np; i++) 
      if (!(np%i)) {
	pfactors.push_back(i);
	ppfactors.push_back(((double)(i*i))/np);
      }
    
    // ideally, Px/Py = I/J
    double ijratio = ((double)(gijk[3]-gijk[0]))/((double)(gijk[4]-gijk[1]));
    
    unsigned int ind = std::lower_bound(ppfactors.begin(), ppfactors.end(), ijratio) - ppfactors.begin();
    if (ind && fabs(ppfactors[ind-1]-ijratio) < fabs(ppfactors[ind]-ijratio)) ind--;
    
    int pi = pfactors[ind];
    int pj = np / pi;
    
    int I = (gijk[3] - gijk[0]), J = (gijk[4] - gijk[1]);
    int iextra = I%pi, jextra = J%pj, i = I/pi, j = J/pj;
    int nri = nr % pi, nrj = nr / pi;
    //
    // pi, pj = # blocks in i, j directions
    // nri, nrj = column, row of this processor's block
    if (lijk) {
      lijk[0] = i*nri + std::min(iextra, nri);
      lijk[3] = lijk[0] + i + (nri < iextra ? 1 : 0);
      lijk[1] = j*nrj + std::min(jextra, nrj);
      lijk[4] = lijk[1] + j + (nrj < jextra ? 1 : 0);
      
      lijk[2] = gijk[2];
      lijk[5] = gijk[5];
    }
	
    if (pip) *pip = pi;
  }  
  return MB_SUCCESS;
}

inline ErrorCode ScdInterface::compute_partition_sqjk(int np, int nr, const int *gijk, int *lijk, int *pjp) 
{
    // square JxK partition

  std::vector<double> pfactors, ppfactors;
  for (int p = 2; p <= np; p++) 
    if (!(np%p)) {
      pfactors.push_back(p);
      ppfactors.push_back(((double)(p*p))/np);
    }
  
    // ideally, Pj/Pk = J/K
  double jkratio = ((double)(gijk[4]-gijk[1]))/((double)(gijk[5]-gijk[2]));

  std::vector<double>::iterator vit  = std::lower_bound(ppfactors.begin(), ppfactors.end(), jkratio);
  if (vit == ppfactors.end()) vit--;
  else if (vit != ppfactors.begin() && fabs(*(vit-1)-jkratio) < fabs((*vit)-jkratio)) vit--;
  int ind = vit - ppfactors.begin();
  
  int pj = pfactors[ind];
  int pk = np / pj;

  int K = (gijk[5] - gijk[2]), J = (gijk[4] - gijk[1]);
  int jextra = J%pj, kextra = K%pk, j = J/pj, k = K/pk;
  int nrj = nr % pj, nrk = nr / pj;
  lijk[1] = j*nrj + std::min(jextra, nrj);
  lijk[4] = lijk[1] + j + (nrj < jextra ? 1 : 0);
  lijk[2] = k*nrk + std::min(kextra, nrk);
  lijk[5] = lijk[2] + k + (nrk < kextra ? 1 : 0);

  lijk[0] = gijk[0];
  lijk[3] = gijk[3];

  if (pjp) *pjp = pj;
  
  return MB_SUCCESS;
}

inline int ScdInterface::gtol(const int *gijk, int i, int j, int k) 
{
  return ((k-gijk[2])*(gijk[3]-gijk[0]+1)*(gijk[4]-gijk[1]+1) + (j-gijk[1])*(gijk[3]-gijk[0]+1) + i-gijk[0]);
}

inline ErrorCode ScdInterface::get_indices(const int *bdy_ind, const int *ldims, int *rdims, int *face_dims, 
                                           std::vector<int> &shared_indices) 
{
  for (int k = face_dims[2]; k <= face_dims[5]; k++)
    for (int j = face_dims[1]; j <= face_dims[4]; j++)
      for (int i = face_dims[0]; i <= face_dims[3]; i++)
        shared_indices.push_back(gtol(ldims, i, j, k));

  if (bdy_ind[0] >= 0) face_dims[0] = face_dims[3] = rdims[3];
  if (bdy_ind[1] >= 0) face_dims[1] = face_dims[4] = rdims[4];
  
  for (int k = face_dims[2]; k <= face_dims[5]; k++)
    for (int j = face_dims[1]; j <= face_dims[4]; j++)
      for (int i = face_dims[0]; i <= face_dims[3]; i++)
        shared_indices.push_back(gtol(rdims, i, j, k));

  return MB_SUCCESS;
}
  
inline ErrorCode ScdInterface::get_neighbor(int pfrom, int np, int part_method,
                                            const int *gdims, const int *ldims,
                                            bool periodic_i, bool periodic_j,
                                            int di, int dj, int dk,
                                            int &pto, int *bdy_ind, int *rdims, int *facedims) 
{
  int dijk[3] = {di, dj, dk};
  
  switch (part_method) {
    case ALLJORKORI:
    case -1:
        return get_neighbor_alljorkori(pfrom, np, gdims, ldims, periodic_i, periodic_j,
                                       dijk, pto, bdy_ind, rdims, facedims);
    case ALLJKBAL:
        return get_neighbor_alljkbal(pfrom, np, gdims, ldims, periodic_i, periodic_j,
                                     dijk, pto, bdy_ind, rdims, facedims);
    case SQIJ:
        return get_neighbor_sqij(pfrom, np, gdims, ldims, periodic_i, periodic_j,
                                 dijk, pto, bdy_ind, rdims, facedims);
    case SQJK:
        return get_neighbor_sqjk(pfrom, np, gdims, ldims, periodic_i, periodic_j,
                                 dijk, pto, bdy_ind, rdims, facedims);
  }

  return MB_FAILURE;
}

inline ScdInterface *ScdBox::sc_impl() const 
{
  return scImpl;
}

inline EntityHandle ScdBox::start_vertex() const
{
  return startVertex;
}
    
inline void ScdBox::start_vertex(EntityHandle startv) 
{
  startVertex = startv;
}
    
inline EntityHandle ScdBox::start_element() const
{
  return startElem;
}
    
inline void ScdBox::start_element(EntityHandle starte) 
{
  startElem = starte;
}
    
inline int ScdBox::num_elements() const
{
  return (!startElem ? 0 : 
          (boxSize[0]- (isPeriodic[0] ? 0 : 1)) * 
          (-1 == boxSize[1] ? 1 : (boxSize[1]-(isPeriodic[1] ? 0 : 1))) * 
          (boxSize[2] == -1 || boxSize[2] == 1 ? 1 : (boxSize[2]-1)));
}
    
inline int ScdBox::num_vertices() const
{
  return boxSize[0] * (!boxSize[1] ? 1 : boxSize[1]) * 
      (!boxSize[2] ? 1 : boxSize[2]);
}
    
inline const int *ScdBox::box_dims() const 
{
  return boxDims;
}

inline const int *ScdBox::global_box_dims() const 
{
  return globalBoxDims;
}

inline void ScdBox::set_global_box_dims(int *box_dims)
{
  for (int i = 0; i < 6; i++) globalBoxDims[i] = box_dims[i];
}

inline int ScdBox::part_method() const 
{
  return partMethod;
}

inline void ScdBox::part_method(int method)
{
  partMethod = method;
}

inline HomCoord ScdBox::box_min() const 
{
  return HomCoord(boxDims, 3);
}

inline HomCoord ScdBox::box_max() const 
{
  return HomCoord(boxDims+3, 3);
}

inline HomCoord ScdBox::box_size() const 
{
  return boxSize;
}

inline void ScdBox::box_size(int *ijk) const 
{
  ijk[0] = boxSize[0];
  ijk[1] = boxSize[1];
  ijk[2] = boxSize[2];
}

inline void ScdBox::box_size(int &i, int &j, int &k) const 
{
  i = boxSize[0];
  j = boxSize[1];
  k = boxSize[2];
}
  
inline EntityHandle ScdBox::get_element(int i, int j, int k) const 
{
  return (!startElem ? MB_ENTITY_NOT_FOUND : 
          startElem + (k-boxDims[2])*boxSizeIJM1 + (j-boxDims[1])*boxSizeIM1 + i-boxDims[0]);
}

inline EntityHandle ScdBox::get_element(HomCoord ijk) const 
{
  return get_element(ijk[0], ijk[1], ijk[2]);
}
  
inline EntityHandle ScdBox::get_vertex(int i, int j, int k) const 
{
  return (vertDat ? startVertex + (k-boxDims[2])*boxSizeIJ + (j-boxDims[1])*boxSize[0] + i-boxDims[0] : 
          get_vertex_from_seq(i, j, k));
}

inline EntityHandle ScdBox::get_vertex(HomCoord ijk) const 
{
  return get_vertex(ijk[0], ijk[1], ijk[2]);
}
  
inline bool ScdBox::contains(const HomCoord ijk) const
{
  return (ijk >= HomCoord(boxDims, 3) && 
          ijk <= HomCoord(boxDims+3, 3));
}

inline bool ScdBox::contains(int i, int j, int k) const
{
  return contains(HomCoord(i, j, k));
}

inline void ScdBox::box_set(EntityHandle this_set) 
{
  boxSet = this_set;
}

inline EntityHandle ScdBox::box_set() 
{
  return boxSet;
}
    
inline ScdVertexData *ScdBox::vert_dat() const
{
  return vertDat;
}
  
inline StructuredElementSeq *ScdBox::elem_seq() const
{
  return elemSeq;
}
  
inline ErrorCode ScdBox::get_params(EntityHandle ent, int &i, int &j, int &k, int &dir) const
{
  HomCoord hc;
  ErrorCode rval = get_params(ent, hc);
  if (MB_SUCCESS == rval) {
    i = hc[0];
    j = hc[1];
    k = hc[2];
    dir = hc[3];
  }
  
  return rval;
}

inline ErrorCode ScdBox::get_params(EntityHandle ent, int &i, int &j, int &k) const 
{
  HomCoord hc;
  ErrorCode rval = get_params(ent, hc);
  if (MB_SUCCESS == rval) {
    i = hc[0];
    j = hc[1];
    k = hc[2];
  }
  
  return rval;
}

inline bool ScdBox::is_periodic_i() const 
{
  return isPeriodic[0];
}

inline bool ScdBox::is_periodic_j() const 
{
  return isPeriodic[1];
}

inline void ScdBox::is_periodic(bool is_periodic_ij[2]) const 
{
  for (int i = 0; i < 2; i++) 
    is_periodic_ij[i] = isPeriodic[i];
}

} // namespace moab
#endif
