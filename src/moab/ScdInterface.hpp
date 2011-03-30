/** \file ScdInterface.hpp
 */
#ifndef SCD_INTERFACE
#define SCD_INTERFACE

#include <vector>    

#include "moab/Interface.hpp"
#include "moab/HomXform.hpp"

namespace moab {

class StructuredElementSeq;
class EntitySequence;
class ScdVertexData;
class EntitySequence;
class ScdBox;
class Core;

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
 * created and tagged with information about the parametric space.  In particular, the BOX_MIN
 * and BOX_MAX tags are used to indicate the lower and upper corners in parametric space (each of these
 * tags is integer size 3).  Structured mesh blocks are also available through ScdBox class objects
 * returned by ScdInterface.  These class objects should be treated only as references to the 
 * structured mesh blocks; that is, the structured mesh referenced by these objects is not deleted
 * when the ScdBox instance is destroyed.  Functions for destroying the actual mesh are are available 
 * on this class, though.
 *
 * Structured mesh blocks are returned in the form of ScdBox class objects.  Each ScdBox instance
 * represents a rectangular block of vertices and possibly elements (edges, quads, or hexes).  The
 * entity handles for a ScdBox are guaranteed to be contiguous, starting at a starting valuable
 * which is also available through the ScdBox class.  However, vertex handles may or may not be
 * contiguous, depending on the construction method.  The start vertex handle is also available from
 * the ScdBox class.
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
     * structured blocks, based on the BOX_MIN tag.  Structured mesh blocks will be stored
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
     * \param coords Coordinates of vertices, in column-major order
     * \param new_box Reference to box of structured mesh
     */
  ErrorCode construct_box(HomCoord low, HomCoord high, double *coords, unsigned int num_coords,
                          ScdBox *& new_box);

    //! Create a structured sequence of vertices, quads, or hexes
    /** Starting handle for the sequence is available from the returned ScdBox.  
     * If creating a structured quad or hex box, subsequent calls must be made to ScdBox::add_vbox, 
     * until all the vertices have been filled in for the box.
     * \param low Lower corner of structured box
     * \param high Higher corner of structured box
     * \param type EntityType, one of MBVERTEX, MBEDGE, MBQUAD, MBHEX
     * \param starting_id Requested start id of entities
     * \param new_box Reference to the newly created box of entities
     */
  ErrorCode create_scd_sequence(HomCoord low, HomCoord high, EntityType type,
                                int starting_id, ScdBox *&new_box);

    //! Return all the structured mesh blocks in this MOAB instance
    /** Return the structured blocks in this MOAB instance.  If these were not searched for
     * at instantiation time, then the search is done now.
     * \param boxes Vector of ScdBox objects representing structured mesh blocks
     */
  ErrorCode find_boxes(std::vector<ScdBox*> &boxes);

    //! Return all the structured mesh blocks in this MOAB instance
    /** Return the structured blocks in this MOAB instance.  If these were not searched for
     * at instantiation time, then the search is done now.
     * \param boxes Range of entity set objects representing structured mesh blocks
     */
  ErrorCode find_boxes(Range &boxes);

    //! Return the tag marking the lower corner of boxes
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_min_tag(bool create_if_missing = true);

    //! Return the tag marking the upper corner of boxes
    /**
     * \param create_if_missing If the tag does not yet exist, create it
     */
  Tag box_max_tag(bool create_if_missing = true);

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
  
private:
    //! Create an entity set for a box, and tag with the parameters
    /** \param low Lower corner parameters for this box
     * \param high Upper corner parameters for this box
     * \param scd_set Entity set created
     */
  ErrorCode create_box_set(const HomCoord low, const HomCoord high,
                           EntityHandle &scd_set);
  
    //! interface instance
  Core *mbImpl;

    //! whether we've searched the database for boxes yet
  bool searchedBoxes;
  
    //! structured mesh blocks; stored as entity set handles since application
    //! controls ScdBox objects
  Range scdBoxes;

    //! tag representing box lower corner
  Tag boxMinTag;

    //! tag representing box upper corner
  Tag boxMaxTag;
  
    //! tag pointing from set to ScdBox
  Tag boxSetTag;
  
};

class ScdBox 
{
  friend class ScdInterface;
  
public:

  ~ScdBox();

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
  int box_dimension() const;
  
    //! Starting vertex handle for this box
  EntityHandle start_vertex() const;
  
    //! Starting entity handle for this box
    /** If this is a vertex box, the start vertex handle is returned.
     */
  EntityHandle start_element() const;

    //! Return the number of elements in the box
    /* Number of elements is (boxSize[0]-1)(boxSize[1]-1)(boxSize[2]-1)
     */
  int num_elements() const;
  
    //! Return the number of vertices in the box
    /* Number of vertices is boxSize[0] * boxSize[1] * boxSize[2]
     */
  int num_vertices() const;
  
    //! Return the lower corner parametric coordinates for this box
  HomCoord box_min() const;
  
    //! Return the lower corner parametric coordinates for this box
    /**
     * \param ijk IJK parameters of lower corner
     */
  void box_min(int *ijk) const;
  
    //! Return the lower corner parametric coordinates for this box
    /**
     * \param i I parameter of lower corner
     * \param j J parameter of lower corner
     * \param k K parameter of lower corner
     */
  void box_min(int &i, int &j, int &k) const;
  
    //! Return the upper corner parametric coordinates for this box
  HomCoord box_max() const;
  
    //! Return the upper corner parametric coordinates for this box
    /**
     * \param ijk IJK parameters of upper corner
     */
  void box_max(int *ijk) const;
  
    //! Return the upper corner parametric coordinates for this box
    /**
     * \param i I parameter of upper corner
     * \param j J parameter of upper corner
     * \param k K parameter of upper corner
     */
  void box_max(int &i, int &j, int &k) const;
  
    //! Return the parameter extents for this box
  HomCoord box_size() const;
  
    //! Return the parametric extents for this box
    /**
     * \param ijk IJK extents of this box
     */
  void box_size(int *ijk) const;
  
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
  
    //! Get the entity(ies) of specified dimension adjacent to parametric element
    /**
     * \param dim Dimension of adjacent entity being requested
     * \param i Parametric coordinates of cell being evaluated
     * \param j Parametric coordinates of cell being evaluated
     * \param k Parametric coordinates of cell being evaluated
     * \param dir Direction (0, 1, or 2), for getting adjacent edges (2d, 3d) or faces (3d) 
     * \param ents EntityHandle vector in which results are returned; should be large enough
     *             to hold results, which could be up to 8 handles
     */
  ErrorCode get_adjancies(int dim, int i, int j, int k, int dir, EntityHandle *ents) const;

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

    //! lower corner
  HomCoord boxMin;
  
    //! upper corner
  HomCoord boxMax;
  
    //! parameter extents
  HomCoord boxSize;

    //! convenience parameters, (boxSize[1]-1)*(boxSize[0]-1) and boxSize[0]-1
  int boxSizeIJ;
  int boxSizeIJM1;
  int boxSizeIM1;
  
};

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
  return (!startElem ? 0 : (boxSize[0]-1) * (-1 == boxSize[1] ? 1 : (boxSize[1]-1)) * 
          (boxSize[2] == -1 ? 1 : (boxSize[2]-1)));
}
    
inline int ScdBox::num_vertices() const
{
  return boxSize[0] * (!boxSize[1] ? 1 : boxSize[1]) * 
      (!boxSize[2] ? 1 : boxSize[2]);
}
    
inline HomCoord ScdBox::box_min() const 
{
  return boxMin;
}

inline void ScdBox::box_min(int *ijk) const 
{
  ijk[0] = boxMin[0];
  ijk[1] = boxMin[1];
  ijk[2] = boxMin[2];
}
  
inline void ScdBox::box_min(int &i, int &j, int &k) const 
{
  i = boxMin[0];
  j = boxMin[1];
  k = boxMin[2];
}
  
inline HomCoord ScdBox::box_max() const 
{
  return boxMax;
}

inline void ScdBox::box_max(int *ijk) const 
{
  ijk[0] = boxMax[0];
  ijk[1] = boxMax[1];
  ijk[2] = boxMax[2];
}
  
inline void ScdBox::box_max(int &i, int &j, int &k) const 
{
  i = boxMax[0];
  j = boxMax[1];
  k = boxMax[2];
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
  
inline EntityHandle ScdBox::get_element(HomCoord ijk) const 
{
  return get_element(ijk[0], ijk[1], ijk[2]);
}
  
inline EntityHandle ScdBox::get_element(int i, int j, int k) const 
{
  return (!startElem ? MB_ENTITY_NOT_FOUND : 
          startElem + (k-boxMin[2])*boxSizeIJM1 + (j-boxMin[1])*boxSizeIM1 + i-boxMin[0]);
}

inline EntityHandle ScdBox::get_vertex(HomCoord ijk) const 
{
  return get_vertex(ijk[0], ijk[1], ijk[2]);
}
  
inline EntityHandle ScdBox::get_vertex(int i, int j, int k) const 
{
  return (!startVertex ? MB_ENTITY_NOT_FOUND : 
          vertDat ? startVertex + (k-boxMin[2])*boxSizeIJ + (j-boxMin[1])*boxSize[0] + i-boxMin[0] : 
          get_vertex_from_seq(i, j, k));
}

inline bool ScdBox::contains(int i, int j, int k) const
{
  return contains(HomCoord(i, j, k));
}

inline bool ScdBox::contains(const HomCoord ijk) const
{
  return (ijk >= boxMin && ijk <= boxMax);
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
} // namespace moab
#endif