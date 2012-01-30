/** 
 * \class moab::Coupler
 * \author Tim Tautges
 *
 * \brief This class couples data between meshes.
 *
 * The coupler interpolates solution data at a set of points.  Data
 * being interpolated resides on a source mesh, in a tag.
 * Applications calling this coupler send in entities, usually points
 * or vertices, and receive back the tag value interpolated at those
 * points.  Entities in the source mesh containing those points 
 * do not have to reside on the same processor.
 *
 * To use, an application should:
 * - instantiate this coupler by calling the constructor collectively
 *   on all processors in the communicator
 * - call locate_points, which locates the points to be interpolated and
 *   (optionally) caches the results in this object
 * - call interpolate, which does the interpolation
 *
 * Multiple interpolations can be done after locating the points.
 *
 */
#ifndef COUPLER_HPP
#define COUPLER_HPP

#include "iBase.h"
#include "moab/Range.hpp"
#include "moab/Interface.hpp"
#include "moab/CartVect.hpp"

namespace moab {

class ParallelComm;

class AdaptiveKDTree;
  
class TupleList;

class Coupler
{
public:

  enum Method {LINEAR_FE, PLAIN_FE};

  enum IntegType {VOLUME};

    /* constructor
     * Constructor, which also optionally initializes the coupler
     * \param pc ParallelComm object to be used with this coupler
     * \param local_elems Local elements in the source mesh
     * \param coupler_id Id of this coupler, should be the same over all procs
     * \param init_tree If true, initializes kdtree inside the constructor
     */
  Coupler(Interface *impl,
            ParallelComm *pc,
            Range &local_elems,
            int coupler_id,
            bool init_tree = true);

    /* destructor
     */
  virtual ~Coupler();
  
    /* \brief Locate points on the source mesh
     * This function finds the element/processor/natural coordinates 
     * containing each point, optionally storing the results locally.
     * \param xyz Point locations (interleaved) being located
     * \param tl Tuple list containing the results, with each tuple
     *           consisting of (p, i), p = proc, i = index on that proc
     * \param store_local If true, stores the tuple list on the Coupler instance
     *
     */
  ErrorCode locate_points(double *xyz, int num_points,
                            TupleList *tl = NULL,
                            bool store_local = true);
  
    /* \brief Locate entities on the source mesh
     * This function finds the element/processor/natural coordinates 
     * containing each entity, optionally storing the results locally.
     * \param ents Entities being located
     * \param tl Tuple list containing the results, with each tuple
     *           consisting of (p, i), p = proc, i = index on that proc
     * \param store_local If true, stores the tuple list on the Coupler instance
     *
     */
  ErrorCode locate_points(Range &ents,
                            TupleList *tl = NULL,
                            bool store_local = true);
  
    /* \brief Interpolate data from the source mesh onto points
     * All entities/points or, if tuple_list is input, only those points
     * are interpolated from the source mesh.  Application should
     * allocate enough memory in interp_vals to hold interpolation results.
     * 
     * If normalization is requested, technique used depends on the coupling
     * method.
     *
     * \param method Interpolation/normalization method
     * \param tag Tag on source mesh holding data to be interpolated
     * \param interp_vals Memory holding interpolated data
     * \param tl Tuple list of points to be interpolated; if NULL, all locations
     *       stored in this object are interpolated
     * \param normalize If true, normalization is done according to method
     */
  ErrorCode interpolate(Coupler::Method method,
                          Tag tag,
                          double *interp_vals,
                          TupleList *tl = NULL,
                          bool normalize = true);

    /* \brief Interpolate data from the source mesh onto points
     * All entities/points or, if tuple_list is input, only those points
     * are interpolated from the source mesh.  Application should
     * allocate enough memory in interp_vals to hold interpolation results.
     * 
     * If normalization is requested, technique used depends on the coupling
     * method.
     *
     * \param method Interpolation/normalization method
     * \param tag_name Name of tag on source mesh holding data to be interpolated
     * \param interp_vals Memory holding interpolated data
     * \param tl Tuple list of points to be interpolated; if NULL, all locations
     *       stored in this object are interpolated
     * \param normalize If true, normalization is done according to method
     */
  ErrorCode interpolate(Coupler::Method method,
                          std::string &tag_name,
                          double *interp_vals,
                          TupleList *tl = NULL,
                          bool normalize = true);

    /* \brief Normalize a field over an entire mesh
     * A field existing on the vertices of elements of a mesh is integrated
     * over all elements in the mesh.  The integrated value is normalized 
     * and the normalization factor is saved to a new tag
     * on the mesh entity set.
     * 
     * \param root_set Entity Set representing the entire mesh
     * \param norm_tag Tag containing field data to integrate
     * \param integ_type Type of integration to perform
     * \param num_integ_pts The number of Gaussian integration points to use in each dimension
     */
  int normalize_mesh(iBase_EntitySetHandle &root_set,
                     const char            *norm_tag,
                     Coupler::IntegType    integ_type,
                     int                   num_integ_pts);

    /* \brief Normalize a field over subsets of entities
     * A field existing on the vertices of elements of a mesh is integrated
     * over subsets of elements identified by the tags and values.  The integrated 
     * values are normalized and the normalization factor is saved to a new tag
     * on the entity sets which contain the elements of a subset.
     * 
     * \param root_set Entity Set from the mesh from which to select subsets
     * \param norm_tag Tag containing field data to integrate
     * \param tag_names Array of tag names used for selecting element subsets
     * \param num_tags Number of tag names
     * \param tag_values Array of tag values passed as strings; the array will be
     *       the same length as that for tag names however some entries may be 
     *       NULL indicating that tag should be matched for existence and not value
     * \param integ_type Type of integration to perform
     * \param num_integ_pts The number of Gaussian integration points to use in each dimension
     */
  int normalize_subset(iBase_EntitySetHandle &root_set,
                       const char            *norm_tag,
                       const char            **tag_names,
                       int                   num_tags,
                       const char            **tag_values,
                       Coupler::IntegType    integ_type,
                       int                   num_integ_pts);

    /* \brief Normalize a field over subsets of entities
     * A field existing on the vertices of elements of a mesh is integrated
     * over subsets of elements identified by the tags and values.  The integrated 
     * values are normalized and the normalization factor is saved to a new tag
     * on the entity sets which contain the elements of a subset.
     * 
     * \param root_set Entity Set from the mesh from which to select subsets
     * \param norm_tag Tag containing field data to integrate
     * \param tag_handles Array of tag handles used for selecting element subsets
     * \param num_tags Number of tag handles
     * \param tag_values Array of tag values passed as strings; the array will be
     *       the same length as that for tag handles however some entries may be 
     *       NULL indicating that tag should be matched for existence and not value
     * \param integ_type Type of integration to perform
     * \param num_integ_pts The number of Gaussian integration points to use in each dimension
     */
  int normalize_subset(iBase_EntitySetHandle &root_set,
                       const char            *norm_tag,
                       iBase_TagHandle       *tag_handles,
                       int                   num_tags,
                       const char            **tag_values,
                       Coupler::IntegType    integ_type,
                       int                   num_integ_pts);

    /* \brief Retrieve groups of entities matching tags and values(if present)
     * Retrieve a vector of vectors of entity handles matching the 
     * tags and values.  The entity set passed is used as the search domain.
     * 
     * \param norm_tag Tag containing field data to integrate
     * \param entity_sets Pointer to vector of vectors of entity set handles
     * \param entity_groups Pointer to vector of vectors of entity handles from each entity set
     * \param integ_type Type of integration to perform
     * \param num_integ_pts The number of Gaussian integration points to use in each dimension
     */
  int do_normalization(const char                                        *norm_tag,
                       std::vector< std::vector<iBase_EntitySetHandle> > &entity_sets,
                       std::vector< std::vector<iBase_EntityHandle> >    &entity_groups,
                       Coupler::IntegType                                integ_type,
                       int                                               num_integ_pts);

    /* \brief Retrieve groups of entities matching tags and values(if present)
     * Retrieve a vector of vectors of entity handles matching the 
     * tags and values.  The entity set passed is used as the search domain.
     * 
     * \param root_set Set from which to search for matching entities
     * \param tag_names Array of tag names used to select entities
     * \param tag_values Array of tag values used to select entities
     * \param num_tags Number of tag names
     * \param entity_sets Pointer to vector of vectors of entity set handles found in the search
     * \param entity_groups Pointer to vector of vectors of entity handles from each entity set
     */
  int get_matching_entities(iBase_EntitySetHandle                             root_set,
                            const char                                        **tag_names,
                            const char                                        **tag_values,
                            int                                               num_tags,
                            std::vector< std::vector<iBase_EntitySetHandle> > *entity_sets,
                            std::vector< std::vector<iBase_EntityHandle> >    *entity_groups);

    /* \brief Retrieve groups of entities matching tags and values(if present)
     * Retrieve a vector of vectors of entity handles matching the 
     * tags and values.  The entity set passed is used as the search domain.
     * 
     * \param root_set Set from which to search for matching entities
     * \param tag_handles Array of tag handles used to select entities
     * \param tag_values Array of tag values used to select entities
     * \param num_tags Number of tag handles
     * \param entity_sets Pointer to vector of vectors of entity set handles found in the search
     * \param entity_groups Pointer to vector of vectors of entity handles from each entity set
     */
  int get_matching_entities(iBase_EntitySetHandle                             root_set,
                            iBase_TagHandle                                   *tag_handles,
                            const char                                        **tag_values,
                            int                                               num_tags,
                            std::vector< std::vector<iBase_EntitySetHandle> > *entity_sets,
                            std::vector< std::vector<iBase_EntityHandle> >    *entity_groups);

    /* \brief Return an array of tuples of tag values for each Entity Set
     * A list of n-tuples will be constructed with 1 n-tuple for each Entity Set.
     * The n-tuple will have an component for each tag given.  It is assumed all
     * of the tags are integer tags.
     * 
     * \param ent_sets Array of Entity Set handles to use for retrieving tag data
     * \param num_sets Number of Entity Sets
     * \param tag_names Array of tag names
     * \param num_tags Number of tag names
     * \param tuples The returned tuple_list structure
     */
  int create_tuples(iBase_EntitySetHandle *ent_sets, 
                    int                   num_sets, 
                    const char            **tag_names, 
                    int                   num_tags,
                    TupleList            **tuples);

    /* \brief Return an array of tuples of tag values for each Entity Set
     * A list of n-tuples will be constructed with 1 n-tuple for each Entity Set.
     * The n-tuple will have an component for each tag given.  It is assumed all
     * of the tags are integer tags.
     * 
     * \param ent_sets Array of Entity Set handles to use for retrieving tag data
     * \param num_sets Number of Entity Sets
     * \param tag_handles Array of tag handles
     * \param num_tags Number of tag handles
     * \param tuples The returned tuple_list structure
     */
  int create_tuples(iBase_EntitySetHandle *ent_sets, 
                    int                   num_sets, 
                    iBase_TagHandle       *tag_handles,
                    int                   num_tags,
                    TupleList            **tuples);

    /* \brief Consolidate an array of n-tuples lists into one n-tuple list with no duplicates
     * An array of list of n-tuples are consolidated into a single list of n-tuples
     * with all duplicates removed.  Only integer columns in the tuple_list are assumed to 
     * be used.
     *
     * \param all_tuples Array of tuple_lists to consolidate to one
     * \param num_tuples Number of tuple_lists
     * \param unique_tuples The consolidated tuple_list with no duplicates
     */
  int consolidate_tuples(TupleList **all_tuples, 
                         int        num_tuples,
                         TupleList **unique_tuples);

    /* \brief Calculate integrated field values for groups of entities
     * An integrated field value, as defined by the field function, 
     * is calculated for each group of entities passed in.
     * 
     * \param groups The vector contains vectors of entity handles, each representing a group
     * \param integ_vals The integrated field values for each group
     * \param norm_tag The tag name of the vertex-based field to be integrated
     * \param num_integ_pts The number of Gaussian integration points to use in each dimension
     * \param integ_type Type of integration to perform
     */
  int get_group_integ_vals(std::vector< std::vector<iBase_EntityHandle> > &groups,
                           std::vector<double> &integ_vals, 
                           const char *norm_tag,
                           int num_integ_pts,
                           Coupler::IntegType integ_type);

    /* \brief Apply a normalization factor to group of entities
     * Multiply a normalization factor with the value of norm_tag for each vertex
     * of each entity in a group.  Save the value back to norm_tag on each vertex.
     *
     * \param entity_sets The vector contains vectors of entity set handles, each containing the members of a group
     * \param norm_factors The normalization factors for each group
     * \param norm_tag The tag to be normalized on each group
     * \param integ_type Type of integration to perform
     */
  int apply_group_norm_factor(std::vector< std::vector<iBase_EntitySetHandle> > &entity_sets,
                              std::vector<double> &norm_factors, 
                              const char *norm_tag,
                              Coupler::IntegType integ_type);

private:

    // given a coordinate position, find all entities containing
    // the point and the natural coords in those ents
  ErrorCode nat_param(double xyz[3], 
                        std::vector<EntityHandle> &entities, 
                        std::vector<CartVect> &nat_coords);
  
  ErrorCode interp_field(EntityHandle elem,
                         CartVect nat_coord, 
                         Tag tag,
                         double &field);

  ErrorCode plain_field_map(EntityHandle elem,
			      Tag tag,
			      double &field);
  
  ErrorCode test_local_box(double *xyz, 
                             int from_proc, int remote_index, int index, 
                             bool &point_located,
                             TupleList *tl = NULL);
  
    /* \brief MOAB instance
     */
  Interface *mbImpl;
  
    /* \brief Initialize the kdtree, locally and across communicator
     */
  ErrorCode initialize_tree();

    /* \brief Kdtree for local mesh
     */
  AdaptiveKDTree *myTree;
  
    /* \brief Local root of the kdtree
     */
  EntityHandle localRoot;

    /* \brief Min/max bounding boxes for all proc tree roots
     */
  std::vector<double> allBoxes;
  
    /* \brief ParallelComm object for this coupler
     */
  ParallelComm *myPc;
  
    /* \brief Id of this coupler
     */
  int myId;
  
    /* \brief Range of locations interpolated onto
     */
  Range myRange;

    /* \brief List of locally mapped tuples
     * Tuples contain the following:
     * n = # mapped points
     * vul[i] = local handle of mapped entity
     * vr[3*i..3*i+2] = natural coordinates in mapped entity
     */
  TupleList *mappedPts;
  
    /* \brief Tuple list of target points and interpolated data
     * Tuples contain the following:
     * n = # target points
     * vi[3*i] = remote proc mapping target point
     * vi[3*i+1] = local index of target point
     * vi[3*i+2] = remote index of target point
     * vr[i] = interpolated data (used by interpolate function)
     */
  TupleList *targetPts;

    /* \brief Locally mapped points
     * Points whose source and target are both local; these
     * points consist of two indices, <target_index, mapped_index>,
     * where target_index is the index in the target points array
     * and mapped_index is the corresponding index into mappedPts
     */
  std::vector<unsigned int> localMappedPts;

    /* \brief Number of iterations of tree building before failing
     *
     */
  int numIts;
};

} // namespace moab

#endif
