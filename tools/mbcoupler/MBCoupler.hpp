/** 
 * \class MBCoupler
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
#ifndef MBCOUPLER_HPP
#define MBCOUPLER_HPP

class MBParallelComm;

#include "MBRange.hpp"
#include "MBInterface.hpp"
#include "MBCartVect.hpp"

extern "C" 
{
  struct tuple_list;
}

class MBAdaptiveKDTree;

class MBCoupler
{
public:

  enum Method {LINEAR_FE, PLAIN_FE};

    /* constructor
     * Constructor, which also optionally initializes the coupler
     * \param pc MBParallelComm object to be used with this coupler
     * \param local_elems Local elements in the source mesh
     * \param coupler_id Id of this coupler, should be the same over all procs
     * \param init_tree If true, initializes kdtree inside the constructor
     */
  MBCoupler(MBInterface *impl,
            MBParallelComm *pc,
            MBRange &local_elems,
            int coupler_id,
            bool init_tree = true);

    /* destructor
     */
  virtual ~MBCoupler();
  
    /* \brief Locate points on the source mesh
     * This function finds the element/processor/natural coordinates 
     * containing each point, optionally storing the results locally.
     * \param xyz Point locations (interleaved) being located
     * \param tl Tuple list containing the results, with each tuple
     *           consisting of (p, i), p = proc, i = index on that proc
     * \param store_local If true, stores the tuple list on the MBCoupler instance
     *
     */
  MBErrorCode locate_points(double *xyz, int num_points,
                            tuple_list *tl = NULL,
                            bool store_local = true);
  
    /* \brief Locate entities on the source mesh
     * This function finds the element/processor/natural coordinates 
     * containing each entity, optionally storing the results locally.
     * \param ents Entities being located
     * \param tl Tuple list containing the results, with each tuple
     *           consisting of (p, i), p = proc, i = index on that proc
     * \param store_local If true, stores the tuple list on the MBCoupler instance
     *
     */
  MBErrorCode locate_points(MBRange &ents,
                            tuple_list *tl = NULL,
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
  MBErrorCode interpolate(MBCoupler::Method method,
                          MBTag tag,
                          double *interp_vals,
                          tuple_list *tl = NULL,
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
  MBErrorCode interpolate(MBCoupler::Method method,
                          std::string &tag_name,
                          double *interp_vals,
                          tuple_list *tl = NULL,
                          bool normalize = true);

private:

    // given a coordinate position, find all entities containing
    // the point and the natural coords in those ents
  MBErrorCode nat_param(double xyz[3], 
                        std::vector<MBEntityHandle> &entities, 
                        std::vector<MBCartVect> &nat_coords);
  
  MBErrorCode interp_field_for_hex(MBEntityHandle elem,
                                   MBCartVect nat_coord, 
                                   MBTag tag,
                                   double &field);

  MBErrorCode plain_field_map(MBEntityHandle elem,
			      MBTag tag,
			      double &field);
  
  MBErrorCode test_local_box(double *xyz, 
                             int from_proc, int remote_index, int index, 
                             bool &point_located,
                             tuple_list *tl = NULL);
  
    /* \brief MOAB instance
     */
  MBInterface *mbImpl;
  
    /* \brief Initialize the kdtree, locally and across communicator
     */
  MBErrorCode initialize_tree();

    /* \brief Kdtree for local mesh
     */
  MBAdaptiveKDTree *myTree;
  
    /* \brief Local root of the kdtree
     */
  MBEntityHandle localRoot;

    /* \brief Min/max bounding boxes for all proc tree roots
     */
  std::vector<double> allBoxes;
  
    /* \brief MBParallelComm object for this coupler
     */
  MBParallelComm *myPc;
  
    /* \brief Id of this coupler
     */
  int myId;
  
    /* \brief Range of locations interpolated onto
     */
  MBRange myRange;

    /* \brief List of locally mapped tuples
     * Tuples contain the following:
     * n = # mapped points
     * vul[i] = local handle of mapped entity
     * vr[3*i..3*i+2] = natural coordinates in mapped entity
     */
  tuple_list *mappedPts;
  
    /* \brief Tuple list of target points and interpolated data
     * Tuples contain the following:
     * n = # target points
     * vi[3*i] = remote proc mapping target point
     * vi[3*i+1] = local index of target point
     * vi[3*i+2] = remote index of target point
     * vr[i] = interpolated data (used by interpolate function)
     */
  tuple_list *targetPts;

    /* \brief Locally mapped points
     * Points whose source and target are both local; these
     * points consist of two indices, <target_index, mapped_index>,
     * where target_index is the index in the target points array
     * and mapped_index is the corresponding index into mappedPts
     */
  std::vector<unsigned int> localMappedPts;

};

#endif
