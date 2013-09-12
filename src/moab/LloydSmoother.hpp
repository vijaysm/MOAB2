/** @class LloydSmoother.cpp \n
 * \brief Perform Lloyd relaxation on a mesh and its dual \n
 *
 * Briefly, Lloyd relaxation is a technique to smooth out a mesh.  The centroid of each cell is computed from its
 * vertex positions, then vertices are placed at the average of their connected cells' centroids, and the process
 * iterates until convergence.
 *
 * In the parallel algorithm, an extra ghost layer of cells is exchanged.  This allows us to compute the centroids
 * for boundary cells on each processor where they appear; this eliminates the need for one round of data exchange
 * (for those centroids) between processors.  New vertex positions must be sent from owning processors to processors
 * sharing those vertices.  Convergence is measured as the maximum distance moved by any vertex.  
 * 
 */

#ifndef LLOYDSMOOTHER_HPP
#define LLOYDSMOOTHER_HPP

#include "moab/Interface.hpp"

namespace moab {

class ParallelComm;
class Error;

class LloydSmoother
{
public:

    /* \brief Constructor
     * Bare constructor, data input to this class through methods.
     * \param impl The MOAB instance for this smoother
     */
  LloydSmoother(Interface *impl);
  
    /* \brief Constructor
     * Convenience constructor, data input directly
     * \param impl The MOAB instance for this smoother
     * \param pcomm The ParallelComm instance by which this mesh is parallel
     * \param elems The mesh to be smoothed
     * \param fixed_tag The tag marking which vertices are fixed
     * \param abs_tol Absolute tolerance measuring convergence
     * \param rel_tol Relative tolerance measuring convergence
     */
  LloydSmoother(Interface *impl, ParallelComm *pc, Range &elems, Tag fixed_tag = 0,
                double abs_tol = -1.0, double rel_tol = 1.0e-6);

    /* \brief Destructor
     */
  ~LloydSmoother();

    /* \brief perform smoothing operation
     */
  ErrorCode perform_smooth();
  
    /* \brief get instance
     */
  Interface *mb_impl() {return mbImpl;}
  
    /* \brief get/set ParallelComm
     */
  ParallelComm *pcomm() {return myPcomm;}
  
    /* \brief get/set ParallelComm
     */
  void pcomm(ParallelComm *pc) {myPcomm = pc;}

    /* \brief get/set elements
     */
  Range &elems() {return myElems;}
  
    /* \brief get/set elements
     */
  const Range &elems() const {return myElems;}
  
    /* \brief get/set fixed tag
     */
  Tag fixed_tag() {return fixedTag;}

    /* \brief get/set fixed tag
     */
  void fixed_tag(Tag fixed) {fixedTag = fixed;}

    /* \brief get/set tolerance
     */
  double abs_tol() {return absTol;}
  
    /* \brief get/set tolerance
     */
  void abs_tol(double tol) {absTol = tol;}
  
    /* \brief get/set tolerance
     */
  double rel_tol() {return relTol;}
  
    /* \brief get/set tolerance
     */
  void rel_tol(double tol) {relTol = tol;}

    /* \brief get/set numIts
     */
  int num_its() {return numIts;}
  void num_its(int num) {numIts = num;}
        
    /* \brief get/set reportIts
     */
  int report_its() {return reportIts;}
  void report_its(int num) {reportIts = num;}
        
  
private:

    //- initialize some things in certain cases
  ErrorCode initialize();
  
    //- MOAB instance
  Interface *mbImpl;

    //- ParallelComm
  ParallelComm *myPcomm;
  
    //- elements to smooth
  Range myElems;
  
    //- tag marking which vertices are fixed, 0 = not fixed, otherwise fixed
  Tag fixedTag;
  
    //- tolerances
  double absTol, relTol;

    //- error handler for this class
  Error *errorHandler;

    //- number of iterations between reporting
  int reportIts;

    //- number of iterations taken during smooth
  int numIts;

    //- keep track of whether I created the fixed tag
  bool iCreatedTag;
};
    
} // namespace moab

#endif
