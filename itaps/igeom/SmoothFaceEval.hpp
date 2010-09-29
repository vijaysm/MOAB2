#ifndef SMOOTH_FACE_EVAL_HPP
#define SMOOTH_FACE_EVAL_HPP

// do we really need iMesh here; maybe go directly to MOAB
//#include "iMesh.h"
#include "moab/Interface.hpp"
#include "moab/Range.hpp"
#include "moab/CartVect.hpp"
#include "MBTagConventions.hpp"
#include "moab/Types.hpp"

#define determ3(p1,q1,p2,q2,p3,q3) ((q3)*((p2)-(p1)) + (q2)*((p1)-(p3)) + (q1)*((p3)-(p2)))
#define sqr(a) ((a)*(a))
#define cube(a) (sqr(a) * (a))
#define quart(a) (sqr(a) * sqr(a))
#define blend(x) (-2.0*(x)*(x)*(x) + 3.0*(x)*(x))

#include <vector>
#include <map>
//#include "MBEntityHandle.hpp"

// work only with CAMAL > = 500
// #if CAMAL_VERSION < 500

// #else

//class RefFace;
#include "moab/GeomTopoTool.hpp"

class SmoothCurveEval;// it is derived from SmoothBase, maybe just need

//! Implement CAMAL geometry callbacks using smooth iMesh 
class SmoothFaceEval // public CMLSurfEval, public SmoothBase
{
public:
   SmoothFaceEval(moab::Interface * mb, moab::EntityHandle surface_set,
         moab::GeomTopoTool * gTool); // entity or entity set

   virtual ~SmoothFaceEval();

   virtual double area();

   virtual void bounding_box(double box_min[3], double box_max[3]);

   virtual void move_to_surface(double& x, double& y, double& z);
   /*
    virtual void move_to_surface(double& x, double& y, double& z,
    double& u_guess, double& v_guess);*/

   virtual bool normal_at(double x, double y, double z, double& nx, double& ny,
         double& nz);

   // initialize normals// they will be stored as tags on nodes
   int init_gradient();

   // some functions for edge evaluations
   moab::ErrorCode evaluate_smooth_edge(moab::EntityHandle eh, double &t,
         moab::CartVect & outv);

   moab::ErrorCode eval_bezier_patch(moab::EntityHandle tri,
         moab::CartVect &areacoord, moab::CartVect &pt);

   void project_to_facet_plane(moab::EntityHandle tri, moab::CartVect &pt,
         moab::CartVect &point_on_plane, double &dist_to_plane);

   void facet_area_coordinate(moab::EntityHandle facet,
         moab::CartVect & pt_on_plane, moab::CartVect & areacoord);

   moab::ErrorCode project_to_facets_main(moab::CartVect &this_point,
         bool trim, bool & outside, moab::CartVect * closest_point_ptr = NULL, // interested only in normal
         moab::CartVect * normal_ptr = NULL); // interested only in closest point

   moab::ErrorCode project_to_facets(
         std::vector<moab::EntityHandle> & facet_list,
         moab::EntityHandle & lastFacet, int interpOrder, double compareTol,
         moab::CartVect &this_point, bool trim, bool & outside,
         moab::CartVect *closest_point_ptr, moab::CartVect * normal_ptr);

   moab::ErrorCode project_to_facet(moab::EntityHandle facet,
         moab::CartVect &pt, moab::CartVect &areacoord,
         moab::CartVect &close_point, bool &outside_facet, double compare_tol);

   bool is_at_vertex(moab::EntityHandle facet, // (IN) facet we are evaluating
         moab::CartVect &pt, // (IN) the point
         moab::CartVect &ac, // (IN) the ac of the point on the facet plane
         double compare_tol, // (IN) return TRUE of closer than this
         moab::CartVect &eval_pt, // (OUT) location at vertex if TRUE
         moab::CartVect *eval_norm_ptr); // (OUT) normal at vertex if TRUE

   moab::ErrorCode project_to_patch(moab::EntityHandle facet, // (IN) the facet where the patch is defined
         moab::CartVect &ac, // (IN) area coordinate initial guess (from linear facet)
         moab::CartVect &pt, // (IN) point we are projecting to patch
         moab::CartVect &eval_pt, // (OUT) The projected point
         moab::CartVect *eval_norm, // (OUT) normal at evaluated point
         bool &outside, // (OUT) the closest point on patch to pt is on an edge
         double compare_tol, // (IN) comparison tolerance
         int edge_id); // (IN) only used if this is to be projected to one
   //      of the edges.  Otherwise, should be -1

   moab::ErrorCode eval_bezier_patch_normal(moab::EntityHandle facet,
         moab::CartVect &areacoord, moab::CartVect &normal);

   // this will be called now from driver...
   moab::ErrorCode compute_tangents_for_each_edge();// they will be used for control points

   moab::ErrorCode get_normals_for_vertices(const moab::EntityHandle * conn2,
         moab::CartVect N[2]);// here we need the gradient tag

   // make this one public, will be called during initialization !!!
   moab::ErrorCode init_edge_control_points(moab::CartVect &P0,
         moab::CartVect &P3, moab::CartVect &N0, moab::CartVect &N3,
         moab::CartVect &T0, moab::CartVect &T3, moab::CartVect * Pi);

   // moved from private, because will be called from PaveDriver
   moab::ErrorCode compute_control_points_on_edges(double min_dot,
         moab::Tag edgeCtrlTag, moab::Tag markTag);

   moab::ErrorCode compute_internal_control_points_on_facets(double min_dot,
         moab::Tag facetCtrlTag, moab::Tag facetEdgeCtrlTag);

   // move from private too
   void DumpModelControlPoints();

   //
   moab::ErrorCode find_loops();

   int eval_counter() {
      return _evaluationsCounter;
   }
private:

   //===========================================================================
   //Function Name: move_ac_inside
   //
   //Member Type:  PRIVATE
   //Description:  find the closest area coordinate to the boundary of the
   //              patch if any of its components are < 0
   //              Return if the ac was modified.
   //===========================================================================
   bool move_ac_inside(moab::CartVect &ac, double tol);

   //===========================================================================
   //Function Name: ac_at_edge
   //
   //Member Type:  PRIVATE
   //Description:  determine the area coordinate of the facet at the edge
   //===========================================================================
   void ac_at_edge(moab::CartVect &fac, // facet area coordinate
         moab::CartVect &eac, // edge area coordinate
         int edge_id); // id of edge

   // some local functions that do not need to be public
   moab::ErrorCode init_bezier_edge(moab::EntityHandle edge, double min_dot);
   //

   moab::ErrorCode find_edges_orientations(moab::EntityHandle edges[3],
         const moab::EntityHandle * conn3, int orient[3]); // maybe we will set it?

   moab::ErrorCode init_facet_control_points(moab::CartVect N[6], // vertex normals (per edge)
         moab::CartVect P[3][5], // edge control points
         moab::CartVect G[6]); // return internal control points

   // those are the bounding box limits;
   // they are adjusted for the control points too in each triangle
   void adjust_bounding_box(moab::CartVect & vect);
   moab::CartVect _minim;
   moab::CartVect _maxim;

   moab::Range _triangles;
   moab::Range _edges;
   moab::Range _nodes;

   //std::vector<double> _fractions;// they are increasing from 0. to 1., do we need these?
   //std::vector<double> _loopLengths;

   // each loop will be actually a vector of moab::EntityHandle, paired with a vector of senses
   // number of loops is decided by the size of _loopEnds
   std::vector<moab::EntityHandle> _loops; // set1, set 3, set 5, ...
   std::vector<char> _senses; // 0 forward, 1 backward: 0, 0, 1, ...
   std::vector<int> _loopEnds;// the first loop starts at 0 always;

   // number of loops is decided by the size of _loopEnds.size()
   // this ref face will be gone, we will replace it with a new call
   //RefFace * _smooth_face;
   //int myDimension;
   //double meshSize;

   // this tag is on edges
   // rval = _mb->tag_create("MARKER", 1, MB_TAG_BIT, _markTag, &value);
   moab::Tag _markTag; // this is a tag used to mark edges when we look for loops

   // this tag is on nodes
   //moab::ErrorCode rval = _mb->tag_create("GRADIENT", 3 * sizeof(double),
   // MB_TAG_DENSE, _gradientTag, &defNormal);
   moab::Tag _gradientTag; // this will be used for normal at nodes

   // this tag is on edges
   //moab::ErrorCode rval = _mb->tag_create("TANGENTS", 6 * sizeof(double),
   // MB_TAG_DENSE, _tangentsTag, &defTangents);
   moab::Tag _tangentsTag; // each edge will have exactly 2 tangents, because there is only
   // one feature edge, and it is periodic
   // the feature edge is exactly on the boundary

   // this tag is on edges
   //moab::ErrorCode rval = _mb->tag_create("CONTROLEDGE", 9 * sizeof(double),
   // MB_TAG_DENSE, _edgeCtrlTag, &defControls);
   moab::Tag _edgeCtrlTag;

   // this tag is on facets (triangles), 6 control points on each facet
   // there are also some 15 points used in evaluation; how to store them?
   //moab::ErrorCode rval = _mb->tag_create("CONTROLFACE", 18 * sizeof(double),
   // MB_TAG_DENSE, _facetCtrlTag, &defControls);
   moab::Tag _facetCtrlTag;

   // these are the 12 points stored for each edge, again
   // it is cheaper this way compared to retrieve the edges every time, determine their orientation, order
   // in triangle, and retrieve the control points from the edge
   // the control points are stored as 12 points, in order : edge 0, 1, and 2, in that order
   //moab::ErrorCode rval = _mb->tag_create("CONTROLEDGEFACE", 36 * sizeof(double),
   // MB_TAG_DENSE, _facetEdgeCtrlTag, &defControls);
   moab::Tag _facetEdgeCtrlTag; //
   // plane of the facet stores as a normal a, b, c and d, distance, for
   // ax+by+cz+d=0
   //moab::ErrorCode rval = _mb->tag_create("PLANE", 4 * sizeof(double),
   // MB_TAG_DENSE, _planeTag, &defPlane);
   moab::Tag _planeTag;

   // counter for calls
   long _evaluationsCounter;

   moab::Interface * _mb;
   moab::EntityHandle _set;
   moab::GeomTopoTool * _my_geomTopoTool;
   moab::EntityHandle _obb_root;
};
// #endif

#endif /* SMOOTH_FACE_EVAL_HPP*/
