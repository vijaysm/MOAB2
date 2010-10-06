/*
 * SmoothCurveEval.hpp
 *
 *  Created on: Jun 9, 2010
 *      Author: iulian
 */

#ifndef SMOOTHCURVEEVAL_HPP_
#define SMOOTHCURVEEVAL_HPP_

#include "moab/Interface.hpp"
#include "moab/Forward.hpp"
#include "moab/CartVect.hpp"

//#include "RefEdge.hpp"
//#include "SmoothFaceEval.hpp"

#include <map>
#include <vector>

class SmoothFaceEval;
//class SmoothVertex;

//class CMLEdgeMesher;
// evaluator for Camal Edge Mesher
// this is really copying what Cubit is doing

class SmoothCurveEval
{
public:
	//SmoothCurveEval(RefEdge * edge, SmoothFaceEval * smoothFaceEval, int loopIndex);
	SmoothCurveEval(moab::Interface * mb, moab::EntityHandle curve); // the new constructor, which will use
	// sense entities to establish the control points on feature edges (geo edges, sets of mesh edges)
	virtual ~SmoothCurveEval();

	virtual double arc_length() ;

	    //! \brief Get the parametric status of the curve.
	    //!
	    //! \return \a true if curve is parametric, \a false otherwise.
	  virtual bool is_parametric() ;

	    //! \brief Get the periodic status of the curve.
	    //!
	    //! \param period The period of the curve if periodic.
	    //!
	    //! \return \a true if curve is periodic, \a false otherwise.
	  virtual bool is_periodic(double& period) ;

	    //! \brief Get the parameter range of the curve.
	    //!
	    //! \param u_start The beginning curve parameter
	    //! \param u_end The ending curve parameter
	    //!
	    //! \note The numerical value of \a u_start may be greater
	    //! than the numerical value of \a u_end.
	  virtual void get_param_range(double& u_start, double& u_end) ;

	    //! Compute the parameter value at a specified distance along the curve.
	    //!
	    //! \param u_root The start parameter from which to compute the distance
	    //! along the curve.
	    //! \param arc_length The distance to move along the curve.
	    //!
	    //! \note For positive values of \a arc_length the distance will be
	    //! computed in the direction of increasing parameter value along the
	    //! curve.  For negative values of \a arc_length the distance will be
	    //! computed in the direction of decreasing parameter value along the
	    //! curve.
	    //!
	    //! \return The parametric coordinate u along the curve
	  virtual double u_from_arc_length(double u_root, double arc_length);


	    //! \brief Evaluate the curve at a specified parameter value.
	    //!
	    //! \param u The parameter at which to evaluate the curve
	    //! \param x The x coordinate of the evaluated point
	    //! \param y The y coordinate of the evaluated point
	    //! \param z The z coordinate of the evaluated point
	  virtual bool position_from_u(double u,
	                               double& x, double& y, double& z ) ;
         
	    //! \brief Move a point near the curve to the closest point on the curve.
	    //!
	    //! \param x The x coordinate of the point
	    //! \param y The y coordinate of the point
	    //! \param z The z coordinate of the point
	  virtual void move_to_curve(double& x, double& y, double& z);

	    //! Get the u parameter value on the curve closest to x,y,z
	    //! and the point on the curve.
	    //!
	    //! \param x The x coordinate of the point
	    //! \param y The y coordinate of the point
	    //! \param z The z coordinate of the point
	    //!
	    //! \return The parametric coordinate u on the curve
	  virtual double u_from_position(double x, double y, double z) ;

	    //! \brief Get the starting point of the curve.
	    //!
	    //! \param x The x coordinate of the start point
	    //! \param y The y coordinate of the start point
	    //! \param z The z coordinate of the start point
	  virtual void start_coordinates(double& x, double& y, double& z) ;

	    //! \brief Get the ending point of the curve.
	    //!
	    //! \param x The x coordinate of the start point
	    //! \param y The y coordinate of the start point
	    //! \param z The z coordinate of the start point
	  virtual void end_coordinates(double& x, double& y, double& z) ;

	  // this will recompute the 2 tangents for each edge, considering the geo edge they are into
	  void compute_tangents_for_each_edge();

	  void compute_control_points_on_boundary_edges(double min_dot,
			  std::map<moab::EntityHandle, SmoothFaceEval*> & mapSurfaces,
			  moab::Tag controlPointsTag, moab::Tag markTag);// min_dot is not used now,
	  // the edges that are computed will be marked with the marker tag, so after that only the interior edges
	  // will be left for control points evaluation

	  // make this public, as we may need it when we call meshing from driver
//	  CMLEdgeMesher * _cmlEdgeMesher;

/*	  void estimate_mesh_count(double curve_mesh_size, int & num_points_out); // this is based on initial settings

	  int get_mesh_count()
		  {return mesh_count;}

	  void set_mesh_count(int iMeshCount)
		  { mesh_count=iMeshCount;}*/

	  moab::ErrorCode evaluate_smooth_edge(moab::EntityHandle eh, double &tt,
	  		moab::CartVect & outv) ;

	  //void create_mesh_edges(std::map<moab::EntityHandle, SmoothVertex*>  & mapVertices);

private:
	  //RefEdge * _ref_edge;
	  //SmoothFaceEval * _smoothFaceEval; // just store the face pointer here,
	  int _loopIndex;

	  int mesh_count;// how many mesh edges will be created on this geo edge set?

	  std::vector<moab::EntityHandle> _entities;// the mesh edges are stored here for fast access
	  moab::EntityHandle startNode, endNode;// those handles are for the nodes in _mb
	  double _leng;
	  std::vector<double> _fractions;// they are increasing from 0. to 1., do we need these?
	  // this will be decided apriori, and eventually reset for paver
	  // fractions will be from 0.0.. to 1.0, they will be decided upon the length of the geo edge

	  moab::Tag _edgeTag;

	  moab::Interface * _mb;
	  moab::EntityHandle _set;



};

#endif /* SMOOTHCURVEEVAL_HPP_ */
