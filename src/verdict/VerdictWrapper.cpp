/*
 * VerdictWrapper.cpp
 *
 *  Created on: Nov 18, 2014
 *
 */

#include "moab/Interface.hpp"
#include "moab/VerdictWrapper.hpp"
#include "Internals.hpp"
#include "moab/verdict.h"

namespace moab
{
VerdictWrapper::VerdictWrapper(Interface * mb) :mbImpl(mb) {
  // TODO Auto-generated constructor stub

}

VerdictWrapper::~VerdictWrapper() {
  // TODO Auto-generated destructor stub
}

static int possibleQuality[MBMAXTYPE][MB_QUALITY_COUNT] = {
    /*
      MB_EDGE_RATIO = 0,  // 0
      |  MB_MAX_EDGE_RATIO , // 1
      |  |  MB_SKEW,            // 2
      |  |  |  MB_TAPER,           // 3
      |  |  |  |  MB_VOLUME,          // 4
      |  |  |  |  |  MB_STRETCH,         // 5
      |  |  |  |  |  |  MB_DIAGONAL,        // 6
      |  |  |  |  |  |  |  MB_DIMENSION,       // 7
      |  |  |  |  |  |  |  |  MB_ODDY,            // 8
      |  |  |  |  |  |  |  |  |  MB_MED_ASPECT_FROBENIUS,// 9
      |  |  |  |  |  |  |  |  |  |  MB_MAX_ASPECT_FROBENIUS, // 10
      |  |  |  |  |  |  |  |  |  |  |  MB_CONDITION,       // 11
      |  |  |  |  |  |  |  |  |  |  |  |  MB_JACOBIAN,        // 12
      |  |  |  |  |  |  |  |  |  |  |  |  |  MB_SCALED_JACOBIAN, // 13
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_SHEAR,           // 14
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_SHAPE,           // 15
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_RELATIVE_SIZE_SQUARED, // 16
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_SHAPE_AND_SIZE,        // 17
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_SHEAR_AND_SIZE,        // 18
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_DISTORTION,            // 19
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_LENGTH,                // 20 only for edge
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_RADIUS_RATIO        // 21 tet
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_ASPECT_BETA      // 22 tet (very similar to 21)
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_ASPECT_RATIO   // 23 tet
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_ASPECT_GAMMA // 24 tet
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_MINIMUM_ANGLE // 25 tet
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  MB_COLLAPSE_RATIO // 26 tet
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
      |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |  |
    */
    /*0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26*/
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBVERTEX
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0},           //  MBEDGE
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBTRI
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBQUAD
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBPOLYGON
     {1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1},           //  MBTET
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBPYRAMID
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBPRISM
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBKNIFE
     {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},           // MBHEX
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},           //  MBPOLYHEDRON
     {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}            //  MBENTITYSET
};

char const * nameQuality [MB_QUALITY_COUNT]
{
   " edge ratio ",          //  MB_EDGE_RATIO = 0,  // 0
   " maximum edge ratio ",  //  MB_MAX_EDGE_RATIO , // 1
   " skew ",                //  MB_SKEW,            // 2
   " taper " ,              //   MB_TAPER,            // 3
   " volume " ,             //   MB_VOLUME,          // 4
   " stretch " ,            //   MB_STRETCH,         // 5
   " diagonal " ,           //   MB_DIAGONAL,        // 6
   " dimension ",           //  MB_DIMENSION,       // 7
   " oddy " ,               //   MB_ODDY,            // 8
   " average Frobenius aspect", //  MB_MED_ASPECT_FROBENIUS,// 9
   " maximum Frobenius aspect ", //   MB_MAX_ASPECT_FROBENIUS, // 10
   " condition number " ,   //   MB_CONDITION,       // 11
   " jacobian " ,           //   MB_JACOBIAN,        // 12
   " scaled jacobian " ,    //   MB_SCALED_JACOBIAN, // 13
   " shear " ,              //   MB_SHEAR,           // 14
   " shape " ,              //  MB_SHAPE,           // 15
   " relative size squared ",//    MB_RELATIVE_SIZE_SQUARED, // 16
   " shape and size " ,      //   MB_SHAPE_AND_SIZE,        // 17
   " shear and size " ,      //   MB_SHEAR_AND_SIZE,        // 18
   " distortion " ,          //   MB_DISTORTION,            // 19
      // next are QuadMetricVals that are not in hex metrics
      // length for edge:
   " length "  ,             //   MB_LENGTH,                // 20 only for edge
   " radius ratio ",         //   MB_RADIUS_RATIO           // 21 tet
   " aspect beta ",          // MB_ASPECT_BETA              // 22 tet
   " aspect ratio ",         // MB_ASPECT_RATIO,            // 23 MBTET
   " aspect gamma ",         // MB_ASPECT_GAMMA             // 24 tet
   " minimum angle ",        //  MB_MINIMUM_ANGLE,         // 25        MBTET
   " collapse ratio "        // MB_COLLAPSE_RATIO,        // 26        MBTET
};

ErrorCode VerdictWrapper::quality_measure(EntityHandle eh, QualityType q, double & quality)
{
  EntityType etype= TYPE_FROM_HANDLE(eh);
  if (possibleQuality[etype][q]==0)
    return MB_NOT_IMPLEMENTED;

  // get coordinates of points
  const EntityHandle * conn = NULL;
  int num_nodes;
  ErrorCode rval = mbImpl->get_connectivity(eh, conn, num_nodes);
  if (rval!=MB_SUCCESS)
    return rval;
  double coordinates[27][3]; // at most 27 nodes per element?
  if (etype!=MBPOLYHEDRON)
  {
    rval = mbImpl->get_coords(conn, num_nodes, &(coordinates[0][0]));
    if (rval!=MB_SUCCESS)
      return rval;
  }
  VerdictFunction func=0;

  if (MBHEX==etype)
  {
    switch(q) {
    case MB_EDGE_RATIO:           func = v_hex_edge_ratio; break;          // 0
    case MB_MAX_EDGE_RATIO:       func = v_hex_max_edge_ratio; break;      // 1
    case MB_SKEW:                 func = v_hex_skew; break;                // 2
    case MB_TAPER:                func = v_hex_taper; break;               // 3
    case MB_VOLUME:               func = v_hex_volume; break;              // 4
    case MB_STRETCH:              func = v_hex_stretch; break;             // 5
    case MB_DIAGONAL:             func = v_hex_diagonal; break;            // 6
    case MB_DIMENSION:            func = v_hex_dimension; break;           // 7
    case MB_ODDY:                 func = v_hex_oddy; break;                // 8
    case MB_MED_ASPECT_FROBENIUS: func = v_hex_med_aspect_frobenius; break;// 9
    case MB_MAX_ASPECT_FROBENIUS: func = v_hex_max_aspect_frobenius; break;// 10
    case MB_CONDITION:            func = v_hex_condition; break;           // 11
    case MB_JACOBIAN:             func = v_hex_jacobian; break;            // 12
    case MB_SCALED_JACOBIAN:      func = v_hex_scaled_jacobian; break;     // 13
    case MB_SHEAR:                func = v_hex_shear; break;               // 14
    case MB_SHAPE:                func = v_hex_shape; break;               // 15
    case MB_RELATIVE_SIZE_SQUARED:func = v_hex_relative_size_squared; break; // 16
    case MB_SHAPE_AND_SIZE:       func = v_hex_shape_and_size; break  ;      // 17
    case MB_SHEAR_AND_SIZE:       func = v_hex_shear_and_size; break;        // 18
    case MB_DISTORTION:           func = v_hex_distortion; break;            // 19
    default :  return MB_FAILURE;
    }
  }
  if (MBEDGE==etype)
  {
    switch (q) {
    case MB_LENGTH: func = v_edge_length; break;                             // 20
    default : return MB_FAILURE;
    }
  }
  if (MBTET == etype)
  {
    switch (q) {
    case MB_EDGE_RATIO:           func = v_tet_edge_ratio; break;           // 0 //! Calculates tet edge ratio metric.
    case MB_RADIUS_RATIO:         func = v_tet_radius_ratio; break;         // 21
    case MB_ASPECT_BETA:          func = v_tet_aspect_beta; break;          // 22
    case MB_ASPECT_RATIO:         func = v_tet_aspect_ratio; break;         // 23
    case MB_ASPECT_GAMMA:         func = v_tet_aspect_gamma; break;         // 24
    case MB_MAX_ASPECT_FROBENIUS: func = v_tet_aspect_frobenius; break;     // 10
    case MB_MINIMUM_ANGLE:        func = v_tet_minimum_angle; break;        // 25
    case MB_COLLAPSE_RATIO:       func = v_tet_collapse_ratio; break;       // 26
    case MB_VOLUME:               func = v_tet_volume; break;               // 4
    case MB_CONDITION:            func = v_tet_condition; break;            // 11
    case MB_JACOBIAN:             func = v_tet_jacobian; break;             // 12
    case MB_SCALED_JACOBIAN:      func = v_tet_scaled_jacobian; break;      // 13
    case MB_SHAPE:                func = v_tet_shape; break;                // 15
    case MB_RELATIVE_SIZE_SQUARED:func = v_tet_relative_size_squared; break;// 16
    case MB_SHAPE_AND_SIZE:       func = v_tet_shape_and_size; break;       // 17
    case MB_DISTORTION:           func = v_tet_distortion; break;           // 19
    default : return MB_FAILURE;
    }
  }
  /*if (MBTRI==etype)
  {
    switch(q) {
    case
    }
  }*/

  if (!func)
    return MB_FAILURE;
  // actual computation happens here
  quality = (*func)(num_nodes, coordinates);

  return MB_SUCCESS;
}
} // end namespace
