/*
 * VerdictWrapper.hpp
 *
 *  Created on: Nov 18, 2014
 *      Author: iulian
 */

#ifndef SRC_VERDICT_MOAB_VERDICTWRAPPER_HPP_
#define SRC_VERDICT_MOAB_VERDICTWRAPPER_HPP_


namespace moab
{

class Interface;

enum QualityType {
  // order exactly from HexMetricVals
  MB_UNDEFINED_QUALITY = -1,
  MB_EDGE_RATIO = 0,        // 0  MBHEX, MBTET
  MB_MAX_EDGE_RATIO ,       // 1  MBHEX,
  MB_SKEW,                  // 2  MBHEX,
  MB_TAPER,                 // 3  MBHEX,
  MB_VOLUME,                // 4  MBHEX, MBTET, MBPRISM, MBKNIFE
  MB_STRETCH,               // 5  MBHEX,
  MB_DIAGONAL,              // 6  MBHEX,
  MB_DIMENSION,             // 7  MBHEX,
  MB_ODDY,                  // 8  MBHEX,
  MB_MED_ASPECT_FROBENIUS,  // 9  MBHEX,
  MB_MAX_ASPECT_FROBENIUS,  // 10 MBHEX, MBTET (aspect_frobenius)
  MB_CONDITION,             // 11 MBHEX, MBTET
  MB_JACOBIAN,              // 12 MBHEX, MBTET
  MB_SCALED_JACOBIAN,       // 13 MBHEX, MBTET
  MB_SHEAR,                 // 14 MBHEX,
  MB_SHAPE,                 // 15 MBHEX, MBTET
  MB_RELATIVE_SIZE_SQUARED, // 16 MBHEX, MBTET
  MB_SHAPE_AND_SIZE,        // 17 MBHEX, MBTET
  MB_SHEAR_AND_SIZE,        // 18 MBHEX,
  MB_DISTORTION,            // 19 MBHEX, MBTET
  // next are QuadMetricVals that are not in hex metrics
  // length for edge:
  MB_LENGTH,                // 20 only for MBEDGE
  MB_RADIUS_RATIO,          // 21        MBTET
  MB_ASPECT_BETA,           // 22        MBTET
  MB_ASPECT_RATIO,          // 23        MBTET
  MB_ASPECT_GAMMA,          // 24        MBTET
  MB_MINIMUM_ANGLE,         // 25        MBTET
  MB_COLLAPSE_RATIO,        // 26        MBTET
  MB_QUALITY_COUNT // used to size the arrays

};

extern char const * nameQuality [MB_QUALITY_COUNT];

class VerdictWrapper {
public:
  VerdictWrapper(Interface * mb);
  virtual ~VerdictWrapper();
  ErrorCode quality_measure(EntityHandle eh, QualityType q, double & quality);
private:
  Interface * mbImpl;

};
} // namespace moab
#endif /* SRC_VERDICT_MOAB_VERDICTWRAPPER_HPP_ */
