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
  MB_EDGE_RATIO = 0,        // 0  MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_MAX_EDGE_RATIO ,       // 1  MBHEX,                           MBQUAD
  MB_SKEW,                  // 2  MBHEX,                           MBQUAD
  MB_TAPER,                 // 3  MBHEX,                           MBQUAD
  MB_VOLUME,                // 4  MBHEX, MBTET, MBPRISM, MBKNIFE
  MB_STRETCH,               // 5  MBHEX,                           MBQUAD
  MB_DIAGONAL,              // 6  MBHEX,
  MB_DIMENSION,             // 7  MBHEX,
  MB_ODDY,                  // 8  MBHEX,                           MBQUAD
  MB_MED_ASPECT_FROBENIUS,  // 9  MBHEX,                           MBQUAD
  MB_MAX_ASPECT_FROBENIUS,  // 10 MBHEX, MBTET (aspect_frobenius)  MBQUAD,  MBTRI (aspect_frobenius)
  MB_CONDITION,             // 11 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_JACOBIAN,              // 12 MBHEX, MBTET,                    MBQUAD
  MB_SCALED_JACOBIAN,       // 13 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_SHEAR,                 // 14 MBHEX,                           MBQUAD,  MBTRI
  MB_SHAPE,                 // 15 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_RELATIVE_SIZE_SQUARED, // 16 MBHEX, MBTET,                    MBQUAD,  MBTRI
  MB_SHAPE_AND_SIZE,        // 17 MBHEX, MBTET,                    MBQUAD
  MB_SHEAR_AND_SIZE,        // 18 MBHEX,                           MBQUAD
  MB_DISTORTION,            // 19 MBHEX, MBTET,                    MBQUAD
  // length for edge:
  MB_LENGTH,                // 20 only for MBEDGE
  // specific to tets
  MB_RADIUS_RATIO,          // 21        MBTET,                    MBQUAD,  MBTRI
  MB_ASPECT_BETA,           // 22        MBTET
  MB_ASPECT_RATIO,          // 23        MBTET,                    MBQUAD,  MBTRI
  MB_ASPECT_GAMMA,          // 24        MBTET
  MB_MINIMUM_ANGLE,         // 25        MBTET,                    MBQUAD,  MBTRI
  MB_COLLAPSE_RATIO,        // 26        MBTET
  // specific to quads
  MB_WARPAGE,               // 27                                  MBQUAD
  MB_AREA,                  // 28                                  MBQUAD,  MBTRI
  MB_MAXIMUM_ANGLE,         // 29                                  MBQUAD,  MBTRI
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
