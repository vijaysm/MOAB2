#include "moab/Core.hpp"
#include "moab/CartVect.hpp"
#include "moab/Range.hpp"
#include "moab/VerdictWrapper.hpp"
#include <iostream>
#include <iomanip>
#include <cstdio>
#include "TestUtil.hpp"

std::string filename = TestDir + "/mbtest1.vtk";

using namespace moab;

void verdict_test1();
void verdict_unit_tests();

int main( int argc, char* argv[] )
{
  if(argc > 1) if (argc > 1) filename = std::string(argv[1]);

  int result = 0;

  result += RUN_TEST(verdict_test1);
  result += RUN_TEST(verdict_unit_tests);

  return result;
}
void verdict_test1()
{
  ErrorCode rval;
  Core moab_core;
  Interface* mb = &moab_core;
  rval = mb->load_mesh( filename.c_str());CHECK_ERR(rval);

  Range entities;
  rval = mb->get_entities_by_handle( 0, entities ); // all entities from the model
  CHECK_ERR(rval);

  VerdictWrapper vw(mb);
  // for size methods/quality, we need a size, to compute relative sizes and stuff
  rval = vw.set_size(1.0);CHECK_ERR(rval);
  for (Range::iterator eit=entities.begin(); eit!=entities.end(); ++eit)
  {
    EntityHandle eh=*eit;
    EntityType etype=TYPE_FROM_HANDLE(eh);
    if (etype==MBVERTEX || etype>MBHEX)
      continue;
    for (int quality=0; quality<MB_QUALITY_COUNT; quality++)
    {
      QualityType q = (QualityType)quality;
      double qm;
      rval = vw.quality_measure(eh, q, qm);
      if (MB_NOT_IMPLEMENTED == rval)
        continue;
      if (MB_FAILURE == rval)
      {
        std::cerr << " failure for entity " << mb->list_entity(eh) << " quality " << vw.quality_name(q) << "\n";
      }
      if (MB_SUCCESS == rval)
      {
        std::cout << "Entity type " << (EntityType)mb->type_from_handle(eh) << " id:" << mb->id_from_handle(eh) << " quality:" <<
            vw.quality_name(q) << " : " << qm << "\n";
      }
    }
    // now compute all qualities for each entity handle
    std::map<QualityType, double> qualities;

    rval = vw.all_quality_measures(eh, qualities);
    if (MB_SUCCESS == rval)
    {
      mb->list_entity(eh);
      for (std::map<QualityType, double>::iterator mit=qualities.begin(); mit!=qualities.end(); ++mit)
      {
        std::cout << "   " << vw.quality_name(mit->first) << " " << mit->second << " \n";
      }
    }

  }

  return;
}


#define MAX_NODES_PER_ELEMENT 27
#define MAX_TESTS_PER_ELEMENT 20

#ifdef VERDICT_USE_FLOAT
#define VERDICT_SIGNIFICANT_FIG 7    // 7 significant figures for floats
#else
#define VERDICT_SIGNIFICANT_FIG 15   // 15 significant figures for doubles
#endif


struct test_case
{
    const char* testname;
    EntityType etype;
    // VerdictFunction function[MAX_TESTS_PER_ELEMENT];
    QualityType function[MAX_TESTS_PER_ELEMENT];
    int num_nodes;
    // note: the 1st dim. of coords must bigger than the maximum num_nodes
    // for any one element being tested
    // double coords[MAX_NODES_PER_ELEMENT][3];
    double coords[MAX_NODES_PER_ELEMENT*3];
    double answer[MAX_TESTS_PER_ELEMENT];
};
using namespace std;

void verdict_unit_tests()
{
// all test cases go here
  test_case testcases[] = {
      {
          "edge calc 1",
          MBEDGE,
          {MB_LENGTH, MB_UNDEFINED_QUALITY},
          2,
          { 0.0,0.0,0.0,
            1.0,1.0,1.0
          },
          { 1.732050807568877, 0.0 }
      },

      {
          "edge calc 2",
          MBEDGE,
          {MB_LENGTH, MB_UNDEFINED_QUALITY},
          2,
          { 0.0,0.0,0.0,
            1.0,0.0,0.0
          },
          { 1.0, 0.0 }
      },

      {
          "edge calc 3",
          MBEDGE,
          {MB_LENGTH, MB_UNDEFINED_QUALITY},
          2,
          { 0.0,0.0,0.0,
            0.0,0.0,0.0
          },
          { 0.0, 0.0 }
      },

      {
          "simple wedge",
          MBPRISM,
          {MB_VOLUME, MB_UNDEFINED_QUALITY},
          6,
          { 0.0,0.0,0.0,
            -1.0,1.0,0.0,
            -1.0,0.0,0.0,
            0.0,0.0,1.0,
            -1.0,1.0,1.0,
            -1.0,0.0,1.0
          },
          { 0.5, 0.0 }
      },

      {
          "singularity wedge",
          MBPRISM,
          {MB_VOLUME, MB_UNDEFINED_QUALITY},
          6,
          { 0.0,0.0,0.0,
            0.0,0.0,0.0,
            0.0,0.0,0.0,
            0.0,0.0,0.0,
            0.0,0.0,0.0,
            0.0,0.0,0.0
          },
          { 0.0, 0.0 }
      },

      {
          "simple tri",
          MBTRI,
          {   MB_AREA, MB_MINIMUM_ANGLE, MB_MAXIMUM_ANGLE, MB_CONDITION,
              MB_SCALED_JACOBIAN, MB_SHAPE, MB_RELATIVE_SIZE_SQUARED, MB_SHAPE_AND_SIZE, MB_DISTORTION, MB_UNDEFINED_QUALITY
          },
          3,
          {
              0.0,0.0,0.0,
              5,0.0,0.0,
              2.5, 4.330127, 0.0,
          },
          { 10.825317499999997,
            5.999999989158383329e+01,
            6.000000021683238316e+01,
            1.000000000000000000e+00,
            9.999999989075266660e-01,
            1.000000000000000000e+00,
            8.533333407912848181e-03,
            8.533333407912848181e-03,
            1.000000000000000000e+00,
            0.0 }
      },

      /*
      {
          "singular tri",
          MBTRI,
          { v_tri_area, v_tri_aspect_ratio, v_tri_condition,
            v_tri_distortion, v_tri_minimum_angle, v_tri_maximum_angle,
            v_tri_relative_size_squared, v_tri_shape, v_tri_shape_and_size,
            MB_UNDEFINED_QUALITY
          },
          3,
          {
            {0,0,0},
            {0.5,0.8660254037,0},
            {1,0,0}
          },
          { 123, 1234, 1234, 1234, 1234, 1234, 1234, 1234, 1234,0}
      },
      */

      /*
      {
          "simple quad",
          MBQUAD,
          { v_quad_skew, MB_UNDEFINED_QUALITY},
          4,
          {
            {0,0,0},
            {1,0,0},
            {1,7,0},
            {0,7,0 }
          },
          { 1.3333333333333333333, 0 }
      },
      */

      {
          "simple quad",
          MBQUAD,
          {   MB_ASPECT_RATIO, MB_SKEW, MB_TAPER, MB_WARPAGE, MB_AREA,
              MB_STRETCH, MB_MAXIMUM_ANGLE, MB_MINIMUM_ANGLE,
              MB_CONDITION, MB_JACOBIAN, MB_SCALED_JACOBIAN, MB_SHEAR,
              MB_SHAPE, MB_RELATIVE_SIZE_SQUARED, MB_SHAPE_AND_SIZE, MB_SHEAR_AND_SIZE,
              MB_DISTORTION, MB_UNDEFINED_QUALITY
          },
          4,
          {
              2.0,0.0,0.0,   //1
              1.0,1.0,2.0,   //2
              0.0,1.0,0.0,   //3
              0.0,0.0,0.0,   //0

          },
          { 1.429964171762364344e+00,
            9.245003270420482089e-02,
            7.453559924999298980e-01,
            7.999999999999998432e-03,
            2.692582403567252314e+00,
            5.773502691896258421e-01,
            9.000000000000000000e+01,
            5.678908923910092454e+01,
            2.307927774486215888e+00,
            1.114172029062311164e+00,
            5.570860145311555822e-01,
            5.570860145311555822e-01,
            4.332891224131210084e-01,
            9.999999999999995559e-01,
            4.332891224131208419e-01,
            5.570860145311553602e-01,
            5.626795729450876360e-01,
            0.0 }
      },

      {
          "tet test",
          MBTET,
          {   /*MB_SHEAR,*/ MB_VOLUME, MB_CONDITION, MB_JACOBIAN,
              MB_SHAPE, MB_RELATIVE_SIZE_SQUARED, MB_SHAPE_AND_SIZE, MB_DISTORTION, MB_UNDEFINED_QUALITY
          },
          4,
          {
              -5.0, -5.0, -5.0,
              -5.0, 5.0, -5.0 ,
              -5.0, -5.0, 5.0 ,
              5.0, -5.0, -5.0 ,

          },
          { 1.666666666666666572e+02,
            1.224744871391589385e+00,
            1.000000000000000000e+03,
            8.399473665965818681e-01,
            3.599999999998462557e-05,
            3.023810519746403209e-05,
            1.000000000000000000e+00,
            0.0
          }
      },

      {
          "hex test",
          MBHEX,
          {   /*MB_ASPECT_RATIO,*/ MB_SKEW, MB_TAPER, MB_VOLUME, MB_STRETCH, MB_DIAGONAL,
              MB_DIMENSION, MB_CONDITION, MB_JACOBIAN, MB_SCALED_JACOBIAN, MB_SHEAR,
              MB_SHAPE, MB_RELATIVE_SIZE_SQUARED, MB_SHEAR_AND_SIZE, MB_SHAPE_AND_SIZE,
              MB_DISTORTION, MB_UNDEFINED_QUALITY
          },
          8,
          {
              -0.2, -0.7, -0.3,  //1
              -0.7, 0.4, -0.6 ,  //2
              -0.5, 0.5, 0.3  ,  //3
              -0.3, -0.5, 0.5 ,  //0

               0.5, -0.8, -0.2,   //5
               0.4, 0.4, -0.6 ,   //6
               0.2, 0.5, 0.2  ,   //7
               0.5, -0.3, 0.8     //4
          },
          {
              2.458897037399689067e-01,
              1.784576525620624188e-01,
              8.130624999999999103e-01,
              6.209702997008308412e-01,
              6.896219298787312768e-01,
              5.245942005845132261e-01,
              1.273059825673506174e+00,
              4.769999999999999241e-01,
              7.778495101180593618e-01,
              7.778495101180593618e-01,
              7.897852810190353345e-01,
              6.738357656250002492e-01,
              5.241428201914338780e-01,
              5.321855695148176579e-01,
              5.847977114834492784e-01,
              0.0
          }
      },

      // keep this one last
      // { 0, {MB_UNDEFINED_QUALITY} , 0, {{0}} , {0} } };
      { 0, MBMAXTYPE, {MB_UNDEFINED_QUALITY}, 0, {0.0, 0.0}, {0.0, 0.0} }

  };


  int i;
  int j = 0;
  double answer_from_lib;
  double tolerance;
//   double norm_answer_from_lib;

#define MAX_STR_LEN 30

  char exponent[MAX_STR_LEN];
  char *base_ptr;
  int base;
  bool passed = true; // have all the tests performed so far passed?

  cout.setf( ios::scientific, ios::floatfield );
  cout.precision(VERDICT_SIGNIFICANT_FIG + 3 );

  ErrorCode merr;
  Interface* iface = new Core();
  VerdictWrapper* vw = new VerdictWrapper(iface);
  EntityHandle dummy=0;

  // loop through each test
  for ( i = 0; testcases[i].testname != 0; i++ )
  {
      cout << endl
           << "[" << i << "]: "
           << "Test case: " << testcases[i].testname
           << endl;

      for ( j = 0;  testcases[i].function[j] != MB_UNDEFINED_QUALITY; j++ )
      {
          /*
          answer_from_lib =
            (testcases[i].function[j])
            (testcases[i].num_nodes, testcases[i].coords);
          */

          cout << "\t #" << j+1 << " TESTING :: " << QualityType_ToString(testcases[i].function[j]) << endl;
          merr = vw->quality_measure(dummy, testcases[i].function[j], answer_from_lib,
                                      testcases[i].num_nodes, testcases[i].etype, testcases[i].coords);MB_CHK_ERR_RET(merr);

          sprintf(exponent, "%e", testcases[i].answer[j]);
          base_ptr = strstr( exponent, "e" );

          base_ptr = &base_ptr[1];

          base = atoi(base_ptr);

          tolerance = pow (10.0, - VERDICT_SIGNIFICANT_FIG ) * pow ( 10.0, base );

          if ( fabs( answer_from_lib - testcases[i].answer[j] ) > tolerance )
          {
              cout << "\t #" << j+1 << " FAILED :: " << QualityType_ToString(testcases[i].function[j]) << endl;

              cout    << "\t\t calculated ( " << answer_from_lib << " ) and "
                      << "expected ( " << testcases[i].answer[j] << ") "
                      << endl;
              passed = false;
          }
          else
          {
              cout << "\t #" << j+1 << " PASSED :: " << QualityType_ToString(testcases[i].function[j]) << endl;
          }
      }
  }
  delete vw;
  delete iface;
  std::cout << endl << "All tests passed ? " << (passed ? 1 : 0) << endl ;
  return;
}
