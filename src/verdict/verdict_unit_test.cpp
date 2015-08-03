/*=========================================================================

  Module:    $RCSfile: verdict_test.cpp,v $

  Copyright (c) 2006 Sandia Corporation.
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/


/*
 *
 * verdict_test.cpp provides routines for testing the quality metrics code
 *
 * This file is part of VERDICT
 *
 */



#define VERDICT_EXPORTS

#include "moab/Core.hpp"
#include "moab/VerdictWrapper.hpp"
// #include "v_vector.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

using namespace std;
using namespace moab;

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


int main( )
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
            { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
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
            { 1.34, 0.30, 0.20, 0.20, 0.23, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
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
            { /*0.0,*/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 }
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
            { /*0.34,*/ 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.34, 0.0, 0.0 }
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
    VerdictWrapper vw(iface);
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
            merr = vw.quality_measure(dummy, testcases[i].function[j], answer_from_lib,
                                        testcases[i].num_nodes, testcases[i].etype, testcases[i].coords);MB_CHK_ERR(merr);

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

    std::cout << endl << "All tests passed ? " << (passed ? 1 : 0) << endl ;

    return (passed ? 0 : 1);
}

