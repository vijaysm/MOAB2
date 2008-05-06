#include <iostream>
#include <assert.h>

#include "MBElemUtil.hpp"

//
// nat_coords_trilinear_hex
// Given an MBEntityHandle defining a MBHEX element defined by a 
// trilinear basis and a set of xyz points in space, it finds the 
// cooresponding natural coordinates and returns them to the pointer
// nat.  The method uses an iterative technique, guessing initial
// coordinates {0, 0, 0} and calculating new coordinates based on the 
// ones just calculated.
// 
void MBElemUtil::nat_coords_trilinear_hex(MBCartVect hex[8], 
                                          MBCartVect xyz,
                                          MBCartVect &ncoords,
                                          double etol)       
{

      MBCartVect nat(0.);
      double A[8]; double B[8]; double C[8]; double D[8];
      double Ax, By, Cz;
      double err = 1.e37;

      double xt, yt, zt;
      double nxi, neta, nmu;
      double pxi, peta, pmu;
      double tmp;

      // Iterative estimate of natural coordinates
      while ( err > etol ) {

            // Estimate the xi-coordinate
            A[0] = (1 - nat[1]) * (1 - nat[2]);
            A[1] = A[0];
            A[2] = (1 + nat[1]) * (1 - nat[2]);        
            A[3] = A[2];    
            A[4] = (1 - nat[1]) * (1 + nat[2]);
            A[5] = A[4]; 
            A[6] = (1 + nat[1]) * (1 + nat[2]);    
            A[7] = A[6];

            Ax = 0;
            for (unsigned int i = 0; i < 8; i++) {
                  Ax = Ax + A[i]*hex[i][0];
            }
            nat[0] = (8*xyz[0] - Ax ) / 
                     (-A[0]*hex[0][0] + A[1]*hex[1][0] + A[2]*hex[2][0] - A[3]*hex[3][0] 
                      -A[4]*hex[4][0] + A[5]*hex[5][0] + A[6]*hex[6][0] - A[7]*hex[7][0]);

            // Estimate the eta-coordinate
            B[0] = (1 - nat[0]) * (1 - nat[2]);
            B[1] = (1 + nat[0]) * (1 - nat[2]);
            B[2] = B[1];
            B[3] = B[0];    
            B[4] = (1 - nat[0]) * (1 + nat[2]);
            B[5] = (1 + nat[0]) * (1 + nat[2]);
            B[6] = B[5];   
            B[7] = B[4];

            By = 0;
            for (unsigned int i = 0; i < 8; i++) {
                  By = By + B[i]*hex[i][1];
            }
            nat[1] = (8*xyz[1] - By ) / 
                     (-B[0]*hex[0][1] - B[1]*hex[1][1] + B[2]*hex[2][1] + B[3]*hex[3][1] 
                      -B[4]*hex[4][1] - B[5]*hex[5][1] + B[6]*hex[6][1] + B[7]*hex[7][1]);

            // Estimate the mu-coordinate
            C[0] = (1 - nat[0]) * (1 - nat[1]);
            C[1] = (1 + nat[0]) * (1 - nat[1]);
            C[2] = (1 + nat[0]) * (1 + nat[1]);
            C[3] = (1 - nat[0]) * (1 + nat[1]);     
            C[4] = C[0];
            C[5] = C[1];
            C[6] = C[2];   
            C[7] = C[3];

            Cz = 0;
            for (unsigned int i = 0; i < 8; i++) {
                  Cz = Cz + C[i]*hex[i][2];
            }
            nat[2] = (8*xyz[2] - Cz ) / 
                     (-C[0]*hex[0][2] - C[1]*hex[1][2] - C[2]*hex[2][2] - C[3]*hex[3][2] 
                      +C[4]*hex[4][2] + C[5]*hex[5][2] + C[6]*hex[6][2] + C[7]*hex[7][2]); 


            // Shortcut variables...
            nxi  = 1 - nat[0];
            neta = 1 - nat[1];
            nmu  = 1 - nat[2];
            pxi  = 1 + nat[0];
            peta = 1 + nat[1];
            pmu  = 1 + nat[2];
	    D[0] = nxi * neta * nmu;
	    D[1] = pxi * neta * nmu;
	    D[2] = pxi * peta * nmu;
	    D[3] = nxi * peta * nmu;
	    D[4] = nxi * neta * pmu;
	    D[5] = pxi * neta * pmu;
	    D[6] = pxi * peta * pmu;
	    D[7] = nxi * peta * pmu;

            // Compute corresponding estimates for x, y, and z to check               
            xt = 0.125 * ( D[0] * hex[0][0] + D[1] * hex[1][0] + 
                           D[2] * hex[2][0] + D[3] * hex[3][0] + 
                           D[4] * hex[4][0] + D[5] * hex[5][0] + 
                           D[6] * hex[6][0] + D[7] * hex[7][0] );

            yt = 0.125 * ( D[0] * hex[0][1] + D[1] * hex[1][1] + 
                           D[2] * hex[2][1] + D[3] * hex[3][1] + 
                           D[4] * hex[4][1] + D[5] * hex[5][1] + 
                           D[6] * hex[6][1] + D[7] * hex[7][1] );

            zt = 0.125 * ( D[0] * hex[0][2] + D[1] * hex[1][2] + 
                           D[2] * hex[2][2] + D[3] * hex[3][2] + 
                           D[4] * hex[4][2] + D[5] * hex[5][2] + 
                           D[6] * hex[6][2] + D[7] * hex[7][2] );

            // Compute error
            err = fabs(xt - xyz[0]);
            tmp = fabs(yt - xyz[1]);
            if (tmp > err) err = tmp;
            tmp = fabs(zt - xyz[2]);
            if (tmp > err) err = tmp;

      }

      ncoords = nat;
      return;

}
