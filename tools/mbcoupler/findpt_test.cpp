#include <iostream>
#include "ElemUtil.hpp"

using namespace std;
using namespace moab;
using namespace ElemUtil;

extern "C"{
#include "errmem.h" //for tmalloc, convenient but not C++
#include "types.h"
}


void test_hex_findpt()
{

    CartVect xyz(.5,.3,.4);
    CartVect rst;
    double dist;

    double *xm[3]; //element coord fields, lex ordering
    const int n=5; //number of nodes per direction (min is 2, for linear element)

    for(int d=0; d<3; d++){
      xm[d]=tmalloc(double, n*n*n);
    }

    double scale = 1./(n-1);
    int node = 0;
    //Stuff xm with sample data
    for(int k=0; k<n; k++){
      for(int j=0; j<n; j++){
        for(int i=0; i<n; i++){

          xm[0][node] = i*scale; 
          xm[1][node] = j*scale;
          xm[2][node] = k*scale;
          
          node++;
        }
      }
    }
        
    hex_findpt(xm, n, xyz, rst, dist);


    cout << "Coords of " << xyz << " are:  "<< rst <<
      " distance: "<< dist << endl;

}



void test_hex_eval()
{
    CartVect rst(.5,.3,.4);
    double value;

    const int n=2; //number of nodes per direction (min is 2, for linear element)
    double *field = tmalloc(double, n*n*n);

    double scale = 1./(n-1);
    int node = 0;

    //Stuff field with sample data
    for(int k=0; k<n; k++){
      for(int j=0; j<n; j++){
        for(int i=0; i<n; i++){
          field[node] = 100*scale*(i); 
          node++;
        }
      }
    }
        
    hex_eval(field, n, rst, value);
    cout << "Value at " << rst << "is: " << value << endl;
}



void test_nat_coords_trilinear_hex2()
{
  CartVect hex[8];
  CartVect xyz(.5,.3,.4);
  CartVect ncoords;;
  double etol;
  
  //Make our sample hex the unit cube [0,1]**3
  hex[0] = CartVect(0,0,0);
  hex[1] = CartVect(1,0,0);
  hex[2] = CartVect(1,1,0);
  hex[3] = CartVect(0,1,0);
  hex[4] = CartVect(0,0,1);
  hex[5] = CartVect(1,0,1);
  hex[6] = CartVect(1,1,1);
  hex[7] = CartVect(0,1,1);

  etol = .1 ; //ignored by nat_coords

  nat_coords_trilinear_hex2(hex, xyz, ncoords, etol);
      
  cout << "Coords of " << xyz << " are:  "<< ncoords << endl;

}



int main(){
  test_nat_coords_trilinear_hex2();
  test_hex_findpt();
  test_hex_eval();
}
