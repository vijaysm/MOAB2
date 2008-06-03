#include <iostream>
#include "MBElemUtil.hpp"

using namespace std;


extern "C"{
#include "errmem.h" //for tmalloc, convenient but not C++
}

int main()
{

    MBCartVect xyz(.5,.3,.4);
    MBCartVect rst;
    double dist;

    MBElemUtil u;

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
        
    u.hex_findpt(xm, n, xyz, rst, dist);


    cout << "Coords of " << xyz << " are:  "<< rst <<
      " distance: "<< dist << endl;

}
