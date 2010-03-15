#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream> 

#include "MBCore.hpp"
#include "MBReadUtilIface.hpp"

using namespace std;

MBErrorCode ReadTriangleOutput( MBCore *mb, std::string fileBase ) {    
  
  //
  // get the read iface from moab
  void* ptr = 0;
  mb->query_interface("MBReadUtilIface", &ptr);
  MBReadUtilIface *iface = reinterpret_cast<MBReadUtilIface*>(ptr);
  //
  string nodeFileName = fileBase+".node";
  ifstream nodeFile (nodeFileName.c_str());
  if (!nodeFile.is_open())
  {
     cout<<"can't open node file .\n";
     return MB_FILE_DOES_NOT_EXIST;
  }
  
  string eleFileName = fileBase+".ele";
  ifstream eleFile (eleFileName.c_str());
  if (!eleFile.is_open())
  {
     cout<<"can't open node file .\n";
     return MB_FILE_DOES_NOT_EXIST;
  }

  string line;
  
  // ignore comment lines that start with #
  
  int num_nodes=0, num_triangles=0;
  while(num_nodes==0)
    {
      getline(nodeFile, line);
      if ('#' == line[0])
	continue;
      stringstream tks(line);
      tks >> num_nodes; // ignore the rest of the line
      cout << "num nodes:" << num_nodes << endl; 
    }
  
  //  allocate a block of vertex handles and read xyz’s into them
  vector<double*> arrays;
  MBEntityHandle startv, *starth;
  MBErrorCode rval = iface->get_node_arrays(2, num_nodes, 0, startv, arrays);
  for (int i = 0; i < num_nodes; i++)
    {
      getline(nodeFile, line);
      if ('#' == line[0])
	continue;
      stringstream tokens(line);
      int nodeId;
      tokens >> nodeId >> arrays[0][i] >> arrays[1][i] ;
    }
  
  
  while(num_triangles==0)
    {
      getline(eleFile, line);
      if ('#' == line[0])
	continue;
      stringstream tks(line);
      tks >> num_triangles; // ignore the rest of the line
      cout << "num triangles:" << num_nodes << endl; 
    }

  MBEntityHandle starte;
  // allocate block of triangle handles and read connectivity into them
  rval = iface->get_element_array(num_triangles, 3, MBTRI, 0, starte, starth);
  
  for (int j = 0; j < num_triangles; j++)
    {
      getline(eleFile, line);
      if ('#' == line[0])
	continue;
      stringstream tokens(line);
      int eleId, node;
      tokens >> eleId;
      for (int k=0; k<3; k++)
	{
	  tokens >> node;
          // vertex starts from 0
          starth[3*j+k] = (MBEntityHandle)(node + (int)startv-1 );
        }
    }

  //       
  return MB_SUCCESS;
}


// .
//  Read Triangle output files 
//  Assume that the format is <filename>.node and <filename>.ele
//   see  http://www.cs.cmu.edu/~quake/triangle.html for details
//
int main(int argc, char **argv) {
  if (3!=argc) {
    cout << "Usage: " << argv[0] << " <filename>  <outFile> " << endl;
    cout << "       <filename>  is the base file name; *.ele and *.node file are read; outFile is a file with an extension recognized by MOAB " << endl;
    return 0;
  }     

  string filename = argv[1];
  char * outfile = argv[2];
  

  // get MOAB instance and read the file                                                                                                  
  MBCore *mb = new MBCore();

   MBErrorCode rval = ReadTriangleOutput(mb, filename);

   if (rval==MB_SUCCESS)
     mb->write_file(outfile); 
   return 0;
}  
