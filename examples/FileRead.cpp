#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream> 

#include "moab/Core.hpp"
#include "moab/ReadUtilIface.hpp"

using namespace std;
using namespace moab;

int comment(string & line)
{
    // if a line starts with '#' is a comment
    // eat white space characters
    size_t found=line.find_first_not_of(" \t");
    if (found==string::npos)
	return 1; // empty line
    if ('#'==line[found])
        return 1; // a comment indeed
    return 0; // a line with some data in it, then

}
ErrorCode ReadTriangleOutput( Interface *mb, string fileBase ) {    
  
  //
  // get the read interface from moab
  ReadUtilIface *iface;
  ErrorCode rval = mb->query_interface(iface);
  //
  if (MB_SUCCESS != rval)
     {
        cout<<"Can't get interface.\n";
        return MB_FAILURE;
     }
  // Triangle default <name>.node
  string nodeFileName = fileBase+".node";
  ifstream nodeFile (nodeFileName.c_str());
  if (!nodeFile.is_open())
  {
     cout<<"can't open node file .\n";
     return MB_FILE_DOES_NOT_EXIST;
  }
  cout << "reading nodes from file " << nodeFileName.c_str() << endl;
  
  string eleFileName = fileBase+".ele";
  ifstream eleFile (eleFileName.c_str());
  if (!eleFile.is_open())
  {
     cout<<"can't open element file .\n";
     return MB_FILE_DOES_NOT_EXIST;
  }
  cout << "reading elements from file " << eleFileName.c_str() << endl;

  string line;
  
  // ignore comment lines that start with #
  
  int num_nodes=0, num_triangles=0;
  while(num_nodes==0)
    {
      getline(nodeFile, line);
      if (comment(line))
	continue;
      stringstream tks(line);
      tks >> num_nodes; // ignore the rest of the first line
                        // maybe will read attributes some other time
      cout << "num nodes:" << num_nodes << endl; 
    }
  
  //  allocate a block of vertex handles and read xyz’s into them
  //  we know the size of the node arrays, and this call will allocate 
  //   needed arrays, coordinate arrays
  //   also, it will return a starting handle for the node sequence
  vector<double*> arrays;
  EntityHandle startv;
  rval = iface->get_node_coords(2, num_nodes, 0, startv, arrays);
  for (int i = 0; i < num_nodes; i++)
    {
      getline(nodeFile, line);
      if (comment(line))
      {
        i--;// read one more line
	continue;
      }
      stringstream tokens(line);
      int nodeId;
      tokens >> nodeId >> arrays[0][i] >> arrays[1][i] ;
    }
  
  // now read the element data from a different file
  // first, find out how many elements are out there
  // first line with data should have it
  while(num_triangles==0)
    {
      getline(eleFile, line);
      if (comment(line))
	continue;
      stringstream tks(line);
      tks >> num_triangles; // ignore the rest of the line
      cout << "num triangles:" << num_triangles << endl; 
    }

  EntityHandle starte;
  EntityHandle *starth; // the connectivity array that will get populated
                          // with triangle data
  // allocate block of triangle handles and read connectivity into them
  rval = iface->get_element_connect(num_triangles, 3, MBTRI, 0, starte, starth);
  
  for (int j = 0; j < num_triangles; j++)
    {
      getline(eleFile, line);
      if (comment(line))
      {
        j--;// read one more line
	continue;
      }
      stringstream tokens(line);
      int eleId;
      unsigned int node;
      tokens >> eleId;
      for (int k=0; k<3; k++)
	{
	  tokens >> node;
          // vertex handles start at startv
            starth[3*j+k] = startv + node - 1;
        }
    }

  mb->release_interface(iface);
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
  Core *mb = new Core();

   ErrorCode rval = ReadTriangleOutput(mb, filename);

   if (rval==MB_SUCCESS)
   {
     cout << "Writing output file " << outfile << endl;
     mb->write_file(outfile); 
   }
   return 0;
}  
