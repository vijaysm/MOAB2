#ifdef WIN32
#pragma warning (disable : 4786)
#endif 

#include <iostream>
#include <assert.h>
#include <vector>


namespace moab{

#define MAX_DEGREE 3
#define MAX_VERTS 100
#define MAX_CHILDRENS 100
#define MAX_LOC 19
#define MAX_CONN 8

class RefTemplates{
public:
  RefTemplates():{}

  ~RefTemplates():{}

  struct refPatterns{
    short int nv_edge; // Number of new vertices on edge
    short int nv_face; // Number of new vertices on face
    short int nv_cell; // Number of new vertices in cell
    short int total_new_verts; // Total number of new vertices per entity
    short int total_new_ents; // Total number of new entities

    int vert_index_bnds[2]; //Lower and upper indices of the new vertices
    double vert_nat_coord[MAX_VERTS][2]; //Natural coordinates of the new vertices

    int ents_conn[MAX_CHILDRENS][MAX_CONN]; //Connectivity of the new entities
    int ents_opphfs[MAX_CHILDRENS][2*MAX_CONN]; // Opposite half-facet map of the new entities
  };

  static const refPatterns refTemplates[2][MAX_DEGREE];
};

}//namesapce moab
  
