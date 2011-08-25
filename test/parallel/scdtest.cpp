#include <string>
#include <iomanip>
#include <iostream>
#include <cassert>

#include "moab/Core.hpp"
#include "moab/ParallelComm.hpp"
#include "moab/HomXform.hpp"
#include "MBParallelConventions.h"
#include "MBTagConventions.hpp"

using namespace std;
using namespace moab;

// Number of cells in each direction:
const int NC = 2;

const int NI = 2;
const int NJ = 2;
const int NK = 1;

// Number of processes:
const int NPROCS = 4;

// Domain size:
const double DSIZE = 10.0;

// MOAB objects:
Interface *mbint = NULL;
ParallelComm *mbpc = NULL;

// Local domain starting and ending hex indexes:
int is, js, ks;
int ie, je, ke;

// Obvious:
int rank;
int size;

void set_local_domain_bounds();
void create_hexes_and_verts();
void resolve_and_exchange();
void error(ErrorCode err);

int main(int argc, char *argv[])
{
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if(size != 4 && size != 2) {
    cerr << "Run this with 2 or 4 processes\n";
    exit(1);
  }
  
  mbint = new Core();
  mbpc  = new ParallelComm(mbint);

  set_local_domain_bounds();
  create_hexes_and_verts();
  resolve_and_exchange();

  MPI_Finalize();
  return 0;
}


void set_local_domain_bounds() 
{
  switch (size) {
    case 2:
        switch(rank) {
          case 0:
              is = 0; ie = NI/2;
              js = 0; je = NJ;
              ks = 0; ke = NK;
              break;


          case 1:
              is = NI/2; ie = NI;
              js = 0; je = NJ;
              ks = 0; ke = NK;
              break;
        }
        break;
        
    case 4:
        switch(rank) {
          case 0:
              is = 0; ie = NI/2;
              js = 0; je = NJ/2;
              ks = 0; ke = NK;
              break;
          case 1:
              is = NI/2; ie = NI;
              js = 0; je = NJ/2;
              ks = 0; ke = NK;
              break;
          case 2:
              is = 0; ie = NI/2;
              js = NJ/2; je = NJ;
              ks = 0; ke = NK;
              break;
          case 3:
              is = NI/2; ie = NI;
              js = NJ/2; je = NJ;
              ks = 0; ke = NK;
              break;
        }
        break;

    default:
        cerr << "Run this with 4 processes\n";
        exit(1);
  }
}


void create_hexes_and_verts()
{
  Core *mbcore = dynamic_cast<Core*>(mbint);
  HomCoord coord_min(0,0,0);
  HomCoord coord_max(ie-is, je-js, ke-ks);
  EntitySequence* vertex_seq = NULL;
  EntitySequence* cell_seq = NULL;
  EntityHandle vs, cs;

  error(mbcore->create_scd_sequence(coord_min, coord_max, MBVERTEX, 1, vs, vertex_seq));
  error(mbcore->create_scd_sequence(coord_min, coord_max, MBHEX, 1, cs, cell_seq));

  HomCoord p1(0,0,0);
  HomCoord p2(1,0,0);
  HomCoord p3(0,1,0);

  error(mbcore->add_vsequence(vertex_seq, cell_seq, p1, p1, p2, p2, p3, p3));

  // Set global id's:
  int gid;
  Tag global_id_tag;
  error(mbint->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, global_id_tag));
  EntityHandle handle = vs;
  int i,j,k;

  ErrorCode err;

  for(k = ks; k < ke + 1; k++)
    for(j = js; j < je + 1; j++)
      for(i = is; i < ie + 1; i++) {
        gid = 1 + i + j*(NI+1) + k*(NI+1)*(NJ+1);
        err = mbint->tag_set_data(global_id_tag, &handle, 1, &gid);
        if(err != MB_SUCCESS) {
          exit(1);
        }
        handle++;
      }

  handle = cs;
  for(k = ks; k < ke; k++)
    for(j = js; j < je; j++)
      for(i = is; i < ie; i++) {
        gid = 1 + i + j*NI + k*NI*NJ;
        error(mbint->tag_set_data(global_id_tag, &handle, 1, &gid));
	handle++;
      }
}


void resolve_and_exchange()
{
  EntityHandle entity_set;

  // Create the entity set:
  error(mbint->create_meshset(MESHSET_SET, entity_set));

  // Get a list of hexes:
  Range range;
  error(mbint->get_entities_by_type(0, MBHEX, range));

  // Add entities to the entity set:
  error(mbint->add_entities(entity_set, range));

  // Add the MATERIAL_SET tag to the entity set:
  Tag tag;
  error(mbint->tag_get_handle(MATERIAL_SET_TAG_NAME, 1, MB_TYPE_INTEGER, tag));
  error(mbint->tag_set_data(tag, &entity_set, 1, &rank));

  // Set up partition sets. This is where MOAB is actually told what
  // entities each process owns:
  error(mbint->get_entities_by_type_and_tag(0, MBENTITYSET,
					    &tag, NULL, 1,
					    mbpc->partition_sets()));

  // Finally, determine which entites are shared and exchange the
  // ghosted entities:
  error(mbpc->resolve_shared_ents(0, -1, -1));
  error(mbpc->exchange_ghost_cells(-1, 0, 1, 0, true));
}

void error(ErrorCode err)
{
  if(err != MB_SUCCESS) {
    cerr << "Error: MOAB function failed\n";
    assert(0);
  }
}
