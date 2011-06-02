#include <string>
#include <iomanip>
#include <iostream>
#include <cassert>

#include <moab/Core.hpp>
#include <moab/Interface.hpp>
#include <moab/ParallelComm.hpp>
#include <moab/HomXform.hpp>
#include <MBParallelConventions.h>
#include <MBTagConventions.hpp>

using namespace std;
using namespace moab;

// Number of cells in each direction:
const int NC = 20;
const int ITERS = 50;

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

Range all_verts;

void set_local_domain_bounds();
void create_hexes_and_verts();
void resolve_and_exchange();
void error(ErrorCode err);
void tag_get_set(Tag tag);

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

  error(mbint->get_entities_by_type(0, MBVERTEX, all_verts));


// Create a tag
  Tag tag;
  error(mbint->tag_get_handle("test_tag", 1, MB_TYPE_DOUBLE, tag, MB_TAG_DENSE|MB_TAG_EXCL));

  Range empty_range;
  tag_get_set(tag);

  int i;
  for(i = 0; i < ITERS; i++) {
    std::cout << i << endl;
    mbpc->exchange_tags(tag, empty_range);
  }

  delete mbpc;
  delete mbint;

  MPI_Finalize();
  return 0;
}

void set_local_domain_bounds() 
{
  switch(rank) {
    case 0:
       
        switch (size) {
          case 2:
              is = 0; ie = NC/2;
              js = 0; je = NC;
              ks = 0; ke = NC;
              break;
       
          case 4:
              is = 0; ie = NC/2;
              js = 0; je = NC/2;
              ks = 0; ke = NC;
              break;
        }
        break;
	
    case 1:
	
        switch(size) {
          case 2:
              is = NC/2; ie = NC;
              js = 0; je = NC;
              ks = 0; ke = NC;
              break;
	
          case 4:
              is = NC/2; ie = NC;
              js = 0; je = NC/2;
              ks = 0; ke = NC;
              break;
        }
        break;
	
    case 2:
        is = 0; ie = NC/2;
        js = NC/2; je = NC;
        ks = 0; ke = NC;
        break;
	
    case 3:
        is = NC/2; ie = NC;
        js = NC/2; je = NC;
        ks = 0; ke = NC;
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
  HomCoord coord_max(NC/2, NC, NC);
  EntitySequence* vertex_seq = NULL;
  EntitySequence* cell_seq = NULL;
  EntityHandle vs, cs;
	
  error(mbcore->create_scd_sequence(coord_min, coord_max, MBVERTEX, 1, vs, vertex_seq));
  error(mbcore->create_scd_sequence(coord_min, coord_max, MBHEX, 1, cs, cell_seq));
	
  HomCoord p1(0,0,0);
  HomCoord p2(NC/2,0,0);
  HomCoord p3(0,NC/2,0);
	
  error(mbcore->add_vsequence(vertex_seq, cell_seq, p1, p1, p2, p2, p3, p3));
	
    // Set global id's:
  int gid;
  Tag global_id_tag;
  error(mbint->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, global_id_tag));
  EntityHandle handle = vs;
  int i,j,k;
	
  ErrorCode err;
	
  for(i = is; i < ie + 1; i++) 
    for(j = js; j < je + 1; j++)
      for(k = ks; k < ke + 1; k++) {   
        gid = k + j*(NC+1) + i*(NC+1)*(NC+1) + 1;
        err = mbint->tag_set_data(global_id_tag, &handle, 1, &gid);
        if(err != MB_SUCCESS) {
          exit(1);
        }
        handle++;
      }
	
  handle = cs;
  for(i = is; i < ie; i++) 
    for(j = js; j < je; j++)
      for(k = ks; k < ke; k++) {       
        gid = k + j*NC + i*NC*NC + 1;
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
	
void tag_get_set(Tag tag)
{
  Range::iterator iter;
  double data;
	
  for(iter = all_verts.begin(); iter != all_verts.end(); iter++) {
    data = 1.0;
    mbint->tag_set_data(tag, &(*iter), 1, &data);
    mbint->tag_get_data(tag, &(*iter), 1, &data);
  }
}
