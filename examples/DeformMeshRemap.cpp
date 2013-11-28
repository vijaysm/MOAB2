/** @example DeformMeshRemap.cpp
 * Description: Account for mesh deformation of a solid due to structural mechanics\n
 * In this example there are two meshes, a "master" and "slave" mesh.  In the master mesh,
 * the solid material is deformed, to mimic what happens when a solid heats up and deforms.
 * The fluid mesh is smoothed to account for those deformations, and tags on the fluid are
 * remapped to those new positions.  Then mesh positions and state variables are transferred
 * to the slave mesh, mimicing another mesh used by some other physics.
 *
 * To run: ./DeformMeshRemap [<master_meshfile> <slave_meshfile>]\n
 * (default values can run if users don't specify the mesh files)
 */

#include "moab/Core.hpp"
#include "moab/Range.hpp"
#include "moab/LloydSmoother.hpp"
#include "moab/ProgOptions.hpp"
#include "MBTagConventions.hpp"

#include <iostream>
#include <assert.h>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

ErrorCode read_file(string &fname, EntityHandle &seth, 
                    Range &solids, Range &solid_elems, Range &fluids, Range &fluid_elems);
void deform_func(double *xold, double *xnew);
ErrorCode deform_master(Range &fluid_elems, Range &solid_elems, Tag &xnew);
ErrorCode smooth_master(int dim, Tag xnew, EntityHandle &master, Range &fluids);
ErrorCode write_to_coords(Range &elems, Tag tagh);

const int SOLID_SETNO = 100, FLUID_SETNO = 200;

Interface *mb;
#define RR(a) if (MB_SUCCESS != rval) {cout << a << endl; return MB_FAILURE;}
    

int main(int argc, char **argv) {

  EntityHandle master, slave;
  ErrorCode rval;

  ProgOptions po("Deformed mesh options");
  po.addOpt<std::string> ("master,m", "Specify the master meshfile name" );
  po.addOpt<std::string> ("slave,s", "Specify the slave meshfile name" );
  po.parseCommandLine(argc, argv);
  std::string foo;
  string masterf, slavef;
  if(!po.getOpt("master", &masterf))
    masterf = string(MESH_DIR) + string("/rodquad.g");
  if(!po.getOpt("slave", &slavef))
    slavef = string(MESH_DIR) + string("/rodtri.g");

  mb = new Core();
  
    // read master/slave files and get fluid/solid material sets
  Range fluids[2], solids[2], solid_elems[2], fluid_elems[2];
  rval = read_file(masterf, master, solids[0], solid_elems[0], fluids[0], fluid_elems[0]); RR("");
  rval = read_file(slavef, slave, solids[1], solid_elems[1], fluids[1], fluid_elems[1]); RR("");

    // deform the master's solid mesh, put results in a new tag
  Tag xnew;
  rval = deform_master(fluid_elems[0], solid_elems[0], xnew); RR("");
  if (debug) write_and_save(solid_elems[0], master, xnew, "deformed.vtk");
  
    // smooth the master mesh
  LloydSmoother *ll = new LloydSmoother(mb, NULL, fluid_elems[0], xnew);
  rval = ll->perform_smooth();
  RR("Failed in lloyd smoothing.");
  cout << "Lloyd smoothing required " << ll->num_its() << " iterations." << endl;
  if (debug) write_and_save(fluid_elems[0], master, xnew, "smoothed.vtk");

    // map new locations to slave
  
  delete ll;
  delete mb;
  
  return MB_SUCCESS;
}

ErrorCode write_and_save(Range &ents, EntithHanlde seth, Tag tagh, const char *filename) 
{
  rval = write_to_coords(ents, tagh); RR("");
  rval = mb->write_file("deformed.vtk", NULL, NULL, &seth, 1); RR("");
  return rval;
}
  
ErrorCode write_to_coords(Range &elems, Tag tagh) 
{
    // write the tag to coordinates
  Range verts;
  ErrorCode rval = mb->get_adjacencies(elems, 0, false, verts, Interface::UNION);
  RR("Failed to get adj vertices.");
  std::vector<double> coords(3*verts.size());
  rval = mb->tag_get_data(tagh, verts, &coords[0]);
  RR("Failed to get tag data.");
  rval = mb->set_coords(verts, &coords[0]);
  RR("Failed to set coordinates.");
  return MB_SUCCESS;
}

void deform_func(double *xold, double *xnew) 
{
  const double RODWIDTH = 0.2, RODHEIGHT = 0.5;
    // function: origin is at middle base of rod, and is .5 high
    // top of rod is (0,.55) on left and (.2,.6) on right
  double delx = 0.5*RODWIDTH;
  
  double xfrac = (xold[0] + .5*RODWIDTH)/RODWIDTH, yfrac = xold[1]/RODHEIGHT;
  xnew[0] = xold[0] + yfrac * delx;
  xnew[1] = xold[1] + yfrac * (1.0 + xfrac) * 0.05;
}
  
ErrorCode deform_master(Range &fluid_elems, Range &solid_elems, Tag &xnew) 
{
    // deform elements with an analytic function

    // create the tag
  ErrorCode rval = mb->tag_get_handle("", 3, MB_TYPE_DOUBLE, xnew, MB_TAG_CREAT|MB_TAG_DENSE);
  RR("Failed to create xnew tag.");
  
    // get all the vertices and coords in the fluid, set xnew to them
  Range verts;
  rval = mb->get_adjacencies(fluid_elems, 0, false, verts, Interface::UNION);
  RR("Failed to get vertices.");
  std::vector<double> coords(3*verts.size(), 0.0);
  rval = mb->get_coords(verts, &coords[0]);
  RR("Failed to get vertex coords.");
  rval = mb->tag_set_data(xnew, verts, &coords[0]);
  RR("Failed to set xnew tag on fluid verts.");
  
    // get all the vertices and coords in the solid
  verts.clear();
  rval = mb->get_adjacencies(solid_elems, 0, false, verts, Interface::UNION);
  RR("Failed to get vertices.");
  coords.resize(3*verts.size(), 0.0);
  rval = mb->get_coords(verts, &coords[0]);
  RR("Failed to get vertex coords.");
  unsigned int num_verts = verts.size();
  for (unsigned int i = 0; i < num_verts; i++)
    deform_func(&coords[3*i], &coords[3*i]);
    
    // set the new tag to those coords
  rval = mb->tag_set_data(xnew, verts, &coords[0]);
  RR("Failed to set tag data.");
  
  return MB_SUCCESS;
}

ErrorCode read_file(string &fname, EntityHandle &seth, 
                    Range &solids, Range &solid_elems, Range &fluids, Range &fluid_elems)
{
    // create meshset
  ErrorCode rval = mb->create_meshset(0, seth);
  RR("Couldn't create master/slave set.");
  rval = mb->load_file(fname.c_str(), &seth);
  RR("Couldn't load master/slave mesh.");

    // get material sets for solid/fluid
  Tag tagh;
  rval = mb->tag_get_handle(MATERIAL_SET_TAG_NAME, tagh); RR("Couldn't get material set tag name.");
  const void *setno_ptr = &SOLID_SETNO;
  rval = mb->get_entities_by_type_and_tag(seth, MBENTITYSET, &tagh, &setno_ptr, 1, solids);
  if (solids.empty()) rval = MB_FAILURE;
  RR("Couldn't get any solid sets.");

    // get solid entities, and dimension
  Range tmp_range;
  for (Range::iterator rit = solids.begin(); rit != solids.end(); rit++) {
    rval = mb->get_entities_by_handle(*rit, tmp_range, true);
    RR("Failed to get entities in solid.");
  }
  int dim = mb->dimension_from_handle(*tmp_range.rbegin());
  assert(dim > 0 && dim < 4);
  
  solid_elems = tmp_range.subset_by_dimension(dim);

  setno_ptr = &FLUID_SETNO;
  rval = mb->get_entities_by_type_and_tag(seth, MBENTITYSET, &tagh, &setno_ptr, 1, fluids);
  if (fluids.empty()) rval = MB_FAILURE;
  RR("Couldn't get any fluid sets.");
  
  for (Range::iterator rit = fluids.begin(); rit != fluids.end(); rit++) {
    rval = mb->get_entities_by_dimension(*rit, dim, fluid_elems, true);
    RR("Failed to get entities in fluid.");
  }
  if (mb->dimension_from_handle(*fluid_elems.begin()) != dim) {
    rval = MB_FAILURE;
    RR("Fluid and solid elements must be same dimension.");
  }
  
  return MB_SUCCESS;
}
