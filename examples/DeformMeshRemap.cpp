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
#include "moab/Skinner.hpp"
#include "moab/LloydSmoother.hpp"
#include "moab/ProgOptions.hpp"
#include "moab/BoundBox.hpp"
#include "moab/SpatialLocator.hpp"
#include "MBTagConventions.hpp"
#include "DataCoupler.hpp"

#ifdef USE_MPI
#include "moab/ParallelComm.hpp"
#endif

#include <iostream>
#include <set>
#include <sstream>
#include <assert.h>

using namespace moab;
using namespace std;

#ifndef MESH_DIR
#define MESH_DIR "."
#endif

ErrorCode read_file(string &fname, EntityHandle &seth, 
                    Range &solids, Range &solid_elems, Range &fluids, Range &fluid_elems);
void deform_func(const BoundBox &bbox, double *xold, double *xnew);
ErrorCode deform_master(Range &fluid_elems, Range &solid_elems, Tag &xnew);
ErrorCode smooth_master(int dim, Tag xnew, EntityHandle &master, Range &fluids);
ErrorCode write_to_coords(Range &elems, Tag tagh);

const int SOLID_SETNO = 100, FLUID_SETNO = 200;

Interface *mb;

const bool debug = true;

class DeformMeshRemap 
{
public:
  //! Enumerator for solid/fluid, master/slave
  enum {MASTER = 0, SLAVE, SOLID, FLUID};

  //! Constructor
  //! If master is NULL, the MOAB part is run in serial;
  //! If slave is NULL but the master isn't, the slave is copied from the master
  //! Create communicators using moab::ParallelComm::get_pcomm
  DeformMeshRemap(Interface *impl, ParallelComm *master = NULL, ParallelComm *slave = NULL);

  //! Destructor
  ~DeformMeshRemap();

  //! Execute the deformed mesh process
  ErrorCode execute();

  //! Add a set number
  ErrorCode add_set_no(int m_or_s, int fluid_or_solid, int set_no);

  //! Remove a set number
  ErrorCode remove_set_no(int m_or_s, int fluid_or_solid, int set_no);

  //! Get the set numbers
  ErrorCode get_set_nos(int m_or_s, int fluid_or_solid, set<int> &set_nos) const;

  //! Get the xNew tag handle
  inline Tag x_new() const { return xNew; }

  //! Get the tag name
  string x_new_name() const { return xNewName; }

  //! Set the tag name
  void x_new_name(const string &name) { xNewName = name; }

  //! Get/set the file name
  string get_file_name(int m_or_s) const;

  //! Get/set the file name
  void set_file_name(int m_or_s, const string &name);

  //! Get/set the x displacement tag names
  string xdisp_name(int idx = 0);
  void xdisp_name(const string &nm, int idx = 0);

private:
  //! Apply a known deformation to the solid elements, putting the results in the xNew tag; also
  //! write current coordinates to the xNew tag for fluid elements
  ErrorCode deform_master(Range &fluid_elems, Range &solid_elems, const char *tag_name = NULL);

  //! Read a file and establish proper ranges
  ErrorCode read_file(int m_or_s, string &fname, EntityHandle &seth);

  //! Write the input tag to the coordinates for the vertices in the input elems
  //! If non-zero tmp_tag is input, save coords to tmp_tag before over-writing with tag value
  ErrorCode write_to_coords(Range &elems, Tag tagh, Tag tmp_tag = 0);

  //! Write the tag to the vertices, then save to the specified file
  //! If restore_coords is true, coords are restored to their initial state after file is written
  ErrorCode write_and_save(Range &ents, EntityHandle seth, Tag tagh, const char *filename,
                           bool restore_coords = false);

  //! Find fluid/solid sets from complement of solid/fluid sets
  ErrorCode find_other_sets(int m_or_s, EntityHandle file_set);

  //! moab interface
  Interface *mbImpl;

#ifdef USE_MPI
  //! ParallelComm for master, slave meshes
  ParallelComm *pcMaster, *pcSlave;
#endif

  //! Material set numbers for fluid materials, for master/slave
  set<int> fluidSetNos[2];

  //! Material set numbers for solid materials, for master/slave
  set<int> solidSetNos[2];

  //! Sets defining master/slave meshes
  EntityHandle masterSet, slaveSet;

  //! Sets in master/slave meshes
  Range fluidSets[2], solidSets[2];

  //! Elements in master/slave meshes
  Range fluidElems[2], solidElems[2];

  //! Filenames for master/slave meshes
  string masterFileName, slaveFileName;

  //! Tag from file, might be 3
  Tag xDisp[3];

  //! Tag used for new positions
  Tag xNew;

  //! Tag name used to read disps from file
  string xDispNames[3];

  //! Tag name used for new positions
  string xNewName;
};

//! Add a set number
inline ErrorCode DeformMeshRemap::add_set_no(int m_or_s, int f_or_s, int set_no) 
{
  set<int> *this_set;
  assert((m_or_s == MASTER || m_or_s == SLAVE) && "m_or_s should be MASTER or SLAVE.");
  if (m_or_s != MASTER && m_or_s != SLAVE) return MB_INDEX_OUT_OF_RANGE;

  switch (f_or_s) {
    case FLUID:
      this_set = &fluidSetNos[m_or_s]; break;
    case SOLID:
      this_set = &solidSetNos[m_or_s]; break;
    default:
      assert(false && "f_or_s should be FLUID or SOLID.");
      return MB_FAILURE;
  }

  this_set->insert(set_no);

  return MB_SUCCESS;
}
  
//! Remove a set number
inline ErrorCode DeformMeshRemap::remove_set_no(int m_or_s, int f_or_s, int set_no) 
{
  set<int> *this_set;
  assert((m_or_s == MASTER || m_or_s == SLAVE) && "m_or_s should be MASTER or SLAVE.");
  if (m_or_s != MASTER && m_or_s != SLAVE) return MB_INDEX_OUT_OF_RANGE;
  switch (f_or_s) {
    case FLUID:
      this_set = &fluidSetNos[m_or_s]; break;
    case SOLID:
      this_set = &solidSetNos[m_or_s]; break;
    default:
      assert(false && "f_or_s should be FLUID or SOLID.");
      return MB_FAILURE;
  }
  set<int>::iterator sit = this_set->find(set_no);
  if (sit != this_set->end()) {
    this_set->erase(*sit);
    return MB_SUCCESS;
  }

  return MB_FAILURE;
}
  
//! Get the set numbers
inline ErrorCode DeformMeshRemap::get_set_nos(int m_or_s, int f_or_s, set<int> &set_nos) const
{
  const set<int> *this_set;
  assert((m_or_s == MASTER || m_or_s == SLAVE) && "m_or_s should be MASTER or SLAVE.");
  if (m_or_s != MASTER && m_or_s != SLAVE) return MB_INDEX_OUT_OF_RANGE;
  switch (f_or_s) {
    case FLUID:
      this_set = &fluidSetNos[m_or_s]; break;
    case SOLID:
      this_set = &solidSetNos[m_or_s]; break;
    default:
      assert(false && "f_or_s should be FLUID or SOLID.");
      return MB_FAILURE;
  }

  set_nos = *this_set;

  return MB_SUCCESS;
}

inline string DeformMeshRemap::xdisp_name(int idx)
{
  return xDispNames[idx];
}

void DeformMeshRemap::xdisp_name(const string &nm, int idx)
{
  xDispNames[idx] = nm;
}

ErrorCode DeformMeshRemap::execute() 
{
  // Read master/slave files and get fluid/solid material sets
  ErrorCode rval = read_file(MASTER, masterFileName, masterSet);CHK_ERR(rval);

  if (solidSetNos[MASTER].empty() || fluidSetNos[MASTER].empty()) {
    rval = find_other_sets(MASTER, masterSet);CHK_SET_ERR(rval, "Failed to find other sets in master mesh");
  }

  bool have_slave = !(slaveFileName == "none");
  if (have_slave) {
    rval = read_file(SLAVE, slaveFileName, slaveSet);CHK_ERR(rval);

    if (solidSetNos[SLAVE].empty() || fluidSetNos[SLAVE].empty()) {
      rval = find_other_sets(SLAVE, slaveSet);CHK_SET_ERR(rval, "Failed to find other sets in slave mesh");
    }
  }

  if (debug) cout << "Constructing data coupler/search tree on master mesh..." << endl;

  Range src_elems = solidElems[MASTER];
  src_elems.merge(fluidElems[MASTER]);

  // Initialize data coupler on source elements
  DataCoupler dc_master(mbImpl, src_elems, 0, NULL);

  Range tgt_verts;
  if (have_slave) {
    // Locate slave vertices in master, orig coords; do this with a data coupler, so you can
    // later interpolate
    Range tmp_range = solidElems[SLAVE];
    tmp_range.merge(fluidElems[SLAVE]);
    rval = mbImpl->get_adjacencies(tmp_range, 0, false, tgt_verts, Interface::UNION);CHK_SET_ERR(rval, "Failed to get target verts");

    // Locate slave vertices, caching results in dc
    if (debug) cout << "Locating slave vertices in master mesh..." << endl;
    rval = dc_master.locate_points(tgt_verts);CHK_SET_ERR(rval, "Point location of tgt verts failed");
    int num_located = dc_master.spatial_locator()->local_num_located();
    if (num_located != (int)tgt_verts.size()) {
      rval = MB_FAILURE;
      cout << "Only " << num_located << " out of " << tgt_verts.size() << " target points successfully located." << endl;
      return rval;
    }
  }

  // Deform the master's solid mesh, put results in a new tag
  if (debug) cout << "Deforming fluid elements in master mesh..." << endl;
  rval = deform_master(fluidElems[MASTER], solidElems[MASTER], "xnew");CHK_ERR(rval);

  { // To isolate the lloyd smoother & delete when done
    if (debug) {
      // Output the skin of smoothed elems, as a check
      // Get the skin; get facets, because we might need to filter on shared entities
      Skinner skinner(mbImpl);
      Range skin;
      rval = skinner.find_skin(0, fluidElems[MASTER], false, skin);CHK_SET_ERR(rval, "Unable to find skin");
      EntityHandle skin_set;
      cout << "Writing skin_mesh.g and fluid_mesh.g." << endl;
      rval = mbImpl->create_meshset(MESHSET_SET, skin_set);CHK_SET_ERR(rval, "Failed to create skin set");
      rval = mbImpl->add_entities(skin_set, skin);CHK_SET_ERR(rval, "Failed to add skin entities to set");
      rval = mbImpl->write_file("skin_mesh.vtk", NULL, NULL, &skin_set, 1);CHK_SET_ERR(rval, "Failure to write skin set");
      rval = mbImpl->remove_entities(skin_set, skin);CHK_SET_ERR(rval, "Failed to remove skin entities from set");
      rval = mbImpl->add_entities(skin_set, fluidElems[MASTER]);CHK_SET_ERR(rval, "Failed to add fluid entities to set");
      rval = mbImpl->write_file("fluid_mesh.vtk", NULL, NULL, &skin_set, 1);CHK_SET_ERR(rval, "Failure to write fluid set");
      rval = mbImpl->delete_entities(&skin_set, 1);CHK_SET_ERR(rval, "Failed to delete skin set");
    }

    // Smooth the master mesh
    if (debug) cout << "Smoothing fluid elements in master mesh..." << endl;
    LloydSmoother ll(mbImpl, NULL, fluidElems[MASTER], xNew);
    rval = ll.perform_smooth();CHK_SET_ERR(rval, "Failed in lloyd smoothing");
    cout << "Lloyd smoothing required " << ll.num_its() << " iterations." << endl;
  }

  // Transfer xNew to coords, for master
  if (debug) cout << "Transferring coords tag to vertex coordinates in master mesh..." << endl;
  rval = write_to_coords(solidElems[MASTER], xNew);CHK_SET_ERR(rval, "Failed writing tag to master fluid verts");
  rval = write_to_coords(fluidElems[MASTER], xNew);CHK_SET_ERR(rval, "Failed writing tag to master fluid verts");

  if (have_slave) {
    // Map new locations to slave
    // Interpolate xNew to slave points
    if (debug) cout << "Interpolating new coordinates to slave vertices..." << endl;
    rval = dc_master.interpolate((int)DataCoupler::VOLUME, "xnew");CHK_SET_ERR(rval, "Failed to interpolate target solution");
    // Transfer xNew to coords, for slave
    if (debug) cout << "Transferring coords tag to vertex coordinates in slave mesh..." << endl;
    rval = write_to_coords(tgt_verts, xNew);CHK_SET_ERR(rval, "Failed writing tag to slave verts");
  }

  if (debug) {
    string str;
#ifdef USE_MPI
    if (pcMaster && pcMaster->size() > 1) 
      str = "PARALLEL=WRITE_PART";
#endif
    if (debug) cout << "Writing smoothed_master.h5m..." << endl;
    rval = mbImpl->write_file("smoothed_master.h5m", NULL, str.c_str(), &masterSet, 1);

    if (have_slave) {
#ifdef USE_MPI
      str.clear();
      if (pcSlave && pcSlave->size() > 1) 
        str = "PARALLEL=WRITE_PART";
#endif
      if (debug) cout << "Writing slave_interp.h5m..." << endl;
      rval = mbImpl->write_file("slave_interp.h5m", NULL, str.c_str(), &slaveSet, 1);
    } // if have_slave
  } // if debug

  if (debug) 
    dc_master.spatial_locator()->get_tree()->tree_stats().print();

  return MB_SUCCESS;
}

string DeformMeshRemap::get_file_name(int m_or_s) const
{
  switch (m_or_s) {
    case MASTER:
      return masterFileName;
    case SLAVE:
      return slaveFileName;
    default:
      assert(false && "m_or_s should be MASTER or SLAVE.");
      return string();
  }
}

void DeformMeshRemap::set_file_name(int m_or_s, const string &name)
{
  switch (m_or_s) {
    case MASTER:
      masterFileName = name; break;
    case SLAVE:
      slaveFileName = name; break;
    default:
      assert(false && "m_or_s should be MASTER or SLAVE.");
  }
}

DeformMeshRemap::DeformMeshRemap(Interface *impl, ParallelComm *master, ParallelComm *slave)  
        : mbImpl(impl), pcMaster(master), pcSlave(slave), masterSet(0), slaveSet(0), xNew(0), xNewName("xnew") 
{
  xDisp[0] = xDisp[1] = xDisp[2] = 0;

  if (!pcSlave && pcMaster)
    pcSlave = pcMaster;
}

DeformMeshRemap::~DeformMeshRemap() 
{
}

int main(int argc, char **argv)
{
  ErrorCode rval;

  ProgOptions po("Deformed mesh options");
  po.addOpt<string> ("master,m", "Specify the master meshfile name" );
  po.addOpt<string> ("worker,w", "Specify the slave/worker meshfile name, or 'none' (no quotes) if master only" );
  po.addOpt<string> ("d1,", "Tag name for displacement x or xyz" );
  po.addOpt<string> ("d2,", "Tag name for displacement y" );
  po.addOpt<string> ("d3,", "Tag name for displacement z" );
  po.addOpt<int> ("fm,", "Specify master fluid material set number(s). If none specified, fluid sets derived from complement of solid sets.");
  po.addOpt<int> ("fs,", "Specify master solid material set number(s). If none specified, solid sets derived from complement of fluid sets.");
  po.addOpt<int> ("sm,", "Specify slave fluid material set number(s). If none specified, fluid sets derived from complement of solid sets.");
  po.addOpt<int> ("ss,", "Specify slave solid material set number(s). If none specified, solid sets derived from complement of fluid sets.");

  po.parseCommandLine(argc, argv);

  mb = new Core();

  DeformMeshRemap *dfr;
#ifdef USE_MPI
  ParallelComm *pc = new ParallelComm(mb, MPI_COMM_WORLD);
  dfr = new DeformMeshRemap(mb, pc);
#else  
  dfr = new DeformMeshRemap(mb);
#endif

  string masterf, slavef;
  if (!po.getOpt("master", &masterf))
    masterf = string(MESH_DIR) + string("/rodquad.g");
  dfr->set_file_name(DeformMeshRemap::MASTER, masterf);

  if (!po.getOpt("worker", &slavef))
    slavef = string(MESH_DIR) + string("/rodtri.g");
  dfr->set_file_name(DeformMeshRemap::SLAVE, slavef);
  if (slavef.empty()) {
    cerr << "Empty slave file name; if no slave, use filename 'none' (no quotes)." << endl;
    return 1;
  }

  vector<int> set_nos;
  po.getOptAllArgs("fm", set_nos);
  for (vector<int>::iterator vit = set_nos.begin(); vit != set_nos.end(); ++vit)
    dfr->add_set_no(DeformMeshRemap::MASTER, DeformMeshRemap::FLUID, *vit);
  set_nos.clear();

  po.getOptAllArgs("fs", set_nos);
  for (vector<int>::iterator vit = set_nos.begin(); vit != set_nos.end(); ++vit)
    dfr->add_set_no(DeformMeshRemap::SLAVE, DeformMeshRemap::FLUID, *vit);
  set_nos.clear();

  po.getOptAllArgs("sm", set_nos);
  for (vector<int>::iterator vit = set_nos.begin(); vit != set_nos.end(); ++vit)
    dfr->add_set_no(DeformMeshRemap::MASTER, DeformMeshRemap::SOLID, *vit);

  po.getOptAllArgs("ss", set_nos);
  for (vector<int>::iterator vit = set_nos.begin(); vit != set_nos.end(); ++vit)
    dfr->add_set_no(DeformMeshRemap::SLAVE, DeformMeshRemap::SOLID, *vit);

  string tnames[3];
  po.getOpt("d1", &tnames[0]);
  po.getOpt("d2", &tnames[1]);
  po.getOpt("d3", &tnames[2]);
  for (int i = 0; i < 3; i++) 
    if (!tnames[i].empty()) dfr->xdisp_name(tnames[i], i);

  rval = dfr->execute();

  delete dfr;
  delete mb;

  return rval;
}

ErrorCode DeformMeshRemap::write_and_save(Range &ents, EntityHandle seth, Tag tagh, const char *filename,
                                          bool restore_coords) 
{
  Tag tmp_tag = 0;
  ErrorCode rval;
  if (restore_coords)
    rval = mbImpl->tag_get_handle("", 3, MB_TYPE_DOUBLE, tmp_tag, MB_TAG_CREAT|MB_TAG_DENSE);

  rval = write_to_coords(ents, tagh, tmp_tag);CHK_ERR(rval);
  rval = mbImpl->write_file(filename, NULL, NULL, &seth, 1);CHK_ERR(rval);
  if (restore_coords) {
    rval = write_to_coords(ents, tmp_tag);CHK_ERR(rval);
    rval = mbImpl->tag_delete(tmp_tag);CHK_ERR(rval);
  }

  return rval;
}
  
ErrorCode DeformMeshRemap::write_to_coords(Range &elems, Tag tagh, Tag tmp_tag) 
{
  // Write the tag to coordinates
  Range verts;
  ErrorCode rval = mbImpl->get_adjacencies(elems, 0, false, verts, Interface::UNION);CHK_SET_ERR(rval, "Failed to get adj vertices");
  vector<double> coords(3*verts.size());

  if (tmp_tag) {
    // Save the coords to tmp_tag first
    rval = mbImpl->get_coords(verts, &coords[0]);CHK_SET_ERR(rval, "Failed to get tmp copy of coords");
    rval = mbImpl->tag_set_data(tmp_tag, verts, &coords[0]);CHK_SET_ERR(rval, "Failed to save tmp copy of coords");
  }

  rval = mbImpl->tag_get_data(tagh, verts, &coords[0]);CHK_SET_ERR(rval, "Failed to get tag data");
  rval = mbImpl->set_coords(verts, &coords[0]);CHK_SET_ERR(rval, "Failed to set coordinates");
  return MB_SUCCESS;
}

void deform_func(const BoundBox &bbox, double *xold, double *xnew) 
{
/*  Deformation function based on max delx and dely at top of rod
    const double RODWIDTH = 0.2, RODHEIGHT = 0.5;
    // function: origin is at middle base of rod, and is .5 high
    // top of rod is (0,.55) on left and (.2,.6) on right
  double delx = 0.5*RODWIDTH;
  
  double xfrac = (xold[0] + .5*RODWIDTH)/RODWIDTH, yfrac = xold[1]/RODHEIGHT;
  xnew[0] = xold[0] + yfrac * delx;
  xnew[1] = xold[1] + yfrac * (1.0 + xfrac) * 0.05;
*/
/* Deformation function based on fraction of bounding box dimension in each direction */
  double frac = 0.01; // taken from approximate relative deformation from LLNL Diablo of XX09 assys
  CartVect *xo = reinterpret_cast<CartVect*>(xold), *xn = reinterpret_cast<CartVect*>(xnew);
  CartVect disp = frac * (*xo - bbox.bMin);
  *xn = *xo + disp;
}

ErrorCode DeformMeshRemap::deform_master(Range &fluid_elems, Range &solid_elems, const char *tag_name) 
{
  // Deform elements with an analytic function
  ErrorCode rval;

  // Get all the vertices and coords in the solid
  Range solid_verts, fluid_verts;
  rval = mbImpl->get_adjacencies(solid_elems, 0, false, solid_verts, Interface::UNION);CHK_SET_ERR(rval, "Failed to get vertices");
  vector<double> coords(3*solid_verts.size()), new_coords(3*solid_verts.size());
  rval = mbImpl->get_coords(solid_verts, &coords[0]);CHK_SET_ERR(rval, "Failed to get vertex coords");
  unsigned int num_verts = solid_verts.size();

  // Get or create the tag
  if (!xDispNames[0].empty() && !xDispNames[1].empty() && !xDispNames[2].empty()) {
    // 3 tags, specifying xyz individual data, integrate into one tag
    rval = mbImpl->tag_get_handle((tag_name ? tag_name : ""), 3, MB_TYPE_DOUBLE,
        xNew, MB_TAG_CREAT | MB_TAG_DENSE);CHK_SET_ERR(rval, "Failed to create xnew tag");
    vector<double> disps(num_verts);
    for (int i = 0; i < 3; i++) {
      rval = mbImpl->tag_get_handle(xDispNames[0].c_str(), 1, MB_TYPE_DOUBLE, xDisp[i]);CHK_SET_ERR(rval, "Failed to get xDisp tag");
      rval = mbImpl->tag_get_data(xDisp[i], solid_verts, &disps[0]);CHK_SET_ERR(rval, "Failed to get xDisp tag values");
      for (unsigned int j = 0; j < num_verts; j++)
        new_coords[3*j+i] = coords[3*j+i] + disps[j];
    }
  }
  else if (!xDispNames[0].empty()) {
    rval = mbImpl->tag_get_handle(xDispNames[0].c_str(), 3, MB_TYPE_DOUBLE, xDisp[0]);CHK_SET_ERR(rval, "Failed to get first xDisp tag");
    xNew = xDisp[0];
    vector<double> disps(3*num_verts);
    rval = mbImpl->tag_get_data(xDisp[0], solid_verts, &disps[0]);
    for (unsigned int j = 0; j < 3*num_verts; j++)
      new_coords[j] = coords[j] + disps[j];
  }
  else {
    // Get the bounding box of the solid mesh
    BoundBox bbox;
    bbox.update(*mbImpl, solid_elems);
  
    for (unsigned int j = 0; j < num_verts; j++)
      deform_func(bbox, &coords[3*j], &new_coords[3*j]);
  }

  if (debug) {
    double len = 0.0;
    for (unsigned int i = 0; i < num_verts; i++) {
      CartVect dx = CartVect(&new_coords[3*i]) - CartVect(&coords[3*i]);
      double tmp_len = dx.length_squared();
      if (tmp_len > len) len = tmp_len;
    }
    Range tmp_elems(fluid_elems);
    tmp_elems.merge(solid_elems);
    BoundBox box;
    box.update(*mbImpl, tmp_elems);
    double max_len = std::max(box.bMax[2] - box.bMin[2], std::max(box.bMax[1] - box.bMin[1], box.bMax[0] - box.bMin[0]));

    cout << "Max displacement = " << len << " (" << 100.0 * len / max_len << "% of max box length)" << endl;
  }

  if (!xNew) {
    rval = mbImpl->tag_get_handle((tag_name ? tag_name : ""), 3, MB_TYPE_DOUBLE, 
                                  xDisp[0], MB_TAG_CREAT|MB_TAG_DENSE);CHK_SET_ERR(rval, "Failed to get xNew tag");
    xNew = xDisp[0];
  }

  // Set the new tag to those coords
  rval = mbImpl->tag_set_data(xNew, solid_verts, &new_coords[0]);CHK_SET_ERR(rval, "Failed to set tag data");

  // Get all the vertices and coords in the fluid, set xnew to them
  rval = mbImpl->get_adjacencies(fluid_elems, 0, false, fluid_verts, Interface::UNION);CHK_SET_ERR(rval, "Failed to get vertices");
  fluid_verts = subtract(fluid_verts, solid_verts);

  if (coords.size() < 3*fluid_verts.size()) coords.resize(3*fluid_verts.size());
  rval = mbImpl->get_coords(fluid_verts, &coords[0]);CHK_SET_ERR(rval, "Failed to get vertex coords");
  rval = mbImpl->tag_set_data(xNew, fluid_verts, &coords[0]);CHK_SET_ERR(rval, "Failed to set xnew tag on fluid verts");

  if (debug) {
    // Save deformed mesh coords to new file for visualizing
    Range tmp_range(fluidElems[MASTER]);
    tmp_range.merge(solidElems[MASTER]);
    rval = write_and_save(tmp_range, masterSet, xNew, "deformed_master.h5m", true);CHK_ERR(rval);
  }

  return MB_SUCCESS;
}

ErrorCode DeformMeshRemap::read_file(int m_or_s, string &fname, EntityHandle &seth)
{
  // Create meshset
  ErrorCode rval = mbImpl->create_meshset(0, seth);CHK_SET_ERR(rval, "Couldn't create master/slave set");
  ostringstream options;
#ifdef USE_MPI
  ParallelComm *pc = (m_or_s == MASTER ? pcMaster : pcSlave);
  if (pc && pc->size() > 1) {
    if (debug) options << "DEBUG_IO=1;CPUTIME;";
    options << "PARALLEL=READ_PART;PARTITION=PARALLEL_PARTITION;PARALLEL_RESOLVE_SHARED_ENTS;"
            << "PARALLEL_GHOSTS=2.0.1;PARALLEL_COMM=" << pc->get_id();
  }
#endif  
  rval = mbImpl->load_file(fname.c_str(), &seth, options.str().c_str());CHK_SET_ERR(rval, "Couldn't load master/slave mesh");

  if (*solidSetNos[m_or_s].begin() == -1 || *fluidSetNos[m_or_s].begin() == -1) return MB_SUCCESS;

  // Get material sets for solid/fluid
  Tag tagh;
  rval = mbImpl->tag_get_handle(MATERIAL_SET_TAG_NAME, tagh);CHK_SET_ERR(rval, "Couldn't get material set tag name");
  for (set<int>::iterator sit = solidSetNos[m_or_s].begin(); sit != solidSetNos[m_or_s].end(); ++sit) {
    Range sets;
    int set_no = *sit;
    const void *setno_ptr = &set_no;
    rval = mbImpl->get_entities_by_type_and_tag(seth, MBENTITYSET, &tagh, &setno_ptr, 1, sets);
    if (MB_SUCCESS != rval || sets.empty()) {
      SET_ERR(MB_FAILURE, "Couldn't find solid set #" << *sit);
    }
    else
      solidSets[m_or_s].merge(sets);
  }

  // Get solid entities, and dimension
  Range tmp_range;
  for (Range::iterator rit = solidSets[m_or_s].begin(); rit != solidSets[m_or_s].end(); ++rit) {
    rval = mbImpl->get_entities_by_handle(*rit, tmp_range, true);CHK_SET_ERR(rval, "Failed to get entities in solid");
  }
  if (!tmp_range.empty()) {
    int dim = mbImpl->dimension_from_handle(*tmp_range.rbegin());
    assert(dim > 0 && dim < 4);
    solidElems[m_or_s] = tmp_range.subset_by_dimension(dim);
  }

  if (debug)
    cout << "Read " << solidElems[m_or_s].size() << " solid elements from " << solidSets[m_or_s].size() <<
    " sets in " << (m_or_s == MASTER ? "master" : "slave") << " mesh." << endl;

  for (set<int>::iterator sit = fluidSetNos[m_or_s].begin(); sit != fluidSetNos[m_or_s].end(); ++sit) {
    Range sets;
    int set_no = *sit;
    const void *setno_ptr = &set_no;
    rval = mbImpl->get_entities_by_type_and_tag(seth, MBENTITYSET, &tagh, &setno_ptr, 1, sets);
    if (MB_SUCCESS != rval || sets.empty()) {
      SET_ERR(MB_FAILURE, "Couldn't find fluid set #" << *sit);
    }
    else
      fluidSets[m_or_s].merge(sets);
  }

  // Get fluid entities, and dimension
  tmp_range.clear();
  for (Range::iterator rit = fluidSets[m_or_s].begin(); rit != fluidSets[m_or_s].end(); ++rit) {
    rval = mbImpl->get_entities_by_handle(*rit, tmp_range, true);CHK_SET_ERR(rval, "Failed to get entities in fluid");
  }
  if (!tmp_range.empty()) {
    int dim = mbImpl->dimension_from_handle(*tmp_range.rbegin());
    assert(dim > 0 && dim < 4);
    fluidElems[m_or_s] = tmp_range.subset_by_dimension(dim);
  }

  if (debug)
    cout << "Read " << fluidElems[m_or_s].size() << " fluid elements from " << fluidSets[m_or_s].size() <<
      " sets in " << (m_or_s == MASTER ? "master" : "slave") << " mesh." << endl;

  return rval;
}

ErrorCode DeformMeshRemap::find_other_sets(int m_or_s, EntityHandle file_set) 
{
  // Solid or fluid sets are missing; find the other
  Range *filled_sets = NULL, *unfilled_sets = NULL, *unfilled_elems = NULL;

  if (fluidSets[m_or_s].empty() && !solidSets[m_or_s].empty()) {
    unfilled_sets = &fluidSets[m_or_s];
    filled_sets = &solidSets[m_or_s];
    unfilled_elems = &fluidElems[m_or_s];
    if (debug)
      cout << "Finding unspecified fluid elements in " << (m_or_s == MASTER ? "master" : "slave") << " mesh...";
  }
  else if (!fluidSets[m_or_s].empty() && solidSets[m_or_s].empty()) {
    filled_sets = &fluidSets[m_or_s];
    unfilled_sets = &solidSets[m_or_s];
    unfilled_elems = &solidElems[m_or_s];
    if (debug)
      cout << "Finding unspecified solid elements in " << (m_or_s == MASTER ? "master" : "slave") << " mesh...";
  }

  // Ok, we know the filled sets, now fill the unfilled sets, and the elems from those
  Tag tagh;
  ErrorCode rval = mbImpl->tag_get_handle(MATERIAL_SET_TAG_NAME, tagh);CHK_SET_ERR(rval, "Couldn't get material set tag name");
  Range matsets;
  rval = mbImpl->get_entities_by_type_and_tag(file_set, MBENTITYSET, &tagh, NULL, 1, matsets);
  if (MB_SUCCESS != rval || matsets.empty()) {
    SET_ERR(MB_FAILURE, "Couldn't get any material sets");
  }
  *unfilled_sets = subtract(matsets, *filled_sets);
  if (unfilled_sets->empty()) {
    SET_ERR(MB_FAILURE, "Failed to find any unfilled material sets");
  }
  Range tmp_range;
  for (Range::iterator rit = unfilled_sets->begin(); rit != unfilled_sets->end(); ++rit) {
    rval = mbImpl->get_entities_by_handle(*rit, tmp_range, true);CHK_SET_ERR(rval, "Failed to get entities in unfilled set");
  }
  int dim = mbImpl->dimension_from_handle(*tmp_range.rbegin());
  assert(dim > 0 && dim < 4);  
  *unfilled_elems = tmp_range.subset_by_dimension(dim);
  if (unfilled_elems->empty()) {
    SET_ERR(MB_FAILURE, "Failed to find any unfilled set entities");
  }

  if (debug) 
    cout << "found " << unfilled_sets->size() << " sets and " << unfilled_elems->size() << " elements." << endl;

  return MB_SUCCESS;
}
