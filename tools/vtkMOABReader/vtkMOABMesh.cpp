#include "vtkMOABMesh.hpp"

#include "vtkFloatArray.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkIntArray.h"
#include "vtkCharArray.h"
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"
#include "vtkUnstructuredGrid.h"
#include "vtkExtractUnstructuredGrid.h"
#include "vtkThreshold.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPointData.h"
#include <sstream>
#include <vector>
#include "assert.h"
#include "MBTagConventions.hpp"

#include "moab/Core.hpp"
#include "moab/CN.hpp"

#define MOABMeshErrorMacro(s) std::cerr s << std::endl
using namespace moab;

const int vtkMOABMesh::vtk_cell_types[] = {
  1, 3, 5, 9, 7, 10, 14, 13, 0, 12, 0, 0, 0};

const bool new_outputs = false;
const bool use_filters = true;

vtkMOABMesh *vtkMOABMesh::instance_ = NULL;

vtkMOABMesh *vtkMOABMesh::instance(Interface *iface, bool create_new) 
{
  if (!instance_ && create_new) 
    instance_ = new vtkMOABMesh(iface);
  
  return instance_;
}

ErrorCode vtkMOABMesh::delete_instance() 
{
  if (instance_) {
    delete instance_;
    instance_ = NULL;
  }
  
  return MB_SUCCESS;
}

vtkMOABMesh::vtkMOABMesh(Interface *iface) 
        : mbImpl(iface), myUG(NULL), iFace(NULL), maxPointId(-1), maxCellId(-1), 
          vtkIdTag(0), outOfDate(true)
{
  if (!mbImpl) mbImpl = new Core();

  mbImpl->query_interface(iFace);
  assert(NULL != iFace);

  int def_val = -1;
  ErrorCode rval = mbImpl->tag_get_handle("__vtkIdTag", 1, MB_TYPE_INTEGER,
                                          vtkIdTag, MB_TAG_DENSE|MB_TAG_CREAT, 
                                          &def_val);
  assert(MB_SUCCESS == rval);
}

ErrorCode vtkMOABMesh::load_file(const char *file_name, const char *options, 
                                 EntityHandle &file_set) 
{
  ErrorCode rval;
  
  rval = mbImpl->create_meshset(MESHSET_SET, file_set);
  if (MB_SUCCESS != rval) return rval;

  fileSets.insert(file_set);
    
  rval = mbImpl->load_file(file_name, &file_set, options);
  if (MB_SUCCESS != rval) return rval;

  outOfDate = true;
  
  fileNames.push_back(std::string(file_name));

  return rval;
}

void vtkMOABMesh::Execute() 
{
  Update();
}

vtkMOABMesh::~vtkMOABMesh()
{
  if (mbImpl) {
    if (iFace)
      mbImpl->release_interface(iFace);

    delete mbImpl;
  }
}

ErrorCode vtkMOABMesh::Update(EntityHandle file_set)
{
  if (!outOfDate) return MB_SUCCESS;
  
    // get the data set & allocate an initial chunk of data
  vtkUnstructuredGrid *ug = this->GetOutput();
  ug->Allocate();
  
  ErrorCode rval  = construct_mesh(file_set);
  if (MB_SUCCESS != rval)
  {
    MOABMeshErrorMacro( << "Failed to construct mesh");
    return rval;
  }
  MOABMeshErrorMacro(<<"Constructed mesh...");

  if (use_filters)
    rval = construct_filters();
  if (MB_SUCCESS != rval)
  {
    MOABMeshErrorMacro( << "Failed to construct filters ");
    return MB_FAILURE;
  }
  MOABMeshErrorMacro(<<"Filters constructed...");
  
    // get all dense tags
  rval = read_tags(file_set);
  MOABMeshErrorMacro(<<"Tags read...");

  outOfDate = false;
  
  MOABMeshErrorMacro(<< "After Update: ug has " << myUG->GetNumberOfPoints()
                     << " points, " << myUG->GetNumberOfCells() << " cells.");

  return rval;
}

ErrorCode vtkMOABMesh::read_tags(EntityHandle file_set) 
{
  ErrorCode rval = read_dense_tags(file_set);
  if (MB_SUCCESS != rval) return rval;

  rval = read_sparse_tags(file_set);
  return rval;
}

ErrorCode vtkMOABMesh::read_dense_tags(EntityHandle file_set) 
{  
    // get all the tags
  std::vector<Tag> tmptags, all_tags;
  ErrorCode rval = mbImpl->tag_get_tags(tmptags);
  if (MB_SUCCESS != rval) return rval;
  
  for (std::vector<Tag>::iterator vit = tmptags.begin(); vit != tmptags.end(); vit++) {
      // skip sparse tags
    TagType ttype;
    rval = mbImpl->tag_get_type(*vit, ttype);
    if (MB_SUCCESS == rval && MB_TAG_DENSE == ttype) all_tags.push_back(*vit);
  }

    // now create field arrays on all 2d and 3d entities
  Range ents2d, ents3d, verts;
  rval = mbImpl->get_entities_by_dimension(file_set, 2, ents2d);
  if (MB_SUCCESS != rval) return rval;

  rval = mbImpl->get_entities_by_dimension(file_set, 3, ents3d);
  if (MB_SUCCESS != rval) return rval;

  rval = mbImpl->get_entities_by_dimension(file_set, 0, verts);
  if (MB_SUCCESS != rval) return rval;

  int *vids;
  void *data;
  if (MB_SUCCESS != rval) return rval;
  std::vector<double> tag_dvals;
  std::vector<int> tag_ivals;
  vtkIntArray *int_array;
  vtkDoubleArray *dbl_array;
  int idef, *idata;
  double ddef, *ddata;
  int min, max;
      
  for (std::vector<Tag>::iterator vit = all_tags.begin(); vit != all_tags.end(); vit++) {
    if (*vit == vtkIdTag) continue;
    
      // create a data array
    DataType dtype;
    rval = mbImpl->tag_get_data_type(*vit, dtype);
    if (MB_SUCCESS != rval) continue;
    std::string tag_name;
    bool has_default = false;
    rval = mbImpl->tag_get_name(*vit, tag_name);
    if (MB_SUCCESS != rval) continue;
    if (MB_TYPE_DOUBLE == dtype) {
      dbl_array = vtkDoubleArray::New();
      dbl_array->SetName(tag_name.c_str());
      if (MB_SUCCESS == mbImpl->tag_get_default_value(*vit, &ddef))
        has_default = true;
    }
    else if (MB_TYPE_INTEGER == dtype) {
      int_array = vtkIntArray::New();
      int_array->SetName(tag_name.c_str());
      if (MB_SUCCESS == mbImpl->tag_get_default_value(*vit, &idef))
        has_default = true;
    }

    if (MB_SUCCESS != rval) continue;

    Range::iterator rit, rit2;
    rit = rit2 = ents2d.begin();
    while (rit != ents2d.end()) {
        // get tag iterator for gids
      rval = mbImpl->tag_iterate(vtkIdTag, rit, ents2d.end(), (void*&)vids);
      if (MB_SUCCESS != rval) continue;
      int count = rit - rit2;
      
      rval = mbImpl->tag_iterate(*vit, rit2, ents2d.end(), data);
      if (MB_SUCCESS != rval) continue;
      
      if (MB_TYPE_DOUBLE == dtype) {
        ddata = (double*)data;
        for (int i = 0; i < count; i++) {
          assert(-1 < vids[i] && vids[i] <= maxCellId);
          if (!has_default || ddata[i] != ddef)
            dbl_array->InsertValue(vids[i], ddata[i]);
        }
      }
      else if (MB_TYPE_INTEGER == dtype) {
        idata = (int*)data;
        for (int i = 0; i < count; i++) {
          assert(-1 < vids[i] && vids[i] <= maxCellId);
          if (!has_default || idata[i] != idef)
            int_array->InsertValue(vids[i], idata[i]);
        }
      }
      
      min = *std::min_element(vids, vids+count);
      max = *std::max_element(vids, vids+count);
      MOABMeshErrorMacro(<< "2d: min = " << min << ", max =  " << max);
    }
    
    rit = rit2 = ents3d.begin();
    while (rit != ents3d.end()) {
        // get tag iterator for vids
      rval = mbImpl->tag_iterate(vtkIdTag, rit, ents3d.end(), (void*&)vids);
      if (MB_SUCCESS != rval) continue;
      int count = rit - rit2;
      
      rval = mbImpl->tag_iterate(*vit, rit2, ents3d.end(), data);
      if (MB_SUCCESS != rval) continue;
      
      if (MB_TYPE_DOUBLE == dtype) {
        ddata = (double*)data;
        for (int i = 0; i < count; i++) {
          assert(-1 < vids[i] && vids[i] <= maxCellId);
          if (!has_default || ddata[i] != ddef)
            dbl_array->InsertValue(vids[i], ddata[i]);
        }
      }
      else if (MB_TYPE_INTEGER == dtype) {
        idata = (int*)data;
        for (int i = 0; i < count; i++) {
          assert(-1 < vids[i] && vids[i] <= maxCellId);
          if (!has_default || idata[i] != idef)
            int_array->InsertValue(vids[i], idata[i]);
        }
      }
      min = *std::min_element(vids, vids+count);
      max = *std::max_element(vids, vids+count);
      MOABMeshErrorMacro(<< "3d: min = " << min << ", max =  " << max);
    }
    
    if (MB_TYPE_DOUBLE == dtype) {
      this->GetOutput()->GetCellData()->AddArray(dbl_array);
      dbl_array->Delete();
      MOABMeshErrorMacro(<< "Read " << dbl_array->GetSize() << " values of dbl tag " << tag_name);
    }
    else if (MB_TYPE_INTEGER == dtype) {
      this->GetOutput()->GetCellData()->AddArray(int_array);
      int_array->Delete();
      MOABMeshErrorMacro(<< "Read " << int_array->GetSize() << " values of int tag " << tag_name);
    }

    rit = rit2 = verts.begin();
    if (MB_TYPE_DOUBLE == dtype) {
      dbl_array = vtkDoubleArray::New();
      dbl_array->SetName(tag_name.c_str());
    }
    else if (MB_TYPE_INTEGER == dtype) {
      int_array = vtkIntArray::New();
      int_array->SetName(tag_name.c_str());
    }
    while (rit != verts.end()) {
        // get tag iterator for vids
      rval = mbImpl->tag_iterate(vtkIdTag, rit, verts.end(), (void*&)vids);
      if (MB_SUCCESS != rval) continue;
      int count = rit - rit2;
      
      rval = mbImpl->tag_iterate(*vit, rit2, verts.end(), data);
      if (MB_SUCCESS != rval) continue;
      
      if (MB_TYPE_DOUBLE == dtype) {
        ddata = (double*)data;
        for (int i = 0; i < count; i++) {
          assert(vids[i] >= 0 && vids[i] <= maxPointId);
          if (!has_default || ddata[i] != ddef)
            dbl_array->InsertValue(vids[i], ddata[i]);
        }
      }
      else if (MB_TYPE_INTEGER == dtype) {
        idata = (int*)data;
        for (int i = 0; i < count; i++) {
          assert(vids[i] >= 0 && vids[i] <= maxPointId);
          if (!has_default || idata[i] != idef)
            int_array->InsertValue(vids[i], idata[i]);
        }
      }
      
      min = *std::min_element(vids, vids+count);
      max = *std::max_element(vids, vids+count);
      MOABMeshErrorMacro(<< "verts: min = " << min << ", max =  " << max);
    }
    
    if (MB_TYPE_DOUBLE == dtype) {
      this->GetOutput()->GetPointData()->AddArray(dbl_array);
      dbl_array->Delete();
      MOABMeshErrorMacro(<< "Read " << dbl_array->GetSize() << " values of dbl tag " << tag_name);
    }
    else if (MB_TYPE_INTEGER == dtype) {
      this->GetOutput()->GetPointData()->AddArray(int_array);
      int_array->Delete();
      MOABMeshErrorMacro(<< "Read " << int_array->GetSize() << " values of int tag " << tag_name);
    }
  }

  return MB_SUCCESS;
}

ErrorCode vtkMOABMesh::read_sparse_tags(EntityHandle file_set) 
{  
    // get all the tags
  std::vector<Tag> tmptags, all_tags;
  ErrorCode rval = mbImpl->tag_get_tags(tmptags);
  if (MB_SUCCESS != rval) return rval;
  
  for (std::vector<Tag>::iterator vit = tmptags.begin(); vit != tmptags.end(); vit++) {
      // skip dense tags
    TagType ttype;
    DataType dtype;
    rval = mbImpl->tag_get_type(*vit, ttype);
    rval = mbImpl->tag_get_data_type(*vit, dtype);
    if (MB_SUCCESS == rval && MB_TAG_SPARSE == ttype && MB_TYPE_INTEGER == dtype) 
      all_tags.push_back(*vit);
  }

    // now create field arrays on all 2d and 3d entities
  Range sets, ents, verts;
  vtkIntArray *int_array;

  Tag gid_tag, gdim_tag;
  rval = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, 1, MB_TYPE_INTEGER, gid_tag);
  if (MB_SUCCESS != rval) return rval;
      
  rval = mbImpl->tag_get_handle(GEOM_DIMENSION_TAG_NAME, 1, MB_TYPE_INTEGER, gdim_tag);
  if (MB_SUCCESS != rval) return rval;
      
  std::vector<int> vids;
  for (std::vector<Tag>::iterator vit = all_tags.begin(); vit != all_tags.end(); vit++) {
    if (*vit == vtkIdTag) continue;

      // if this is a geometry tag, loop
    int lmax = (*vit == gdim_tag ? 3 : 0);
    static int lvals[] = {0, 1, 2, 3};
    static const char *lnames[] = {"GeomVertex", "GeomCurve", "GeomSurface", "GeomVolume"};
    
    for (int l = 1; l <= lmax; l++) {
      sets.clear();
      int *lval = lvals+l;
      rval = mbImpl->get_entities_by_type_and_tag(file_set, MBENTITYSET, &(*vit), 
                                                  (const void* const*)(lmax ? &lval : NULL), 1, sets);
      if (MB_SUCCESS != rval || sets.empty()) continue;
      
        // create a data array
      std::string tag_name;
      bool has_default = false;
      if (lmax) tag_name = std::string(lnames[l]);
      else {
        rval = mbImpl->tag_get_name(*vit, tag_name);
        if (MB_SUCCESS != rval) continue;
      }
      if (MB_SUCCESS != rval) continue;
      int_array = vtkIntArray::New();
      int_array->SetName(tag_name.c_str());
      bool had_ents = false;
      
        // loop over sets then entities
      for (Range::iterator rit = sets.begin(); rit != sets.end(); rit++) {
          // get the tag value
        int this_val;
        rval = mbImpl->tag_get_data((lmax ? gid_tag : *vit), &(*rit), 1, &this_val);
        if (MB_SUCCESS != rval) continue;

          // get the non-vertex entities, and their vtk ids
        ents.clear();
        for (int d = 1; d <= 3; d++) {
          rval = mbImpl->get_entities_by_dimension(*rit, d, ents, true);
          if (MB_SUCCESS != rval) continue;
        }
        if (ents.empty()) continue;
        vids.resize(ents.size());
        rval = mbImpl->tag_get_data(vtkIdTag, ents, &vids[0]);
        if (MB_SUCCESS != rval || ents.empty()) continue;

        std::cout << "Tag " << tag_name << ", value " << this_val << ", entities:" << std::endl;
        ents.print("   ");
        
        for (unsigned int e = 0; e < vids.size(); e++) {
          assert(-1 != vids[e]);
          int_array->InsertValue(vids[e], this_val);
        }

        had_ents = true;
      }

        // add the data array to the output
      if (had_ents) this->GetOutput()->GetCellData()->AddArray(int_array);
      int_array->Delete();
    }
  }

  return MB_SUCCESS;
}

ErrorCode vtkMOABMesh::construct_mesh(EntityHandle file_set) 
{
    // construct the vtk representation of the mesh
  
    // get all the hexes and quads
  Range all_elems;
  ErrorCode result = MB_SUCCESS, tmp_result;
  for (int dim = 0; dim <= 3; dim++) 
  {
    tmp_result = mbImpl->get_entities_by_dimension(file_set, dim, all_elems);
    if (tmp_result != MB_SUCCESS) result = tmp_result;
  }
  if (MB_SUCCESS != result)
    {
    MOABMeshErrorMacro( << "Failure getting hexes from mesh. " );
    return result;
    }

  MOABMeshErrorMacro(<< "Read " << all_elems.size() << " entities from MOAB.");

    // create the elements
  int success = this->create_elements(file_set);
  if (MB_SUCCESS != result)
    {
    MOABMeshErrorMacro( << "Problem filling in quad data. " );
    return result;
    }

  return MB_SUCCESS;
  
}

ErrorCode vtkMOABMesh::create_points_vertices(EntityHandle file_set, Range &verts) 
{
    // get the global id tag
  ErrorCode result;

  result = mbImpl->get_entities_by_type(file_set, MBVERTEX, verts);
  if (MB_SUCCESS != result)
  {
    MOABMeshErrorMacro( << "Couldn't gather vertices. " );
    return result;
  }

  MOABMeshErrorMacro(<< "Gathered " << verts.size() << " vertices from MOAB.");
  
    // assign ids to the vertices
  std::vector<int> vids(verts.size());
  for (unsigned int i = 0; i < verts.size(); i++)
    vids[i] = ++maxPointId;

  result = mbImpl->tag_set_data(vtkIdTag, verts, &vids[0]);
  if (MB_SUCCESS != result)
  {
    MOABMeshErrorMacro( << "Couldn't set ids on vertices. " );
    return result;
  }
  
    // allocate and fill in coordinate arrays
  std::vector<double*> coords(3);
  coords[0] = new double[verts.size()];
  coords[1] = new double[verts.size()];
  coords[2] = new double[verts.size()];
  result = iFace->get_node_coords(3, verts.size(), verts,
                                  0, 0, coords);
  if (MB_SUCCESS != result)
  {
    MOABMeshErrorMacro( << "Couldn't get nodal coordinates. " );
    return result;
  }

    // put these data into a point array
  vtkPoints *points = vtkPoints::New();
  int dum;
  points->SetNumberOfPoints(verts.size());
  assert(MB_SUCCESS == result);
  unsigned int i = 0;
  for (Range::const_iterator rit = verts.begin(); rit != verts.end(); rit++, i++)
  {
    points->SetPoint(vids[i], coords[0][i], coords[1][i], coords[2][i]);
  }
  myUG->SetPoints(points);
  points->Delete();

  return MB_SUCCESS;
}

ErrorCode vtkMOABMesh::create_elements(EntityHandle file_set)
{
    // get the vertices
  Range verts;
  ErrorCode result;

    // create points/vertices in vtk database
  result = create_points_vertices(file_set, verts);
  if (MB_SUCCESS != result)
  {
    MOABMeshErrorMacro( << "Couldn't create points/vertices. " );
    return result;
  }

  MOABMeshErrorMacro(<< "After create_points_vertices: ug has " << myUG->GetNumberOfPoints()
                     << " points, " << myUG->GetNumberOfCells() << " cells.");
  
    // for the remaining elements, add them individually
  int ids[CN::MAX_NODES_PER_ELEMENT];
  vtkIdType vtkids[CN::MAX_NODES_PER_ELEMENT];
  const EntityHandle *connect;
  int num_connect;
  bool first = true;

  for (EntityType this_type = MBEDGE; this_type != MBENTITYSET; this_type++) {

      // don't try to represent elements vtk doesn't understand
    if (vtk_cell_types[this_type] == 0) continue;
    
    Range elems;
    result = mbImpl->get_entities_by_type(file_set, this_type, elems);
    if (MB_SUCCESS != result)
    {
      MOABMeshErrorMacro( << "Couldn't get elements. " );
      return result;
    }

    std::vector<int> eids(elems.size());
    result = mbImpl->tag_get_data(vtkIdTag, elems, &eids[0]);
    if (MB_SUCCESS != result)
    {
      MOABMeshErrorMacro( << "Couldn't get elements vtkIdTag. " );
      return result;
    }
    
    int e = 0;
    bool changed = false;
    for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++, e++) {
      if (-1 != eids[e]) continue;
      
      changed = true;
      
        // get the connectivity of these elements
      result = mbImpl->get_connectivity(*rit, connect, num_connect, true);
      if (MB_SUCCESS != result)
      {
        MOABMeshErrorMacro( << "Couldn't get element connectivity. " );
        return result;
      }

        // get the id tag for these vertices
      result = mbImpl->tag_get_data(vtkIdTag, connect, num_connect, ids);
      if (MB_SUCCESS != result)
      {
        MOABMeshErrorMacro( << "Couldn't get vertex ids for element. " );
        return result;
      }

        // ok, now insert this cell
      for (int i = 0; i < num_connect; i++) vtkids[i] = ids[i];
      int last_id = myUG->InsertNextCell(vtk_cell_types[this_type], num_connect, vtkids);
      eids[e] = last_id;
      maxCellId = std::max(last_id, maxCellId);
    }

    if (changed) {
      result = mbImpl->tag_set_data(vtkIdTag, elems, &eids[0]);
      if (MB_SUCCESS != result)
      {
        MOABMeshErrorMacro( << "Couldn't save element ids. " );
        return result;
      }
    }
  }
  
  MOABMeshErrorMacro(<< "After creating cells: ug has " << myUG->GetNumberOfPoints()
                     << " points, " << myUG->GetNumberOfCells() << " cells.");
  
  return MB_SUCCESS;
}

moab::ErrorCode vtkMOABMesh::construct_filters() 
{
    // apply threshold and type filters to the output to get multiple actors
    // corresponding to dual surfaces and curves, then group the dual actors
    // together using a group filter

  vtkUnstructuredGrid *ug = this->GetOutput();
/*
    // first, get the non-dual mesh
  vtkExtractUnstructuredGrid *primal = vtkExtractUnstructuredGrid::New();
  primal->SetInput(this->GetOutput());
  primal->SetCellMinimum(0);
  primal->SetCellMaximum(this->MaxPrimalId);

    // set merging on so points aren't duplicated
  primal->SetMerging(1);

    // now do dual surfaces; do threshold-based extraction for now
  MBTag vtkIdTag;
  MBErrorCode result = mbImpl->tag_get_handle(GLOBAL_ID_TAG_NAME, vtkIdTag);
  assert(MB_SUCCESS == result && 0 != vtkIdTag);
  
  int ds_id;
  for (ds_id = 0; ds_id < this->NumberOfDualSurfaces; ds_id++) {
    vtkThreshold *ds_filter = vtkThreshold::New();
    ds_filter->SelectInputScalars(DUAL_SURF_ATTRIBUTE_NAME);
    ds_filter->SetAttributeModeToUseCellData();
    ds_filter->ThresholdBetween(((double)ds_id-0.5), ((double)ds_id+0.5));
    ds_filter->SetInput(ug);
    this->add_name(ds_filter->GetOutput(), "dual_surf_", ds_id);
  }
  
    // same for dual curves
  int dc_id;
  for (dc_id = 0; dc_id < this->NumberOfDualCurves; dc_id++) {
    vtkThreshold *dc_filter = vtkThreshold::New();
    dc_filter->SelectInputScalars(DUAL_CURVE_ATTRIBUTE_NAME);
    dc_filter->SetAttributeModeToUseCellData();
    dc_filter->ThresholdBetween(((double)dc_id-0.5), ((double)dc_id+0.5));
    dc_filter->SetInput(ug);
    this->add_name(dc_filter->GetOutput(), "dual_curve_", dc_id);
  }

    // lastly, get the dual vertices and put those in a group
    // first, get the non-dual mesh
  vtkExtractUnstructuredGrid *dual_verts = vtkExtractUnstructuredGrid::New();
  dual_verts->SetCellMinimum(this->DualVertexIdOffset);
  dual_verts->SetCellMaximum(this->DualVertexIdOffset+this->NumberOfDualVertices-1);
  this->add_name(dual_verts->GetOutput(), "dual_verts", 0);
*/
  return MB_SUCCESS;
}

void vtkMOABMesh::add_name(vtkUnstructuredGrid *output, const char *prefix,
                           const int id) 
{
  vtkCharArray* nmArray =  vtkCharArray::New();
  nmArray->SetName("Name");
  vtkstd::ostringstream name;
  name << prefix << id << "\0";
  nmArray->SetNumberOfTuples(static_cast<vtkIdType>(name.str().length()));
  char* copy = nmArray->GetPointer(0);
  memcpy(copy, name.str().c_str(), name.str().length());
  output->GetFieldData()->AddArray(nmArray);
  nmArray->Delete();
}
