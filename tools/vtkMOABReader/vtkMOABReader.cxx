/*=========================================================================

  Program:   Visualization Toolkit
  Module:    $RCSfile$

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkMOABReader.h"

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

#define IS_BUILDING_MB
#include "moab/Core.hpp"
#undef IS_BUILDING_MB
#include "moab/CN.hpp"

using namespace moab;

vtkCxxRevisionMacro(vtkMOABReader, "$Revision$")
vtkStandardNewMacro(vtkMOABReader)

moab::Interface *gMB = NULL;

const int vtkMOABReader::vtk_cell_types[] = {
  1, 3, 5, 9, 7, 10, 14, 13, 0, 12, 0, 0, 0};

const bool new_outputs = false;
const bool use_filters = true;

vtkMOABReader::vtkMOABReader()
{
  this->FileName = NULL;
  VtkOffsetIdTag = 0;
  
}

vtkMOABReader::~vtkMOABReader()
{
  if (NULL != gMB) {
    delete gMB;
    gMB = NULL;
  }
}

void vtkMOABReader::Execute()
{
  this->DebugOn();

    // initialize MOAB & read file
  if (NULL == gMB) gMB = new Core();
  
    // make an offset id tag
  moab::ErrorCode result;
  int success;
  if (0 == VtkOffsetIdTag) {
    
    result = gMB->tag_get_handle("__vtk_offset_id_tag", VtkOffsetIdTag);
    if (MB_SUCCESS != result) {
      int dum = -1;
      result = gMB->tag_create("__vtk_offset_id_tag", sizeof(int), MB_TAG_SPARSE, VtkOffsetIdTag, &dum);
    }
  }
  
  assert (MB_SUCCESS == result && VtkOffsetIdTag != 0);
  result = gMB->load_mesh(this->GetFileName());
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to open file " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Read MOAB file...");

    // get the data set & allocate an initial chunk of data
  vtkUnstructuredGrid *ug = this->GetOutput();
  ug->Allocate();
  
  success = construct_mesh();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct mesh from " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Constructed mesh...");

  if (use_filters)
    success = construct_filters();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct filters ");
    return;
    }
  vtkDebugMacro(<<"Filters constructed...");
  
    // get all dense tags
  success = read_tags();
  vtkDebugMacro(<<"Tags read...");
  
}

int vtkMOABReader::read_tags() 
{
    // get all the tags
  std::vector<Tag> tmptags, all_tags;
  moab::ErrorCode rval = gMB->tag_get_tags(tmptags);
  if (MB_SUCCESS != rval) return rval;
  
  for (std::vector<Tag>::iterator vit = tmptags.begin(); vit != tmptags.end(); vit++) {
      // skip sparse tags
    TagType ttype;
    rval = gMB->tag_get_type(*vit, ttype);
    if (MB_SUCCESS == rval && MB_TAG_DENSE == ttype) all_tags.push_back(*vit);
  }

    // now create field arrays on all 2d and 3d entities
  Range ents2d, ents3d, verts;
  rval = gMB->get_entities_by_dimension(0, 2, ents2d);
  if (MB_SUCCESS != rval) return rval;

  rval = gMB->get_entities_by_dimension(0, 3, ents3d);
  if (MB_SUCCESS != rval) return rval;

  rval = gMB->get_entities_by_dimension(0, 0, verts);
  if (MB_SUCCESS != rval) return rval;

  int *gids;
  void *data;
  moab::Tag gid_tag;
  rval = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != rval) return rval;
  std::vector<double> tag_dvals;
  std::vector<int> tag_ivals;
  vtkIntArray *int_array;
  vtkDoubleArray *dbl_array;
  int idef, *idata;
  double ddef, *ddata;
  int min, max;
      
  for (std::vector<Tag>::iterator vit = all_tags.begin(); vit != all_tags.end(); vit++) {
    if (*vit == gid_tag) continue;
    
      // create a data array
    DataType dtype;
    rval = gMB->tag_get_data_type(*vit, dtype);
    if (MB_SUCCESS != rval) continue;
    std::string tag_name;
    bool has_default = false;
    rval = gMB->tag_get_name(*vit, tag_name);
    if (MB_SUCCESS != rval) continue;
    if (MB_TYPE_DOUBLE == dtype) {
      dbl_array = vtkDoubleArray::New();
      dbl_array->SetName(tag_name.c_str());
      if (MB_SUCCESS == gMB->tag_get_default_value(*vit, &ddef))
        has_default = true;
    }
    else if (MB_TYPE_INTEGER == dtype) {
      int_array = vtkIntArray::New();
      int_array->SetName(tag_name.c_str());
      if (MB_SUCCESS == gMB->tag_get_default_value(*vit, &idef))
        has_default = true;
    }

    if (MB_SUCCESS != rval) continue;

    Range::iterator rit, rit2;
    rit = rit2 = ents2d.begin();
    while (rit != ents2d.end()) {
        // get tag iterator for gids
      rval = gMB->tag_iterate(gid_tag, rit, ents2d.end(), (void*&)gids);
      if (MB_SUCCESS != rval) continue;
      int count = rit - rit2;
      
      rval = gMB->tag_iterate(*vit, rit2, ents2d.end(), data);
      if (MB_SUCCESS != rval) continue;
      
      if (MB_TYPE_DOUBLE == dtype) {
        ddata = (double*)data;
        for (int i = 0; i < count; i++)
          if (!has_default || ddata[i] != ddef)
            dbl_array->InsertValue(gids[i], ddata[i]);
      }
      else if (MB_TYPE_INTEGER == dtype) {
        idata = (int*)data;
        for (int i = 0; i < count; i++)
          if (!has_default || idata[i] != idef)
            int_array->InsertValue(gids[i], idata[i]);
      }
      
      min = *std::min_element(gids, gids+count);
      max = *std::max_element(gids, gids+count);
      vtkDebugMacro(<< "2d: min = " << min << ", max =  " << max);
    }
    
    rit = rit2 = ents3d.begin();
    while (rit != ents3d.end()) {
        // get tag iterator for gids
      rval = gMB->tag_iterate(gid_tag, rit, ents3d.end(), (void*&)gids);
      if (MB_SUCCESS != rval) continue;
      int count = rit - rit2;
      
      int *gid_data = (int*)data;
      rval = gMB->tag_iterate(*vit, rit2, ents3d.end(), data);
      if (MB_SUCCESS != rval) continue;
      
      if (MB_TYPE_DOUBLE == dtype)
        for (int i = 0; i < count; i++)
          if (!has_default || ddata[i] != ddef)
            dbl_array->InsertValue(gids[i], ((double*)data)[i]);
      else if (MB_TYPE_INTEGER == dtype)
        for (int i = 0; i < count; i++)
          if (!has_default || idata[i] != idef)
            int_array->InsertValue(gids[i], ((int*)data)[i]);

      min = *std::min_element(gids, gids+count);
      max = *std::max_element(gids, gids+count);
      vtkDebugMacro(<< "3d: min = " << min << ", max =  " << max);
    }
    
    if (MB_TYPE_DOUBLE == dtype) {
      this->GetOutput()->GetCellData()->AddArray(dbl_array);
      dbl_array->Delete();
      vtkDebugMacro(<< "Read " << dbl_array->GetSize() << " values of dbl tag " << tag_name);
    }
    else if (MB_TYPE_INTEGER == dtype) {
      this->GetOutput()->GetCellData()->AddArray(int_array);
      int_array->Delete();
      vtkDebugMacro(<< "Read " << int_array->GetSize() << " values of int tag " << tag_name);
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
        // get tag iterator for gids
      rval = gMB->tag_iterate(gid_tag, rit, verts.end(), (void*&)gids);
      if (MB_SUCCESS != rval) continue;
      int count = rit - rit2;
      
      int *gid_data = (int*)data;
      rval = gMB->tag_iterate(*vit, rit2, verts.end(), data);
      if (MB_SUCCESS != rval) continue;
      
      if (MB_TYPE_DOUBLE == dtype)
        for (int i = 0; i < count; i++) {
          assert(gids[i] == i);
          if (!has_default || ddata[i] != ddef)
            dbl_array->InsertValue(gids[i], ((double*)data)[i]);
        }
      else if (MB_TYPE_INTEGER == dtype)
        for (int i = 0; i < count; i++) {
          assert(gids[i] == i);
          if (!has_default || idata[i] != idef)
            int_array->InsertValue(gids[i], ((int*)data)[i]);
        }

      min = *std::min_element(gids, gids+count);
      max = *std::max_element(gids, gids+count);
      vtkDebugMacro(<< "verts: min = " << min << ", max =  " << max);
    }
    
    if (MB_TYPE_DOUBLE == dtype) {
      this->GetOutput()->GetPointData()->AddArray(dbl_array);
      dbl_array->Delete();
      vtkDebugMacro(<< "Read " << dbl_array->GetSize() << " values of dbl tag " << tag_name);
    }
    else if (MB_TYPE_INTEGER == dtype) {
      this->GetOutput()->GetPointData()->AddArray(int_array);
      int_array->Delete();
      vtkDebugMacro(<< "Read " << int_array->GetSize() << " values of int tag " << tag_name);
    }

  }

  return 0;
}

int vtkMOABReader::construct_mesh() 
{
    // construct the vtk representation of the mesh; just do hexes and quads 
    // for now
  moab::WriteUtilIface *iface = NULL;
  gMB->query_interface("WriteUtilIface", reinterpret_cast<void**>(&iface));
  assert(NULL != iface);
  
    // get all the hexes and quads
  moab::Range all_elems;
  moab::ErrorCode result = MB_SUCCESS, tmp_result;
  for (int dim = 0; dim <= 3; dim++) 
  {
    tmp_result = gMB->get_entities_by_dimension(0, dim, all_elems);
    if (tmp_result != MB_SUCCESS) result = tmp_result;
  }
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failure getting hexes from mesh. " );
    return result;
    }

  vtkDebugMacro(<< "Read " << all_elems.size() << " entities from MOAB.");

    // get the data set
  vtkUnstructuredGrid *ug = this->GetOutput();
  
    // create the elements
  int success = this->create_elements(iface, ug);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Problem filling in quad data. " );
    return result;
    }

  return MB_SUCCESS;
  
}

int vtkMOABReader::create_points_vertices(WriteUtilIface *iface,
                                                vtkUnstructuredGrid *&ug,
                                                const Range &verts) 
{
    // get the global id tag
  moab::Tag gid_tag;
  moab::ErrorCode result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) {
    int dum = -1;
    result = gMB->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                             gid_tag, &dum);
    if (MB_SUCCESS != result)
    {
      vtkErrorMacro( << "Couldn't create global id tag. " );
      return result;
    }
  }
  assert(0 != gid_tag);

    // allocate and fill in coordinate arrays
  std::vector<double*> coords(3);
  coords[0] = new double[verts.size()];
  coords[1] = new double[verts.size()];
  coords[2] = new double[verts.size()];
  result = iface->get_node_coords(3, verts.size(), verts,
                                  gid_tag, 0, coords);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't get nodal coordinates. " );
    return result;
  }

    // put these data into a point array
  vtkPoints *points = vtkPoints::New();
  int i = 0;
  int dum;
  points->SetNumberOfPoints(verts.size());
  for (moab::Range::const_iterator rit = verts.begin(); rit != verts.end(); rit++)
  {
    result =  gMB->tag_get_data(gid_tag, &(*rit), 1, &dum);
    assert(MB_SUCCESS == result && dum == i);
    points->SetPoint(i, coords[0][i], coords[1][i], coords[2][i]);
    i++;
  }
  ug->SetPoints(points);
  points->Delete();
  this->MaxPointId = verts.size() - 1;

  return MB_SUCCESS;
}

int vtkMOABReader::create_elements(WriteUtilIface *iface,
                                         vtkUnstructuredGrid *&ug)
{
    // get the global id tag
  moab::Tag gid_tag;
  moab::ErrorCode result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) {
    int dum = -1;
    result = gMB->tag_create(GLOBAL_ID_TAG_NAME, sizeof(int), MB_TAG_DENSE,
                             gid_tag, &dum);
    if (MB_SUCCESS != result)
    {
      vtkErrorMacro( << "Couldn't create global id tag. " );
      return result;
    }
  }
  assert(0 != gid_tag);

    // get the vertices
  Range verts;
  result = gMB->get_entities_by_type(0, MBVERTEX, verts);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't gather vertices. " );
    return result;
  }

  vtkDebugMacro(<< "Gathered " << verts.size() << " vertices from MOAB.");
  
    // assign ids to the vertices
  result = iface->assign_ids(verts, 0, 0);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't assign vertex ids. " );
    return result;
  }

    // create points/vertices in vtk database
  int success = this->create_points_vertices(iface, ug, verts);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't create points/vertices. " );
    return result;
  }

  vtkDebugMacro(<< "After create_points_vertices: ug has " << ug->GetNumberOfPoints()
                << " points, " << ug->GetNumberOfCells() << " cells.");
  
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
    result = gMB->get_entities_by_type(0, this_type, elems);
    if (MB_SUCCESS != result)
    {
      vtkErrorMacro( << "Couldn't get elements. " );
      return result;
    }

    std::vector<int> eids;
    for (Range::iterator rit = elems.begin(); rit != elems.end(); rit++) {
      
      // get the connectivity of these elements
      result = gMB->get_connectivity(*rit, connect, num_connect, true);
      if (MB_SUCCESS != result)
      {
        vtkErrorMacro( << "Couldn't get element connectivity. " );
        return result;
      }

        // get the id tag for these vertices
      result = gMB->tag_get_data(gid_tag, connect, num_connect, ids);
      if (MB_SUCCESS != result)
      {
        vtkErrorMacro( << "Couldn't get vertex ids for element. " );
        return result;
      }

        // ok, now insert this cell
      for (int i = 0; i < num_connect; i++) vtkids[i] = ids[i];
      int last_id = ug->InsertNextCell(vtk_cell_types[this_type], num_connect, vtkids);
      eids.push_back(last_id);
      if (!new_outputs) {
        assert(first || last_id == this->MaxCellId+1);
        this->MaxCellId = last_id;
      }
    }
    
    result = gMB->tag_set_data(gid_tag, elems, &eids[0]);
    if (MB_SUCCESS != result)
    {
      vtkErrorMacro( << "Couldn't save element ids. " );
      return result;
    }
    
  }
  
  return MB_SUCCESS;
}

void vtkMOABReader::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "Max point id: " << MaxPointId << std::endl;
  os << indent << "Max cell id: " << MaxCellId << std::endl;
  
  this->Superclass::PrintSelf(os,indent);
}

int vtkMOABReader::construct_filters() 
{
    // apply threshold and type filters to the output to get multiple actors
    // corresponding to dual surfaces and curves, then group the dual actors
    // together using a group filter

  vtkUnstructuredGrid *ug = this->GetOutput(0);
/*
    // first, get the non-dual mesh
  vtkExtractUnstructuredGrid *primal = vtkExtractUnstructuredGrid::New();
  primal->SetInput(this->GetOutput());
  primal->SetCellMinimum(0);
  primal->SetCellMaximum(this->MaxPrimalId);

    // set merging on so points aren't duplicated
  primal->SetMerging(1);

    // now do dual surfaces; do threshold-based extraction for now
  MBTag gid_tag;
  MBErrorCode result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  assert(MB_SUCCESS == result && 0 != gid_tag);
  
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

void vtkMOABReader::add_name(vtkUnstructuredGrid *output, const char *prefix,
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

void vtkMOABReader::Update()
{
  vtkDebugMacro("In update");;
  int i;
  
  this->UpdateInformation();
  this->UpdateData(0);
  
  for (i = 0; i < this->GetNumberOfOutputs(); i++)
    {
    if ( this->GetOutput(i) )
      {
      this->GetOutput(i)->DataHasBeenGenerated();
      }
    }
}
