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
#include <sstream>

#define IS_BUILDING_MB
#include "MBCore.hpp"
#undef IS_BUILDING_MB

#define DUAL_SURF_ATTRIBUTE_NAME "DualSurfaceId"
#define DUAL_CURVE_ATTRIBUTE_NAME "DualCurveId"

vtkCxxRevisionMacro(vtkMOABReader, "$Revision$");
vtkStandardNewMacro(vtkMOABReader);

MBInterface *gMB = NULL;

const int vtkMOABReader::vtk_cell_types[] = {
  1, 3, 5, 9, 7, 10, 14, 13, 0, 12, 0, 0, 0};

const bool new_outputs = false;

vtkMOABReader::vtkMOABReader()
{
  this->NumberOfDualSurfaces = 0;
  this->NumberOfDualCurves = 0;
  this->FileName = NULL;
  VtkOffsetIdTag = 0;
  
}

vtkMOABReader::~vtkMOABReader()
{
}

void vtkMOABReader::Execute()
{
  this->DebugOn();

    // initialize MOAB & read file
  if (NULL == gMB) gMB = new MBCore();
  
    // make an offset id tag
  MBErrorCode result;
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
  
  result = construct_mesh();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct mesh from " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Constructed mesh...");

  DualTool dt(gMB);
  result = construct_dual(dt);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct dual from " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Constructed dual...");

//  result = modify_pipeline(dt);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to modify pipeline with threshold and group filters ");
    return;
    }
  vtkDebugMacro(<<"Modified pipeline...");
  
}

MBErrorCode vtkMOABReader::construct_mesh() 
{
    // construct the vtk representation of the mesh; just do hexes and quads 
    // for now
  MBWriteUtilIface *iface = NULL;
  gMB->query_interface("MBWriteUtilIface", reinterpret_cast<void**>(&iface));
  assert(NULL != iface);
  
    // get all the hexes and quads
  MBRange all_elems;
  MBErrorCode result = MB_SUCCESS, tmp_result;
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
  result = this->create_elements(iface, ug);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Problem filling in quad data. " );
    return result;
    }

  this->MaxPrimalId = this->MaxCellId;

  return MB_SUCCESS;
  
}

MBErrorCode vtkMOABReader::create_points_vertices(MBWriteUtilIface *iface,
                                                  vtkUnstructuredGrid *&ug,
                                                  const MBRange &verts) 
{
    // get the global id tag
  MBTag gid_tag;
  MBErrorCode result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
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
  result = iface->get_node_arrays(3, verts.size(), verts,
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
  for (MBRange::const_iterator rit = verts.begin(); rit != verts.end(); rit++)
  {
    result =  gMB->tag_get_data(gid_tag, &(*rit), 1, &dum);
    assert(MB_SUCCESS == result && dum == i);
    points->SetPoint(i, coords[0][i], coords[1][i], coords[2][i]);
    i++;
  }
  ug->SetPoints(points);
  points->Delete();
  this->MaxPointId = verts.size() - 1;

    // create vertices at the nodes; assumes nodes ordered in increasing ids
    // (which was checked in an assert during point creation)
  int last_pt;
  for (i = 0; i < verts.size(); i++) {
    last_pt = ug->InsertNextCell(vtk_cell_types[MBVERTEX], 1, &i);
  }

  if (!new_outputs) {
    assert(last_pt == verts.size()-1 && ug->GetNumberOfCells() == verts.size());
    this->MaxCellId = verts.size()-1;
  }

  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::create_elements(MBWriteUtilIface *iface,
                                           vtkUnstructuredGrid *&ug)
{
    // get the global id tag
  MBTag gid_tag;
  MBErrorCode result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
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
  MBRange verts;
  result = gMB->get_entities_by_type(0, MBVERTEX, verts);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't gather vertices. " );
    return result;
  }

  vtkDebugMacro(<< "Gathered " << verts.size() << " vertices from MOAB.");
  
    // assign ids to the vertices; keep track of how many so we can set point ids
    // for dual vertices later
  result = iface->assign_ids(verts, 0, 0);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't assign vertex ids. " );
    return result;
  }

    // create points/vertices in vtk database
  result = this->create_points_vertices(iface, ug, verts);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't create points/vertices. " );
    return result;
  }

  vtkDebugMacro(<< "After create_points_vertices: ug has " << ug->GetNumberOfPoints()
                << " points, " << ug->GetNumberOfCells() << " cells.");
  
    // for the remaining elements, add them individually
  int ids[31];
  const MBEntityHandle *connect;
  int num_connect;

  for (MBEntityType this_type = MBEDGE; this_type != MBENTITYSET; this_type++) {

      // don't try to represent elements vtk doesn't understand
    if (vtk_cell_types[this_type] == 0) continue;
    
    MBRange elems;
    result = gMB->get_entities_by_type(0, this_type, elems);
    if (MB_SUCCESS != result)
    {
      vtkErrorMacro( << "Couldn't get elements. " );
      return result;
    }
  
    for (MBRange::iterator rit = elems.begin(); rit != elems.end(); rit++) {
      
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
      int last_id = ug->InsertNextCell(vtk_cell_types[this_type], num_connect, ids);
      if (!new_outputs) {
        assert(last_id == this->MaxCellId+1);
        this->MaxCellId = last_id;
      }
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::construct_dual(DualTool &dt) 
{
    // construct dual
  MBErrorCode result = dt.construct_hex_dual();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct dual. ");
    assert(false);
    return result;
    }

    // find # dual surfaces & curves
  MBTag ds_tag = dt.dualSurface_tag();
  MBRange ds_range;
  result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                             &ds_tag, NULL, 1, ds_range);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to get number of dual surfaces. ");
    return result;
    }
  this->NumberOfDualSurfaces = ds_range.size();
  
  MBTag dc_tag = dt.dualCurve_tag();
  MBRange dc_range;
  result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                             &dc_tag, NULL, 1, dc_range);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to get number of dual curves. ");
    return result;
    }
  this->NumberOfDualCurves = dc_range.size();
  
  MBEntityHandle dual_ent;

    // get the data set
  vtkUnstructuredGrid *ug = this->GetOutput();
  
    // gather the points defining the whole dual
  gather_points(dt, ug);

    // get the vertices of the dual and put in the ug
  get_vertex_polys(dt, ug);

    // get the global id tag
  int ds_id;
  MBTag gid_tag;
  result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) 
  {
    vtkErrorMacro( << "Failed to get global id tag from MOAB. ");
    return result;
  }

  vtkIntArray *ds_idarray;
  if (!new_outputs) {
      // set up a field data array corresponding to the dual surface id;
      // will set that on each polygon to use as scalar data later
    ds_idarray = vtkIntArray::New();
    ds_idarray->SetName(DUAL_SURF_ATTRIBUTE_NAME);
    for (int i = this->MaxCellId; i >= 0; i--)
      ds_idarray->InsertValue(i, -1);
  }

    // get the polygons for each dual surface and put in ug
  for (MBRange::iterator rit = ds_range.begin(); rit != ds_range.end(); rit++) {
    result = gMB->tag_get_data(gid_tag, &(*rit), 1, &ds_id);

    vtkUnstructuredGrid *surf_ug;
    
    if (new_outputs) {
        // make a new ug for this dual surf
      surf_ug = vtkUnstructuredGrid::New();
      this->AddOutput(surf_ug);
      this->add_name(surf_ug, "dual_surf_", ds_id);
        // make the polygons which represent this dual surf; returns the id of the first polygon
        // set up a field data array corresponding to the dual surface id;
        // will set that on each polygon to use as scalar data later
      ds_idarray = vtkIntArray::New();
      ds_idarray->SetName(DUAL_SURF_ATTRIBUTE_NAME);

    }
    else 
      surf_ug = ug;

    int first = get_dual_surf_polys(dt, *rit, ds_id, ds_idarray, surf_ug);
    assert(-1 != first);

      // assign the starting polygon id in the dual surface set - id map
    result = gMB->tag_set_data(VtkOffsetIdTag, &(*rit), 1, &first);
    assert(MB_SUCCESS == result);

    if (new_outputs) {
      surf_ug->SetPoints(ug->GetPoints());
      surf_ug->GetCellData()->AddArray(ds_idarray);
      ds_idarray->Delete();
      surf_ug->Delete();
    }
  }

    // get the polygons for each dual curve and put in 
    // 2..NumberOfDualSurfaces+1..2+NumberOfDualSurfaces+NumberOfDualCurves 
    // polydata instances
  int dc_id;
  int max_ds_id = this->MaxCellId;
  vtkIntArray *dc_idarray;
  if (!new_outputs) {
      // set up a field data array corresponding to the dual curve id;
      // will set that on each polyline to use as scalar data later
    dc_idarray = vtkIntArray::New();
    dc_idarray->SetName(DUAL_CURVE_ATTRIBUTE_NAME);
      // set dual curve ids to default value for dual surf polygons
    for (int i = this->MaxCellId; i >=0; i--)
      dc_idarray->InsertValue(i, -1);
  }
  
  for (MBRange::iterator rit = dc_range.begin(); rit != dc_range.end(); rit++) {
    result = gMB->tag_get_data(gid_tag, &(*rit), 1, &dc_id);

    vtkUnstructuredGrid *curve_ug;
    if (new_outputs) {
      
        // make a new ug for this dual curve
      curve_ug = vtkUnstructuredGrid::New();
      this->AddOutput(curve_ug);
      this->add_name(curve_ug, "dual_curve_", dc_id);

        // set up a field data array corresponding to the dual curve id;
        // will set that on each polyline to use as scalar data later
      dc_idarray = vtkIntArray::New();
      dc_idarray->SetName(DUAL_CURVE_ATTRIBUTE_NAME);
    }
    else 
      curve_ug = ug;
      
    
      // make the polylines which represent this dual curve; returns the id of the first polyline
    int first = get_dual_curve_polys(dt, *rit, dc_id, dc_idarray, curve_ug);
    assert(-1 != first);

      // assign the starting polygon id in the dual surface set - id map
    result = gMB->tag_set_data(VtkOffsetIdTag, &(*rit), 1, &first);

    if (new_outputs) {
      curve_ug->SetPoints(ug->GetPoints());
      curve_ug->GetCellData()->AddArray(dc_idarray);
      dc_idarray->Delete();
      curve_ug->Delete();
    }
  }

  if (!new_outputs) {
      // set ds_idarray to default value for dual curve polylines
    for (int i = this->MaxCellId; i > max_ds_id; i--)
      ds_idarray->InsertValue(i, -1);

    ug->GetCellData()->AddArray(dc_idarray);
    ug->GetCellData()->AddArray(ds_idarray);
    ds_idarray->Delete();
    dc_idarray->Delete();
  }

  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::get_vertex_polys(DualTool &dt,
                                            vtkUnstructuredGrid *&ug)
{
  MBTag dcell_tag = dt.isDualCell_tag();
  assert(0 != dcell_tag);

    // get the vertices
  int dum_tag = 0x1;
  int *dum_tag_ptr = &dum_tag;
  DualVertexRange.clear();
  MBErrorCode result = gMB->get_entities_by_type_and_tag(0, MBVERTEX, &dcell_tag,
                                                         (const void **)&dum_tag_ptr, 1, DualVertexRange);
  if (MB_SUCCESS != result) return result;

    // get the graphics points on the vertices
  std::vector<DualTool::GraphicsPoint> points;
  result = dt.get_graphics_points(DualVertexRange, points);
  if (MB_SUCCESS != result) return result;

    // make another output for the vertices
  vtkUnstructuredGrid *vug;
  if (new_outputs) {
    vug = vtkUnstructuredGrid::New();
    this->AddOutput(vug);
    this->add_name(vug, "dual_vertices", 0);
  }
  else
    vug = ug;
  
    // need to keep the offset between 1st MOAB vertex and corresponding vtk element id
  std::vector<DualTool::GraphicsPoint>::const_iterator vit = points.begin(); 
  this->DualVertexIdOffset = vug->InsertNextCell(vtk_cell_types[MBVERTEX], 1, const_cast<int*>(&(vit->id)));
  if (!new_outputs) {
    assert(this->MaxCellId+1 == this->DualVertexIdOffset);
    this->MaxCellId = this->DualVertexIdOffset;
  }
  
  this->NumberOfDualVertices = points.size();

  int last;
  
    // add each graphics point to the cell array
  for (vit++; vit != points.end(); vit++) {
    last = vug->InsertNextCell(vtk_cell_types[MBVERTEX], 1, const_cast<int*>(&(vit->id)));
    if (!new_outputs) {
      assert(this->MaxCellId+1 == last);
      this->MaxCellId = last;
    }
  }

  if (new_outputs) {
      // borrow the points array from the reader output
    vug->SetPoints(ug->GetPoints());
    vug->Delete();
  }
  
  return MB_SUCCESS;
}

int vtkMOABReader::get_dual_surf_polys(DualTool &dt,
                                       MBEntityHandle dual_ent,
                                       const int ds_id,
                                       vtkIntArray *&ds_idarray,
                                       vtkUnstructuredGrid *&ug) 
{
    // should be an entity set
  assert(gMB->type_from_handle(dual_ent) == MBENTITYSET);

  MBRange two_cells;
  MBErrorCode result = gMB->get_entities_by_handle(dual_ent, two_cells);
  if (MB_SUCCESS != result) return result;
  
  std::vector<DualTool::GraphicsPoint> points;
  std::vector<int> npts;

  int first = -1, tmp_first;
  std::vector<int> pt_ids;
  
  for (MBRange::iterator rit = two_cells.begin(); rit != two_cells.end(); rit++) {
      // get the cells for this 2cell
    points.clear();
    npts.clear();
    result = dt.get_graphics_points(*rit, npts, points);
    if (MB_SUCCESS != result) return result;
  
    // add each polygon to the cell array
    std::vector<DualTool::GraphicsPoint>::iterator vit = points.begin(); 
    for (std::vector<int>::iterator nit = npts.begin(); nit != npts.end(); 
         nit++) {
        // first need to make a list of point ids
      pt_ids.reserve(*nit);
      for (int i = 0; i < *nit; i++)
        pt_ids[i] = vit++->id;

#ifdef NDEBUG
      for (int i = 0; i < *nit; i++)
        assert(pt_ids[i] < this->MaxPointId && pt_ids[i] >= 0);
#endif
      
        // now make the actual polygon
      tmp_first = ug->InsertNextCell(vtk_cell_types[MBPOLYGON], *nit, &pt_ids[0]);
      if (-1 == first) first = tmp_first;

        // put the dual surface id as a scalar on the data array
      ds_idarray->InsertValue(tmp_first, ds_id);
      if (!new_outputs) {
        assert(this->MaxCellId+1 == tmp_first);
        this->MaxCellId = tmp_first;
      }
    }
  }

  return first;
}

int vtkMOABReader::get_dual_curve_polys(DualTool &dt,
                                        MBEntityHandle dual_ent,
                                        const int dc_id,
                                        vtkIntArray *&dc_idarray,
                                        vtkUnstructuredGrid *&ug) 
{
    // should be an entity set
  assert(gMB->type_from_handle(dual_ent) == MBENTITYSET);

  MBRange one_cells;
  MBErrorCode result = gMB->get_entities_by_handle(dual_ent, one_cells);
  if (MB_SUCCESS != result) return result;
  
  std::vector<DualTool::GraphicsPoint> points;
  std::vector<int> npts;

  int first = -1, tmp_first;
  std::vector<int> pt_ids;
  
  for (MBRange::iterator rit = one_cells.begin(); rit != one_cells.end(); rit++) {
      // get the cells for this 2cell
    points.clear();
    npts.clear();
    result = dt.get_graphics_points(*rit, npts, points);
    if (MB_SUCCESS != result) return result;
  
    // add each polyline to the cell array
    std::vector<DualTool::GraphicsPoint>::iterator vit = points.begin(); 
    for (std::vector<int>::iterator nit = npts.begin(); nit != npts.end(); 
         nit++) {
        // first need to make a list of point ids
      pt_ids.reserve(*nit);
      for (int i = 0; i < *nit; i++)
        pt_ids[i] = vit++->id;

#ifdef NDEBUG
      for (int i = 0; i < *nit; i++)
        assert(pt_ids[i] < this->MaxPointId && pt_ids[i] >= 0);
#endif
        // now make the actual polyline; hardwire the polyline type
      tmp_first = ug->InsertNextCell(VTK_POLY_LINE, *nit, &pt_ids[0]);
      if (-1 == first) first = tmp_first;

        // put the dual curve id as a scalar on the data array
      dc_idarray->InsertValue(tmp_first, dc_id);
      if (!new_outputs) {
        assert(this->MaxCellId+1 == tmp_first);
        this->MaxCellId = tmp_first;
      }
    }
  }
  
  return first;
}

MBErrorCode vtkMOABReader::gather_points(DualTool &dt, vtkUnstructuredGrid *&ug)
{
  int i;
  MBRange dum_range;
  MBEntityHandle dual_ent;
  std::vector<DualTool::GraphicsPoint> gps;

    // get the dual surface sets
  MBTag ds_tag = dt.dualSurface_tag();
  MBErrorCode result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                                         &ds_tag, NULL, 1, dum_range);
  if (MB_SUCCESS != result) return result;

    // get the points defined on all the entities in these sets, re-assigning
    // ids at the same time
  result = dt.get_graphics_points(dum_range, gps, true, MaxPointId+1);
  if (MB_SUCCESS != result) return result;

    // get the point array
  vtkPoints *points = ug->GetPoints();

    // increase point array size to hold these
  assert(points->GetNumberOfPoints() == MaxPointId+1);
  
  std::vector<DualTool::GraphicsPoint>::reverse_iterator pit;
    // iterate from reverse, so we only re-allocate on the point list once
  for (pit = gps.rbegin(); pit != gps.rend(); pit++) {
    points->InsertPoint(pit->id, pit->xyz);
    if (pit->id > this->MaxPointId) this->MaxPointId = pit->id;
  }

  assert(points->GetNumberOfPoints() == MaxPointId+1);
}

void vtkMOABReader::PrintSelf(ostream& os, vtkIndent indent)
{
  os << indent << "Number of dual surfaces: " << NumberOfDualSurfaces << std::endl;
  os << indent << "Number of dual curves: " << NumberOfDualCurves << std::endl;
  
  this->Superclass::PrintSelf(os,indent);
}

int vtkMOABReader::GetNumberOfDualSurfaces()
{
  return NumberOfDualSurfaces;
}

int vtkMOABReader::GetNumberOfDualCurves() 
{
  return NumberOfDualCurves;
}

MBErrorCode vtkMOABReader::modify_pipeline(DualTool &dt) 
{
    // apply threshold and type filters to the output to get multiple actors
    // corresponding to dual surfaces and curves, then group the dual actors
    // together using a group filter

  this->SetNumberOfOutputs(3 + this->NumberOfDualSurfaces + this->NumberOfDualCurves);

  vtkUnstructuredGrid *ug = this->GetOutput(0);

    // first, get the non-dual mesh
  vtkExtractUnstructuredGrid *primal = vtkExtractUnstructuredGrid::New();
  primal->SetInput(this->GetOutput());
  this->SetNthOutput(1, primal->GetOutput());
  primal->SetCellMinimum(0);
  primal->SetCellMaximum(this->MaxPrimalId);

    // add the primal as an output
    //primal->Delete();
  
    // now do dual surfaces; do threshold-based extraction for now
  MBTag gid_tag;
  MBErrorCode result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  assert(MB_SUCCESS == result && 0 != gid_tag);
  
  int ds_id;
  MBTag ds_tag = dt.dualSurface_tag();
  MBRange ds_range;
  result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                             &ds_tag, NULL, 1, ds_range);
  for (MBRange::iterator rit = ds_range.begin(); rit != ds_range.end(); rit++) {
    result = gMB->tag_get_data(gid_tag, &(*rit), 1, &ds_id);

    vtkThreshold *ds_filter = vtkThreshold::New();
    ds_filter->SelectInputScalars(DUAL_SURF_ATTRIBUTE_NAME);
    ds_filter->SetAttributeModeToUseCellData();
    ds_filter->ThresholdBetween(((double)ds_id-0.5), ((double)ds_id+0.5));
    ds_filter->SetInput(ug);
    this->add_name(ds_filter->GetOutput(), "dual_surf_", ds_id);
    this->SetNthOutput(ds_id+2, ds_filter->GetOutput());
      //ds_filter->Delete();
  }
  
    // same for dual curves
  int dc_id;
  MBTag dc_tag = dt.dualSurface_tag();
  MBRange dc_range;
  result = gMB->get_entities_by_type_and_tag(0, MBENTITYSET, 
                                             &dc_tag, NULL, 1, dc_range);
  for (MBRange::iterator rit = dc_range.begin(); rit != dc_range.end(); rit++) {
    result = gMB->tag_get_data(gid_tag, &(*rit), 1, &dc_id);

    vtkThreshold *dc_filter = vtkThreshold::New();
    dc_filter->SelectInputScalars(DUAL_CURVE_ATTRIBUTE_NAME);
    dc_filter->SetAttributeModeToUseCellData();
    dc_filter->ThresholdBetween(((double)dc_id-0.5), ((double)dc_id+0.5));
    dc_filter->SetInput(ug);
    this->add_name(dc_filter->GetOutput(), "dual_curve_", dc_id);
    this->SetNthOutput(this->NumberOfDualSurfaces + 2 + dc_id, 
                       dc_filter->GetOutput());
      //dc_filter->Delete();
  }

    // lastly, get the dual vertices and put those in a group
    // first, get the non-dual mesh
  vtkExtractUnstructuredGrid *dual_verts = vtkExtractUnstructuredGrid::New();
  dual_verts->SetCellMinimum(this->DualVertexIdOffset);
  dual_verts->SetCellMaximum(this->DualVertexIdOffset+this->NumberOfDualVertices-1);
  this->add_name(dual_verts->GetOutput(), "dual_verts", 0);
  this->SetNthOutput(this->NumberOfDualSurfaces + this->NumberOfDualCurves +
                     2, dual_verts->GetOutput());
    //dual_verts->Delete();

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
