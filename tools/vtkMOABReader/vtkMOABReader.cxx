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
#include "vtkPolyData.h"
#include "vtkObjectFactory.h"
#include "vtkUnstructuredGrid.h"

#define IS_BUILDING_MB
#include "MBCore.hpp"
#undef IS_BUILDING_MB

vtkCxxRevisionMacro(vtkMOABReader, "$Revision$");
vtkStandardNewMacro(vtkMOABReader);

MBInterface *gMB = NULL;

const int vtkMOABReader::vtk_cell_types[] = {
  1, 3, 5, 9, 7, 10, 14, 13, 0, 12, 0, 0, 0};

vtkMOABReader::vtkMOABReader()
{
  this->vtkSource::SetNthOutput(0, vtkUnstructuredGrid::New());
  // Releasing data for pipeline parallism.
  // Filters will know it is empty. 
  this->Outputs[0]->ReleaseData();
  this->Outputs[0]->Delete();
  this->NumberOfDualSurfaces = 0;
  this->NumberOfDualCurves = 0;
}

vtkMOABReader::~vtkMOABReader()
{
}

//----------------------------------------------------------------------------
vtkUnstructuredGrid *vtkMOABReader::GetOutput()
{
  if (this->NumberOfOutputs < 1)
    {
    return NULL;
    }
  
  return (vtkUnstructuredGrid *)(this->Outputs[0]);
}

//----------------------------------------------------------------------------
void vtkMOABReader::SetOutput(vtkUnstructuredGrid *output)
{
  this->vtkSource::SetNthOutput(0, output);
}


//----------------------------------------------------------------------------
// I do not think this should be here, but I do not want to remove it now.
void vtkMOABReader::ComputeInputUpdateExtents(vtkDataObject *)
{
}

void vtkMOABReader::Execute()
{

    // initialize MOAB & read file
  if (NULL == gMB) gMB = new MBCore();
  
  MBErrorCode result = gMB->load_mesh(this->GetFileName());
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to open file " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Read MOAB file...");

  result = construct_mesh();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct mesh from " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Constructed mesh...");

  result = construct_dual();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct dual from " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Constructed dual...");

}

MBErrorCode vtkMOABReader::construct_mesh() 
{
    // construct the vtk representation of the mesh; just do hexes and quads 
    // for now
  MBWriteUtilIface *iface = NULL;
  gMB->query_interface("MBWriteUtilIface", reinterpret_cast<void**>(&iface));
  assert(NULL != iface);
  
    // get all the hexes and quads
  MBRange all_elems, quads, hexes;
  MBErrorCode result = MB_SUCCESS;
//  result = gMB->get_entities_by_type(0, MBQUAD, quads);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failure getting quads from mesh. " );
    return result;
    }
  result = gMB->get_entities_by_type(0, MBHEX, hexes);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failure getting hexes from mesh. " );
    return result;
    }
//  std::copy(quads.begin(), quads.end(), mb_range_inserter(all_elems));
  std::copy(hexes.begin(), hexes.end(), mb_range_inserter(all_elems));

    // create the data set
  vtkUnstructuredGrid *ug = vtkUnstructuredGrid::New();
  
    // create the point and vertex data for the nodes
  result = this->create_points_vertices(iface, ug, all_elems);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Problem filling in point and vertex data. " );
    return result;
    }
  
    // create the quads
//  result = this->create_elements(iface, ug, quads);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Problem filling in quad data. " );
    return result;
    }
  
    // create the hexes
  result = this->create_elements(iface, ug, hexes);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Problem filling in quad data. " );
    return result;
    }

  this->AddOutput(ug);
  ug->Delete();

  return MB_SUCCESS;
  
}

MBErrorCode vtkMOABReader::create_points_vertices(MBWriteUtilIface *iface,
                                                  vtkUnstructuredGrid *&ug,
                                                  const MBRange &all_elems) 
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

    // gather the nodes used in the elements
  MBRange all_nodes;
  result = iface->gather_nodes_from_elements(all_elems, 0, all_nodes);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't gather nodes from elements. " );
    return result;
  }
  
    // allocate and fill in coordinate arrays
  std::vector<double*> coords(3);
  coords[0] = new double[all_nodes.size()];
  coords[1] = new double[all_nodes.size()];
  coords[2] = new double[all_nodes.size()];
  result = iface->get_node_arrays(3, all_nodes.size(), all_nodes,
                                  gid_tag, 0, coords);
  if (MB_SUCCESS != result)
  {
    vtkErrorMacro( << "Couldn't get nodal coordinates. " );
    return result;
  }

    // put these data into a point array
  vtkPoints *points = vtkPoints::New();
  vtkFloatArray* pcoords = vtkFloatArray::New();
  pcoords->SetNumberOfComponents(3);
  pcoords->SetNumberOfTuples(all_nodes.size());
  int i = 0;
  int dum;
  for (MBRange::iterator rit = all_nodes.begin(); rit != all_nodes.end(); rit++)
  {
    assert(MB_SUCCESS == gMB->tag_get_data(gid_tag, &(*rit), 1, &dum) &&
           dum == i);
    pcoords->SetTuple3(i, coords[0][i], coords[1][i], coords[2][i]);
    i++;
  }
  points->SetData(pcoords);
  ug->SetPoints(points);
  points->Delete();
  pcoords->Delete();

    // create vertices at the nodes; assumes nodes ordered in increasing ids
    // (which was checked in an assert during point creation)
  vtkCellArray *cells = vtkCellArray::New();
  int *idArray = cells->WritePointer(all_nodes.size(), 
                                     2*all_nodes.size());
  for (i = 0; i < all_nodes.size(); i++)
    idArray[i] = i;

    // set the cell types
  int *type_array = new int[all_nodes.size()];
  for (i = 0; i < all_nodes.size(); i++)
    type_array[i] = vtk_cell_types[MBVERTEX];
  
  ug->SetCells(type_array, cells);

  cells->Delete();
  delete [] type_array;
  
  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::create_elements(MBWriteUtilIface *iface,
                                           vtkUnstructuredGrid *&ug,
                                           const MBRange &elems) 
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

    // get the connectivity of these elements
  MBEntityType this_type = gMB->type_from_handle(*(elems.begin()));
  assert(this_type != MBMAXTYPE);
  int num_verts = MBCN::VerticesPerEntity(this_type);
  
  vtkCellArray *cells = vtkCellArray::New();
  int *idArray = cells->WritePointer(elems.size(), 
                                     (num_verts+1)*elems.size());

  result = iface->get_element_array(elems.size(), num_verts, gid_tag,
                                    elems, gid_tag, 0, idArray);
  if (MB_SUCCESS != result) return result;

    // need to shift vertex ids over one, then put
    // # vertices in at start, for each element
  int i, j;
    // set index to last element of array
  int index = (num_verts+1) * elems.size() - 1;
  for (i = 0; i < elems.size(); i++) {
      // shift vertex ids up
    for (j = 0; j < num_verts; j++) {
      idArray[index] = idArray[index-1];
      index--;
    }
    idArray[index] = num_verts;
    index--;
  }
  assert(index == -1);

    // set the cell types
  int *type_array = new int[elems.size()];
  for (i = 0; i < elems.size(); i++)
    type_array[i] = vtk_cell_types[this_type];
  
  ug->SetCells(type_array, cells);

  cells->Delete();
  delete [] type_array;

  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::construct_dual() 
{
    // construct dual
  DualTool dt(gMB);
  MBErrorCode result = dt.construct_hex_dual();
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to construct dual. ");
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
    // gather the points defining the whole dual
  vtkPoints *points = vtkPoints::New();
  gather_points(dt, points);

    // points go into the 1st output polydata
  vtkPolyData *pd = vtkPolyData::New();
  pd->SetPoints(points);

    // get the vertices of the dual and put in that same polydata
  vtkCellArray *vpolys = vtkCellArray::New();
  get_vertex_polys(dt, *vpolys);
  pd->SetVerts(vpolys);

    // now add that polydata as an output
  this->AddOutput(pd);
  vpolys->Delete();
  pd->Delete();

    // get the polygons for each dual surface and put in 
    // 2..1+NumberOfDualSurfaces polydata instances
  int ds_id;
  MBTag gid_tag;
  result = gMB->tag_get_handle(GLOBAL_ID_TAG_NAME, gid_tag);
  if (MB_SUCCESS != result) 
  {
    vtkErrorMacro( << "Failed to get global id tag from MOAB. ");
    return result;
  }
  for (MBRange::iterator rit = ds_range.begin(); rit != ds_range.end(); rit++) {
      // I *think* each dual surf is a polydata, not sure
    vtkPolyData *pd = vtkPolyData::New();
    result = gMB->tag_get_data(gid_tag, &(*rit), 1, &ds_id);

      // get the polygons which represent this dual surf
    vtkCellArray *polys = vtkCellArray::New();
    get_dual_surf_polys(dt, *rit, *polys);
      
    pd->SetPoints(points);
    pd->SetPolys(polys);
      
    this->AddOutput(pd);
    polys->Delete();
    pd->Delete();
  }

    // get the polygons for each dual curve and put in 
    // 2..NumberOfDualSurfaces+1..2+NumberOfDualSurfaces+NumberOfDualCurves 
    // polydata instances
  int dc_id;
  for (MBRange::iterator rit = dc_range.begin(); rit != dc_range.end(); rit++) {
      // I *think* each dual curve is a polydata, not sure
    vtkPolyData *pd = vtkPolyData::New();
    result = gMB->tag_get_data(gid_tag, &(*rit), 1, &dc_id);

      // get the polygons which represent this dual surf
    vtkCellArray *polys = vtkCellArray::New();
    get_dual_curve_polys(dt, *rit, *polys);
      
    pd->SetPoints(points);
    pd->SetLines(polys);
      
    this->AddOutput(pd);
    polys->Delete();
    pd->Delete();
  }

  points->Delete();

  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::get_vertex_polys(DualTool &dt,
                                             vtkCellArray &vpolys) 
{
    // make a single polydata for all dual vertices
  MBTag dcell_tag = dt.isDualCell_tag();
  assert(0 != dcell_tag);

    // get the vertices
  int dum_tag = 0x1;
  int *dum_tag_ptr = &dum_tag;
  MBRange dum_range;
  MBErrorCode result = gMB->get_entities_by_type_and_tag(0, MBVERTEX, &dcell_tag,
                                                         (const void **)&dum_tag_ptr, 1, dum_range);
  if (MB_SUCCESS != result) return result;

    // get the graphics points on the vertices
  std::vector<DualTool::GraphicsPoint> points;
  result = dt.get_graphics_points(dum_range, points);
  if (MB_SUCCESS != result) return result;
  
    // add each graphics point to the cell array
  for (std::vector<DualTool::GraphicsPoint>::iterator vit = points.begin(); 
       vit != points.end(); vit++)
    vpolys.InsertNextCell(1, &(vit->id));

  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::get_dual_surf_polys(DualTool &dt,
                                                MBEntityHandle dual_ent,
                                                vtkCellArray &polys) 
{
    // should be an entity set
  assert(gMB->type_from_handle(dual_ent) == MBENTITYSET);

  MBRange two_cells;
  MBErrorCode result = gMB->get_entities_by_handle(dual_ent, two_cells);
  if (MB_SUCCESS != result) return result;
  
  std::vector<DualTool::GraphicsPoint> points;
  std::vector<int> npts;
  
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
      polys.InsertNextCell(*nit);
      for (int i = 0; i < *nit; i++)
        polys.InsertCellPoint(vit++->id);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::get_dual_curve_polys(DualTool &dt,
                                                MBEntityHandle dual_ent,
                                                vtkCellArray &polys) 
{
    // should be an entity set
  assert(gMB->type_from_handle(dual_ent) == MBENTITYSET);

  MBRange one_cells;
  MBErrorCode result = gMB->get_entities_by_handle(dual_ent, one_cells);
  if (MB_SUCCESS != result) return result;
  
  std::vector<DualTool::GraphicsPoint> points;
  std::vector<int> npts;
  
  for (MBRange::iterator rit = one_cells.begin(); rit != one_cells.end(); rit++) {
      // get the cells for this 1cell
    points.clear();
    npts.clear();
    result = dt.get_graphics_points(*rit, npts, points);
    if (MB_SUCCESS != result) return result;
  
    // add each polyline to the cell array
    std::vector<DualTool::GraphicsPoint>::iterator vit = points.begin(); 
    for (std::vector<int>::iterator nit = npts.begin(); nit != npts.end(); 
         nit++) {
      polys.InsertNextCell(*nit);
      for (int i = 0; i < *nit; i++)
        polys.InsertCellPoint(vit++->id);
    }
  }
  
  return MB_SUCCESS;
}

MBErrorCode vtkMOABReader::gather_points(DualTool &dt, vtkPoints *&points)
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
  result = dt.get_graphics_points(dum_range, gps, true);
  if (MB_SUCCESS != result) return result;

    // put these points: put into a vtkPoints
    // 
  vtkFloatArray* pcoords = vtkFloatArray::New();
  pcoords->SetNumberOfComponents(3);
  pcoords->SetNumberOfTuples(gps.size());
  std::vector<DualTool::GraphicsPoint>::iterator pit;
  i = 0;
  for (pit = gps.begin(); pit != gps.end(); pit++)
  {
    assert(gps[i].id == i);
    pcoords->SetTuple3(i, gps[i].xyz[0], gps[i].xyz[1], gps[i].xyz[2]);
    i++;
  }
  points->SetData(pcoords);
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

