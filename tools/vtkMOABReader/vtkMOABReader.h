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
// .NAME vtkMOABReader - read vtk unstructured grid data file
// .SECTION Description
// vtkMOABReader is a source object that reads ASCII or binary 
// unstructured grid data files in vtk format. (see text for format details).
// The output of this reader is a single vtkUnstructuredGrid data object.
// The superclass of this class, vtkDataReader, provides many methods for
// controlling the reading of the data file, see vtkDataReader for more
// information.
// .SECTION Caveats
// Binary files written on one system may not be readable on other systems.
// .SECTION See Also
// vtkUnstructuredGrid vtkDataReader

#ifndef __vtkMOABReader_h
#define __vtkMOABReader_h

#include "vtkSetGet.h"
#include "vtkPoints.h"
#include "vtkUnstructuredGridSource.h"

#include "MBInterface.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBRange.hpp"
#include "DualTool.hpp"

class vtkIntArray;

#include <map>

class VTK_EXPORT vtkMOABReader : public vtkUnstructuredGridSource
{
public:
  static vtkMOABReader *New();
  vtkTypeRevisionMacro(vtkMOABReader,vtkUnstructuredGridSource);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Specify file name of the Exodus file.
  vtkSetStringMacro(FileName);
  vtkGetStringMacro(FileName);

  int GetNumberOfDualSurfaces();
  int GetNumberOfDualCurves();
  
  virtual void Update();
  
protected:
  vtkMOABReader();
  ~vtkMOABReader();

  void Execute();

private:
  vtkMOABReader(const vtkMOABReader&);  // Not implemented.
  void operator=(const vtkMOABReader&);  // Not implemented.


  static const int vtk_cell_types[];
  int NumberOfDualSurfaces;
  int NumberOfDualCurves;
  int NumberOfDualVertices;

  int MaxPointId;
  int MaxCellId;
  int MaxPrimalId;
  int DualVertexIdOffset;

  MBTag VtkOffsetIdTag;
  
  MBRange DualVertexRange;

  char *FileName;
  
  MBErrorCode construct_mesh();

  MBErrorCode create_points_vertices(MBWriteUtilIface *iface,
                                     vtkUnstructuredGrid *&ug,
                                     const MBRange &all_elems);
  
  MBErrorCode create_elements(MBWriteUtilIface *iface,
                              vtkUnstructuredGrid *&ug);
  MBErrorCode construct_dual(DualTool &dt);
  
  MBErrorCode get_vertex_polys(DualTool &dt,
                               vtkUnstructuredGrid *&ug);
  
  int get_dual_surf_polys(DualTool &dt,
                          MBEntityHandle dual_ent,
                          const int ds_id,
                          vtkIntArray *&ds_idarray,
                          vtkUnstructuredGrid *&ug);
  
  int get_dual_curve_polys(DualTool &dt,
                           MBEntityHandle dual_ent,
                           const int dc_id,
                           vtkIntArray *&dc_idarray,
                           vtkUnstructuredGrid *&ug);
  
  MBErrorCode gather_points(DualTool &dt, vtkUnstructuredGrid *&ug);

  MBErrorCode modify_pipeline(DualTool &dt);

  void add_name(vtkUnstructuredGrid *output, const char *prefix,
                const int id);
  
};

#endif


