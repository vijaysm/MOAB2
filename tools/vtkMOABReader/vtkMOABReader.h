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

#include "vtkDataReader.h"
#include "vtkSetGet.h"
#include "vtkMOABReaderConfigure.h"
#include "vtkUnstructuredGrid.h"

#include "MBInterface.hpp"
#include "MBWriteUtilIface.hpp"
#include "MBRange.hpp"
#include "DualTool.hpp"

class VTK_vtkMOABReader_EXPORT vtkMOABReader : public vtkDataReader
{
public:
  static vtkMOABReader *New();
  vtkTypeRevisionMacro(vtkMOABReader,vtkDataReader);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Get the output of this reader.
  vtkUnstructuredGrid *GetOutput();
  vtkUnstructuredGrid *GetOutput(int idx)
    {return (vtkUnstructuredGrid *) this->vtkSource::GetOutput(idx); };
  void SetOutput(vtkUnstructuredGrid *output);
  
  int GetNumberOfDualSurfaces();
  int GetNumberOfDualCurves();
  
protected:
  vtkMOABReader();
  ~vtkMOABReader();

  void Execute();

  // Since the Outputs[0] has the same UpdateExtent format
  // as the generic DataObject we can copy the UpdateExtent
  // as a default behavior.
  void ComputeInputUpdateExtents(vtkDataObject *output);
  
private:
  vtkMOABReader(const vtkMOABReader&);  // Not implemented.
  void operator=(const vtkMOABReader&);  // Not implemented.


  static const int vtk_cell_types[];
  int NumberOfDualSurfaces;
  int NumberOfDualCurves;
  MBErrorCode construct_mesh();
  
  MBErrorCode create_points_vertices(MBWriteUtilIface *iface,
                                     vtkUnstructuredGrid *&ug,
                                     const MBRange &all_elems);
  
  MBErrorCode create_elements(MBWriteUtilIface *iface,
                              vtkUnstructuredGrid *&ug,
                              const MBRange &elems);
  
  MBErrorCode construct_dual();
  
  MBErrorCode get_vertex_polys(DualTool &dt,
                               vtkCellArray &vpolys);
  
  MBErrorCode get_dual_surf_polys(DualTool &dt,
                                  MBEntityHandle dual_ent,
                                  vtkCellArray &polys);
  
  MBErrorCode get_dual_curve_polys(DualTool &dt,
                                   MBEntityHandle dual_ent,
                                   vtkCellArray &polys);
  
  MBErrorCode gather_points(DualTool &dt, vtkPoints *&points);
  
};

#endif


