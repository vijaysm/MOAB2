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

#include "moab/Interface.hpp"
#include "moab/WriteUtilIface.hpp"
#include "moab/Range.hpp"

class vtkIntArray;

using namespace moab;

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

  virtual void Update();
  
protected:
  vtkMOABReader();
  ~vtkMOABReader();

  void Execute();

private:
  vtkMOABReader(const vtkMOABReader&);  // Not implemented.
  void operator=(const vtkMOABReader&);  // Not implemented.


  static const int vtk_cell_types[];

  int MaxPointId;
  int MaxCellId;

  Tag VtkOffsetIdTag;
  
  char *FileName;
  
  int construct_mesh();

  int create_points_vertices(WriteUtilIface *iface,
                                   vtkUnstructuredGrid *&ug,
                                   const Range &all_elems);
  
  int create_elements(WriteUtilIface *iface,
                                  vtkUnstructuredGrid *&ug);
  
  int construct_filters();

  int read_tags();
  

  void add_name(vtkUnstructuredGrid *output, const char *prefix,
                const int id);
  
};

#endif


