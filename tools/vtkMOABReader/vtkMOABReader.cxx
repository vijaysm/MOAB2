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

#include "vtkUnstructuredGrid.h"
#include "vtkMOABMesh.hpp"
#include "vtkObjectFactory.h"

using namespace moab;

vtkCxxRevisionMacro(vtkMOABReader, "$Revision$")
vtkStandardNewMacro(vtkMOABReader)

const bool new_outputs = false;
const bool use_filters = true;

vtkMOABReader::vtkMOABReader()
{
  this->FileName = NULL;
  iConstructedMOAB = false;
}

vtkMOABReader::~vtkMOABReader()
{
  if (iConstructedMOAB && vtkMOABMesh::instance(NULL, false))
    vtkMOABMesh::instance()->delete_instance();
}

void vtkMOABReader::Execute()
{
  this->DebugOn();

    // initialize MOAB & read file
  vtkMOABMesh *vkm = vtkMOABMesh::instance(NULL, false);
  if (NULL == vkm) {
    vkm = vtkMOABMesh::instance();
    iConstructedMOAB = true;
    vkm->SetOutput(this->GetOutput());
  }
  
  moab::ErrorCode result = MB_SUCCESS;
  if (!vkm->file_loaded(this->GetFileName())) {

    EntityHandle file_set;
    result = vkm->load_file(this->GetFileName(), NULL, file_set);
    if (MB_SUCCESS != result)
    {
      vtkErrorMacro( << "Failed to open file " << this->GetFileName() );
      return;
    }
    vtkDebugMacro(<<"Read MOAB file...");

    vkm->Update(file_set);
  }

  vtkDebugMacro(<< "After Execute: ug has " << this->GetOutput()->GetNumberOfPoints()
                << " points, " << this->GetOutput()->GetNumberOfCells() << " cells.");
}

void vtkMOABReader::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

void vtkMOABReader::Update()
{
}
