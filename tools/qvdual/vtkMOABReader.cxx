/*
 * Copyright 2004 Sandia Corporation.
 * Under the terms of Contract DE-AC04-94AL85000, there is a non-exclusive
 * license for use of this work by or on behalf of the
 * U.S. Government. Redistribution and use in source and binary forms, with
 * or without modification, are permitted provided that this Notice and any
 * statement of authorship are reproduced on all copies.
 */
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

#include "vtkObjectFactory.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMOABUtils.h"
#include "vtkPolyDataMapper.h"
#include "vtkDataSetMapper.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkProperty.h"
#include "vtkExtractEdges.h"
#include "vtkTubeFilter.h"
#include "vtkRenderer.h"
#include "vtkExtractGeometry.h"
#include "vtkPlane.h"

#include "assert.h"

vtkCxxRevisionMacro(vtkMOABReader, "$Revision$");
vtkStandardNewMacro(vtkMOABReader);

vtkMOABReader::vtkMOABReader()
{
  this->FileName = NULL;
}

vtkMOABReader::~vtkMOABReader()
{
}

void vtkMOABReader::Execute()
{
    //this->DebugOn();

    // assert that MOAB has been initialized
  MBErrorCode result;
  assert(NULL != vtkMOABUtils::mbImpl);
  
  result = vtkMOABUtils::mbImpl->load_mesh(this->GetFileName());
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to open file " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Read MOAB file...");

    // initialize the vtk data on top of MOAB
  vtkUnstructuredGrid *ug = this->GetOutput();
  
  result = vtkMOABUtils::make_vertex_points(ug);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to make vertex points for " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Vertex points constructed.");

    // now make the cells
  result = vtkMOABUtils::make_cells(ug);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to make cells for " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Cells constructed.");

  bool tubes = true;

  vtkMOABUtils::myUG = ug;
  
    // make sure there's a mapper, actor for the whole mesh in ug, put in renderer
  vtkPolyDataMapper *poly_mapper;
  vtkDataSetMapper *set_mapper;
  vtkTubeFilter *tube_filter;
  vtkExtractEdges *edge_filter;

  vtkActor *mesh_actor = vtkActor::New();
  
  if (tubes) {
      // extract edges and build a tube filter for them
    edge_filter = vtkExtractEdges::New();
    edge_filter->SetInput(ug);
    tube_filter = vtkTubeFilter::New();
    tube_filter->SetInput(edge_filter->GetOutput());
    tube_filter->SetNumberOfSides(6);
    tube_filter->SetRadius(0.005);
    poly_mapper =  vtkPolyDataMapper::New();
    poly_mapper->SetInput(tube_filter->GetOutput());
    mesh_actor->SetMapper(poly_mapper);
    poly_mapper->ImmediateModeRenderingOn();
  }
  else {
    set_mapper = vtkDataSetMapper::New();
    set_mapper->SetInput(ug);
    mesh_actor->SetMapper(set_mapper);
    set_mapper->ImmediateModeRenderingOn();
  }
  
  vtkMOABUtils::myRen->AddActor(mesh_actor);
  result = vtkMOABUtils::mbImpl->tag_set_data(vtkMOABUtils::vtkSetActorTag, 
                                              NULL, 0, &mesh_actor);

  if (MB_SUCCESS != result) {
    vtkErrorMacro(<< "Failed to set actor for mesh in vtkMOABReader::Execute().");
    return;
  }
    
    // now turn around and set a different property for the mesh, because we want the tubes
    // to be shaded in red
  vtkMOABUtils::actorProperties[mesh_actor] = NULL;
  vtkProperty *this_prop = vtkMOABUtils::get_property(mesh_actor, true);
  this_prop->SetRepresentationToSurface();
  this_prop->SetColor(0.0, 1.0, 0.0);
  this_prop->SetEdgeColor(0.0, 1.0, 0.0);
//  mesh_actor->VisibilityOff();
  

    /*
    // center the camera on the center of the ug
  vtkMOABUtils::myRen->GetActiveCamera()->SetFocalPoint(ug->GetPoint(1));
  vtkMOABUtils::myRen->GetActiveCamera()->SetPosition(0, 0, 50.0);
  vtkMOABUtils::myRen->GetActiveCamera()->SetViewUp(0, 1.0, 0.0);
  
  std::cout << "Set focal point to " 
            << ug->GetPoint(1)[0] 
            << ", "
            << ug->GetPoint(1)[1] 
            << ", "
            << ug->GetPoint(1)[2] 
            << std::endl;
    */

  mesh_actor->Delete();
  if (tubes) {
    tube_filter->Delete();
    edge_filter->Delete();
    poly_mapper->Delete();
  }
  else {
    set_mapper->Delete();
  }
      
    // construct actors and prop assemblies for the sets
  result = vtkMOABUtils::update_all_actors(0, ug, false);
  if (MB_SUCCESS != result)
    {
    vtkErrorMacro( << "Failed to update " << this->GetFileName() );
    return;
    }
  vtkDebugMacro(<<"Set actors updated.");
  
}
