#include <vtkActor.h>
#include <vtkCompositeDataGeometryFilter.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataWriter.h>

#include "vtkMOABModelReader.h"
#include "vtkMOABReader.h"
#include "moab/Core.hpp"

#include "vtkXMLCompositeDataWriter.h"
#include "vtkXMLMultiBlockDataWriter.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkGeometryFilter.h"
#include "vtkCellData.h"
#include "vtkPointData.h"

#include "vtkDataSetSurfaceFilter.h"

#include "vtkConeSource.h"

int main(int argc, char* argv[])
{
  vtkSmartPointer<vtkRenderer> ren =
    vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renWin =
    vtkSmartPointer<vtkRenderWindow>::New();
  vtkSmartPointer<vtkRenderWindowInteractor> iren =
    vtkSmartPointer<vtkRenderWindowInteractor>::New();
  vtkSmartPointer<vtkInteractorStyleTrackballCamera> style =
    vtkSmartPointer<vtkInteractorStyleTrackballCamera>::New();

  vtkSmartPointer<vtkMOABReader> reader =
    vtkSmartPointer<vtkMOABReader>::New();
//  reader->SetFileName("/home/sankhesh/Projects/SiMBA/bld/meshkit/src/meshkit/rgg/cgd.h5m");
 reader->SetFileName("/home/sankhesh/Projects/SiMBA/bld/meshkit/src/meshkit/rgg/twoassm_out.h5m");
//  reader->SetFileName("/home/sankhesh/Projects/SiMBA/bld/meshkit/src/meshkit/rgg/test/vhtr-hexvertex/vhtr-hexvertex.h5m");
//  reader->Update();
//  vtkSmartPointer<vtkMultiBlockDataSet> mb = reader->GetOutput();

//  vtkSmartPointer<vtkDataSetSurfaceFilter> dsf =
//      vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
//  dsf->SetInputData(mb->GetBlock(0));

//  vtkSmartPointer<vtkGeometryFilter> geom =
//      vtkSmartPointer<vtkGeometryFilter>::New();

//  geom->SetInputData(mb->GetBlock(0));
//  geom->Update();

//  vtkSmartPointer<vtkPolyData> pd = geom->GetOutput();
//  std::cout << "After passing through geometry filter: Block 0 has " <<
//               pd->GetNumberOfPoints() << " points, " << pd->GetNumberOfCells()
//            << " cells." << std::endl;
//  std::cout << "Point Data for AKDTree has " << pd->GetPointData()->GetArray("AKDTree")->GetNumberOfTuples() << " tuples." << std::endl;
//  std::cout << "Cell Data for AKDTree has " << pd->GetCellData()->GetArray("AKDTree")->GetNumberOfTuples() << " tuples." << std::endl;

//  vtkSmartPointer<vtkPolyDataWriter> pdWriter =
//    vtkSmartPointer<vtkPolyDataWriter>::New();
//  pdWriter->SetFileName("DatasetPDcgd.vtk");
//    //  pdWriter->SetFileName("BlockPD.vtk");
//  pdWriter->SetInputConnection(dsf->GetOutputPort());
//  pdWriter->Update();

//  vtkSmartPointer<vtkXMLCompositeDataWriter> writer =
//      vtkSmartPointer<vtkXMLMultiBlockDataWriter>::New();
//  writer->SetFileName("Multiblock.vtm");
//  writer->SetInputConnection(reader->GetOutputPort());
//  writer->Write();


//  reader->Update();
//  std::cout << vtkPolyData::SafeDownCast(reader->GetOutput())->GetPoint(0)[0] << " " <<
//               vtkPolyData::SafeDownCast(reader->GetOutput())->GetPoint(0)[1] << " " <<
//               vtkPolyData::SafeDownCast(reader->GetOutput())->GetPoint(0)[2] << std::endl;
  vtkSmartPointer<vtkCompositeDataGeometryFilter> compositeGeomFilter =
    vtkSmartPointer<vtkCompositeDataGeometryFilter>::New();
  compositeGeomFilter->SetInputConnection(reader->GetOutputPort());


  vtkSmartPointer<vtkPolyDataMapper> mapper =
    vtkSmartPointer<vtkPolyDataMapper>::New();
//  mapper->SetInputConnection(reader->GetOutputPort());
  mapper->SetInputConnection(compositeGeomFilter->GetOutputPort());
//  mapper->SetInputConnection(dsf->GetOutputPort());
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->SetMapper(mapper);
//  actor->GetProperty()->SetRepresentationToSurface();
//  actor->GetProperty()->EdgeVisibilityOn();
//  actor->GetProperty()->SetPointSize(200.0);


  ren->AddActor(actor);
  renWin->AddRenderer(ren);
  iren->SetRenderWindow(renWin);
  iren->SetInteractorStyle(style);
  renWin->Render();
  iren->Initialize();
  iren->Start();

  return 0;
}
