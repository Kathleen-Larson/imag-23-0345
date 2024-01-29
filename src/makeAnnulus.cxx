#include <iostream>
#include <fstream>
#include <armadillo>

#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkIdList.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkCellLocator.h>



const int nDims = 3;



vtkSmartPointer<vtkPolyData> readPolyData(std::string filename) {
  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();
  return reader->GetOutput();
}



void writePolyData(vtkSmartPointer<vtkPolyData> polyData, std::string filename) {
  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polyData);
  writer->SetFileName(filename.c_str());
  writer->Write();
}



bool isUnique(vtkIdType p, std::vector<vtkIdType> vec) {
  for(unsigned int i = 0; i < static_cast<int>(vec.size()); i++) {
    if(p == vec.at(i)) {
      return false;
    }
  }
  return true;
}



bool isInsideSurface(vtkSmartPointer<vtkPolyData> surface, float x, float y, float z, float tol) {

  vtkSmartPointer<vtkPoints> pt = vtkSmartPointer<vtkPoints>::New();
  double p0[nDims] = {x,y,z};
  pt->InsertNextPoint(p0);

  vtkSmartPointer<vtkPolyData> ptPolyData = vtkSmartPointer<vtkPolyData>::New();
  ptPolyData->SetPoints(pt);

  vtkSmartPointer<vtkSelectEnclosedPoints> selectEnclosedPointsFilter = vtkSmartPointer<vtkSelectEnclosedPoints>::New();
  selectEnclosedPointsFilter->SetSurfaceData(surface);
  selectEnclosedPointsFilter->SetInputData(ptPolyData);
  selectEnclosedPointsFilter->Update();

  if(selectEnclosedPointsFilter->IsInside(0))
    return true;
  return false;
}




int main(int argc, char * argv[]) {
  
  if(argc < 4) {
    std::cout << "Usage: surfaceOuterFileName surfaceInnerFileName annulusFileName" << std::endl;
    return -1;
  }

  // Parse args
  std::string surfaceOuterFileName = argv[1];
  std::string surfaceInnerFileName = argv[2];
  std::string annulusFileName = argv[3];


  // Read bounding surfaces
  vtkSmartPointer<vtkPolyData> surfaceOuter = readPolyData(surfaceOuterFileName);
  vtkSmartPointer<vtkPolyData> surfaceInner = readPolyData(surfaceInnerFileName);


  // Set outer labels to 1
  vtkSmartPointer<vtkFloatArray> surfaceOuterLabels = vtkSmartPointer<vtkFloatArray>::New();
  surfaceOuterLabels->SetNumberOfComponents(1);
  surfaceOuterLabels->SetNumberOfValues(surfaceOuter->GetPoints()->GetNumberOfPoints());
  surfaceOuterLabels->SetName("Surface Ids");
  for(unsigned int p = 0; p < surfaceOuter->GetPoints()->GetNumberOfPoints(); p++) {
    surfaceOuterLabels->SetValue(p,1);
  }
  surfaceOuter->GetPointData()->AddArray(surfaceOuterLabels);


  // Set inner labels to 2
  vtkSmartPointer<vtkFloatArray> surfaceInnerLabels = vtkSmartPointer<vtkFloatArray>::New();
  surfaceInnerLabels->SetNumberOfComponents(1);
  surfaceInnerLabels->SetNumberOfValues(surfaceInner->GetPoints()->GetNumberOfPoints());
  surfaceInnerLabels->SetName("Surface Ids");
  for(unsigned int p = 0; p < surfaceInner->GetPoints()->GetNumberOfPoints(); p++) {
    surfaceInnerLabels->SetValue(p,2);
  }
  surfaceInner->GetPointData()->AddArray(surfaceInnerLabels);


  // Make sure everything has normals
  vtkSmartPointer<vtkPolyDataNormals> surfaceOuterNormalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
  surfaceOuterNormalsFilter->SetInputData(surfaceOuter);
  surfaceOuterNormalsFilter->ComputePointNormalsOn();
  surfaceOuterNormalsFilter->ComputeCellNormalsOn();
  surfaceOuterNormalsFilter->FlipNormalsOff();
  surfaceOuterNormalsFilter->Update();

  vtkSmartPointer<vtkPolyDataNormals> surfaceInnerNormalsFilter = vtkSmartPointer<vtkPolyDataNormals>::New();
  surfaceInnerNormalsFilter->SetInputData(surfaceInner);
  surfaceInnerNormalsFilter->ComputePointNormalsOn();
  surfaceInnerNormalsFilter->ComputeCellNormalsOn();
  surfaceInnerNormalsFilter->FlipNormalsOff();
  surfaceInnerNormalsFilter->Update();
  

  // Combine surfaces into single polydata
  vtkSmartPointer<vtkAppendPolyData> appendPolyDataFilter = vtkSmartPointer<vtkAppendPolyData>::New();
  appendPolyDataFilter->AddInputData(surfaceOuterNormalsFilter->GetOutput());
  appendPolyDataFilter->AddInputData(surfaceInnerNormalsFilter->GetOutput());
  appendPolyDataFilter->Update();

  vtkSmartPointer<vtkCleanPolyData> cleanPolyDataFilter = vtkSmartPointer<vtkCleanPolyData>::New();
  cleanPolyDataFilter->SetInputData(appendPolyDataFilter->GetOutput());
  cleanPolyDataFilter->Update();

  vtkSmartPointer<vtkPolyData> annulus = cleanPolyDataFilter->GetOutput();
  annulus ->BuildLinks();
  writePolyData(annulus,annulusFileName);
  

  return 0;
}
