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

#include <vtkTriangleFilter.h>
#include <vtkDecimatePro.h>


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



int main(int argc, char * argv[]) {
  
  if(argc < 4) {
    std::cout << "Usage: inputMeshFileName outputMeshFileName decimationFactor" << std::endl;
    return -1;
  }

  // Read inputs
  std::string inputMeshFileName = argv[1];
  std::string outputMeshFileName = argv[2];
  float decimationFactor = atof(argv[3]);

  vtkSmartPointer<vtkPolyData> inputMesh = readPolyData(inputMeshFileName);


  // Decimate
  vtkSmartPointer<vtkDecimatePro> decimateFilter = vtkSmartPointer<vtkDecimatePro>::New();
  decimateFilter->SetInputData(inputMesh);
  decimateFilter->SetTargetReduction(decimationFactor);
  decimateFilter->PreserveTopologyOn();
  decimateFilter->Update();

  vtkSmartPointer<vtkPolyData> decimatedMesh = decimateFilter->GetOutput();
  

  // Clean polydata
  vtkSmartPointer<vtkCleanPolyData> cleanPolyDataFilter = vtkSmartPointer<vtkCleanPolyData>::New();
  cleanPolyDataFilter->SetInputData(decimatedMesh);
  cleanPolyDataFilter->Update();

  std::cout << inputMesh->GetPoints()->GetNumberOfPoints() << " points reduced to " << cleanPolyDataFilter->GetOutput()->GetPoints()->GetNumberOfPoints() << std::endl;

  writePolyData(cleanPolyDataFilter->GetOutput(),outputMeshFileName);
  

  return 0;
}
