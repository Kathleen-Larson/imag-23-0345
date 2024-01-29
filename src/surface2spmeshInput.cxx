#include <iostream>
#include <fstream>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkCell.h>
#include <vtkTriangle.h>
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
#include <vtkPLYWriter.h>
#include <vtkSTLWriter.h>



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
    std::cout << "Usage: surfaceOuterFileName spmeshFile_base spmeshFile_descriptor" << std::endl;
    return -1;
  }

  // Parse args
  std::string surfaceOuterFileName = argv[1];
  std::string spmeshFileName = argv[2];
  std::string spmeshFile_descriptor = argv[3];

    
  // Read surface
  vtkSmartPointer<vtkPolyData> surfaceOuter0 = readPolyData(surfaceOuterFileName);
  vtkSmartPointer<vtkCleanPolyData> cleanPolyDataFilter = vtkSmartPointer<vtkCleanPolyData>::New();
  cleanPolyDataFilter->SetInputData(surfaceOuter0);
  cleanPolyDataFilter->Update();
  vtkSmartPointer<vtkPolyData> surfaceOuter = cleanPolyDataFilter->GetOutput();
  
  unsigned int nBoundaryPts = surfaceOuter->GetPoints()->GetNumberOfPoints();
  unsigned int nCellsTotal = surfaceOuter->GetNumberOfCells();

  double bounds[nDims*2];
  surfaceOuter->GetPoints()->GetBounds(bounds);
  
  unsigned int nFacets = 0;
  for(unsigned int n = 0; n < nCellsTotal; n++) {
    vtkSmartPointer<vtkCell> cell = surfaceOuter->GetCell(n);
    if(cell->GetCellType() == 5) {
      nFacets++;
    }
  }
  

  // Create .nml file
  std::ofstream spmeshFile(spmeshFileName.c_str());

  int nChars = spmeshFile_descriptor.length();
  std::string bufferLine(nChars,42);
  spmeshFile << "#REM " << bufferLine << "\n";
  spmeshFile << "#REM " << spmeshFile_descriptor << "\n";
  spmeshFile << "#REM " << bufferLine << "\n\n";

  spmeshFile << "#MODEL_NAME MESH_DOMAIN\n\n";

  spmeshFile << "#REM (1: BDY DATA/ 2: MESH DATA/ 3:BOTH)\n";
  spmeshFile << "#INPUT_TYPE 1\n\n";

  spmeshFile << "#REM (0: NO MEMBER/ 1: WITH NUMBER)\n";
  spmeshFile << "#NUMBERING 1\n\n";

  spmeshFile << "#REM (1: THERMO/ 2: ELASTIC/ 3:FLUID)\n";
  spmeshFile << "#ANALY_TYPE 1\n\n";
  
  spmeshFile << "#DIMENSION " << nDims << "\n\n"; // # of dimensions

  spmeshFile << "#BOUND_BOX"; // bounding box of model
  for(unsigned int d = 0; d < nDims*2; d+=2) {
    spmeshFile << " " << floor(bounds[d]) - 0.0005;
  }
  for(unsigned int d = 1; d < nDims*2; d+=2) {
    spmeshFile << " " << ceil(bounds[d]) + 0.0005;
  }
  spmeshFile << "\n\n";
    
    
  spmeshFile << "#BDY_NODES " << nBoundaryPts << "\n";
  spmeshFile << "#BDY_ELEM " << nCellsTotal << "\n\n";

  spmeshFile << "#REM TITLE PARAMS: NUM_NODE, DIMENSION\n";
  spmeshFile << "#BDY_NODE_COOR " << nBoundaryPts << " " << nDims << "\n";
  for(unsigned int p = 0; p < nBoundaryPts; p++) {
    double p0[nDims];
    surfaceOuter->GetPoint(p,p0);
    spmeshFile << p+1 << " " << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
  }
  spmeshFile << "\n";

  spmeshFile << "#BDY_ELEM_LIST " << nCellsTotal << " 13\n";
  for(unsigned int f = 0; f < nFacets; f++) {
    vtkSmartPointer<vtkIdList> cellPts = vtkSmartPointer<vtkIdList>::New();
    surfaceOuter->GetCellPoints(f,cellPts);
    if(cellPts->GetNumberOfIds()==3) {
      spmeshFile << f+1;
      for(unsigned int p = 0; p < cellPts->GetNumberOfIds(); p++) {
	spmeshFile << " " << cellPts->GetId(p)+1;
      }
      spmeshFile << " 1 0\n";
    } else {
      std::cout << "cell not a triangle" << std::endl;
    }
  }
  spmeshFile << "\n";

  spmeshFile << "#END_OF_FILE\n";
  spmeshFile.close();
  
  return 0;
}
