#include <iostream>
#include <fstream>
#include <armadillo>

#include <vtkMath.h>
#include <vtkFloatArray.h>
#include <vtkIdList.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellLocator.h>
#include <vtkTetra.h>
#include <vtkLine.h>


const int nDims = 3;



void writeUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> ugrid, std::string filename) {
  vtkSmartPointer<vtkUnstructuredGridWriter> uGridWriter = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  uGridWriter->SetInputData(ugrid);
  uGridWriter->SetFileName(filename.c_str());
  uGridWriter->Write();
}



vtkSmartPointer<vtkPolyData> readPolyData(std::string filename) {

  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  return reader->GetOutput();
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
  selectEnclosedPointsFilter->SetTolerance(tol);
  selectEnclosedPointsFilter->Update();

  if(selectEnclosedPointsFilter->IsInside(0))
    return true;
  return false;
}



bool isUnique(vtkIdType p, std::vector<vtkIdType> vec) {
  for(unsigned int i = 0; i < static_cast<int>(vec.size()); i++) {
    if(p == vec.at(i)) {
      return false;
    }
  }
  return true;
}



int main(int argc, char * argv[]) {

  // Parse args
  if(argc < 7) {
    std::cout << "Usage: spmeshFileName_base annulusFileName outerSurfaceFileName innerSurfaceFileName voroFileName outputUGridFileName" << std::endl;
    return -1;
  }

  std::string spmeshFileName_base = argv[1];
  std::string annulusFileName = argv[2];
  std::string outerSurfaceFileName = argv[3];
  std::string innerSurfaceFileName = argv[4];
  std::string voroFileName = argv[5];
  std::string outputUGridFileName = argv[6];

  vtkSmartPointer<vtkPolyData> annulus = readPolyData(annulusFileName);
  unsigned int nAnnulusPts = annulus->GetPoints()->GetNumberOfPoints();
  vtkSmartPointer<vtkFloatArray> annulusIDs = vtkFloatArray::SafeDownCast(annulus->GetPointData()->GetArray("Surface Ids"));

  vtkSmartPointer<vtkPolyData> outerSurface = readPolyData(outerSurfaceFileName);
  vtkSmartPointer<vtkPolyData> innerSurface = readPolyData(innerSurfaceFileName);
  
  std::string spmeshFileName_node = spmeshFileName_base + ".nod";
  std::string spmeshFileName_elem = spmeshFileName_base + ".elm";

  
  // Read spmesh node file and set points
  vtkSmartPointer<vtkPoints> spmeshPts = vtkSmartPointer<vtkPoints>::New();
  std::fstream spmeshFile_node(spmeshFileName_node.c_str());
  std::string line;
  unsigned int line_count = 0;

  while(std::getline(spmeshFile_node,line)){
    std::istringstream line_stream(line);
    if(line_count > 0) {
      unsigned int n;
      float x, y, z;
      line_stream >> n >> x >> y >> z;

      double p0[nDims] = {x,y,z};
      spmeshPts->InsertPoint(n-1,p0);
    }
    
    line_count++;
  }  
  unsigned int nSpmeshPts = spmeshPts->GetNumberOfPoints();

  
  // Remove incorrect points
  vtkSmartPointer<vtkCellLocator> annulusCellLocator = vtkSmartPointer<vtkCellLocator>::New();
  annulusCellLocator->SetDataSet(annulus);
  annulusCellLocator->BuildLocator();
  
  vtkSmartPointer<vtkCellLocator> outerSurfaceCellLocator = vtkSmartPointer<vtkCellLocator>::New();
  outerSurfaceCellLocator->SetDataSet(outerSurface);
  outerSurfaceCellLocator->BuildLocator();

  vtkSmartPointer<vtkCellLocator> innerSurfaceCellLocator = vtkSmartPointer<vtkCellLocator>::New();
  innerSurfaceCellLocator->SetDataSet(innerSurface);
  innerSurfaceCellLocator->BuildLocator();
  
  vtkSmartPointer<vtkPoints> outputUGridPts = vtkSmartPointer<vtkPoints>::New();
  arma::vec spmesh2ugrid_ptCorrespondance(nSpmeshPts);
  unsigned int UGridPtCount = 0;
  unsigned int nancount = 0;
  
  for(unsigned int p = 0; p < nSpmeshPts; p++) {
    double p0[nDims];
    spmeshPts->GetPoint(p,p0);
    bool isValidPt = false;
    
    // Find closest cell on annulus (input to spmesh)
    double closestPoint[nDims], closestPointDist2;
    vtkIdType cellId;
    int subId;
    annulusCellLocator->FindClosestPoint(p0,closestPoint,cellId,subId,closestPointDist2);

    if(sqrt(closestPointDist2) < 0.1) { // if right next to an annulus cell (aka on the boundary)
      isValidPt = true;
    }
    else { // figure out if it's inside the annular ring or not
      double closestPoint_outer[nDims], closestPointDist2_outer;
      vtkIdType cellId_outer;
      outerSurfaceCellLocator->FindClosestPoint(p0,closestPoint_outer,cellId_outer,subId,closestPointDist2_outer);
    
      double closestPoint_inner[nDims], closestPointDist2_inner;
      vtkIdType cellId_inner;
      innerSurfaceCellLocator->FindClosestPoint(p0,closestPoint_inner,cellId_inner,subId,closestPointDist2_inner);

      arma::vec vec2outer(nDims), vec2inner(nDims);
      for(unsigned int d = 0; d < nDims; d++) {
	vec2outer(d) = p0[d] - closestPoint_outer[d];
	vec2inner(d) = p0[d] - closestPoint_inner[d];
      }

      // if the point is inside annulus ring, vectors should be pointing away from each other
      double angle = acos(arma::dot(vec2outer,vec2inner) / (sqrt(arma::sum(arma::square(vec2outer))) * sqrt(arma::sum(arma::square(vec2inner)))));
      if(abs(angle) > arma::datum::pi/4) {
	isValidPt = true;
      }
    }
    //isValidPt = true;
    if(isValidPt) {
      outputUGridPts->InsertNextPoint(p0);
      spmesh2ugrid_ptCorrespondance(p) = UGridPtCount;
      UGridPtCount++;
    } else {
      spmesh2ugrid_ptCorrespondance(p) = arma::datum::nan;
      nancount++;
    }
  }
  unsigned int nUGridPts = outputUGridPts->GetNumberOfPoints();
  vtkSmartPointer<vtkUnstructuredGrid> outputUGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  outputUGrid->SetPoints(outputUGridPts);

  
  // Now set boundary labels
  vtkSmartPointer<vtkFloatArray> boundaryIDs = vtkSmartPointer<vtkFloatArray>::New();
  boundaryIDs->SetNumberOfComponents(1);
  boundaryIDs->SetNumberOfTuples(nUGridPts);
  boundaryIDs->SetName("Boundary Ids");
    
  for(unsigned int p = 0; p < nUGridPts; p++) {
     double p0[nDims];
    outputUGridPts->GetPoint(p,p0);
    double closestPoint[nDims], closestPointDist2;
    vtkIdType cellId;
    int subId;
    annulusCellLocator->FindClosestPoint(p0,closestPoint,cellId,subId,closestPointDist2);

    if(sqrt(closestPointDist2) < 0.01) {
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
      annulus->GetCellPoints(cellId,ptIds);

      unsigned int label = annulusIDs->GetValue(ptIds->GetId(0));
      boundaryIDs->SetValue(p,label);
    } else {
      boundaryIDs->SetValue(p,0);
    }
  }
  outputUGrid->GetPointData()->AddArray(boundaryIDs);
  outputUGrid->BuildLinks();
  

  // Read spmesh element file and set connecting lines
  std::fstream spmeshFile_elem(spmeshFileName_elem.c_str());
  line_count = 0;

  std::vector<std::vector<vtkIdType> > neighbors;
  for(unsigned int n = 0; n < nUGridPts; n++) {
    std::vector<vtkIdType> neighborList;
    neighborList.push_back(n);
    neighbors.push_back(neighborList);
  }
  
  while(std::getline(spmeshFile_elem,line)){
    std::istringstream line_stream(line);
    if(line_count > 0) {
      unsigned int elm, n1, n2, n3, n4;
      line_stream >> elm >> n1 >> n2 >> n3 >> n4;

      arma::urowvec tetIds(4);
      tetIds(0) = n1 - 1;
      tetIds(1) = n2 - 1;
      tetIds(2) = n3 - 1;
      tetIds(3) = n4 - 1;

      // Make tets
      arma::rowvec tetIdCorrespondances(4);
      for(unsigned int n = 0; n < tetIds.n_elem; n++) {
	tetIdCorrespondances(n) = spmesh2ugrid_ptCorrespondance(tetIds(n));
      }
      if(!tetIdCorrespondances.has_nan()) {
	vtkSmartPointer<vtkIdList> tetIdsNew = vtkSmartPointer<vtkIdList>::New();
	tetIdsNew->SetNumberOfIds(4);
	for(unsigned int n = 0; n < tetIdCorrespondances.n_elem; n++) {
	  tetIdsNew->SetId(n,tetIdCorrespondances(n));
	}
	outputUGrid->InsertNextCell(VTK_TETRA,tetIdsNew);
      }
    }
    line_count++;
  }

  outputUGrid->BuildLinks();
  writeUnstructuredGrid(outputUGrid,outputUGridFileName);


  // Get Bounds
  double bounds[nDims*2];
  outputUGrid->GetPoints()->GetBounds(bounds);

  for(unsigned int d = 0; d < nDims*2; d+=2) {
    std::cout << floor(bounds[d]) << " " << ceil(bounds[d+1]) << " " << std::flush;
  }
  std::cout << std::endl;
  

  // Export voro file points
  std::ofstream voroFile(voroFileName.c_str());
  for(unsigned int p = 0; p < nUGridPts; p++) {
    double p0[nDims];
    outputUGrid->GetPoint(p,p0);
    voroFile << p << " " << p0[0] << " " << p0[1] << " " << p0[2] << "\n";
  }
  voroFile.close();

  
  return 0;
}
