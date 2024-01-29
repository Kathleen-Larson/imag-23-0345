#include <iostream>
#include <string>
#include <stdlib.h>
#include <armadillo>
#include <limits>
#include <sstream>
#include <cmath>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkLine.h>
#include <vtkTriangle.h>
#include <vtkPolygon.h>
#include <vtkMath.h>

#include <vtkPolyData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataNormals.h>
#include <vtkAppendPolyData.h>
#include <vtkCleanPolyData.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGridWriter.h>

#include <vtkSphereSource.h>
#include <vtkSphere.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkCellLocator.h>
#include <vtkPointLocator.h>

#include <vtkOBBTree.h>
#include <vtkKdTreePointLocator.h>



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



vtkSmartPointer<vtkUnstructuredGrid> readUnstructuredGrid(std::string filename) {

  vtkSmartPointer<vtkUnstructuredGridReader> uGridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
  uGridReader->SetFileName(filename.c_str());
  uGridReader->Update();

  return uGridReader->GetOutput();
}



void writeUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> uGrid, std::string filename) {

  vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  writer->SetInputData(uGrid);
  writer->SetFileName(filename.c_str());
  writer->Write();
}



bool isValid(float label, std::vector<int> bounds) {
  unsigned int nBounds = static_cast<int>(bounds.size());
  bool valid = false;
  for(unsigned int b = 0; b < nBounds; b+=2) {
    if(label >= bounds.at(b) && label <= bounds.at(b+1)) {
      valid = true;
    }
  }

  if(valid) {
    return true;
  } else {
    return false;
  }
}



vtkSmartPointer<vtkDataArray> GetNormals(vtkSmartPointer<vtkPolyData> polydata) {
  vtkSmartPointer<vtkPolyDataNormals> normalsGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
  normalsGenerator->SetInputData(polydata);
  normalsGenerator->ComputePointNormalsOn();
  normalsGenerator->ComputeCellNormalsOff();
  normalsGenerator->Update();

  vtkSmartPointer<vtkDataArray> normals = normalsGenerator->GetOutput()->GetPointData()->GetNormals();
  return normals;
}




void visualizePoint(float x, float y, float z, std::string filename) {

  vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
  sphereSource->SetCenter(x,y,z);
  sphereSource->SetRadius(0.25);
  sphereSource->SetPhiResolution(5);
  sphereSource->SetThetaResolution(5);
  sphereSource->Update();

  writePolyData(sphereSource->GetOutput(),filename);
}



float distance_pt2polydata(vtkSmartPointer<vtkPolyData> surface, float x, float y, float z) {
  
  vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(surface);
  cellLocator->BuildLocator();

  // Find cell closest to point
  double p0[nDims] = {x, y, z};
  double closestPoint[nDims], closestPointDist2;
  vtkIdType cellId;
  int subId;
  cellLocator->FindClosestPoint(p0,closestPoint,cellId,subId,closestPointDist2);

  return sqrt(closestPointDist2);
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

  std::cout << "Is inside?  " << selectEnclosedPointsFilter->IsInside(0) << std::endl;

  if(selectEnclosedPointsFilter->IsInside(0))
    return true;
  return false;
  
}



arma::mat findIntersections(vtkSmartPointer<vtkPolyData> surface, float x0, float y0, float z0, float x1, float y1, float z1, float tol) {
  double p0[nDims] = {x0,y0,z0};
  double p1[nDims] = {x1,y1,z1};

  vtkSmartPointer<vtkOBBTree> OBBTree = vtkSmartPointer<vtkOBBTree>::New();
  OBBTree->SetDataSet(surface);
  OBBTree->BuildLocator();
  
  vtkSmartPointer<vtkPoints> intersectionPoints = vtkSmartPointer<vtkPoints>::New();
  OBBTree->IntersectWithLine(p0,p1,intersectionPoints,NULL);

  unsigned int nInts = intersectionPoints->GetNumberOfPoints();
  arma::mat intersections(nInts,nDims);
  
  if(nInts > 0) {
    for(unsigned int p = 0; p < nInts; p++) {
      double pIntersection[nDims];
      intersectionPoints->GetPoint(p,pIntersection);
      for(unsigned int d = 0; d < nDims; d++) {
	intersections(p,d) = pIntersection[d];
      }
    }
  } else {
    intersections.clear();
  }

  return intersections;
}


float manualInterpolation(vtkSmartPointer<vtkUnstructuredGrid> uGrid, float x, float y, float z, unsigned int n, vtkSmartPointer<vtkFloatArray> labelArray) {
  vtkSmartPointer<vtkKdTreePointLocator> pointLocator = vtkSmartPointer<vtkKdTreePointLocator>::New();
  pointLocator->SetDataSet(uGrid);
  pointLocator->BuildLocator();

  double p0[nDims] = {x, y, z};
  vtkSmartPointer<vtkIdList> nClosestPtIds = vtkSmartPointer<vtkIdList>::New();
  pointLocator->FindClosestNPoints(n,p0,nClosestPtIds);


  // Interpolate selected points
  arma::vec ptVals(n,arma::fill::zeros), ptDistances(n,arma::fill::zeros);
  for(unsigned int i = 0; i < n; i++) {
    ptVals(i) = labelArray->GetValue(nClosestPtIds->GetId(i));
    double p1[nDims];
    uGrid->GetPoints()->GetPoint(nClosestPtIds->GetId(i),p1);
    ptDistances(i) = std::sqrt(vtkMath::Distance2BetweenPoints(p0,p1));
  }
  arma::vec ptDistancesNorm = ptDistances / arma::sum(ptDistances);
  float interpVal = arma::sum(ptDistancesNorm % ptVals);
  
  return interpVal;
}



arma::mat manualInterpolation2(vtkSmartPointer<vtkUnstructuredGrid> uGrid, float x, float y, float z, unsigned int n) {

  // Get point IDs
  vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet(uGrid);
  pointLocator->BuildLocator();

  double p0[nDims] = {x, y, z};
  vtkSmartPointer<vtkIdList> nClosestPtIds = vtkSmartPointer<vtkIdList>::New();
  pointLocator->FindClosestNPoints(n,p0,nClosestPtIds);


  // Get interpolation weighting
  arma::mat interpData(n,2,arma::fill::zeros);
  for(unsigned int i = 0; i < n; i++) {
    unsigned int id = nClosestPtIds->GetId(i);
    interpData(i,0) = id;
    
    double p1[nDims];
    uGrid->GetPoints()->GetPoint(id,p1);
    interpData(i,1) = std::sqrt(vtkMath::Distance2BetweenPoints(p0,p1));  
  }

  interpData.col(1) = interpData.col(1) / arma::sum(interpData.col(1));

  return interpData;
}


bool isFloatArray(vtkSmartPointer<vtkPolyData> polydata, std::string arrayName) {
  unsigned int nPts = polydata->GetPoints()->GetNumberOfPoints();

  vtkSmartPointer<vtkFloatArray> testArray = vtkFloatArray::SafeDownCast(polydata->GetPointData()->GetArray(arrayName.c_str()));
  
  if(testArray) {
    std::cout << "Array found" << std::endl;
    return true;
  }
  return false;
}


int main(int argc, char * argv[]) {

  // Parse inputs
  if(argc < 7) {
    std::cout << "Usage: FDMmeshFileName inputSurfaceOuterFileName inputSurfaceInnerFileName outputSurfaceFileName stepSize integrationDirection(1: out->in, -1:in->out) boundsArrayName bounds(must be in pairs of lower/upper bounds)" << std::endl;
    return -1;
  }
  
  std::string FDMmeshFileName = argv[1];
  std::string inputSurfaceOuterFileName = argv[2];
  std::string inputSurfaceInnerFileName = argv[3];
  std::string outputSurfaceFileName = argv[4];
  float stepSize = -1*atof(argv[5]);
  int integrationDirection = atoi(argv[6]);
  int argc_boundsInfoStart = 7;

  

  // Read data
  vtkSmartPointer<vtkUnstructuredGrid> inputUGrid = readUnstructuredGrid(FDMmeshFileName);
  vtkSmartPointer<vtkFloatArray> xGradient = vtkFloatArray::SafeDownCast(inputUGrid->GetPointData()->GetArray("Solution gradientX"));
  vtkSmartPointer<vtkFloatArray> yGradient = vtkFloatArray::SafeDownCast(inputUGrid->GetPointData()->GetArray("Solution gradientY"));
  vtkSmartPointer<vtkFloatArray> zGradient = vtkFloatArray::SafeDownCast(inputUGrid->GetPointData()->GetArray("Solution gradientZ"));
  unsigned int nUGridPts = inputUGrid->GetPoints()->GetNumberOfPoints();

  vtkSmartPointer<vtkPolyData> surfaceStart = vtkSmartPointer<vtkPolyData>::New();
  vtkSmartPointer<vtkPolyData> surfaceEnd = vtkSmartPointer<vtkPolyData>::New();

  if(integrationDirection == 1) { // outer-->inner
    surfaceStart = readPolyData(inputSurfaceOuterFileName);
    surfaceEnd = readPolyData(inputSurfaceInnerFileName);
  }
  else if(integrationDirection == -1) { // inner-->outer
    surfaceStart = readPolyData(inputSurfaceInnerFileName);
    surfaceEnd = readPolyData(inputSurfaceOuterFileName);
  }
  else {
    std::cout << "Invalid integration direction (1: out->in, -1:in->out)" << std::endl;
    return -1;
  }

  vtkSmartPointer<vtkPolyData> outputSurface = vtkSmartPointer<vtkPolyData>::New();
  outputSurface->DeepCopy(surfaceStart);
  vtkSmartPointer<vtkDataArray> outputSurface_normals = GetNormals(outputSurface);
  unsigned int outputSurface_nPts = outputSurface->GetPoints()->GetNumberOfPoints();

  vtkSmartPointer<vtkFloatArray> boundsArray = vtkSmartPointer<vtkFloatArray>::New();
  unsigned int nBounds = argc - (argc_boundsInfoStart);

  bool isBounds;
  std::vector<int> bounds;
  if(nBounds == 0) {
    isBounds = false;
    bounds.clear();
  }
  else {   
    if(nBounds % 2 != 0) {
      std::cout << "Must provide even number of bounds for parcellation" << std::endl;
      return -1;
    }
    else {
      isBounds = true;
      std::string boundsArrayName = argv[argc_boundsInfoStart];
      boundsArray = vtkFloatArray::SafeDownCast(surfaceStart->GetPointData()->GetArray(boundsArrayName.c_str()));
      for(unsigned int n = 0; n < nBounds; n++) {
	bounds.push_back(atoi(argv[argc - nBounds + n]));
      }
    }
  }
  

  // Convert polydata (bounding surfaces) to ugrids
  unsigned int surfaceStart_nPts = surfaceStart->GetPoints()->GetNumberOfPoints();
  vtkSmartPointer<vtkPoints> surfaceStart_UGridPts = vtkSmartPointer<vtkPoints>::New();
  for(unsigned int p = 0; p < surfaceStart_nPts; p++) {
    double p0[nDims];
    surfaceStart->GetPoint(p,p0);
    surfaceStart_UGridPts->InsertNextPoint(p0);
  }

  vtkSmartPointer<vtkUnstructuredGrid> surfaceStart_UGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  surfaceStart_UGrid->SetPoints(surfaceStart_UGridPts);
  for(unsigned int n = 0; n < surfaceStart->GetNumberOfCells(); n++) {
    surfaceStart->GetCell(n);
    vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
    surfaceStart->GetCellPoints(n,cellPtIds);

    if(cellPtIds->GetNumberOfIds() == 3) {
      surfaceStart_UGrid->InsertNextCell(VTK_TRIANGLE,cellPtIds);
    }
  }
  
  unsigned int surfaceEnd_nPts = surfaceEnd->GetPoints()->GetNumberOfPoints();
  vtkSmartPointer<vtkPoints> surfaceEnd_UGridPts = vtkSmartPointer<vtkPoints>::New();
  for(unsigned int p = 0; p < surfaceEnd_nPts; p++) {
    double p0[nDims];
    surfaceEnd->GetPoint(p,p0);
    surfaceEnd_UGridPts->InsertNextPoint(p0);
  }
  
  vtkSmartPointer<vtkUnstructuredGrid> surfaceEnd_UGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  surfaceEnd_UGrid->SetPoints(surfaceEnd_UGridPts);
  for(unsigned int n = 0; n < surfaceEnd->GetNumberOfCells(); n++) {
    surfaceEnd->GetCell(n);
    vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
    surfaceEnd->GetCellPoints(n,cellPtIds);

    if(cellPtIds->GetNumberOfIds() == 3) {
      surfaceEnd_UGrid->InsertNextCell(VTK_TRIANGLE,cellPtIds);
    }
  }
  

  // Initialize label arrays
  std::string thickness_array_name = "SBL-Thickness";
  std::string completed_array_name = "Completed";
  vtkSmartPointer<vtkFloatArray> laplacianThickness_array = vtkFloatArray::SafeDownCast(surfaceStart->GetPointData()->GetArray(thickness_array_name.c_str()));
  vtkSmartPointer<vtkFloatArray> completed_array = vtkFloatArray::SafeDownCast(surfaceStart->GetPointData()->GetArray(completed_array_name.c_str()));
    
  if(!laplacianThickness_array || !completed_array) {
    laplacianThickness_array = vtkSmartPointer<vtkFloatArray>::New();
    laplacianThickness_array->SetNumberOfComponents(1);
    laplacianThickness_array->SetNumberOfValues(surfaceStart_nPts);
    laplacianThickness_array->SetName(thickness_array_name.c_str());

    completed_array = vtkSmartPointer<vtkFloatArray>::New();
    completed_array->SetNumberOfComponents(1);
    completed_array->SetNumberOfValues(surfaceStart_nPts);
    completed_array->SetName(completed_array_name.c_str());
    
    for(unsigned int p = 0; p < outputSurface_nPts; p++) {
      laplacianThickness_array->SetValue(p,0);
      completed_array->SetValue(p,0);
    }
    
    outputSurface->GetPointData()->AddArray(laplacianThickness_array);
    outputSurface->GetPointData()->AddArray(completed_array);
    outputSurface->BuildLinks();
  }



  // Integration
  float print_update_thresh = 0.05;
  float print_update_thresh_init = print_update_thresh;

  float maxLen = 10.0;

  for(unsigned int p = 0; p < outputSurface_nPts; p++) {
    
    // Do we proceed?
    bool proceed = true;
    int completed = static_cast<int>(completed_array->GetValue(p));

    if(isBounds) {
      float label = boundsArray->GetValue(p);
      if(!isValid(label,bounds) || completed == 1) {
	proceed = false;
      }
    }

    // Ok we do
    if(proceed) {
      double p0[nDims], p1[nDims];
      surfaceStart->GetPoint(p,p0);
      
      // Initialize the line integration
      float len = 0;
      unsigned int it = 0;
      bool isIntersection = false;
      std::vector<int> integrationPtIds;

      double normal[nDims];
      arma::vec direction(nDims,arma::fill::zeros);
      outputSurface_normals->GetTuple(p,normal);
      direction(0) = normal[0];
      direction(1) = normal[1];
      direction(2) = normal[2];
      direction = integrationDirection * stepSize * (direction / sqrt(arma::sum(arma::square(direction))));
      
      for(unsigned int d = 0; d < nDims; d++) {
	p1[d] = p0[d] + direction(d);
      }
      len += sqrt(arma::sum(arma::square(direction)));
      
      while(!isIntersection && len < maxLen) {
	// Update
	for(unsigned int d = 0; d < nDims; d++) {
	  p0[d] = p1[d];
	}
	direction.zeros();

	arma::mat ids_weighting = manualInterpolation2(inputUGrid,p0[0],p0[1],p0[2],6);
	for(unsigned int n = 0; n < ids_weighting.n_rows; n++) {
          vtkIdType id = ids_weighting(n,0);
          direction(0) += xGradient->GetValue(id) / ids_weighting(n,1);
          direction(1) += yGradient->GetValue(id) / ids_weighting(n,1);
          direction(2) += zGradient->GetValue(id) / ids_weighting(n,1);
        }

	direction = integrationDirection * stepSize * (direction / sqrt(arma::sum(arma::square(direction))));
	
	// Extend along search direction
	for(unsigned int d = 0; d < nDims; d++) {
	  p1[d] = p0[d] + direction(d);
	}
      
	// Check for intersection
	arma::mat pIntersections = findIntersections(surfaceEnd,p0[0],p0[1],p0[2],p1[0],p1[1],p1[2],0.001);
	
	if(pIntersections.n_rows == 0) {
	  len += sqrt(arma::sum(arma::square(direction)));
	}
	else if(pIntersections.n_rows == 1) {
	  double pIntersection[nDims] = {pIntersections(0),pIntersections(1),pIntersections(2)};
	  len += sqrt(vtkMath::Distance2BetweenPoints(p0,pIntersection));	  
	  isIntersection = true;
	}
	else if(pIntersections.n_rows > 1) {
	  arma::vec dIntersections(pIntersections.n_rows);
	  for(unsigned int q = 0; q < dIntersections.n_elem; q++) {
	    double pIntersection[nDims] = {pIntersections.row(q)(0),pIntersections.row(q)(1),pIntersections.row(q)(2)};
	    dIntersections(q) = sqrt(vtkMath::Distance2BetweenPoints(p0,pIntersection));
	  }
	  arma::uword qmin = index_min(dIntersections);
	  len += dIntersections(qmin);
	  	
	  double pEnd[nDims] = {pIntersections.row(qmin)(0),pIntersections.row(qmin)(1),pIntersections.row(qmin)(2)};
	  isIntersection = true;
	}
	it++;
      }
      laplacianThickness_array->SetValue(p,len);
      completed_array->SetValue(p,1);

      // Show progress and output current iteration 
      float percent_done = static_cast<float>(p)/static_cast<float>(outputSurface_nPts);
      if(percent_done > print_update_thresh) {
	int int_percent_done = static_cast<int>(100*percent_done);
	std::cout << int_percent_done << "%....... " << std::endl;
	print_update_thresh += print_update_thresh_init;

	outputSurface->GetPointData()->AddArray(laplacianThickness_array);
	outputSurface->GetPointData()->AddArray(completed_array);
	outputSurface->BuildLinks();
	writePolyData(outputSurface,outputSurfaceFileName);
      }
    }
  }
  std::cout << "Done!" << std::endl;


  // Output data
  outputSurface->GetPointData()->AddArray(laplacianThickness_array);
  outputSurface->GetPointData()->AddArray(completed_array);
  outputSurface->BuildLinks();
  writePolyData(outputSurface,outputSurfaceFileName);
  
  
  return 0;
}
