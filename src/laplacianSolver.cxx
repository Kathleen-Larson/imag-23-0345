#include <iostream>
#include <string>
#include <stdlib.h>
#include <armadillo>
#include <limits>
#include <sstream>
#include <cmath>
#include <time.h>

#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkPointData.h>
#include <vtkPointLocator.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
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
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkImplicitBoolean.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkDelaunay3D.h>
#include <vtkGeometryFilter.h>
#include <vtkCellLocator.h>
#include <vtkExtractEdges.h>
#include <vtkImplicitPolyDataDistance.h>
#include <vtkDistancePolyDataFilter.h>



const int nDims = 3;



vtkSmartPointer<vtkPolyData> readPolyData(std::string filename) {

  vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(filename.c_str());
  reader->Update();

  return reader->GetOutput();
}



vtkSmartPointer<vtkUnstructuredGrid> readUnstructuredGrid(std::string filename) {

  vtkSmartPointer<vtkUnstructuredGridReader> uGridReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
  uGridReader->SetFileName(filename.c_str());
  uGridReader->Update();

  return uGridReader->GetOutput();
}



void writePolyData(vtkSmartPointer<vtkPolyData> polyData, std::string filename) {

  vtkSmartPointer<vtkXMLPolyDataWriter> writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetInputData(polyData);
  writer->SetFileName(filename.c_str());
  writer->Write();
}



void writeUnstructuredGrid(vtkSmartPointer<vtkUnstructuredGrid> uGrid, std::string filename) {

  vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
  writer->SetInputData(uGrid);
  writer->SetFileName(filename.c_str());
  writer->Write();
}



float polygonArea(arma::mat vertexCoords) {
  unsigned int nVerts = vertexCoords.n_rows;
  arma::rowvec area_vec(nDims);
  arma::rowvec v0 = vertexCoords.row(0);

  for(unsigned int k = 1; k < nVerts - 1; k++) {
    arma::rowvec v0vk = vertexCoords.row(k) - vertexCoords.row(0);
    arma::rowvec v0vk1 = vertexCoords.row(k+1) - vertexCoords.row(0);

    arma::rowvec crossProd(nDims);
    crossProd(0) = v0vk(1)*v0vk1(2) - v0vk(2)*v0vk1(1);
    crossProd(1) = v0vk(0)*v0vk1(2) - v0vk(2)*v0vk1(0);
    crossProd(2) = v0vk(0)*v0vk1(1) - v0vk(1)*v0vk1(0);

    area_vec += 0.5*crossProd;
  }

  arma::vec magn = arma::sqrt(arma::sum(arma::square(area_vec),1));
  return magn(0);
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



arma::mat findIntersections(vtkSmartPointer<vtkUnstructuredGrid> uGrid, float x0, float y0, float z0, float x1, float y1, float z1, float tol) {
  double p0[nDims] = {x0,y0,z0};
  double p1[nDims] = {x1,y1,z1};

  vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
  vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
  cellLocator->SetDataSet(uGrid);
  cellLocator->BuildLocator();
  cellLocator->FindCellsAlongLine(p0,p1,tol,cellIds);

  unsigned int nInts = cellIds->GetNumberOfIds();
  arma::mat intersections(nInts,nDims);
  for(unsigned int i = 0; i < nInts; i++) {
    vtkSmartPointer<vtkIdList> cellPtIds = vtkSmartPointer<vtkIdList>::New();
    uGrid->GetCellPoints(cellIds->GetId(i),cellPtIds);

    vtkSmartPointer<vtkPoints> polygonPoints = vtkSmartPointer<vtkPoints>::New();
    for(unsigned int j = 0; j < 3; j++) {
      polygonPoints->InsertNextPoint(uGrid->GetPoint(cellPtIds->GetId(j)));
    }

    vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New();
    polygon->GetPoints()->DeepCopy(polygonPoints);
    polygon->GetPointIds()->SetNumberOfIds(3);
    for(unsigned int j = 0; j < 3; j++) {
      polygon->GetPointIds()->SetId(j,j);
    }

    double pIntersection[nDims], pCoords[nDims];
    double t;
    int subId;
    vtkIdType id = polygon->IntersectWithLine(p0,p1,tol,t,pIntersection,pCoords,subId);
    for(unsigned int d = 0; d < nDims; d++) {
      intersections(i,d) = pIntersection[d];
    }
  }

  return intersections;
}



float manualInterpolation(vtkSmartPointer<vtkUnstructuredGrid> ugrid, vtkIdType p, unsigned int n, arma::vec interpolationData) {

  vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
  pointLocator->SetDataSet(ugrid);
  pointLocator->BuildLocator();
  
  double p0[nDims];
  ugrid->GetPoint(p,p0);

  vtkSmartPointer<vtkIdList> nClosestPtIds = vtkSmartPointer<vtkIdList>::New();
  pointLocator->FindClosestNPoints(n,p0,nClosestPtIds);


  // Interpolate selected points
  arma::vec ptVals(n,arma::fill::zeros), ptDistances(n,arma::fill::zeros);
  std::vector<int> nanValues;

  for(unsigned int i = 1; i < n; i++) {
    unsigned int ptId = nClosestPtIds->GetId(i);

    arma::vec interpolationValue(1);
    interpolationValue(0) = interpolationData(nClosestPtIds->GetId(i));

    if(interpolationValue.has_nan() || ptId==p) {
      ptVals(i) = arma::datum::nan;
      ptDistances(i) = arma::datum::nan;
      nanValues.push_back(i);
    }
    else {
      ptVals(i) = interpolationData(nClosestPtIds->GetId(i));
      double p1[nDims];
      ugrid->GetPoints()->GetPoint(nClosestPtIds->GetId(i),p1);
      ptDistances(i) = std::sqrt(vtkMath::Distance2BetweenPoints(p0,p1));
    }
  }

  int n_nanValues = static_cast<int>(nanValues.size());
  arma::uvec nanRows(n_nanValues);
  for(unsigned int j = 0; j < n_nanValues; j++) {
    nanRows(j) = nanValues.at(j);
  }
  ptDistances.shed_rows(nanRows);
  ptVals.shed_rows(nanRows);

  arma::vec ptDistancesNorm = ptDistances / arma::sum(ptDistances);
  float interpVal = arma::sum(ptDistancesNorm % ptVals);
  return interpVal;
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

  // Parse inputs
  if(argc < 7) {
    std::cout << "Usage: inputFDMuGridFileName annulusFileName outerSurfaceFileName innerSurfaceFileName voronoiFileName outputFDMuGridFileName" << std::endl;
    return -1;
  }

  std::string inputUGridFileName = argv[1];
  std::string annulusFileName = argv[2];
  std::string outerSurfaceFileName = argv[3];
  std::string innerSurfaceFileName = argv[4];
  std::string voronoiFileName = argv[5];
  std::string outputFDMuGridFileName = argv[6];


  // Read input data
  vtkSmartPointer<vtkUnstructuredGrid> inputUGrid = readUnstructuredGrid(inputUGridFileName);
  unsigned int nUGridPts = inputUGrid->GetPoints()->GetNumberOfPoints();

  vtkSmartPointer<vtkPolyData> annulus = readPolyData(annulusFileName);
  vtkSmartPointer<vtkDataArray> annulusPtNormals = annulus->GetPointData()->GetNormals();
  vtkSmartPointer<vtkFloatArray> inputBoundaryLabels = vtkFloatArray::SafeDownCast(inputUGrid->GetPointData()->GetArray("Boundary Ids"));

  vtkSmartPointer<vtkPolyData> outerSurface = readPolyData(outerSurfaceFileName);
  unsigned int outerSurface_nPts = outerSurface->GetPoints()->GetNumberOfPoints();

  vtkSmartPointer<vtkPolyData> innerSurface = readPolyData(innerSurfaceFileName);
  unsigned int innerSurface_nPts = innerSurface->GetPoints()->GetNumberOfPoints();


  // Read Voronoi diagram file
  std::fstream voronoiFile(voronoiFileName.c_str());
  std::string line;
  arma::field<arma::ivec> adjacentPts_voro(nUGridPts);
  arma::field<arma::vec> distances_voro(nUGridPts), vFaceAreas_voro(nUGridPts);

  while(std::getline(voronoiFile,line)) {
    std::istringstream line_stream(line);
    unsigned int p, nFaces;
    line_stream >> p >> nFaces;
    
    double p0[nDims];
    inputUGrid->GetPoint(p,p0);

    arma::ivec adjacentPts_p(nFaces);
    arma::vec distances_p(nFaces);
    for(unsigned int n = 0; n < nFaces; n++) {
      // Get neighbor node IDs
      int q;
      line_stream >> q;
      adjacentPts_p(n) = q;

      // Get distance between points
      if(q >= 0 ) {
	double q0[nDims];
	inputUGrid->GetPoint(q,q0);
	distances_p(n) = std::sqrt(vtkMath::Distance2BetweenPoints(p0,q0));
      } else {
	distances_p(n) = arma::datum::inf;
      }
    }  

    // Get Voronoi face areas
    arma::vec vFaceAreas_p(nFaces);
    for(unsigned int n = 0; n < nFaces; n++) {
      float vFaceArea;
      line_stream >> vFaceArea;
      vFaceAreas_p(n) = vFaceArea;
    }

    // Remove invalid neighbor IDs
    arma::uvec invalidIds = arma::find(adjacentPts_p < 0);
    adjacentPts_p.shed_rows(invalidIds);
    distances_p.shed_rows(invalidIds);
    vFaceAreas_p.shed_rows(invalidIds);
    
    // Build adjacent point field data
    if(p < nUGridPts) {
      adjacentPts_voro(p) = adjacentPts_p;
      distances_voro(p) = distances_p;
      vFaceAreas_voro(p) = vFaceAreas_p;
    } else {
      std::cout << "hmmm" << std::endl;
    }
  }



  /* Build FDM polydata
     - Get neighbors from Voro and compare to SPMESH
  */

  // Set points based on SPMESH
  vtkSmartPointer<vtkUnstructuredGrid> FDMuGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
  vtkSmartPointer<vtkPoints> FDMuGridPoints = vtkSmartPointer<vtkPoints>::New();
  vtkSmartPointer<vtkFloatArray> boundaryLabels = vtkSmartPointer<vtkFloatArray>::New();
  boundaryLabels->SetNumberOfComponents(1);
  boundaryLabels->SetNumberOfTuples(nUGridPts);
  boundaryLabels->SetName("Boundary Ids");

  for(unsigned int p = 0; p < nUGridPts; p++) {
    double p0[nDims];
    inputUGrid->GetPoints()->GetPoint(p,p0);
    FDMuGridPoints->InsertNextPoint(p0);
    boundaryLabels->SetValue(p,inputBoundaryLabels->GetValue(p));
  }
  FDMuGrid->SetPoints(FDMuGridPoints);
  FDMuGrid->GetPointData()->AddArray(boundaryLabels);


  // Identify neighbors
  arma::field<arma::ivec> adjacentPts(nUGridPts);
  arma::field<arma::vec> distances(nUGridPts), vFaceAreas(nUGridPts);
  for(unsigned int p = 0; p < nUGridPts; p++) {

    // Get all point neighbors from SPMESH ugrid
    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    inputUGrid->GetPointCells(p,cellIds);
    std::vector<vtkIdType> neighbors_SPMESH;

    for(unsigned int i = 0; i < cellIds->GetNumberOfIds(); i++) {
      unsigned int cell = cellIds->GetId(i);
      vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
      inputUGrid->GetCellPoints(cell,ptIds);    
      
      for(unsigned int j = 0; j < ptIds->GetNumberOfIds(); j++) {
	int q = ptIds->GetId(j);
	if(q!=p && isUnique(q,neighbors_SPMESH)) {
	  neighbors_SPMESH.push_back(q);
	}
      }
    }

    // Compare voro output to SPMESH neighbors
    arma::ivec adjacentPts_spmesh_p(0,arma::fill::none);
    arma::vec distances_spmesh_p(0,arma::fill::none), vFaceAreas_spmesh_p(0,arma::fill::none);

    arma::ivec adjacentPts_voro_p = adjacentPts_voro(p);
    arma::vec distances_voro_p = distances_voro(p);
    arma::vec vFaceAreas_voro_p = vFaceAreas_voro(p);

    for(unsigned int n = 0; n < adjacentPts_voro_p.n_elem; n++) {
      int q = adjacentPts_voro_p(n);
      if(!isUnique(q,neighbors_SPMESH)) {
	arma::irowvec q_rowvec(1);
	arma::vec distances_q_rowvec(1), vFaceAreas_q_rowvec(1);
	q_rowvec(0) = q;
	distances_q_rowvec(0) = distances_voro_p(n);
	vFaceAreas_q_rowvec(0) = vFaceAreas_voro_p(n);
	
	adjacentPts_spmesh_p.insert_rows(adjacentPts_spmesh_p.n_elem,q_rowvec);
	distances_spmesh_p.insert_rows(distances_spmesh_p.n_elem,distances_q_rowvec);
	vFaceAreas_spmesh_p.insert_rows(vFaceAreas_spmesh_p.n_elem,vFaceAreas_q_rowvec);

	if(p < q) {
          vtkSmartPointer<vtkIdList> lineIds = vtkSmartPointer<vtkIdList>::New();
          lineIds->SetNumberOfIds(2);
          lineIds->SetId(0,p);
          lineIds->SetId(1,q);
          FDMuGrid->InsertNextCell(VTK_LINE,lineIds);
        }
      }
    }
    adjacentPts(p) = adjacentPts_spmesh_p;
    distances(p) = distances_spmesh_p;
    vFaceAreas(p) = vFaceAreas_spmesh_p;
  }
  FDMuGrid->BuildLinks();
 

  // Double check adjacency lists
  for(unsigned int p = 0; p < nUGridPts; p++) {
    arma::ivec adjacentPts_p = adjacentPts(p);
    arma::vec distances_p = distances(p);
    arma::vec vFaceAreas_p = vFaceAreas(p);

    for(unsigned int n = 0; n < adjacentPts_p.n_elem; n++) {
      int q = adjacentPts_p(n);
      arma::ivec adjacentPts_q = adjacentPts(q);
      arma::uvec inds = find(adjacentPts_q==p);
      
      if(inds.n_elem==0) {
	adjacentPts_q.insert_rows(0,1);
	adjacentPts_q(0) = p;

	arma::vec distances_q = distances(q);
	distances_q.insert_rows(0,1);
	distances_q(0) = distances_p(n);

	arma::vec vFaceAreas_q = vFaceAreas(q);
	vFaceAreas_q.insert_rows(0,1);
	vFaceAreas_q(0) = vFaceAreas_p(n);

      } else if(inds.n_elem > 1) {
	std::cout << "huh??" << std::endl;
      }
    }
  }

  

  /* Solve for gradient of Laplacian equation over the FDM mesh
     - Using methods from (Sukumar and Bolander, 2003)
     - Calculate derivative at each point by estimated a local uniform grid
  */

  // Initialize gradient
  arma::mat UGradient0(nUGridPts,nDims,arma::fill::zeros), UGradient1(nUGridPts,nDims,arma::fill::zeros);

  for(unsigned int p = 0; p < nUGridPts; p++) {
    double p0[nDims];
    FDMuGridPoints->GetPoint(p,p0);    
    
    if(boundaryLabels->GetValue(p) == 0) { // Interior node --> Voronoi cell method
      UGradient0(p,0) = 0;
      UGradient0(p,1) = 0;
      UGradient0(p,2) = 0;
      
    } else { // Boundary node
      // Get normal from annulus
      vtkSmartPointer<vtkPointLocator> pointLocator = vtkSmartPointer<vtkPointLocator>::New();
      pointLocator->SetDataSet(annulus);
      pointLocator->BuildLocator();
      double closestPoint[nDims], pMesh[nDims];
      vtkIdType pMeshId = pointLocator->FindClosestPoint(p0);
      annulus->GetPoint(pMeshId,pMesh);

      float dist2vertex = std::sqrt(vtkMath::Distance2BetweenPoints(p0,pMesh));
      double normal[nDims];
      
      if(dist2vertex < 0.01) { 
	// If shared vertex with annulus
	annulusPtNormals->GetTuple(pMeshId,normal);
	UGradient0(p,0) = normal[0];
	UGradient0(p,1) = normal[1];
	UGradient0(p,2) = normal[2];
	UGradient0.row(p) = UGradient0.row(p);

      } else {
	// If vertex lies in the middle of an annulus cell
	vtkSmartPointer<vtkCellLocator> cellLocator = vtkSmartPointer<vtkCellLocator>::New();
	cellLocator->SetDataSet(annulus);
	cellLocator->BuildLocator();
	double closestPoint[nDims], closestPointDist2;
	vtkIdType cellId;
	int subId;
	cellLocator->FindClosestPoint(p0,closestPoint,cellId,subId,closestPointDist2);

	// Interpolate normal based on closest annulus vertices
	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
	annulus->GetCellPoints(cellId,ptIds);
	unsigned int nIds = ptIds->GetNumberOfIds();
	arma::vec dists2p(nIds,arma::fill::zeros);
	arma::uvec has_nan_p(nIds,arma::fill::zeros);
	for(unsigned int i = 0; i < nIds; i++) {
	  double q0[nDims];
	  FDMuGridPoints->GetPoint(ptIds->GetId(i),q0);
	  dists2p(i) = std::sqrt(vtkMath::Distance2BetweenPoints(p0,q0));
	  if(dists2p.has_nan()) {
	    has_nan_p(i) = 1;
	  }
	}
	
	arma::uvec has_nan_p_i = find(has_nan_p == 1);
	dists2p.shed_rows(has_nan_p_i);
	for(unsigned int i = 0; i < dists2p.n_elem; i++) {
	  annulusPtNormals->GetTuple(ptIds->GetId(i),normal);	  
	  UGradient0(p,0) += normal[0] * dists2p(i) / arma::sum(dists2p);
	  UGradient0(p,1) += normal[1] * dists2p(i) / arma::sum(dists2p);
	  UGradient0(p,2) += normal[2] * dists2p(i) / arma::sum(dists2p);
	}
      }
    }
  }
  UGradient1 = UGradient0;


  // Iteratively solve for the gradient
  unsigned int it = 0;
  float L_inf = 1, tol = 0.001;
  
  while(L_inf > tol && it < 100) {
    for(unsigned int p = 0; p < nUGridPts; p++) {
      if(boundaryLabels->GetValue(p) == 0) {
	arma::ivec adjacentPts_p = adjacentPts(p);
	arma::vec distances_p = distances(p);
	arma::vec vFaceAreas_p = vFaceAreas(p);
	
	float alphaI = 0;
        for(unsigned int n = 0; n < adjacentPts_p.n_elem; n++) {
          alphaI += vFaceAreas_p(n) / distances_p(n);
        }


	for(unsigned int d = 0; d < nDims; d++) {
	  float sum_J = 0;
	  for(unsigned int j = 0; j < adjacentPts_p.n_elem; j++) {
	    unsigned int q = adjacentPts_p(j);

	    float vFaceAreas_pj = vFaceAreas_p(j);
	    float distances_pj = distances_p(j);
	    float dU_j = UGradient0(q,d);

	    sum_J += (vFaceAreas_pj * dU_j / distances_pj) / alphaI;
	  }
	  UGradient1(p,d) = sum_J;
	}
      }
    }
    arma::mat diff = arma::abs(UGradient0 - UGradient1);
    L_inf = diff.max();
    UGradient0 = UGradient1;
    it++;
  }
  

  // Normalize
  arma::mat UGradientNorm(nUGridPts,nDims,arma::fill::zeros);
  for(unsigned int p = 0; p < nUGridPts; p++) {
    UGradientNorm.row(p) = UGradient1.row(p) / sqrt(arma::sum(arma::square(UGradient1.row(p))));
  }

  
  // Add gradient components as scalar arrays
  // Gradient_x
  std::string xGradientArrayName = "Solution gradientX";
  vtkSmartPointer<vtkFloatArray> xGradient = vtkSmartPointer<vtkFloatArray>::New();
  xGradient->SetNumberOfComponents(1);
  xGradient->SetNumberOfTuples(nUGridPts);
  xGradient->SetName(xGradientArrayName.c_str());

  // Gradient_y
  std::string yGradientArrayName = "Solution gradientY";
  vtkSmartPointer<vtkFloatArray> yGradient = vtkSmartPointer<vtkFloatArray>::New();
  yGradient->SetNumberOfComponents(1);
  yGradient->SetNumberOfTuples(nUGridPts);
  yGradient->SetName(yGradientArrayName.c_str());

  // Gradient_z
  std::string zGradientArrayName = "Solution gradientZ";
  vtkSmartPointer<vtkFloatArray> zGradient = vtkSmartPointer<vtkFloatArray>::New();
  zGradient->SetNumberOfComponents(1);
  zGradient->SetNumberOfTuples(nUGridPts);
  zGradient->SetName(zGradientArrayName.c_str());


  // Assign
  for(unsigned int p = 0; p < nUGridPts; p++) {
    // Check if nan
    if(boundaryLabels->GetValue(p)==0 && UGradientNorm.row(p).has_nan()) {
      for(unsigned int d = 0; d < nDims; d++) {
	if(UGradientNorm.row(p).col(d).has_nan()) {
	  UGradientNorm(p,d) = manualInterpolation(FDMuGrid,p,4,UGradientNorm.col(d));
	}
      }
    }
    xGradient->SetValue(p,UGradientNorm(p,0));
    yGradient->SetValue(p,UGradientNorm(p,1));
    zGradient->SetValue(p,UGradientNorm(p,2));
  }
  
  FDMuGrid->GetPointData()->AddArray(xGradient);
  FDMuGrid->GetPointData()->AddArray(yGradient);
  FDMuGrid->GetPointData()->AddArray(zGradient);
  FDMuGrid->BuildLinks();
  writeUnstructuredGrid(FDMuGrid,outputFDMuGridFileName);

  return 0;
}
