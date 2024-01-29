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



bool isUnique(vtkIdType p, std::vector<vtkIdType> vec) {
  for(unsigned int i = 0; i < static_cast<int>(vec.size()); i++) {
    if(p == vec.at(i)) {
      return false;
    }
  }
  return true;
}



int main(int argc, char * argv[]) {
  
  if(argc < 1) {
    std::cout << "Usage: meshFileName stat(Min, Max, Mean, Median)" << std::endl;
    return -1;
  }


  // Parse args
  std::string meshFileName = argv[1];
  std::string stat = argv[2];

  vtkSmartPointer<vtkPolyData> mesh = readPolyData(meshFileName);
  unsigned int nPts = mesh->GetPoints()->GetNumberOfPoints();

  // Get average size edge
  std::vector<std::vector<vtkIdType> > neighborPts;
  for(unsigned int p = 0; p < nPts; p++) {
    std::vector<vtkIdType> neighborList;
    neighborList.push_back(p);
    neighborPts.push_back(neighborList);
  }

  std::vector<float> distances_vec;
  for(unsigned int p = 0; p < nPts; p++) {
    double p0[nDims];
    mesh->GetPoint(p,p0);

    vtkSmartPointer<vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();
    mesh->GetPointCells(p,cellIds);
    
    for(unsigned int i = 0; i < cellIds->GetNumberOfIds(); i++) {
      unsigned int cellId = cellIds->GetId(i);

      int cellType = mesh->GetCellType(cellId);
      if(cellType == 3 || cellType == 5 || cellType==10) { // this only works with lines, triangles, or tetrahedra
	vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
	mesh->GetCellPoints(cellId,ptIds);
	
	for(unsigned int j = 0; j < ptIds->GetNumberOfIds(); j++) {
	  unsigned int q = ptIds->GetId(j);
	  
	  if(isUnique(q,neighborPts.at(p))) {
	    double q0[nDims];
	    mesh->GetPoint(q,q0);
	    distances_vec.push_back(std::sqrt(vtkMath::Distance2BetweenPoints(p0,q0)));
	    neighborPts.at(p).push_back(q);
	    neighborPts.at(q).push_back(p);
	  }
	}
      }
    }
  }


  // Calculate desired statistic
  arma::vec distances(static_cast<int>(distances_vec.size()));
  for(unsigned int n = 0; n < distances.n_elem; n++) {
    distances(n) = distances_vec.at(n);
  }

  if(stat == "Min") {
    std::cout << distances.min() << std::endl;
  } else if(stat == "Max") {
    std::cout << distances.max() << std::endl;
  } else if(stat == "Mean") {
    std::cout << arma::mean(distances) << std::endl;
  } else if(stat=="Median") {
    std::cout << arma::median(distances) << std::endl;
  } else {
    std::cout << "Invalid statistic (stat = Min, Max, Mean, Median)" << std::endl;
  }

  return 0;
}
