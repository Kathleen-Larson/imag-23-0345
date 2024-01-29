#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPointData.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkFloatArray.h>
#include <vtkXMLPolyDataWriter.h>
#include <vtkIdTypeArray.h>
#include <vtkIdList.h>
#include <vtkTriangleFilter.h>
#include <vtkCellArray.h>
#include <vtkSelectEnclosedPoints.h>
#include <vtkClipPolyData.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkFillHolesFilter.h>
#include <vtkCleanPolyData.h>
#include <vtkFeatureEdges.h>
#include <vtkStripper.h>
#include <vtkAppendPolyData.h>
#include <vtkCenterOfMass.h>
#include <vtkDelaunay3D.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkDataArray.h>

#include <itkImage.h>
#include <itkImageFileReader.h>
#include <itkImageFileWriter.h>
#include <itkChangeInformationImageFilter.h>



const unsigned int nDims = 3;

typedef itk::Image<float,nDims> FloatImageType;
typedef itk::ImageFileReader<FloatImageType> FloatImageFileReaderType;
typedef itk::ImageFileWriter<FloatImageType> FloatImageFileWriterType;



FloatImageType::Pointer readFloatImage(std::string filename) {
  FloatImageFileReaderType::Pointer imageReader = FloatImageFileReaderType::New();
  imageReader->SetFileName(filename.c_str());
  imageReader->Update();

  FloatImageType::Pointer image = imageReader->GetOutput();
  return image;
}



vtkSmartPointer<vtkPolyData> readXMLPolyData(std::string filename) {
  vtkSmartPointer<vtkXMLPolyDataReader> meshreader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
  meshreader->SetFileName(filename.c_str());
  meshreader->Update();
  vtkSmartPointer<vtkPolyData> mesh = meshreader->GetOutput();
  return mesh;
}



void writeXMLPolyData(vtkSmartPointer<vtkPolyData> mesh, std::string filename) {
  vtkSmartPointer<vtkXMLPolyDataWriter> meshWriter = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  meshWriter->SetInputData(mesh);
  meshWriter->SetFileName(filename.c_str());
  meshWriter->Write();
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



unsigned int MajorityVote(std::vector<unsigned int> votes, std::vector<int> bounds) {
  double vtkPnt[nDims];
  unsigned int nVotes = votes.size();
  std::vector< unsigned int> candidates;
  std::vector< unsigned int> count;

  unsigned int nCandidates = 0;
  count.resize(nVotes,0);
  for(unsigned int i = 0; i < nVotes; i++) {
    unsigned int current = votes[i];
    if(!isValid(current,bounds)) {
      continue;
    }
    bool found = false;
    for(unsigned int j = 0; j < nCandidates; j++) {
      if(current == candidates[j]) {
        found = true;
        count[j]++;
      }
    }
    if(!found) {
      candidates.push_back(current);
      count[nCandidates]++;
      nCandidates++;
    }
  }

  unsigned int maxVote = 0, choice = 0;
  for(unsigned int i = 0; i < nCandidates; i++) {
    if(count[i] > maxVote) {
      maxVote = count[i];
      choice = candidates[i];
    }
  }
  return choice ;
}



vtkSmartPointer<vtkFloatArray> AssignLabels(vtkSmartPointer<vtkPolyData> mesh, FloatImageType::Pointer image, std::string labelArrayName, std::vector<int> bounds) {

  unsigned int nPts = mesh->GetNumberOfPoints();
  
  vtkSmartPointer<vtkPolyDataNormals> normalsGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
  normalsGenerator->SetInputData(mesh);
  normalsGenerator->ComputePointNormalsOn();
  normalsGenerator->ComputeCellNormalsOff();
  normalsGenerator->FlipNormalsOff();
  normalsGenerator->Update();
  vtkSmartPointer<vtkDataArray> meshNormals = normalsGenerator->GetOutput()->GetPointData()->GetNormals();
  
  vtkSmartPointer<vtkFloatArray> labels = vtkSmartPointer<vtkFloatArray>::New() ;
  labels->SetNumberOfComponents(1);
  labels->SetNumberOfValues(nPts);
  labels->SetName(labelArrayName.c_str());

  for(unsigned int p = 0; p < nPts; p++) {
    double vtkPnt[3];
    mesh->GetPoint(p,vtkPnt);

    itk::ContinuousIndex<float,3> current;
    FloatImageType::IndexType currentIndex;
    current[0] = -vtkPnt[0];
    current[1] = -vtkPnt[1];
    current[2] = vtkPnt[2];

    bool inside = image->TransformPhysicalPointToIndex(current,currentIndex);
    if(inside) {
      FloatImageType::PixelType pixel = image->GetPixel(currentIndex);

      int nTry = 0;
      while(nTry < 3) {
	//std::cout << pixel << std::endl;

	if(isValid(pixel,bounds)) {
	  labels->SetValue(p,pixel);
	  nTry = 3;
	}
	else {
	  double vtkPntNew[nDims], normal[nDims];
	  meshNormals->GetTuple(p,normal);
	  
	  for(unsigned int d = 0; d < nDims; d++) {
	    vtkPntNew[d] = vtkPnt[d] - normal[d];
	  }
	  current[0] = -vtkPntNew[0];
	  current[1] = -vtkPntNew[1];
	  current[2] = vtkPntNew[2];
	  
	  image->TransformPhysicalPointToIndex(current,currentIndex);
	  pixel = image->GetPixel(currentIndex);
	  nTry++;
	}
	//std::cout << nTry << std::endl;
      }
    }
    else {
      std::cout << "not inside. huh. " << std::endl;
    }
  }

  return labels;
}



unsigned int FillHoles(vtkSmartPointer<vtkPolyData> mesh, vtkSmartPointer<vtkFloatArray> labels, std::vector<int> bounds) {
  mesh->BuildLinks();
  unsigned int nPts = labels->GetNumberOfTuples();

  vtkSmartPointer <vtkIdList> cellIds = vtkSmartPointer<vtkIdList>::New();

  unsigned int nFail = 0;
  for(unsigned int p = 0; p < nPts; p++) {
    std::vector<unsigned int> candidates, neighbors, uniqueNeighbors;
    float pixel = labels->GetValue(p);
    if(!isValid(pixel,bounds)) {
      neighbors.empty();
      mesh->GetPointCells(p,cellIds);
      unsigned int nCells = cellIds->GetNumberOfIds();
      for(unsigned int i = 0; i < nCells; i++) {
        vtkSmartPointer<vtkIdList> ptIds = vtkSmartPointer<vtkIdList>::New();
        mesh->GetCellPoints(cellIds->GetId(i),ptIds);
        vtkIdType nPts = ptIds->GetNumberOfIds();
        for(unsigned j = 0; j < nPts; j++) {
          if(ptIds->GetId(j) == p) continue;
          neighbors.push_back(ptIds->GetId(j));
        }
      }

      unsigned int nNeighbors = 0;
      for(unsigned int i = 0; i < neighbors.size(); i++) {
        bool found = false;
        for(unsigned int j = 0; j < nNeighbors; j++) {
          if(uniqueNeighbors[j] == neighbors[i]) {
            found = true;
            continue;
          }
        }
        if(!found) {
          uniqueNeighbors.push_back(neighbors[i]);
          nNeighbors++;
        }
      }
      candidates.resize(nNeighbors,0);
      for(unsigned int n = 0; n < nNeighbors; n++) {
        candidates[n] = labels->GetValue(uniqueNeighbors[n]);
      }
      pixel = MajorityVote(candidates,bounds);
      labels->SetValue(p, pixel);
      if(!isValid(pixel,bounds))
        nFail++;
    }
  }

  return nFail;
}



int main(int argc, char *argv[]) {

  if(argc < 3) {
    std::cout << "Usage: " << std::endl;
    std::cout << "inputMeshFileName inputLabelMapFileName outputMeshFileName labelArrayName bounds" << std::endl;
    return -1;
  }

  std::string inputMeshFileName = argv[1];
  std::string labelFileName = argv[2];
  std::string outputMeshFileName = argv[3];
  std::string labelArrayName = argv[4];

  vtkSmartPointer<vtkPolyData> inputMesh = readXMLPolyData(inputMeshFileName);
  vtkSmartPointer<vtkPolyData> parcellatedMesh = vtkSmartPointer<vtkPolyData>::New();
  parcellatedMesh->DeepCopy(inputMesh);

  FloatImageType::Pointer image = readFloatImage(labelFileName);

  unsigned int nBounds = argc-5;
  if(nBounds == 0) {
    std::cout << "Must provide at least 1 set of bounds for parcellation" << std::endl;
    return -1;
  }
  if(nBounds % 2 != 0) {
    std::cout << "Must provide even number of bounds for parcellation" << std::endl;
    return -1;
  }

  std::vector<int> bounds;
  for(unsigned int n = 0; n < nBounds; n++) {
    bounds.push_back(atoi(argv[argc - nBounds + n]));
  }

  vtkSmartPointer<vtkFloatArray> labels = AssignLabels(parcellatedMesh,image,labelArrayName,bounds);

  unsigned int it = 0;
  unsigned int nHoles = 1;
  while(nHoles = FillHoles(parcellatedMesh,labels,bounds)) {
    it++;
    if( it > 100 ) {
      break;
    }
  }
  std::cout << nHoles << " remaining after " << it << " iterations" << std::endl;

  parcellatedMesh->GetPointData()->AddArray(labels);
  writeXMLPolyData(parcellatedMesh,outputMeshFileName);

  return 0;
}
