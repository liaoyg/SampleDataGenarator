#include <vtkGenericDataObjectReader.h>
#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointData.h>
#include <vtkCellArray.h>
#include <vtkDataArray.h>
#include <vtkPointSource.h>
#include <vtkOctreePointLocator.h>
#include <vtkDataSetCollection.h>
#include <vtkMath.h>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

int main ( int argc, char *argv[] )
{
  // Ensure a filename was specified
  if(argc != 2)
    {
    std::cerr << "Usage: " << argv[0] << " InputFilename" << endl;
    return EXIT_FAILURE;
    }

  // Get the filename from the command line
  std::string inputFilename = argv[1];
  int fileNum = 1;

  for (int i = 0; i < fileNum; i++) {
	  stringstream ss;
	  ss << inputFilename << "." << i+1 << ".vtk" << endl;
	  string fileName;
	  ss >> fileName;
	  // Get all data from the file
	  vtkSmartPointer<vtkGenericDataObjectReader> reader =
		  vtkSmartPointer<vtkGenericDataObjectReader>::New();
	  reader->SetFileName(fileName.c_str());
	  reader->SetReadAllScalars(true);
	  reader->SetReadAllVectors(true);
	  reader->Update();

	  //output parameter
	  int dim = 64;

	  if (reader->IsFileUnstructuredGrid())
	  {
		  vtkUnstructuredGrid* output = reader->GetUnstructuredGridOutput();

		  //output->GetAttributes();
		  vtkFieldData* fd = output->GetAttributesAsFieldData(0);
		  vtkPointData* pd = output->GetPointData();
		  //Create the tree
		  vtkSmartPointer<vtkOctreePointLocator> octree =
			  vtkSmartPointer<vtkOctreePointLocator>::New();
		  octree->SetDataSet(output);
		  octree->BuildLocator();
		  std::cout << "Number of points in tree: " << octree->GetDataSet()->GetNumberOfPoints() << std::endl;

		  //Get vectors from data point
		  
		  //ouput to raw data
		  //string vectorName = "vorticity";
		  int pd_array_num = pd->GetNumberOfArrays();
		  for (int aryId = 0; aryId < pd_array_num; aryId++)
		  {
			  //t cd_array_num = c
			  
			  const char* vectorName = pd->GetArrayName(aryId);

			  if (vectorName == NULL)
				  std::cout << "No such arry named: " << vectorName << std::endl;
			  std::cout << "Array Name: " << vectorName << std::endl;
			  ss.clear();
			  ss << inputFilename << "_" << dim << "_" << i << "_" << vectorName << ".raw" << endl;
			  std::string outputFilename;
			  ss >> outputFilename;
			  ofstream out(outputFilename, ios::binary);
			  std::cout << "output is a unstructured data: " << outputFilename << std::endl;

			  if (!pd->HasArray(vectorName))
				  std::cout << "No such arry named: " << vectorName << std::endl;
			  vtkDataArray* vectorsflds = pd->GetArray(vectorName);
			  if (vectorsflds->GetNumberOfComponents() != 3)
				  continue;
			  double* bound = output->GetBounds();
			  vector<double> datasetbox(bound, bound + 6);

			  for (int i = 0; i < dim; i++) { //z
				  double z = double((datasetbox[5] - datasetbox[4])*(i + 0.5) / dim + datasetbox[4]);
				  for (int j = 0; j < dim; j++) { //y
					  double y = double((datasetbox[3] - datasetbox[2])*(j + 0.5) / dim + datasetbox[2]);
					  for (int k = 0; k < dim; k++) //x
					  {
						  double x = double((datasetbox[1] - datasetbox[0])*(k + 0.5) / dim + datasetbox[0]);
						  double point[3] = { x, y, z };
						  vtkSmartPointer<vtkIdList> result =
							  vtkSmartPointer<vtkIdList>::New();
						  double radius = (datasetbox[1] - datasetbox[0]) / dim;
						  //octree->FindClosestNPoints(4, point, result); //find 8 closest point
						  octree->FindPointsWithinRadius(radius, point, result);
						  //int closestId = octree->FindClosestPoint(point);
						  float vectorres[3] = { 0, 0, 0 };
						  double accudist = 0;
						  for (int n = 0; n < result->GetNumberOfIds(); n++) {
							  double* referpt = output->GetPoint(result->GetId(n));
							  double dis = vtkMath::Distance2BetweenPoints(referpt, point);
							  double* refervctor = vectorsflds->GetTuple3(result->GetId(n));
							  for (int coord = 0; coord < 3; coord++)
								  vectorres[coord] += refervctor[coord] * dis;
							  accudist += dis;
						  }
						  if (accudist != 0) {
							  for (int coord = 0; coord < 3; coord++)
								  vectorres[coord] /= accudist;
						  }
						  //double* vectorclose = output->GetPoint(closestId);
						  vtkMath::Normalize(vectorres);
						  for (int coord = 0; coord < 3; coord++) {
							  //vectorres[coord] = vectorclose[coord];
							  out.write(reinterpret_cast<char*>(&vectorres[coord]), sizeof(vectorres[coord]));
						  }
						  //out2 << vectorres[0] << " " << vectorres[1] << " " << vectorres[2] << " " << endl;
					  }
				  }
			  }
			  out.close();
			  //out2.close();
		  }
	  }
  }
  return EXIT_SUCCESS;
}
