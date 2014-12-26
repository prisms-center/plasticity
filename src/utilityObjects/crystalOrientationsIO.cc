//class to read/write crystal orientations from/to external text files
#include <fstream>
#include <iostream>
#include <sstream>

template <int dim>
class crystalOrientationsIO{
public:
  crystalOrientationsIO();
  void loadOrientations(std::string _voxelFileName, 
			unsigned int headerLines,
			std::string _orientationFileName,
			unsigned int _numPts[], 
			double _stencil[]);
  void loadEulerAngles(std::string _eulerFileName);
  unsigned int getMaterialID(double _coords[]);
  void addToOutputOrientations(std::vector<double>& _orientationsInfo);
  void writeOutputOrientations();
  std::map<unsigned int, std::vector<double> > eulerAngles;
  std::vector<std::vector<double> > outputOrientations;
private:
  std::map<double,std::map<double, std::map<double, unsigned int> > > inputVoxelData;
  ConditionalOStream  pcout;  
};

//constructor
template <int dim>
crystalOrientationsIO<dim>::crystalOrientationsIO():
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{}

//addToOutputOrientations adds data to be written out to output oreintations file
template <int dim>
void crystalOrientationsIO<dim>::addToOutputOrientations(std::vector<double>& _orientationsInfo){
  outputOrientations.push_back(_orientationsInfo);
}

//writeOutputOreintations writes outputOrientations to file
template <int dim>
void crystalOrientationsIO<dim>::writeOutputOrientations(){
  pcout << "writing orientations data to file\n";
  std::string fileName("orientationsOutputProc");
  fileName += boost::to_string(Utilities::MPI::this_mpi_process(MPI_COMM_WORLD));
  std::ofstream file(fileName.c_str());
  char buffer[200];
  if (file.is_open()){
    for (std::vector<std::vector<double> >::iterator it = outputOrientations.begin() ; it != outputOrientations.end(); ++it){
      for (std::vector<double>::iterator it2 = it->begin() ; it2 != it->end(); ++it2){
	sprintf(buffer, "%8.2e ",*it2);
	file << buffer;
      }
      file << std::endl;
    }
    file.close();
  }
  else {
    pcout << "Unable to open file for writing orientations\n";
    exit(1);
  }
  
  //join files from all processors into a single file on processor 0 
  //and delete individual processor files 
  MPI_Barrier(MPI_COMM_WORLD);
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0){
     std::string fileName2("orientationsOutput");
      std::ofstream file2(fileName2.c_str());
      for (unsigned int proc=0; proc<Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD); proc++){
	std::string fileName3("orientationsOutputProc");
	fileName3 += boost::to_string(proc);
	std::ofstream file3(fileName3.c_str(), std::ofstream::in);
	file2 << file3.rdbuf();
	//delete file from processor proc
	remove(fileName3.c_str());
      }
      file2.close();
  }
}

//loadEulerAngles reads the orientation euler angles file   
template <int dim>
void crystalOrientationsIO<dim>::loadEulerAngles(std::string _eulerFileName){
 //check if dim==3
  if (dim!=3) {
    pcout << "loadEulerAngles only implemented for dim==3\n";
    exit(1);
  }

  //open data file
  std::ifstream eulerDataFile(_eulerFileName.c_str());
  //read data 
  std::string line; 
  double value; 
  unsigned int id;
  if (eulerDataFile.is_open()){
    pcout << "reading orientation euler angles file\n";
    //skip header lines
    for (unsigned int i=0; i<1; i++) std::getline (eulerDataFile,line);
    //read data
    while (getline (eulerDataFile,line)){
      std::stringstream ss(line);
      unsigned int id; 
      ss >> id; 
      //double temp;
      //ss >> temp;
      eulerAngles[id]=std::vector<double>(3);
      ss >> eulerAngles[id][0];
      ss >> eulerAngles[id][1];
      ss >> eulerAngles[id][2];
      //pcout << id << " " << eulerAngles[id][0] << " " << eulerAngles[id][1] << " " << eulerAngles[id][2] << std::endl;
    }
  }
  else{
    pcout << "Unable to open eulerDataFile\n"; 
    exit(1);
  }
}

//loadOrientations reads the voxel data file and orientations file   
template <int dim>
void crystalOrientationsIO<dim>::loadOrientations(std::string _voxelFileName, 
						  unsigned int headerLines,
						  std::string _orientationFileName,
						  unsigned int _numPts[], 
						  double _stencil[]){
  //check if dim==3
  if (dim!=3) {
    pcout << "voxelDataFile read only implemented for dim==3\n";
    exit(1);
  }

  //open voxel data file
  std::ifstream voxelDataFile(_voxelFileName.c_str());
  //read voxel data 
  std::string line; 
  double value; 
  unsigned int id;
  if (voxelDataFile.is_open()){
    pcout << "reading voxel data file\n";
    //skip header lines
    for (unsigned int i=0; i<headerLines; i++) std::getline (voxelDataFile,line);
    //read data
    for (unsigned int x=0; x<_numPts[0]; x++){
      double xVal=x*_stencil[0];
      if (inputVoxelData.count(xVal)==0) inputVoxelData[xVal]=std::map<double, std::map<double, unsigned int> >();
      for (unsigned int y=0; y<_numPts[1]; y++){
	double yVal=y*_stencil[1];
	if (inputVoxelData[xVal].count(yVal)==0) inputVoxelData[xVal][yVal]=std::map<double, unsigned int>();
	std::getline (voxelDataFile,line);
	std::stringstream ss(line);
	for (unsigned int z=0; z<_numPts[2]; z++){
	  double zVal=z*_stencil[2];
	  ss >> inputVoxelData[xVal][yVal][zVal];
	  //pcout <<  inputVoxelData[xVal][yVal][zVal] << " ";
	}
	//pcout << "\n";
      }
    }
  }
  else {
    pcout << "Unable to open file voxelDataFile\n"; 
    exit(1);
  }
}

//return materialID closest to given (x,y,z)
template <int dim>
unsigned int crystalOrientationsIO<dim>::getMaterialID(double _coords[]){
  if (inputVoxelData.size()==0){
    pcout << "inputVoxelData not initialized\n";
    exit(1);
  }

  //find nearest point
  //iterator to nearest x slice
  std::map<double,std::map<double, std::map<double, unsigned int> > >::iterator itx=inputVoxelData.lower_bound(_coords[0]); 
  if(itx != inputVoxelData.end()) --itx;
  //iterator to nearest y slice
  std::map<double, std::map<double, unsigned int> >::iterator ity=itx->second.lower_bound(_coords[1]); 
  if(ity == itx->second.end()) --ity;
  //iterator to nearest z slice
  std::map<double, unsigned int>::iterator itz=ity->second.lower_bound(_coords[2]); 
  if(itz != ity->second.end()) --itz;
  return itz->second; 
}
