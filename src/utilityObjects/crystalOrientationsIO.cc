#include "../../include/crystalOrientationsIO.h"

//constructor
template <int dim>
crystalOrientationsIO<dim>::crystalOrientationsIO():
pcout (std::cout, dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{}

  //loadOrientationVector reads the orientation euler angles file
  template <int dim>
  void crystalOrientationsIO<dim>::loadOrientationVector(std::string _eulerFileName, bool _enableMultiphase, unsigned int _numVoxelData){
    //check if dim==3
    if (dim!=3) {
      pcout << "loadOrientationVector only implemented for dim==3\n";
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
      unsigned int counter=0;
      if (_enableMultiphase){
        counter=1;
      }
      unsigned int numberOfAdditionalVoxelInfo=_numVoxelData;
      while (getline (eulerDataFile,line)){
        std::stringstream ss(line);

        ss >> id;
        //double temp;
        //ss >> temp;
        eulerAngles[id]=std::vector<double>(dim+counter+numberOfAdditionalVoxelInfo);
        ss >> eulerAngles[id][0];
        ss >> eulerAngles[id][1];
        ss >> eulerAngles[id][2];
        if (_enableMultiphase){
          ss >> eulerAngles[id][dim];
        }
        if (numberOfAdditionalVoxelInfo>0){
          for (unsigned int j=0;j<numberOfAdditionalVoxelInfo;j++){
            ss >> eulerAngles[id][dim+counter+j];
          }
        }

      }
    }
    else{
      pcout << "Unable to open eulerDataFile\n";
      exit(1);
    }

///////////Number Of Grains In The Microstructure/////////
    numberOfGrains=id;
  }

  //loadOrientations reads the voxel data file and orientations file
  template <int dim>
  void crystalOrientationsIO<dim>::loadOrientations(std::string _voxelFileName,
    unsigned int headerLines,
    std::string _orientationFileName,
    std::vector<unsigned int> _numPts,
    std::vector<double> _span){
      //check if dim==3
      if (dim!=3) {
        pcout << "voxelDataFile read only implemented for dim==3\n";
        exit(1);
      }

      double _stencil[3]={_span[0]/(_numPts[0]), _span[1]/(_numPts[1]), _span[2]/(_numPts[2])}; // Dimensions of voxel

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
          double xVal=x*_stencil[0]+_stencil[0]/2;
          if (inputVoxelData.count(xVal)==0) inputVoxelData[xVal]=std::map<double, std::map<double, unsigned int> >();
          for (unsigned int y=0; y<_numPts[1]; y++){
            double yVal=y*_stencil[1]+_stencil[1]/2;
            if (inputVoxelData[xVal].count(yVal)==0) inputVoxelData[xVal][yVal]=std::map<double, unsigned int>();
            std::getline (voxelDataFile,line);
            std::stringstream ss(line);
            for (unsigned int z=0; z<_numPts[2]; z++){
              double zVal=z*_stencil[2]+_stencil[2]/2;
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

  double dist_to_lower,dist_to_upper;
  //find nearest point
  //iterator to nearest x slice
  std::map<double,std::map<double, std::map<double, unsigned int> > >::iterator itx;
  std::map<double,std::map<double, std::map<double, unsigned int> > >::iterator itx_upper=inputVoxelData.lower_bound(_coords[0]);
  std::map<double,std::map<double, std::map<double, unsigned int> > >::iterator itx_lower = itx_upper; itx_lower--;
  std::map<double,std::map<double, std::map<double, unsigned int> > >::iterator itx_begin=inputVoxelData.begin();

  if (itx_upper == inputVoxelData.end()) {
    itx=itx_lower;
  }
  else if (itx_upper==itx_begin) {
    itx=itx_upper;
  }
  else{
    dist_to_lower = std::abs(itx_lower->first - _coords[0]);
    dist_to_upper = std::abs(itx_upper->first - _coords[0]);
    if (dist_to_upper >= dist_to_lower) {
      itx=itx_lower;
    }
    else {
      itx=itx_upper;
    }
  }
  if(itx == inputVoxelData.end()) --itx;

  //iterator to nearest y slice
  std::map<double, std::map<double, unsigned int> >::iterator ity;
  std::map<double, std::map<double, unsigned int> >::iterator ity_upper=itx->second.lower_bound(_coords[1]);
  std::map<double, std::map<double, unsigned int> >::iterator ity_lower=ity_upper; ity_lower--;
  std::map<double, std::map<double, unsigned int> >::iterator ity_begin=itx->second.begin();

  if (ity_upper == itx->second.end()) {
    ity=ity_lower;
  }
  else if (ity_upper==ity_begin) {
    ity=ity_upper;
  }
  else{
    dist_to_lower = std::abs(ity_lower->first - _coords[1]);
    dist_to_upper = std::abs(ity_upper->first - _coords[1]);
    if (dist_to_upper >= dist_to_lower) {
      ity=ity_lower;
    }
    else {
      ity=ity_upper;
    }
  }

  if(ity == itx->second.end()) --ity;
  //iterator to nearest z slice
  std::map<double, unsigned int>::iterator itz;
  std::map<double, unsigned int>::iterator itz_upper=ity->second.lower_bound(_coords[2]);
  std::map<double, unsigned int>::iterator itz_lower=itz_upper; itz_lower--;
  std::map<double, unsigned int>::iterator itz_begin=ity->second.begin();
  if (itz_upper == ity->second.end()) {
    itz=itz_lower;
  }
  else if (itz_upper==itz_begin) {
    itz=itz_upper;
  }
  else{
    dist_to_lower = std::abs(itz_lower->first - _coords[2]);
    dist_to_upper = std::abs(itz_upper->first - _coords[2]);
    if (dist_to_upper >= dist_to_lower) {
      itz=itz_lower;
    }
    else {
      itz=itz_upper;
    }
  }
  if(itz == ity->second.end()) --itz;
  return itz->second;
}

    #include "../../include/crystalOrientationsIO_template_instantiations.h"
