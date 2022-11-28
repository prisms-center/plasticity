//class to read/write crystal orientations from/to external text files
#include <fstream>
#include <iostream>
#include <sstream>
#include "dealIIheaders.h"


template <int dim>
class crystalOrientationsIO{
public:
  crystalOrientationsIO();
  void loadOrientations(std::string _voxelFileName,
			unsigned int headerLines,
			std::string _orientationFileName,
			std::vector<unsigned int> _numPts,
			std::vector<double> _span);
  void loadOrientationVector(std::string _eulerFileName, bool _enableMultiphase, unsigned int _numVoxelData);
  unsigned int getMaterialID(double _coords[]);
  std::map<unsigned int, std::vector<double> > eulerAngles;
  unsigned int numberOfGrains;
private:
  std::map<double,std::map<double, std::map<double, unsigned int> > > inputVoxelData;
  dealii::ConditionalOStream  pcout;
};
