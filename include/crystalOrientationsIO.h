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
  void loadOrientationVector(std::string _eulerFileName);
  unsigned int getMaterialID(double _coords[]);
  void addToOutputOrientations(std::vector<double>& _orientationsInfo);
  void writeOutputOrientations(bool writeOutput,std::string outputDirectory);
  std::map<unsigned int, std::vector<double> > eulerAngles;
  std::vector<std::vector<double> > outputOrientations;
private:
  std::map<double,std::map<double, std::map<double, unsigned int> > > inputVoxelData;
  dealii::ConditionalOStream  pcout;
};
