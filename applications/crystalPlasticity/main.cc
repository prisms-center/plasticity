//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#include "../../include/crystalPlasticity.h"

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);

      ParameterHandler parameter_handler;

      std::list<std::string> args;
      for (int i=1; i<argc; ++i) args.push_back (argv[i]);

      if (args.size() == 0){
        std::cerr<<"Provide name of a parameter file."<<std::endl;
        exit (1);
      }

      const std::string parameter_file = args.front ();
      userInputParameters userInputs(parameter_file,parameter_handler);

      crystalPlasticity<3> problem(userInputs);

      //reading materials atlas files
      if(!userInputs.readExternalMesh){
        problem.orientations.loadOrientations(userInputs.grainIDFile,
  					    userInputs.headerLinesGrainIDFile,
  					    userInputs.grainOrientationsFile,
  					    userInputs.numPts,
  					    userInputs.span);}
      problem.orientations.loadOrientationVector(userInputs.grainOrientationsFile, userInputs.enableMultiphase, userInputs.additionalVoxelInfo);

      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Exception on processing: " << std::endl
		<< exc.what() << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      std::cerr << "Unknown exception!" << std::endl
		<< "Aborting!" << std::endl
		<< "----------------------------------------------------"
		<< std::endl;
      return 1;
    }

  return 0;
}
