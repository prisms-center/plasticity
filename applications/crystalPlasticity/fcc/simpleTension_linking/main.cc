//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#include "../../../../include/crystalPlasticity.h"

//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);

      ParameterHandler parameter_handler;

      userInputParameters userInputs("prm.in",parameter_handler);

      crystalPlasticity<3> problem(userInputs);

      //reading materials atlas files
      problem.orientations.loadOrientations(userInputs.grainIDFile,
					    userInputs.headerLinesGrainIDFile,
					    userInputs.grainOrientationsFile,
					    userInputs.numPts,
					    userInputs.span);
      problem.orientations.loadOrientationVector(userInputs.grainOrientationsFile);
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
