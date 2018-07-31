//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

//parameters
#include "parameters.h"

//FCC model header
#include "../../../../src/materialModels/crystalPlasticity/fcc/model.h"


//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);
      crystalPlasticity<3> problem;



      //reading materials atlas files
      double stencil[3]={spanX/(numPts[0]-1), spanY/(numPts[1]-1), spanZ/(numPts[2]-1)}; // Dimensions of voxel
      problem.orientations.loadOrientations(grainIDFile,
					    headerLinesGrainIDFile,
					    grainOrientationsFile,
					    numPts,
					    stencil);
      problem.orientations.loadOrientationVector("orientations.txt");
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
