//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;
 
//parameters
#include "parameters.h"

//FCC model header
#include "../../../../src/materialModels/crystalPlasticity/hcp/model.h"
 
//Specify Dirichlet boundary conditions 
template <int dim>
void crystalPlasticity<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
  //back boundary:   u_x=0
  if (node[0] == 0.0){
    if (dof==0) {flag=true; value=0.0;}
  }
  //front boundary:  u_x=0.001
  if (node[0] == spanX){
    if (dof==0) {flag=true; value=0.001;}
  }
  //left boundary:   u_y=0
  if (node[1] == 0.0){
    if (dof==1) {flag=true; value=0.0;}
  }
  //bottom boundary: u_z=0
  if (node[2] == 0.0){
    if (dof==2) {flag=true; value=0.0;}
  }
}
  
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
      problem.orientations.loadOrientationVector(grainOrientationsFile);
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

