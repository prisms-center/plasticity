//tension BVP
//general headers
#include <fstream>
#include <sstream>

//parameters file
#include "parameters.h"

//dealIIheaders
#include "../../../src/materialModels/continuumPlasticity/continuumPlasticity.h"

//Specify Dirichlet boundary conditions 
template <int dim>
void continuumPlasticity<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
  //back boundary:   u_y=u_z=0
  if (node[0] == 0.0){
    if (dof!=0) {flag=true; value=0.0;}
  }
  //front boundary:  u_y=u_z=0
  if (node[0] == 5.0){
    if (dof!=0) {flag=true; value=0.0;}
  }
  //left boundary:   u_y=u_z=0
  if (node[1] == 0.0){
    if (dof!=0) {flag=true; value=0.0;}
  }
  //right boundary:   u_y=u_z=0
  if (node[1] == 1.0){
    if (dof!=0) {flag=true; value=0.0;}
  }
  //bottom boundary: u=0
  if (node[2] == 0.0){
    {flag=true; value=0.0;}
  }
  //top boundary: u_x=g, u_y=u_z=0
  if (node[2] == 1.0){
    if (dof==0) {flag=true; value=totalDisplacement/totalNumIncrements;}
    else {flag=true; value=0.0;}
  }
}


//main
int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);
      continuumPlasticity<3> problem;

      //Read material parameters
      problem.properties.lambda = lame_lambda;
      problem.properties.mu = lame_mu;
      problem.properties.tau_y = yield_stress;
      problem.properties.K = strain_hardening;
      problem.properties.H = kinematic_hardening;

      //Read pfunction names for strain energy density and yield functions
      problem.properties.strainEnergyModel = strain_energy_function;
      problem.properties.yieldModel = yield_function;
      problem.properties.isoHardeningModel = iso_hardening_function;

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

