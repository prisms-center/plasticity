//tension BVP
//general headers
#include <fstream>
#include <sstream>
#include <iostream>
using namespace std;

#include "../../include/continuumPlasticity.h"

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
      continuumPlasticity<3> problem(userInputs);

      //Read material parameters
      problem.properties.lambda = userInputs.lame_lambda;
      problem.properties.mu = userInputs.lame_mu;
      problem.properties.tau_y = userInputs.yield_stress;
      problem.properties.K = userInputs.strain_hardening;
      problem.properties.H = userInputs.kinematic_hardening;

      //Read pfunction names for strain energy density and yield functions
      problem.properties.strainEnergyModel = userInputs.strain_energy_function;
      problem.properties.yieldModel = userInputs.yield_function;
      problem.properties.isoHardeningModel = userInputs.iso_hardening_function;

      //Copy to properties the convergence criteria
      problem.properties.stopOnConvergenceFailure = userInputs.stopOnConvergenceFailure;
      problem.run();
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
