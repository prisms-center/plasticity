//methods to allow for pre/post iteration level updates
#include "../../include/ellipticBVP.h"

//method called before each iteration
template <int dim>
void ellipticBVP<dim>::updateBeforeIteration(){
  //default method does nothing
}

//method called after each iteration
template <int dim>
void ellipticBVP<dim>::updateAfterIteration(){
  //default method does nothing
}

//method called after each iteration
template <int dim>
bool ellipticBVP<dim>::testConvergenceAfterIteration(){
  //default method resets solution to previously converged solution if resetIncrement flagis true
  if (resetIncrement){
    solution=oldSolution;
    solutionWithGhosts=oldSolution;

    resetIncrement=false;
    char buffer[100];
    sprintf(buffer,
	    "current increment reset by model. Restarting increment with loadFactorSetByModel: %12.6e\n",
	    loadFactorSetByModel);
    pcout << buffer;
    return false;
  }
  return true;
}
#include "../../include/ellipticBVP_template_instantiations.h"
