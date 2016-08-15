//methods to allow for pre/post iteration level updates

#ifndef UPDATEITERATION_H
#define UPDATEITERATION_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

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

#endif
