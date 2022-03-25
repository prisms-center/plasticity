//solve non linear system of equations for ellipticBVP class
#include "../../include/ellipticBVP.h"

//solve non-linear system of equations
template <int dim>
bool ellipticBVP<dim>::solveNonLinearSystem(){
  //residuals
  double relNorm=1.0, initialNorm=1.0e-16, currentNorm=0.0;

  //non linear iterations
  char buffer[200];
  currentIteration=0;
  while (currentIteration < userInputs.maxNonLinearIterations){
    //call updateBeforeIteration, if any
    updateBeforeIteration();

    //Calling assemble
    computing_timer.enter_section("assembly");
    assemble();    
    computing_timer.exit_section("assembly");

    if (!resetIncrement){
      //Calculate residual norms and check for convergence
      currentNorm=residual.l2_norm();
      initialNorm=std::max(initialNorm, currentNorm);
      relNorm=currentNorm/initialNorm;
      //print iteration information
      sprintf(buffer,
        "nonlinear iteration %3u [current residual: %8.2e, initial residual: %8.2e, relative residual: %8.2e]\n",
        currentIteration,
        currentNorm,
        initialNorm,
        relNorm);
        pcout << buffer;
        //if not converged, solveLinearSystem Ax=b
        computing_timer.enter_section("solve");
        solveLinearSystem(constraints, jacobian, residual, solution, solutionWithGhosts, solutionIncWithGhosts);
        computing_timer.exit_section("solve");
        currentIteration++;
      }
      //call updateAfterIteration, if any
      updateAfterIteration();
    }

    //check if maxNonLinearIterations reached

    //update old solution to new converged solution
    oldSolution=solution;
    return true;
  }
  #include "../../include/ellipticBVP_template_instantiations.h"
