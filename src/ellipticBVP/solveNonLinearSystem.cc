//solve non linear system of equations for ellipticBVP class

#ifndef SOLVENONLINEAR_ELLIPTICBVP_H
#define SOLVENONLINEAR_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//solve non-linear system of equations
template <int dim>
bool ellipticBVP<dim>::solveNonLinearSystem(){
  //residuals
  double relNorm=1.0, initialNorm=1.0e-16, currentNorm=0.0;

  //non linear iterations
  char buffer[200];
  currentIteration=0;
  while (currentIteration < maxNonLinearIterations){
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
      
      //check for convergence in abs tolerance
      if (currentNorm<absNonLinearTolerance){
	pcout << "nonlinear iterations converged in absolute norm\n";
	break; 
      }
      //check for convergence in relative tolerance
      else if(relNorm<relNonLinearTolerance){
	pcout << "nonlinear iterations converged in relative norm\n";
	break; 
      }
      
      //if not converged, solveLinearSystem Ax=b
      computing_timer.enter_section("solve");
      solveLinearSystem(constraints, jacobian, residual, solution, solutionWithGhosts, solutionIncWithGhosts);
      computing_timer.exit_section("solve");
      currentIteration++;
    }
    
    //convergence test after iteration
    bool convFlag=testConvergenceAfterIteration();
    if (!convFlag) {return false;}
    
    //call updateAfterIteration, if any
    updateAfterIteration();
  }
  
  //check if maxNonLinearIterations reached
  if (currentIteration >= maxNonLinearIterations){
    pcout <<  "nonlinear iterations did not converge in maxNonLinearIterations\n";
    if (stopOnConvergenceFailure) {exit (1);}
    else {pcout << "stopOnConvergenceFailure==false, so marching ahead\n";}
  }

  //update old solution to new converged solution
  oldSolution=solution;
  return true;
}
#endif
