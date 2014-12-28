//solve Ax=b 

#ifndef SOLVELINEAR_ELLIPTICBVP_H
#define SOLVELINEAR_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//solve linear system of equations AX=b using iterative solver
template <int dim>
void ellipticBVP<dim>::solveLinearSystem(){
  vectorType completely_distributed_solutionInc (locally_owned_dofs, mpi_communicator);
#ifdef linearSolverType
  SolverControl solver_control(maxLinearSolverIterations, relLinearSolverTolerance*residual.l2_norm());
  linearSolverType solver(solver_control, mpi_communicator);
  PETScWrappers::PreconditionJacobi preconditioner(jacobian);
#else
  pcout << "\nError: solverType not defined. This is required for ELLIPTIC BVP.\n\n";
  exit (-1);
#endif
  //solve Ax=b
  try{
    solver.solve (jacobian, completely_distributed_solutionInc, residual, preconditioner);
    char buffer[200];
    sprintf(buffer, 
	    "linear system solved in %3u iterations\n",
	    solver_control.last_step());
    pcout << buffer;
  }
  catch (...) {
    pcout << "\nWarning: solver did not converge in "
	  << solver_control.last_step()
	  << " iterations as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";     
  }
  constraints.distribute (completely_distributed_solutionInc);
	solutionIncWithGhosts=completely_distributed_solutionInc;
  solution+=completely_distributed_solutionInc; 
  solutionWithGhosts=solution;
}

#endif
