//solve Ax=b 

#ifndef SOLVELINEAR_ELLIPTICBVP_H
#define SOLVELINEAR_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//solve linear system of equations AX=b using iterative solver
template <int dim>
void ellipticBVP<dim>::solveLinearSystem(){
  vectorType completely_distributed_solution (locally_owned_dofs, mpi_communicator);
#ifdef linearSolverType
  SolverControl solver_control(maxLinearSolverIterations, relLinearSolverTolerance*residual.l2_norm());
  linearSolverType solver(solver_control, mpi_communicator);
  PETScWrappers::PreconditionBlockJacobi preconditioner(jacobian);
#else
  pcout << "\nError: solverType not defined. This is required for ELLIPTIC BVP.\n\n";
  exit (-1);
#endif
  //solve Ax=b
  try{
    solver.solve (jacobian, completely_distributed_solution, residual, preconditioner);
    char buffer[200];
    sprintf(buffer, 
	    "linear system solved in %3u iterations [final norm: %8.2e, initial norm: %8.2e, tol criterion: %8.2e]\n",
	    solver_control.last_step(),
	    solver_control.last_value(),
	    solver_control.initial_value(),
	    relLinearSolverTolerance);
    pcout << buffer;
  }
  catch (...) {
    pcout << "\nWarning: solver did not converge in "
	  << solver_control.last_step()
	  << " iterations as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";     
  }
  constraints.distribute (completely_distributed_solution);
  solution+=completely_distributed_solution;
  solutionWithGhosts=solution;
}

#endif
