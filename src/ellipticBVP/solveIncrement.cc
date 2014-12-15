//solve increment method for ellipticBVP class

#ifndef SOLVEINC_ELLIPTICBVP_H
#define SOLVEINC_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//solve linear system of equations AX=b using iterative solver
template <int dim>
void ellipticBVP<dim>::solveIncrement(){
  vectorType completely_distributed_solution (locally_owned_dofs, mpi_communicator);
#ifdef solverType
  SolverControl solver_control(maxSolverIterations, relSolverTolerance*residual.l2_norm());
  solverType solver(solver_control, mpi_communicator);
  PETScWrappers::PreconditionBlockJacobi preconditioner(jacobian);
  try{
    solver.solve (jacobian, completely_distributed_solution, residual, preconditioner);
    char buffer[200];
    sprintf(buffer, 
	    "solved in %3u iterations (final residual: %10.4e, initial residual: %10.4e, tol: %10.4e)\n",
	    solver_control.last_step(),
	    solver_control.last_value(),
	    solver_control.initial_value(),
	    solver_control.tolerance());
    pcout << buffer;
  }
  catch (...) {
    pcout << "\nWarning: solver did not converge in "
	  << solver_control.last_step()
	  << " iterations as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";     
  }
  constraints.distribute (completely_distributed_solution);
  solution+=completely_distributed_solution;
  solutionLocal=solution;

#else
  pcout << "\nError: solverType not defined. This is required for ELLIPTIC BVP.\n\n";
  exit (-1);
#endif
}

#endif
