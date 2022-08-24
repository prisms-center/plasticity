//solve Ax=b
#include "../../include/ellipticBVP.h"

//solve linear system of equations AX=b using iterative solver
template <int dim>

#if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
void ellipticBVP<dim>::solveLinearSystem(ConstraintMatrix& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts){
#else
void ellipticBVP<dim>::solveLinearSystem(AffineConstraints<double>& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts){
#endif

  vectorType completely_distributed_solutionInc (locally_owned_dofs, mpi_communicator);
  SolverControl solver_control(userInputs.maxLinearSolverIterations, userInputs.relLinearSolverTolerance*b.l2_norm());
  PETScWrappers::SolverBiCG solver(solver_control, mpi_communicator);
  PETScWrappers::PreconditionJacobi preconditioner(A);

  //solve Ax=b
  try{
    solver.solve (A, completely_distributed_solutionInc, b, preconditioner);
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
  constraintmatrix.distribute (completely_distributed_solutionInc);
  dxGhosts=completely_distributed_solutionInc;
  x+=completely_distributed_solutionInc;
  xGhosts=x;
}

//solve linear system of equations AX=b using iterative solver
//This method is for solving the linear system of equations that arise in
//the projection of scalar post-processing fields
template <int dim>
#if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
void ellipticBVP<dim>::solveLinearSystem2(ConstraintMatrix& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts){
#else
void ellipticBVP<dim>::solveLinearSystem2(AffineConstraints<double>& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts){
#endif
  vectorType completely_distributed_solutionInc (locally_owned_dofs_Scalar, mpi_communicator);
  SolverControl solver_control(userInputs.maxLinearSolverIterations, userInputs.relLinearSolverTolerance*b.l2_norm());
  PETScWrappers::SolverCG solver(solver_control, mpi_communicator);
  PETScWrappers::PreconditionJacobi preconditioner(A);

  //solve Ax=b
  try{
    solver.solve (A, completely_distributed_solutionInc, b, preconditioner);
  }
  catch (...) {
    pcout << "\nWarning: solver did not converge in "
	  << solver_control.last_step()
	  << " iterations as per set tolerances. consider increasing maxSolverIterations or decreasing relSolverTolerance.\n";
  }
  constraintmatrix.distribute (completely_distributed_solutionInc);
  dxGhosts=completely_distributed_solutionInc;
  x+=completely_distributed_solutionInc;
  xGhosts=x;
}
#include "../../include/ellipticBVP_template_instantiations.h"
