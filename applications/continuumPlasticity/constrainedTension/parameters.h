/**
 *Lagrangian basis function order
 */
#define feOrder   1
/**
 *Gaussian quadrature point order
 */
#define quadOrder 2 
/**
 *Exponential value used to define mesh refinement
 */
#define meshRefineFactor 2
/**
 *Total applied displacement
 */
#define totalDisplacement 0.05
/**
 *Flag to write output files
 */
#define writeOutput true
/**
 *Solve type for linear solves
 */
#define linearSolverType PETScWrappers::SolverCG
/**
 *Number of increments (i.e. loads steps, pseudo-time steps)
 */
#define totalNumIncrements 200
/**
 *Maximum iterations within the linear solve
 */
#define maxLinearSolverIterations 5000
/**
 *Relative tolerance for the linear solver
 */
#define relLinearSolverTolerance  1.0e-14
/**
 *Maximum iterations within the nonlinear solve (Newton-Raphson iterations)
 */
#define maxNonLinearIterations 30
/**
 *Absolute tolerance for nonlinear solver
 */
#define absNonLinearTolerance 1.0e-15
/**
 *Relative tolerance for the nonlinear solver
 */
#define relNonLinearTolerance 1.0e-10
/**
 *Flag to stop the problem if convergence fails
 */
#define stopOnConvergenceFailure true

/**
 *Lame' material parameter, lambda
 */
#define lame_lambda 100.6582e9
/**
 *Lame' material parameter, mu
 */
#define lame_mu 45.6473e9
/**
 *Value for yield stress (Kirchhoff)
 */
#define yield_stress 33.014025e6
/**
 *Linear isotropic strain hardening coefficient
 */
#define strain_hardening 1000
/**
 *Kinematic strain hardening coefficient
 */
#define kinematic_hardening 0.0
/**
 *Strain energy density function ("quadlog", "stvenkir", or "neohook")
 */
#define strain_energy_function "quadlog"
/**
 *Yield function (currently only "von_mises")
 */
#define yield_function "von_mises"
