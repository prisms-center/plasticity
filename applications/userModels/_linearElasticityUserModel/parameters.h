/**
 *Lagrangian basis function order
 */
#define feOrder   1
/**
 *Gaussian quadrature point order
 */
#define quadOrder 2 
/**
 *Set the length of the domain in the x-direction, from zero to the specified length
 */
#define spanX 5.0
/**
 *Set the length of the domain in the y-direction, from zero to the specified length
 */
#define spanY 1.0
/**
 *Set the length of the domain in the z-direction, from zero to the specified length
 */
#define spanZ 1.0
/**
 *The number of elements in the x-direction
 *The number of elements in each direction is 2^(refineFactor) * subdivisions
 *For optimal performance, use meshRefineFactor primarily to determine the element size
 */
#define subdivisionsX 5
/**
 *The number of elements in the y-direction
 *The number of elements in each direction is 2^(refineFactor) * subdivisions
 *For optimal performance, use meshRefineFactor primarily to determine the element size
 */
#define subdivisionsY 1
/**
 *The number of elements in the z-direction
 *The number of elements in each direction is 2^(refineFactor) * subdivisions
 *For optimal performance, use meshRefineFactor primarily to determine the element size
 */
#define subdivisionsZ 1
/**
 *Exponential value used to define mesh refinement (2^n*2^n*2^n elements(3->8*8*8 =512 elements) )
 */
#define meshRefineFactor 2
/**
 *Write .eps file of the mesh (Only written for serial runs and if number of elements < 10000)
 */
#define writeMeshToEPS  true
/**
 *Flag to write output vtu and pvtu files
 */
#define writeOutput true
/**
 *Define directory to ouput results
 */
#define outputDirectory "."
/**
 *Specify how frequently to write output files (i.e. output every n steps; using 0 or 1 will output every step)
 */
#define skipOutputSteps 0
/**
 *Flag to output the equivalent plastic strain field
 */
#define output_alpha true
/**
 *Flag to output the von Mises (Kirchhoff) stress field
 */
#define output_tau_vm true
/**
 *Total applied displacement
 */
#define totalDisplacement 0.5
/**
*Solve type for linear solves
 */
#define linearSolverType PETScWrappers::SolverCG
/**
 *Number of increments (i.e. loads steps, pseudo-time steps)
 */
#define totalNumIncrements 10
/**
 *Maximum iterations within the linear solve
 */
#define maxLinearSolverIterations 8000
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
 *User model parameters. '#define enableUserModel' is required to 
 *enable the user model functionality
 */
#define enableUs1erModel
/**
 * Material properties (Lame' parameters)
 */
#define lame_lambda 100.6582e9
#define lame_mu 45.6473e9
