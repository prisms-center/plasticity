/**
 *Lagrangian basis function order
 */
const int feOrder=1;
/**
 *Gaussian quadrature point order
 */
const int quadOrder = 2;
/**
 *Set the length of the domain in the x-direction, from zero to the specified length
 */
const float spanX =5.0;
/**
 *Set the length of the domain in the y-direction, from zero to the specified length
 */
const float spanY =1.0;
/**
 *Set the length of the domain in the z-direction, from zero to the specified length
 */
const float spanZ =1.0;
/**
 *The number of elements in the x-direction
 *The number of elements in each direction is 2^(refineFactor) * subdivisions
 *For optimal performance, use meshRefineFactor primarily to determine the element size
 */
const int subdivisionsX= 5;
/**
 *The number of elements in the y-direction
 *The number of elements in each direction is 2^(refineFactor) * subdivisions
 *For optimal performance, use meshRefineFactor primarily to determine the element size
 */
const int subdivisionsY= 1;
/**
 *The number of elements in the z-direction
 *The number of elements in each direction is 2^(refineFactor) * subdivisions
 *For optimal performance, use meshRefineFactor primarily to determine the element size
 */
const int subdivisionsZ= 1;
/**
 *Exponential value used to define mesh refinement (2^n*2^n*2^n elements(3->8*8*8 =512 elements) )
 */
const int meshRefineFactor =2;
/**
 *Write .eps file of the mesh (Only written for serial runs and if number of elements < 10000)
 */
const bool writeMeshToEPS = true;
/**
 *Flag to write output vtu and pvtu files
 */
const bool writeOutput = true;
/**
 *Define directory to ouput results
 */
const char outputDirectory[2] = ".";
/**
 *Specify how frequently to write output files (i.e. output every n steps; using 0 or 1 will output every step)
 */
const int skipOutputSteps= 0;
/**
 *Flag to output the equivalent plastic strain field
 */
const bool output_alpha =true;
/**
 *Flag to output the von Mises (Kirchhoff) stress field
 */
const bool output_tau_vm= true;
/**
 *Total applied displacement
 */
const float totalDisplacement =0.05;
/**
*Solve type for linear solves
 */
#define linearSolverType PETScWrappers::SolverCG
/**
 *Number of increments (i.e. loads steps, pseudo-time steps)
 */
const int totalNumIncrements =100;
/**
 *Maximum iterations within the linear solve
 */
const int maxLinearSolverIterations =8000;
/**
 *Relative tolerance for the linear solver
 */
const float relLinearSolverTolerance = 1.0e-14;
/**
 *Maximum iterations within the nonlinear solve (Newton-Raphson iterations)
 */
const int maxNonLinearIterations= 4;
/**
 *Absolute tolerance for nonlinear solver
 */
const float absNonLinearTolerance =1.0e-15;
/**
 *Relative tolerance for the nonlinear solver
 */
const float relNonLinearTolerance =1.0e-10;
/**
 *Flag to stop the problem if convergence fails
 */
const bool stopOnConvergenceFailure = false;

/**
 *Lame' material parameter, lambda
 */
const float lame_lambda =100.6582e9;
/**
 *Lame' material parameter, mu
 */
const float lame_mu =45.6473e9;
/**
 *Value for yield stress (Kirchhoff)
 */
const float yield_stress =33.014025e3;
/**
 *Linear isotropic strain hardening coefficient
 */
const float strain_hardening =2.0259e9;
/**
 *Kinematic strain hardening coefficient
 */
const float  kinematic_hardening =1e6;
/**
 *Strain energy density function ("quadlog", "stvenkir", or "neohook")
 */
#define strain_energy_function "quadlog"
/**
 *Yield function (currently only "von_mises")
 */
#define yield_function "von_mises"
/**
 *Isotropic hardening function
 */
#define iso_hardening_function "linear_hardening"
