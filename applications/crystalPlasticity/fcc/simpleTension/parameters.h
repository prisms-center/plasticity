/*FE parameters*/
#define feOrder   1 // Basis function interpolation order (1-linear)
#define quadOrder 2 // Quadrature point order n^3 (2->8 quadrature points)

/*Mesh parameters*/
//Set the length of the domain in all three dimensions
//Each axes spans from zero to the specified length
#define spanX 1.0
#define spanY 1.0
#define spanZ 1.0
// The number of elements in each direction is 2^(refineFactor) * subdivisions
// For optimal performance, use meshRefineFactor primarily to determine the element size
#define subdivisionsX 1
#define subdivisionsY 1
#define subdivisionsZ 1
#define meshRefineFactor 3 // 2^n*2^n*2^n elements(3->8*8*8 =512 elements)
#define writeMeshToEPS  true //Only written for serial runs and if number of elements < 10000

/*Solution output parameters*/
#define writeOutput true // flag to write output vtu and pvtu files
#define outputDirectory "."
#define skipOutputSteps 0
#define output_Eqv_strain true
#define output_Eqv_stress true
#define output_Grain_ID   true

/*Solver parameters*/
#define linearSolverType PETScWrappers::SolverCG // Type of linear solver
#define totalNumIncrements 100 // No. of increments
#define maxLinearSolverIterations 50000 // Maximum iterations for linear solver
#define relLinearSolverTolerance  1.0e-10 // Relative linear solver tolerance
#define maxNonLinearIterations 4 // Maximum no. of non-linear iterations
#define absNonLinearTolerance 1.0e-18 // Non-linear solver tolerance
#define relNonLinearTolerance 1.0e-3 // Relative non-linear solver tolerance
#define stopOnConvergenceFailure false // Flag to stop problem if convergence fails

/*Adaptive time-stepping parameters*/
#define enableAdaptiveTimeStepping false //Flag to enable adaptive time steps
#define adaptiveLoadStepFactor 0.5 // Load step factor
#define adaptiveLoadIncreaseFactor 1.25
#define succesiveIncForIncreasingTimeStep 10

//Elastic Parameters
double elasticStiffness[6][6]={{170.0e3, 124.0e3, 124.0e3, 0, 0, 0},
				   {124.0e3, 170.0e3, 124.0e3, 0, 0, 0},
				   {124.0e3, 124.0e3, 170.0e3, 0, 0, 0},
				   {0, 0, 0, 75.0e3, 0, 0},
				   {0, 0, 0, 0, 75.0e3, 0},
				   {0, 0, 0, 0, 0, 75.0e3}}; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)



//Crystal Plasticity parameters
#define numSlipSystems 12 // generally 12 for FCC
#define latentHardeningRatio 1.4  //q1

double initialSlipResistance[numSlipSystems]= {16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0, 16.0}; //CRSS of the slip sytems
double initialHardeningModulus[numSlipSystems]= {180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0, 180.0}; //Hardening moduli of slip systems
double powerLawExponent[numSlipSystems]= {2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25, 2.25}; // Power law coefficient
double saturationStress[numSlipSystems]= {148.0, 148.0, 148.0, 148.0, 148.0, 148.0, 148.0, 148.0, 148.0, 148.0, 148.0, 148.0}; // Saturation stress


//Backstress factor

#define backstressFactor 0.0 //(Ratio between backstress and CRSS during load reversal)


//Slip systems files
#define slipDirectionsFile "slipDirections.txt" // Slip Directions File
#define slipNormalsFile "slipNormals.txt" // Slip Normals File


// Crystal Plasticity Constitutive model tolerances (for advanced users)
#define modelStressTolerance 1.0e-6 // Stress tolerance for the yield surface (MPa)
#define modelMaxSlipSearchIterations 1 // Maximum no. of active slip search iterations
#define modelMaxSolverIterations 4 // Maximum no. of iterations to achieve non-linear convergence
#define modelMaxPlasticSlipL2Norm 0.8 // L2-Norm of plastic slip strain-used for load-step adaptivity


//Read Input Microstructure
unsigned int numPts[3]={20, 20, 22}; // No. of voxels in x,y and z directions
#define grainIDFile "grainID.txt" // Grain ID File
#define headerLinesGrainIDFile 5 // No. of header Lines
#define grainOrientationsFile "orientations.txt" // Slip Normals File
