//FE parameters
#define feOrder   1 // Basis function interpolation order (1-linear)
#define quadOrder 2 // Quadrature point order n^3 (2->8 quadrature points)
#define meshRefineFactor 3 // 2^n*2^n*2^n elements(3->8*8*8 =512 elements)

//Solution output parameters
#define writeOutput true // flag to write output vtu and pvtu files
#define outputDirectory "."
#define skipOutputSteps 0
#define output_Eqv_strain true
#define output_Eqv_stress true
#define output_Grain_ID   true

#define linearSolverType PETScWrappers::SolverCG // Type of linear solver
#define totalNumIncrements 100 // No. of increments
#define maxLinearSolverIterations 50000 // Maximum iterations for linear solver
#define relLinearSolverTolerance  1.0e-10 // Relative linear solver tolerance
#define maxNonLinearIterations 5 // Maximum no. of non-linear iterations
#define absNonLinearTolerance 1.0e-18 // Non-linear solver tolerance
#define relNonLinearTolerance 1.0e-6 // Relative non-linear solver tolerance
#define stopOnConvergenceFailure false // Flag to stop problem if convergence fails


//Elastic Parameters

#define c11  97.7e3 // C11 (MPa)
#define c12 82.7e3 // C12 (MPa)
#define c44 37.5e3  // C44 (MPa)

double elasticStiffness[6][6]={{97.7e3, 82.7e3, 82.7e3, 0, 0, 0},
				   {82.7e3, 97.7e3, 82.7e3, 0, 0, 0},
				   {82.7e3, 82.7e3, 97.7e3, 0, 0, 0},
				   {0, 0, 0, 37.5e3, 0, 0},
				   {0, 0, 0, 0, 37.5e3, 0}, 
				   {0, 0, 0, 0, 0, 37.5e3}}; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)


//Crystal Plasticity parameters

#define numSlipSystems 12 //
#define latentHardeningRatio 1.4  //q1

double initialSlipResistance[numSlipSystems]= {200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};
double initialHardeningModulus[numSlipSystems]= {1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0};
double powerLawExponent[numSlipSystems]= {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
double saturationStress[numSlipSystems]= {500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0};

// Crystal Plasticity Constitutive model parameters

#define modelStressTolerance 1.0e-6 // Stress tolerance for the yield surface (MPa)
#define modelMaxSlipSearchIterations 20 // Maximum no. of active slip search iterations
#define modelMaxSolverIterations 10 // Maximum no. of iterations to achieve non-linear convergence
#define modelMaxPlasticSlipL2Norm 0.8 // L2-Norm of plastic slip strain-used for load-step adaptivity
#define adaptiveLoadStepFactor 0.5 // Load step factor

// Read Input Microstructure

unsigned int numPts[3]={20, 20, 22}; // No. of voxels in x,y and z directions






