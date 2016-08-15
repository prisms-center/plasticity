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
#define maxNonLinearIterations 5 // Maximum no. of non-linear iterations
#define absNonLinearTolerance 1.0e-18 // Non-linear solver tolerance
#define relNonLinearTolerance 1.0e-6 // Relative non-linear solver tolerance
#define stopOnConvergenceFailure false // Flag to stop problem if convergence fails

//Elastic Parameters
double elasticStiffness[6][6]={{59.3e3, 25.7e3, 21.4e3, 0, 0, 0},
				   {25.7e3, 59.3e3, 21.4e3, 0, 0, 0},
				   {21.4e3, 21.4e3, 61.5e3, 0, 0, 0},
				   {0, 0, 0, 16.4e3, 0, 0},
				   {0, 0, 0, 0, 16.4e3, 0}, 
				   {0, 0, 0, 0, 0, 16.8e3}}; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)


//Crystal Plasticity 
//slip parameters
#define numSlipSystems 18 // Total No. of slip systems (slip)
#define latentHardeningRatio 1.4  //q1

double initialSlipResistance[numSlipSystems]= {25.0, 25.0, 25.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0, 68.0};
double initialHardeningModulus[numSlipSystems]= {100.0, 100.0, 100.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0, 130.0};
double powerLawExponent[numSlipSystems]= {1.1, 1.1, 1.1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8};
double saturationStress[numSlipSystems]= {70.0, 70.0, 70.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0, 210.0};

//Twin parameters

#define numTwinSystems 6 // No. of twin systems

double initialSlipResistanceTwin[numTwinSystems]= {40.0, 40.0, 40.0, 40.0, 40.0, 40.0};
double initialHardeningModulusTwin[numTwinSystems]= {50.0, 50.0, 50.0, 50.0, 50.0, 50.0};
double powerLawExponentTwin[numTwinSystems]= {1.1, 1.1, 1.1, 1.1, 1.1, 1.1};
double saturationStressTwin[numTwinSystems]= {50.0, 50.0, 50.0, 50.0, 50.0, 50.0};

#define twinThresholdFraction 0.25 //
#define twinSaturationFactor 0.25 //
#define twinShear 0.129 // characteristic twin shear


// Crystal Plasticity Constitutive model parameters

#define modelStressTolerance 1.0e-6 // Stress tolerance for the yield surface (MPa)
#define modelMaxSlipSearchIterations 20 // Maximum no. of active slip search iterations
#define modelMaxSolverIterations 10 // Maximum no. of iterations to achieve non-linear convergence
#define modelMaxPlasticSlipL2Norm 0.8 // L2-Norm of plastic slip strain-used for load-step adaptivity
#define adaptiveLoadStepFactor 0.5 // Load step factor



// Read Input Microstructure

unsigned int numPts[3]={20, 20, 22}; // No. of voxels in x,y and z directions






