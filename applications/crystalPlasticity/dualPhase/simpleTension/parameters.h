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
#define output_Phase_ID true

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


#define multiplePhase true

//Phase1


double elasticStiffness1[6][6]={{97.7e3, 82.7e3, 82.7e3, 0, 0, 0},
    {82.7e3, 97.7e3, 82.7e3, 0, 0, 0},
    {82.7e3, 82.7e3, 97.7e3, 0, 0, 0},
    {0, 0, 0, 37.5e3, 0, 0},
    {0, 0, 0, 0, 37.5e3, 0},
    {0, 0, 0, 0, 0, 37.5e3}}; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)


//Crystal Plasticity parameters

#define numSlipSystems1 12 //
#define latentHardeningRatio1 1.4  //q1

double initialSlipResistance1[numSlipSystems1]= {200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0};  //CRSS of slip sytems
double initialHardeningModulus1[numSlipSystems1]= {1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0, 1500.0}; //Hardening moduli of slip systems
double powerLawExponent1[numSlipSystems1]= {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; // Power law coefficient
double saturationStress1[numSlipSystems1]= {500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0, 500.0}; // Saturation stress


//Slip systems files
#define slipDirectionsFile1 "slipDirections1.txt" // Slip Directions File
#define slipNormalsFile1 "slipNormals1.txt" // Slip Normals File



//Phase2

//Elastic Parameters
double elasticStiffness2[6][6]={{162e3, 92.0e3, 69.0e3, 0, 0, 0},
    {92.0e3, 162e3, 69.0e3, 0, 0, 0},
    {69.0e3, 69.0e3, 180e3, 0, 0, 0},
    {0, 0, 0, 46.7e3, 0, 0},
    {0, 0, 0, 0, 46.7e3, 0},
    {0, 0, 0, 0, 0, 35e3}}; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)


//Crystal Plasticity
//slip parameters
#define numSlipSystems2 18 // Total No. of slip systems (slip)
#define latentHardeningRatio2 1.4  //q1

double initialSlipResistance2[numSlipSystems2]= {49.0, 49.0, 49.0, 37.0, 37.0, 37.0, 197.0, 197.0, 197.0, 197.0, 197.0, 197.0, 200.0, 200.0, 200.0, 200.0, 200.0, 200.0}; //CRSS of slip sytems
double initialHardeningModulus2[numSlipSystems2]= {150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0, 150.0}; //Hardening moduli of slip systems
double powerLawExponent2[numSlipSystems2]= {1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5}; // Power law coefficient
double saturationStress2[numSlipSystems2]= {200.0, 200.0, 200.0, 200.0, 200.0, 200.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0, 550.0}; // Saturation stress

//Twin parameters

#define numTwinSystems 6 // No. of twin systems

double initialSlipResistanceTwin[numTwinSystems]= {213.0, 213.0, 213.0, 213.0, 213.0, 213.0}; //CRSS of twin sytems
double initialHardeningModulusTwin[numTwinSystems]= {1000.0, 1000.0, 1000.0, 1000.0, 1000.0, 1000.0}; //Hardening moduli of twin systems
double powerLawExponentTwin[numTwinSystems]= {1.5, 1.5, 1.5, 1.5, 1.5, 1.5};// Power law coefficient
double saturationStressTwin[numTwinSystems]= {700.0, 700.0, 700.0, 700.0, 700.0, 700.0}; // Saturation stress

#define twinThresholdFraction 0.25 // threshold fraction of characteristic twin shear (<1)
#define twinSaturationFactor 0.25 // twin growth saturation factor  (<(1-twinThresholdFraction))
#define twinShear 0.169 // characteristic twin shear


//Slip systems files
#define slipDirectionsFile2 "slipDirections2.txt" // Slip Directions File
#define slipNormalsFile2 "slipNormals2.txt" // Slip Normals File

//Twin systems files
#define twinDirectionsFile "twinDirections.txt" // Slip Directions File
#define twinNormalsFile "twinNormals.txt" // Slip Normals File





// Crystal Plasticity Constitutive model tolerances (for advanced users)
#define modelStressTolerance 1.0e-6 // Stress tolerance for the yield surface (MPa)
#define modelMaxSlipSearchIterations 1 // Maximum no. of active slip search iterations
#define modelMaxSolverIterations 3 // Maximum no. of iterations to achieve non-linear convergence
#define modelMaxPlasticSlipL2Norm 0.8 // L2-Norm of plastic slip strain-used for load-step adaptivity


//Read Input Microstructure
unsigned int numPts[3]={20, 20, 22}; // No. of voxels in x,y and z directions
#define grainIDFile "grainID.txt" // Grain ID File
#define headerLinesGrainIDFile 5 // No. of header Lines
#define grainOrientationsFile "orientations1.txt" // Slip Normals File






