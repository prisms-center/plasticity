#define feOrder   1 // Basis function interpolation order (1-linear)
#define quadOrder 2 // Quadrature point order n^3 (2->8 quadrature points)
#define meshRefineFactor 3 // 2^n*2^n*2^n elements(3->8*8*8 =512 elements)
#define writeOutput true // flag to write output vtu and pvtu files)
#define linearSolverType PETScWrappers::SolverCG // Type of linear solver
#define totalNumIncrements 100 // No. of increments
#define maxLinearSolverIterations 50000 // Maximum iterations for linear solver
#define relLinearSolverTolerance  1.0e-10 // Relative linear solver tolerance
#define maxNonLinearIterations 5 // Maximum no. of non-linear iterations
#define absNonLinearTolerance 1.0e-18 // Non-linear solver tolerance
#define relNonLinearTolerance 1.0e-6 // Relative non-linear solver tolerance
#define stopOnConvergenceFailure false // Flag to stop problem if convergence fails

//Elastic Parameters

#define c11 59.3e3 // C11 (MPa)
#define c12 25.7e3 // C12 (MPa)
#define c13 21.4e3  // C44 (MPa)
#define c33 61.5e3 // C33 (MPa)
#define c44 16.4e3 // C44 (MPa)



//Crystal Plasticity parameters

#define numSlipSystems 24 // Total No. of slip systems (slip+twin)
#define numTwinSystems 6 // No. of twin systems
#define latentHardeningRatio 1.4  //q1
#define powerLawExponent1 1.1  //a_basal
#define powerLawExponent2 0.8  //a_prismatic
#define powerLawExponent3 0.8  //a_pyramidal<a>
#define powerLawExponent4 0.8  //a_pyramidal<c+a>
#define powerLawExponent5 1.1  //a_twin<c+a>
#define initialSlipResistance1 25.0 // CRSS s0_basal(MPa)
#define initialSlipResistance2 68.0 // CRSS s0_prismatic(MPa)
#define initialSlipResistance3 68.0 // CRSS s0_pyramidal<a>(MPa)
#define initialSlipResistance4 68.0 // CRSS s0_pyramidal<c+a>(MPa)
#define initialSlipResistance5 40.0 // CRSS s0_twin<c+a>(MPa)
#define saturationStress1 70.0 //s_s_basal(MPa)
#define saturationStress2 210.0 //s_s_prismatic(MPa)
#define saturationStress3 210.0 //s_s_pyramidal<a>(MPa)
#define saturationStress4 210.0 //s_s_pyramidal<c+a>(MPa)
#define saturationStress5 50.0 //s_s_twin<c+a>(MPa)
#define initialHardeningModulus1 100.0 //h0_basal(MPa)
#define initialHardeningModulus2 130.0 //h0_prismatic(MPa)
#define initialHardeningModulus3 130.0 //h0_pyramidal<a>(MPa)
#define initialHardeningModulus4 130.0 //h0_pyramidal<c+a>(MPa)
#define initialHardeningModulus5 50.0 //h0_twin<c+a>(MPa)

// Read Input Microstructure

unsigned int numPts[3]={20, 20, 22}; // No. of voxels in x,y and z directions






