// Class to load in the user input from parameters.h and the variable definition
// part of equations.h

#ifndef INCLUDE_USERINPUTPARAMETERS_H_
#define INCLUDE_USERINPUTPARAMETERS_H_

#include "dealIIheaders.h"
#include "model_variables.h"
#include "varBCs.h"
#include "inputFileReader.h"
#include "variableAttributeLoader.h"
#include "nucleationParameters.h"

template <int dim>
class userInputParameters
{

public:
  /*FE parameters*/
  unsigned int feOrder; // Basis function interpolation order (1-linear)
  unsigned int quadOrder;// Quadrature point order n^3 (2->8 quadrature points)

  /*Mesh parameters*/
  //Set the length of the domain in all three dimensions
  //Each axes spans from zero to the specified length
  std::vector<double> span;

  // The number of elements in each direction is 2^(refineFactor) * subdivisions
  // For optimal performance, use meshRefineFactor primarily to determine the element size
  std::vector<double> subdivisions;

  unsigned int meshRefineFactor; // 2^n*2^n*2^n elements(3->8*8*8 =512 elements)
  bool writeMeshToEPS; //Only written for serial runs and if number of elements < 10000

  /*Solution output parameters*/
  bool writeOutput; // flag to write output vtu and pvtu files
  std::string outputDirectory;
  unsigned int skipOutputSteps;
  bool output_Eqv_strain;
  bool output_Eqv_stress;
  bool output_Grain_ID;

  /*Solver parameters*/
  std::string linearSolverType; // Type of linear solver
  unsigned int totalNumIncrements; // No. of increments
  unsigned int maxLinearSolverIterations; // Maximum iterations for linear solver
  unsigned int maxNonLinearIterations; // Maximum no. of non-linear iterations
  double relLinearSolverTolerance; // Relative linear solver tolerance
  double absNonLinearTolerance; // Non-linear solver tolerance
  double relNonLinearTolerance; // Relative non-linear solver tolerance
  bool stopOnConvergenceFailure; // Flag to stop problem if convergence fails

  /*Adaptive time-stepping parameters*/
  bool enableAdaptiveTimeStepping; //Flag to enable adaptive time steps
  double adaptiveLoadStepFactor; // Load step factor
  double adaptiveLoadIncreaseFactor;
  double succesiveIncForIncreasingTimeStep;

  //Elastic Parameters
  FullMatrix<double> elasticStiffness; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems; // generally 12 for FCC
  double latentHardeningRatio; //q1

  std::vector<double> initialSlipResistance; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent; // Power law coefficient
  std::vector<double> saturationStress; // Saturation stress

  //Backstress factor
  double backstressFactor; //(Ratio between backstress and CRSS during load reversal)

  //Slip systems files
  std::string slipDirectionsFile; // Slip Directions File
  std::string slipNormalsFile; // Slip Normals File

  // Crystal Plasticity Constitutive model tolerances (for advanced users)
  double modelStressTolerance; // Stress tolerance for the yield surface (MPa)
  unsigned int modelMaxSlipSearchIterations; // Maximum no. of active slip search iterations
  unsigned int modelMaxSolverIterations; // Maximum no. of iterations to achieve non-linear convergence
  double modelMaxPlasticSlipL2Norm; // L2-Norm of plastic slip strain-used for load-step adaptivity

  //Read Input Microstructure
  std::vector<unsigned int> numPts; // No. of voxels in x,y and z directions
  std::string grainIDFile; // Grain ID File
  unsigned int headerLinesGrainIDFile; // No. of header Lines
  std::string grainOrientationsFile; // Slip Normals File

private:
};

#endif /* INCLUDE_USERINPUTPARAMETERS_H_ */
