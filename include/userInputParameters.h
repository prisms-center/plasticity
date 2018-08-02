// Class to load in the user input from parameters.h and the variable definition
// part of equations.h

#ifndef INCLUDE_USERINPUTPARAMETERS_H_
#define INCLUDE_USERINPUTPARAMETERS_H_

#include "dealIIheaders.h"

class userInputParameters
{

public:

  userInputParameters(std::string input_file_name, dealii::ParameterHandler & parameter_handler);

  void declare_parameters(dealii::ParameterHandler & parameter_handler);

  unsigned int dim;

  /*FE parameters*/
  unsigned int feOrder; // Basis function interpolation order (1-linear)
  unsigned int quadOrder;// Quadrature point order n^3 (2->8 quadrature points)

  /*Mesh parameters*/
  //Set the length of the domain in all three dimensions
  //Each axes spans from zero to the specified length
  std::vector<double> span;

  // The number of elements in each direction is 2^(refineFactor) * subdivisions
  // For optimal performance, use meshRefineFactor primarily to determine the element size
  std::vector<unsigned int> subdivisions;

  unsigned int meshRefineFactor; // 2^n*2^n*2^n elements(3->8*8*8 =512 elements)

  bool writeMeshToEPS; //Only written for serial runs and if number of elements < 10000

  bool readExternalMesh;
  std::string externalMeshFileName;

  double delT; // Time increment
  double totalTime; // Total simulation time

  std::string BCfilename; // Boundary conditions file
  unsigned int BCheaderLines; // No. of header Lines in BC file
  unsigned int NumberofBCs; // No. of boundary conditions

  /*Solution output parameters*/
  bool writeOutput; // flag to write output vtu and pvtu files
  std::string outputDirectory;
  unsigned int skipOutputSteps;
  bool output_Eqv_strain;
  bool output_Eqv_stress;
  bool output_Grain_ID;
  bool output_Twin;

  /*Solver parameters*/
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

  std::string crystalStructure;

  //Elastic Parameters
  std::vector<std::vector<double>> elasticStiffness; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems;
  double latentHardeningRatio; //q1
  std::vector<double> initialSlipResistance; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent; // Power law coefficient
  std::vector<double> saturationStress; // Saturation stress
  std::string slipDirectionsFile; // Slip Directions File
  std::string slipNormalsFile; // Slip Normals File

  bool enableTwinning;
  unsigned int numTwinSystems;
  std::vector<double> initialSlipResistanceTwin; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulusTwin; //Hardening moduli of slip systems
  std::vector<double> powerLawExponentTwin; // Power law coefficient
  std::vector<double> saturationStressTwin; // Saturation stress
  std::string twinDirectionsFile; // Twin Directions File
  std::string twinNormalsFile; // Twin Normals File

  double backstressFactor; //Ratio between backstress and CRSS during load reversal
  double twinThresholdFraction; // threshold fraction of characteristic twin shear (<1)
  double twinSaturationFactor; // twin growth saturation factor  (<(1-twinThresholdFraction))
  double twinShear; // characteristic twin shear

  // Crystal Plasticity Constitutive model tolerances (for advanced users)
  double modelStressTolerance; // Stress tolerance for the yield surface (MPa)
  unsigned int modelMaxSlipSearchIterations; // Maximum no. of active slip search iterations
  unsigned int modelMaxSolverIterations; // Maximum no. of iterations to achieve non-linear convergence
  double modelMaxPlasticSlipL2Norm; // L2-Norm of plastic slip strain-used for load-step adaptivity

  //Read Input Microstructure
  std::string grainIDFile; // Grain ID File
  std::vector<unsigned int> numPts; // No. of voxels in x,y and z directions

  std::string grainOrientationsFile; // Grain orientations file
  unsigned int headerLinesGrainIDFile; // No. of header Lines in grain orientations file

private:
  dealii::ConditionalOStream  pcout;
};

#endif /* INCLUDE_USERINPUTPARAMETERS_H_ */
