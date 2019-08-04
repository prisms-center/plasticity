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
  double externalMeshParameter;
  double delT; // Time increment
  double totalTime; // Total simulation time

  std::string BCfilename; // Boundary conditions file
  unsigned int BCheaderLines; // No. of header Lines in BC file
  unsigned int NumberofBCs; // No. of boundary conditions
  bool useVelocityGrad; // Specify whether to use velocity gradient for BC
  std::vector<std::vector<double>> targetVelGrad; // 	Velocity gradient required for L based boundary conditions
  bool enableCyclicLoading;
  unsigned int cyclicLoadingFace;
  unsigned int cyclicLoadingDOF;
  double quarterCycleTime;

  /*Solution output parameters*/
  bool writeOutput, writeQuadratureOutput; // flag to write output vtu and pvtu files
  std::string outputDirectory;
  unsigned int skipOutputSteps, skipQuadratureOutputSteps;
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
  unsigned int additionalVoxelInfo; // Additional Voxel info in addition to three orientation components
  bool enableMultiphase; //Flag to indicate if Multiphase is enabled
  unsigned int numberofPhases; // Number of phases

  bool enableUserMaterialModel; //Flag to indicate if User material Model is enabled

    bool enableUserMaterialModel1; //Flag to indicate if User material Model is enabled for Phase 1
  unsigned int numberofUserMatConstants1; // Number of User Material Constants
  unsigned int numberofUserMatStateVar1; // Number of User Material State Variables
  std::vector<double> UserMatConstants1; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  std::vector<double> UserMatStateVar1; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)


  //Elastic Parameters
  std::vector<std::vector<double>> elasticStiffness1; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems1;
  double latentHardeningRatio1; //q1
  std::vector<double> initialSlipResistance1; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus1; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent1; // Power law coefficient
  std::vector<double> saturationStress1; // Saturation stress
  std::vector<double> C_1_slip1; // Kinematic hardening constants
  std::vector<double> C_2_slip1; // Kinematic hardening constants
  std::string slipDirectionsFile1; // Slip Directions File
  std::string slipNormalsFile1; // Slip Normals File
  bool enableKinematicHardening1;
  bool enableTwinning1;
  unsigned int numTwinSystems1;
  std::vector<double> initialSlipResistanceTwin1; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulusTwin1; //Hardening moduli of slip systems
  std::vector<double> powerLawExponentTwin1; // Power law coefficient
  std::vector<double> saturationStressTwin1; // Saturation stress
  std::vector<double> C_1_twin1; // Kinematic hardening constants
  std::vector<double> C_2_twin1; // Kinematic hardening constants
  std::string twinDirectionsFile1; // Twin Directions File
  std::string twinNormalsFile1; // Twin Normals File

  double twinThresholdFraction1; // threshold fraction of characteristic twin shear (<1)
  double twinSaturationFactor1; // twin growth saturation factor  (<(1-twinThresholdFraction))
  double twinShear1; // characteristic twin shear

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

  bool enableUserMaterialModel2; //Flag to indicate if User material Model is enabled
  unsigned int numberofUserMatConstants2; // Number of User Material Constants
  unsigned int numberofUserMatStateVar2; // Number of User Material State Variables
  std::vector<double> UserMatConstants2; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  std::vector<double> UserMatStateVar2; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  //Elastic Parameters
  std::vector<std::vector<double>> elasticStiffness2; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems2;
  double latentHardeningRatio2; //q1
  std::vector<double> initialSlipResistance2; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus2; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent2; // Power law coefficient
  std::vector<double> saturationStress2; // Saturation stress
  std::vector<double> C_1_slip2; // Kinematic hardening constants
  std::vector<double> C_2_slip2; // Kinematic hardening constants
  std::string slipDirectionsFile2; // Slip Directions File
  std::string slipNormalsFile2; // Slip Normals File
  bool enableKinematicHardening2;
  bool enableTwinning2;
  unsigned int numTwinSystems2;
  std::vector<double> initialSlipResistanceTwin2; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulusTwin2; //Hardening moduli of slip systems
  std::vector<double> powerLawExponentTwin2; // Power law coefficient
  std::vector<double> saturationStressTwin2; // Saturation stress
  std::vector<double> C_1_twin2; // Kinematic hardening constants
  std::vector<double> C_2_twin2; // Kinematic hardening constants
  std::string twinDirectionsFile2; // Twin Directions File
  std::string twinNormalsFile2; // Twin Normals File

  double twinThresholdFraction2; // threshold fraction of characteristic twin shear (<1)
  double twinSaturationFactor2; // twin growth saturation factor  (<(1-twinThresholdFraction))
  double twinShear2; // characteristic twin shear

  bool enableUserMaterialModel3; //Flag to indicate if User material Model is enabled
  unsigned int numberofUserMatConstants3; // Number of User Material Constants
  unsigned int numberofUserMatStateVar3; // Number of User Material State Variables
  std::vector<double> UserMatConstants3; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  std::vector<double> UserMatStateVar3; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  //Elastic Parameters
  std::vector<std::vector<double>> elasticStiffness3; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems3;
  double latentHardeningRatio3; //q1
  std::vector<double> initialSlipResistance3; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus3; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent3; // Power law coefficient
  std::vector<double> saturationStress3; // Saturation stress
  std::vector<double> C_1_slip3; // Kinematic hardening constants
  std::vector<double> C_2_slip3; // Kinematic hardening constants
  std::string slipDirectionsFile3; // Slip Directions File
  std::string slipNormalsFile3; // Slip Normals File
  bool enableKinematicHardening3;
  bool enableTwinning3;
  unsigned int numTwinSystems3;
  std::vector<double> initialSlipResistanceTwin3; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulusTwin3; //Hardening moduli of slip systems
  std::vector<double> powerLawExponentTwin3; // Power law coefficient
  std::vector<double> saturationStressTwin3; // Saturation stress
  std::vector<double> C_1_twin3; // Kinematic hardening constants
  std::vector<double> C_2_twin3; // Kinematic hardening constants
  std::string twinDirectionsFile3; // Twin Directions File
  std::string twinNormalsFile3; // Twin Normals File

  double twinThresholdFraction3; // threshold fraction of characteristic twin shear (<1)
  double twinSaturationFactor3; // twin growth saturation factor  (<(1-twinThresholdFraction))
  double twinShear3; // characteristic twin shear

  bool enableUserMaterialModel4; //Flag to indicate if User material Model is enabled
  unsigned int numberofUserMatConstants4; // Number of User Material Constants
  unsigned int numberofUserMatStateVar4; // Number of User Material State Variables
  std::vector<double> UserMatConstants4; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  std::vector<double> UserMatStateVar4; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  //Elastic Parameters
  std::vector<std::vector<double>> elasticStiffness4; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems4;
  double latentHardeningRatio4; //q1
  std::vector<double> initialSlipResistance4; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus4; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent4; // Power law coefficient
  std::vector<double> saturationStress4; // Saturation stress
  std::vector<double> C_1_slip4; // Kinematic hardening constants
  std::vector<double> C_2_slip4; // Kinematic hardening constants
  std::string slipDirectionsFile4; // Slip Directions File
  std::string slipNormalsFile4; // Slip Normals File
  bool enableKinematicHardening4;
  bool enableTwinning4;
  unsigned int numTwinSystems4;
  std::vector<double> initialSlipResistanceTwin4; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulusTwin4; //Hardening moduli of slip systems
  std::vector<double> powerLawExponentTwin4; // Power law coefficient
  std::vector<double> saturationStressTwin4; // Saturation stress
  std::vector<double> C_1_twin4; // Kinematic hardening constants
  std::vector<double> C_2_twin4; // Kinematic hardening constants
  std::string twinDirectionsFile4; // Twin Directions File
  std::string twinNormalsFile4; // Twin Normals File

  double twinThresholdFraction4; // threshold fraction of characteristic twin shear (<1)
  double twinSaturationFactor4; // twin growth saturation factor  (<(1-twinThresholdFraction))
  double twinShear4; // characteristic twin shear

private:
  dealii::ConditionalOStream  pcout;
};

#endif /* INCLUDE_USERINPUTPARAMETERS_H_ */
