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

  //Using Taylor Model//
  bool flagTaylorModel;
  unsigned int numberTaylorSubsteps;
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
  double criticalDeltaFCriteria;  //Critical DeltaF Criteria
  double totalTime; // Total simulation time

  bool enableSimpleBCs; // Specify whether to use Simple (Basic) BCs
  std::string BCfilename; // Boundary conditions file
  unsigned int BCheaderLines; // No. of header Lines in BC file
  unsigned int NumberofBCs; // No. of boundary conditions
  bool useVelocityGrad; // Specify whether to use velocity gradient for BC
  std::vector<std::vector<double>> targetVelGrad; // 	Velocity gradient required for L based boundary conditions

  bool enableTabularBCs; // Specify whether to use Tabular BCs
  std::string Tabular_BCfilename; // Tabular BCs file
  unsigned int tabularBCs_InputStepNumber; //Number of Input data for Tabular BCs
  unsigned int tabularNumberofBCs; //Number of Tabular BCs
  std::vector<double> tabularTimeInput; //Table for Time intervals of Tabular BCs
  std::vector<double> tabularDispInput; //Table for Displacements of Tabular BCs

  bool enableNeumannBCs; // Specify whether to use Tabular Neumann BCs
  std::string Tabular_NeumannBCfilename; // Tabular Neumann BCs file
  unsigned int tabularNeumannBCs_InputStepNumber; //Number of Input data for Tabular Neumann BCs
  unsigned int neumannBCsNumber; //Number of Tabular Neumann BCs
  std::vector<double> tabularTimeNeumannBCs; //Table for Time intervals of Tabular Neumann BCs
  std::vector<int> neumannBCsBoundaryID,dofNeumannBCs; //Vector including the BoundaryID and dof that Neumann BCs are applied

  bool enableTorsionBCs; // Specify whether to use Tabular Torsion BCs
  unsigned int torsionAxis; //Axis of Torsion (x=0; y=1; z=2)
  std::vector<double> tabularTimeInputTorsion; //Table for Time intervals of Tabular Torsion BCs
  std::vector<double> tabularTorsionBCsInput; //Table for Anglular velocity (radians) of Tabular Torsion BCs
  std::vector<double> centerTorsion; //Center point for Torsion (x,y) for torsion axis=z or (y,z) for torsion axis=x or (z,x) for torsion axis=y

  bool enableDICpipeline; // Specify whether to use DIC pipeline
  unsigned int DIC_InputStepNumber, X_dic, Y_dic, Z_dic; //Number of Input data for DIC experiment
  std::vector<double> timeInputDIC; //Table for Time intervals of DIC experiment input
  std::string DIC_BCfilename1,DIC_BCfilename2,DIC_BCfilename3,DIC_BCfilename4; // DIC Boundary conditions file


  bool enableCyclicLoading;
  unsigned int cyclicLoadingFace;
  unsigned int cyclicLoadingDOF;
  double quarterCycleTime;

  bool enableNodalDisplacementBCs;
  unsigned int numberOfNodalBCs;
  double nodalDisplacementBCsTolerance;
  std::string nodalDisplacement_BCfilename;

  bool enablePeriodicBCs,enableTabularPeriodicBCs; // Specify whether to use Periodic BCs or Tabular Periodic BCs
  double periodicTabularTime;
  std::string Periodic_BCfilename; // PEriodic BCs Constraints file
  unsigned int numberVerticesConstraint; //Number of Vertices Constraints for Periodic BCs
  unsigned int numberEdgesConstraint; //Number of Edges Constraints for Periodic BCs
  unsigned int numberFacesConstraint; //Number of Faces Constraints for Periodic BCs
  std::vector<std::vector<int>> periodicBCsInput,tabularPeriodicCoef; // 	Periodic BCs Input
  std::vector<std::vector<double>> periodicBCsInput2,tabularPeriodicTimeInput; // 	Periodic BCs Input

  //INDENTATION
  bool enableIndentationBCs; // Specify whether to use Indentation BCs
  unsigned int indentationKeyFrames; // How many key frames in indentation
  unsigned int freezeActiveSetSteps;
  bool freezeActiveSetNewIteration;
  std::vector<int> indentationKeyFrameTimeIncs; // the increment number of each key frame
  std::string Indentation_BCfilename; // Indentation BCs Constraints file
  std::vector<double> refinementCenter; // the center point of the refinement zone
  std::vector<double> refinementZoneSize; // Radius of refinement zone from center
  std::vector<int> refinementFactor; // Number of cell refinements for plastic zone
  double activeSetCriterionCoefficient; // The number to multiply the stiffness by for the active set criterion
  double activeSetLambdaTolerance; // The value to add to the criterion to account for errors in lambda/mass
  bool debugIndentation; // If true, debugging information will be written about the initial active set for each increment


  //Continuum
  bool continuum_Isotropic; // Specify whether to use Isotropic Continuum Elasto-Plasticity model
  double lame_lambda, lame_mu, yield_stress, strain_hardening, kinematic_hardening;
  std::string strain_energy_function, yield_function, iso_hardening_function;


  /*Solution output parameters*/
  bool writeOutput; // flag to write output vtu and pvtu files
  bool writeQuadratureOutput; // flag to write quadrature points output
  bool writeGrainAveragedOutput; // flag to write Grain Average output
  bool tabularOutput; // flag to write outputs based on a time table
  std::vector<double> tabularTimeOutput; // 	Periodic BCs Input
  std::string outputDirectory;
  unsigned int skipOutputSteps, skipQuadratureOutputSteps, skipGrainAveragedOutputSteps, numberOfGrainAverageDataOutput;
  bool output_Eqv_strain;
  bool output_Eqv_stress;
  bool output_Time;
  bool output_alpha;
  bool output_Indenter_Load;
  bool output_Grain_ID;
  bool output_Twin;

  bool flagUserDefinedAverageOutput;
  unsigned int numberUserDefinedAverageOutput;

  bool output_Var1,output_Var2,output_Var3,output_Var4,output_Var5,output_Var6,output_Var7,output_Var8,output_Var9,output_Var10;
  bool output_Var11,output_Var12,output_Var13,output_Var14,output_Var15,output_Var16,output_Var17,output_Var18,output_Var19,output_Var20;
  bool output_Var21,output_Var22,output_Var23,output_Var24;

  bool outputCellCenters_Var1,outputCellCenters_Var2,outputCellCenters_Var3,outputCellCenters_Var4,outputCellCenters_Var5,outputCellCenters_Var6,outputCellCenters_Var7,outputCellCenters_Var8,outputCellCenters_Var9,outputCellCenters_Var10;
  bool outputCellCenters_Var11,outputCellCenters_Var12,outputCellCenters_Var13,outputCellCenters_Var14,outputCellCenters_Var15,outputCellCenters_Var16,outputCellCenters_Var17,outputCellCenters_Var18,outputCellCenters_Var19,outputCellCenters_Var20;
  bool outputCellCenters_Var21,outputCellCenters_Var22,outputCellCenters_Var23,outputCellCenters_Var24,outputCellCenters_Var25,outputCellCenters_Var26,outputCellCenters_Var27,outputCellCenters_Var28;

  bool flagBufferLayer; //flag to exclude the data inside the buffer layers and only output the results from the real specimen.
  // This is specifically good for HEDM characterization when you're adding buffer layers to exclude the BCs effects, but one doesn't want to include the data inside the buffer layers.
  unsigned int dimBufferLayer; // The dimension in which the buffer layer is applied. x=0, y=1, and z=2.
  double lowerBufferLayer,upperBufferLayer; //Upper and lower bounds for Buffer layer.

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
  std::string latentHardeningRatioFileName1; //q1
  std::vector<double> initialSlipResistance1; //CRSS of the slip sytems
  std::vector<double> initialHardeningModulus1; //Hardening moduli of slip systems
  std::vector<double> powerLawExponent1; // Power law coefficient
  std::vector<double> saturationStress1; // Saturation stress
  std::vector<double> C_1_slip1; // Kinematic hardening constants
  std::vector<double> C_2_slip1; // Kinematic hardening constants
  std::string slipDirectionsFile1; // Slip Directions File
  std::string slipNormalsFile1; // Slip Normals File
  bool enableKinematicHardening1,enableAdvRateDepModel;
  bool enableTwinning1;
  bool enableAdvancedTwinModel; // Flag to indicate if Advanced Twinning Model enabled
  bool enableOneTwinSys_Reorien; // Flag to indicate if one twin system reorientation is allowed
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
  double criteriaTwinVisual; //In the case of Advanced twin model, the integration point with Twin volumes larger than this Critical Value is considered twined during visualization

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

  bool enableElementDeletion; //Flag to indicate if element deletion is enabled
  unsigned int deletionGrainID; //The grainID for which the element are deleted

  bool enableUserMaterialModel2; //Flag to indicate if User material Model is enabled
  unsigned int numberofUserMatConstants2; // Number of User Material Constants
  unsigned int numberofUserMatStateVar2; // Number of User Material State Variables
  std::vector<double> UserMatConstants2; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  std::vector<double> UserMatStateVar2; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)
  //Elastic Parameters
  std::vector<std::vector<double>> elasticStiffness2; // 	Elastic Stiffness Matrix -Voigt Notation (MPa)

  //Crystal Plasticity parameters
  unsigned int numSlipSystems2;
  std::string latentHardeningRatioFileName2; //q2
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
  std::string latentHardeningRatioFileName3; //q3
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
  std::string latentHardeningRatioFileName4; //q4
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
