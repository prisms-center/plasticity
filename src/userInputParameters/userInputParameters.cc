// Functions for userInputParameters class

#include "../../include/userInputParameters.h"

userInputParameters::userInputParameters(std::string inputfile, dealii::ParameterHandler & parameter_handler):
pcout (std::cout, dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
{

  declare_parameters(parameter_handler);

  #if (DEAL_II_VERSION_MAJOR < 9 && DEAL_II_VERSION_MINOR < 5)
  parameter_handler.read_input(inputfile);
  #else
  parameter_handler.parse_input(inputfile);
  #endif

  flagTaylorModel = parameter_handler.get_bool("Flag To Use Taylor Model");
  dim=parameter_handler.get_integer("Number of dimensions");

  feOrder=parameter_handler.get_integer("Order of finite elements");
  quadOrder=parameter_handler.get_integer("Order of quadrature");

  span.push_back(parameter_handler.get_double("Domain size X"));
  if (dim > 1){
    span.push_back(parameter_handler.get_double("Domain size Y"));
    if (dim > 2){
      span.push_back(parameter_handler.get_double("Domain size Z"));
    }
  }

  subdivisions.push_back(parameter_handler.get_integer("Subdivisions X"));
  if (dim > 1){
    subdivisions.push_back(parameter_handler.get_integer("Subdivisions Y"));
    if (dim > 2){
      subdivisions.push_back(parameter_handler.get_integer("Subdivisions Z"));
    }
  }

  meshRefineFactor = parameter_handler.get_integer("Refine factor");

  writeMeshToEPS = parameter_handler.get_bool("Write Mesh To EPS");

  //External mesh parameters
  readExternalMesh = parameter_handler.get_bool("Use external mesh");
  externalMeshFileName = parameter_handler.get("Name of file containing external mesh");
  externalMeshParameter = parameter_handler.get_double("External mesh parameter");

  //output parameters
  writeOutput = parameter_handler.get_bool("Write Output");
  outputDirectory = parameter_handler.get("Output Directory");

  tabularOutput = parameter_handler.get_bool("Tabular Output");
  tabularTimeOutput=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Time Output Table")));

  skipOutputSteps=parameter_handler.get_integer("Skip Output Steps");
  writeQuadratureOutput = parameter_handler.get_bool("Write Quadrature Output");
  writeGrainAveragedOutput = parameter_handler.get_bool("Write Grain Averaged Output");
  skipQuadratureOutputSteps=parameter_handler.get_integer("Skip Quadrature Output Steps");
  skipGrainAveragedOutputSteps=parameter_handler.get_integer("Skip Grain Averaged Output Steps");
  numberOfGrainAverageDataOutput=parameter_handler.get_integer("Number of Grain Averaged Output Variables");

  if(skipOutputSteps<=0)
  skipOutputSteps=1;

  if(skipQuadratureOutputSteps<=0)
  skipQuadratureOutputSteps=1;

  delT=parameter_handler.get_double("Time increments");
  criticalDeltaFCriteria=parameter_handler.get_double("critical DeltaF Criteria");
  numberTaylorSubsteps=parameter_handler.get_integer("Number of Taylor Substeps");
  totalTime=parameter_handler.get_double("Total time");

  enableSimpleBCs=parameter_handler.get_bool("Use Simple BCs");
  BCfilename=parameter_handler.get("Boundary condition filename");
  BCheaderLines=parameter_handler.get_integer("BC file number of header lines");
  NumberofBCs=parameter_handler.get_integer("Number of boundary conditions");

  useVelocityGrad=parameter_handler.get_bool("Use velocity gradient BC");

  targetVelGrad.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Velocity gradient row 1"))));
  targetVelGrad.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Velocity gradient row 2"))));
  targetVelGrad.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Velocity gradient row 3"))));

  enableTabularBCs=parameter_handler.get_bool("Use Tabular BCs");
  Tabular_BCfilename=parameter_handler.get("Tabular Boundary condition filename");
  tabularBCs_InputStepNumber=parameter_handler.get_integer("Number of time data for Tabular BCs");
  tabularNumberofBCs=parameter_handler.get_integer("Number of tabular boundary conditions");
  tabularTimeInput=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Time Table")));

  enableNeumannBCs=parameter_handler.get_bool("Use Neumann BCs");
  Tabular_NeumannBCfilename=parameter_handler.get("Tabular Neumann Boundary condition filename");
  tabularNeumannBCs_InputStepNumber=parameter_handler.get_integer("Number of time data for Tabular Neumann BCs");
  neumannBCsNumber=parameter_handler.get_integer("Number of tabular Neumann boundary conditions");
  neumannBCsBoundaryID=dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("Boundary IDs of Neumann BCs")));
  dofNeumannBCs=dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("dof of Neumann BCs")));
  tabularTimeNeumannBCs=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Time Neumann BCs")));


  enableTorsionBCs=parameter_handler.get_bool("Use Torsion BCs");
  torsionAxis=parameter_handler.get_integer("Torsion Axis  x 0  y 1  z 2");
  tabularTimeInputTorsion=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Time Table for Torsion")));
  tabularTorsionBCsInput=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Torsion BCs")));
  centerTorsion=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Center point for Torsion")));

  enableDICpipeline=parameter_handler.get_bool("Use DIC pipeline");
  DIC_InputStepNumber=parameter_handler.get_integer("Number of Input data for DIC experiment");
  timeInputDIC=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("DIC Time Table")));
  X_dic=parameter_handler.get_integer("Number of Points in DIC input in X direction");
  Y_dic=parameter_handler.get_integer("Number of Points in DIC input in Y direction");
  Z_dic=parameter_handler.get_integer("Number of Points in DIC input in Z direction");
  DIC_BCfilename1=parameter_handler.get("DIC Boundary condition filename 1");
  DIC_BCfilename2=parameter_handler.get("DIC Boundary condition filename 2");
  DIC_BCfilename3=parameter_handler.get("DIC Boundary condition filename 3");
  DIC_BCfilename4=parameter_handler.get("DIC Boundary condition filename 4");

  enableCyclicLoading=parameter_handler.get_bool("Enable cyclic loading");
  cyclicLoadingFace=parameter_handler.get_integer("Cyclic loading face");
  cyclicLoadingDOF=parameter_handler.get_integer("Cyclic loading direction");
  quarterCycleTime=parameter_handler.get_double("Quarter cycle time");

  enableNodalDisplacementBCs=parameter_handler.get_bool("Enable Nodal Displacement BCs");
  numberOfNodalBCs=parameter_handler.get_integer("Number of Nodal Displacement BCs");
  nodalDisplacementBCsTolerance = parameter_handler.get_double("Tolerance for Nodal Displacement BCs");
  nodalDisplacement_BCfilename=parameter_handler.get("Nodal Displacement BCs filename");

  enablePeriodicBCs=parameter_handler.get_bool("Use Periodic BCs");
  Periodic_BCfilename=parameter_handler.get("Periodic Boundary condition Constraint filename");
  numberVerticesConstraint=parameter_handler.get_integer("Number of Vertices Constraints");
  numberEdgesConstraint=parameter_handler.get_integer("Number of Edges Constraints");
  numberFacesConstraint=parameter_handler.get_integer("Number of Faces Constraints");
  periodicBCsInput.push_back(dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("Vertices Periodic BCs row 1"))));
  periodicBCsInput2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Vertices Periodic BCs row 2"))));
  enableTabularPeriodicBCs=parameter_handler.get_bool("Use Tabular Periodic BCs");
  periodicTabularTime=parameter_handler.get_double("periodic Tabular time");
  tabularPeriodicTimeInput.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Periodic Time Table"))));
  tabularPeriodicCoef.push_back(dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("Tabular Periodic Time Table Coefficient"))));
  enableIndentationBCs=parameter_handler.get_bool("Use Indentation BCs");
  indentationKeyFrames=parameter_handler.get_integer("Indentation Key Frame Number");
  freezeActiveSetSteps=parameter_handler.get_integer("Indentation Active Set Freeze Iterations");
  freezeActiveSetNewIteration=parameter_handler.get_bool("Indentation Active Set Freeze on Iteration 1");
  activeSetCriterionCoefficient=parameter_handler.get_double("Active Set Criterion Coefficient");
  activeSetLambdaTolerance=parameter_handler.get_double("Active Set Lambda Tolerance");
  indentationKeyFrameTimeIncs=dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("Indentation Key Frame Increment Numbers")));
  Indentation_BCfilename=parameter_handler.get("Indentation Boundary condition Constraint filename");
  refinementCenter=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Refinement Zone Center")));
  refinementZoneSize=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Refinement Zone Size")));
  refinementFactor=dealii::Utilities::string_to_int(dealii::Utilities::split_string_list(parameter_handler.get("Refinement Factor")));
  debugIndentation=parameter_handler.get_bool("Debug Indentation Active Set");
  continuum_Isotropic = parameter_handler.get_bool("Continuum Isotropic");
  if(continuum_Isotropic){
      pcout<<"\n Isotropic Plasticity enabled (not crystal plasticity) \n";
      lame_lambda = parameter_handler.get_double("lame_lambda");
      lame_mu = parameter_handler.get_double("lame_mu");
      yield_stress = parameter_handler.get_double("yield_stress");
      strain_hardening = parameter_handler.get_double("strain_hardening");
      kinematic_hardening = parameter_handler.get_double("kinematic_hardening");
      strain_energy_function=parameter_handler.get("strain_energy_function");
      yield_function=parameter_handler.get("yield_function");
      iso_hardening_function=parameter_handler.get("iso_hardening_function");
  }
  output_Eqv_strain = parameter_handler.get_bool("Output Equivalent strain");
  output_Eqv_stress = parameter_handler.get_bool("Output Equivalent stress");
  output_Time = parameter_handler.get_bool("Output Time");
  output_alpha = parameter_handler.get_bool("Output Equivalent plastic strain alpha");
  output_Indenter_Load = parameter_handler.get_bool("Output Indenter Load");
  output_Grain_ID = parameter_handler.get_bool("Output Grain ID");
  output_Twin = parameter_handler.get_bool("Output Twin fractions");

  flagUserDefinedAverageOutput = parameter_handler.get_bool("Output Userdefined Average Variable");
  numberUserDefinedAverageOutput = parameter_handler.get_integer("Number of Output Userdefined Average Variable");

  output_Var1 = parameter_handler.get_bool("Output Variable 1");
  output_Var2 = parameter_handler.get_bool("Output Variable 2");
  output_Var3 = parameter_handler.get_bool("Output Variable 3");
  output_Var4 = parameter_handler.get_bool("Output Variable 4");
  output_Var5 = parameter_handler.get_bool("Output Variable 5");
  output_Var6 = parameter_handler.get_bool("Output Variable 6");
  output_Var7 = parameter_handler.get_bool("Output Variable 7");
  output_Var8 = parameter_handler.get_bool("Output Variable 8");
  output_Var9 = parameter_handler.get_bool("Output Variable 9");
  output_Var10 = parameter_handler.get_bool("Output Variable 10");
  output_Var11 = parameter_handler.get_bool("Output Variable 11");
  output_Var12 = parameter_handler.get_bool("Output Variable 12");
  output_Var13 = parameter_handler.get_bool("Output Variable 13");
  output_Var14 = parameter_handler.get_bool("Output Variable 14");
  output_Var15 = parameter_handler.get_bool("Output Variable 15");
  output_Var16 = parameter_handler.get_bool("Output Variable 16");
  output_Var17 = parameter_handler.get_bool("Output Variable 17");
  output_Var18 = parameter_handler.get_bool("Output Variable 18");
  output_Var19 = parameter_handler.get_bool("Output Variable 19");
  output_Var20 = parameter_handler.get_bool("Output Variable 20");
  output_Var21 = parameter_handler.get_bool("Output Variable 21");
  output_Var22 = parameter_handler.get_bool("Output Variable 22");
  output_Var23 = parameter_handler.get_bool("Output Variable 23");
  output_Var24 = parameter_handler.get_bool("Output Variable 24");

  outputCellCenters_Var1 = parameter_handler.get_bool("Output CellCenters Variable 1");
  outputCellCenters_Var2 = parameter_handler.get_bool("Output CellCenters Variable 2");
  outputCellCenters_Var3 = parameter_handler.get_bool("Output CellCenters Variable 3");
  outputCellCenters_Var4 = parameter_handler.get_bool("Output CellCenters Variable 4");
  outputCellCenters_Var5 = parameter_handler.get_bool("Output CellCenters Variable 5");
  outputCellCenters_Var6 = parameter_handler.get_bool("Output CellCenters Variable 6");
  outputCellCenters_Var7 = parameter_handler.get_bool("Output CellCenters Variable 7");
  outputCellCenters_Var8 = parameter_handler.get_bool("Output CellCenters Variable 8");
  outputCellCenters_Var9 = parameter_handler.get_bool("Output CellCenters Variable 9");
  outputCellCenters_Var10 = parameter_handler.get_bool("Output CellCenters Variable 10");
  outputCellCenters_Var11 = parameter_handler.get_bool("Output CellCenters Variable 11");
  outputCellCenters_Var12 = parameter_handler.get_bool("Output CellCenters Variable 12");
  outputCellCenters_Var13 = parameter_handler.get_bool("Output CellCenters Variable 13");
  outputCellCenters_Var14 = parameter_handler.get_bool("Output CellCenters Variable 14");
  outputCellCenters_Var15 = parameter_handler.get_bool("Output CellCenters Variable 15");
  outputCellCenters_Var16 = parameter_handler.get_bool("Output CellCenters Variable 16");
  outputCellCenters_Var17 = parameter_handler.get_bool("Output CellCenters Variable 17");
  outputCellCenters_Var18 = parameter_handler.get_bool("Output CellCenters Variable 18");
  outputCellCenters_Var19 = parameter_handler.get_bool("Output CellCenters Variable 19");
  outputCellCenters_Var20 = parameter_handler.get_bool("Output CellCenters Variable 20");
  outputCellCenters_Var21 = parameter_handler.get_bool("Output CellCenters Variable 21");
  outputCellCenters_Var22 = parameter_handler.get_bool("Output CellCenters Variable 22");
  outputCellCenters_Var23 = parameter_handler.get_bool("Output CellCenters Variable 23");
  outputCellCenters_Var24 = parameter_handler.get_bool("Output CellCenters Variable 24");
  outputCellCenters_Var25 = parameter_handler.get_bool("Output CellCenters Variable 25");
  outputCellCenters_Var26 = parameter_handler.get_bool("Output CellCenters Variable 26");
  outputCellCenters_Var27 = parameter_handler.get_bool("Output CellCenters Variable 27");
  outputCellCenters_Var28 = parameter_handler.get_bool("Output CellCenters Variable 28");

  flagBufferLayer = parameter_handler.get_bool("Output Buffer Layer Removal Feature");
  dimBufferLayer=parameter_handler.get_integer("Output Buffer Layer Removal dimension x0 y1 z2");
  lowerBufferLayer=parameter_handler.get_double("Output Buffer Layer Removal Lower Bound");
  upperBufferLayer=parameter_handler.get_double("Output Buffer Layer Removal Upper Bound");

  maxLinearSolverIterations=parameter_handler.get_integer("Maximum linear solver iterations");
  maxNonLinearIterations=parameter_handler.get_integer("Maximum non linear iterations");
  relLinearSolverTolerance=parameter_handler.get_double("Relative linear solver tolerance");
  absNonLinearTolerance=parameter_handler.get_double("Absolute nonLinear solver tolerance");
  relNonLinearTolerance=parameter_handler.get_double("Relative nonLinear solver tolerance");
  stopOnConvergenceFailure = parameter_handler.get_bool("Stop on convergence failure");
  enableAdaptiveTimeStepping = parameter_handler.get_bool("Enable adaptive Time stepping");
  adaptiveLoadStepFactor=parameter_handler.get_double("Adaptive load step factor");
  adaptiveLoadIncreaseFactor=parameter_handler.get_double("Adaptive load increase Factor");
  succesiveIncForIncreasingTimeStep=parameter_handler.get_double("Succesive increment for increasing time step");


  additionalVoxelInfo=parameter_handler.get_integer("Additional Voxel info");
  enableMultiphase = parameter_handler.get_bool("Enable Multiphase");
  numberofPhases=parameter_handler.get_integer("Number of Phases");

  enableUserMaterialModel = parameter_handler.get_bool("Enable User Material Model");

  enableUserMaterialModel1 = parameter_handler.get_bool("Enable User Material Model 1");
  numberofUserMatConstants1=parameter_handler.get_integer("Number of User Material Constants 1");
  numberofUserMatStateVar1=parameter_handler.get_integer("Number of User Material State Variables 1");
  UserMatConstants1=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material Constants 1")));
  UserMatStateVar1=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material State Variables Initial Values 1")));

  //elasticStiffness.reinit(6,6);
  elasticStiffness1.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 1"))));
  elasticStiffness1.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 2"))));
  elasticStiffness1.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 3"))));
  elasticStiffness1.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 4"))));
  elasticStiffness1.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 5"))));
  elasticStiffness1.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 6"))));

  enableAdvRateDepModel = parameter_handler.get_bool("Advanced Rate Dependent Model enabled");
  numSlipSystems1=parameter_handler.get_integer("Number of Slip Systems");

  latentHardeningRatioFileName1=parameter_handler.get("Latent Hardening Ratio filename");

  initialSlipResistance1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance")));
  initialHardeningModulus1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus")));
  powerLawExponent1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent")));
  saturationStress1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress")));
  slipDirectionsFile1 = parameter_handler.get("Slip Directions File");
  slipNormalsFile1 = parameter_handler.get("Slip Normals File");
  enableKinematicHardening1 = parameter_handler.get_bool("Enable Kinematic Hardening");

  if(enableKinematicHardening1){
    C_1_slip1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Slip Kinematic Hardening")));
    C_2_slip1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Slip Kinematic Hardening")));
  }


  enableTwinning1 = parameter_handler.get_bool("Twinning enabled");
  enableAdvancedTwinModel = parameter_handler.get_bool("Advanced Twinning Model enabled");
  enableOneTwinSys_Reorien = parameter_handler.get_bool("One twin system Reorientation enabled");
  if(enableTwinning1){
    numTwinSystems1=parameter_handler.get_integer("Number of Twin Systems");
    initialSlipResistanceTwin1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance Twin")));
    initialHardeningModulusTwin1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus Twin")));
    powerLawExponentTwin1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent Twin")));
    saturationStressTwin1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress Twin")));
    twinDirectionsFile1 = parameter_handler.get("Twin Directions File");
    twinNormalsFile1 = parameter_handler.get("Twin Normals File");
    if(enableKinematicHardening1){
      C_1_twin1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Twin Kinematic Hardening")));
      C_2_twin1 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Twin Kinematic Hardening")));
    }


  }
  else{
    pcout<<"Twinning is not enabled \n";

    numTwinSystems1=1;
    initialSlipResistanceTwin1.push_back(10e14);
    initialHardeningModulusTwin1.push_back(1);
    powerLawExponentTwin1.push_back(0);
    saturationStressTwin1.push_back(10e15);
    C_1_twin1.push_back(0);
    C_2_twin1.push_back(0);

  }


  twinThresholdFraction1=parameter_handler.get_double("Twin Threshold Fraction");
  twinSaturationFactor1=parameter_handler.get_double("Twin Saturation Factor");
  twinShear1=parameter_handler.get_double("Characteristic Twin Shear");
  criteriaTwinVisual=parameter_handler.get_double("Critical Value for Twin Visualization");


  modelStressTolerance=parameter_handler.get_double("Stress Tolerance");
  modelMaxSlipSearchIterations=parameter_handler.get_integer("Max Slip Search Iterations");
  modelMaxSolverIterations=parameter_handler.get_integer("Max Solver Iterations");
  modelMaxPlasticSlipL2Norm=parameter_handler.get_double("Max Plastic Slip L2 Norm");

  grainIDFile = parameter_handler.get("Grain ID file name");
  numPts.push_back(parameter_handler.get_double("Voxels in X direction"));
  if (dim > 1){
    numPts.push_back(parameter_handler.get_double("Voxels in Y direction"));
    if (dim > 2){
      numPts.push_back(parameter_handler.get_double("Voxels in Z direction"));
    }
  }

  grainOrientationsFile = parameter_handler.get("Orientations file name");
  headerLinesGrainIDFile=parameter_handler.get_integer("Header Lines GrainID File");
  enableElementDeletion=parameter_handler.get_bool("Use Element Deletion");
  deletionGrainID=parameter_handler.get_integer("GrainID of Deleted Elements");



  if (enableMultiphase&&(numberofPhases>=2)){
    enableUserMaterialModel2 = parameter_handler.get_bool("Enable User Material Model 2");
    numberofUserMatConstants2=parameter_handler.get_integer("Number of User Material Constants 2");
    numberofUserMatStateVar2=parameter_handler.get_integer("Number of User Material State Variables 2");
    UserMatConstants2=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material Constants 2")));
    UserMatStateVar2=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material State Variables Initial Values 2")));

    //elasticStiffness.reinit(6,6);
    elasticStiffness2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 2 row 1"))));
    elasticStiffness2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 2 row 2"))));
    elasticStiffness2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 2 row 3"))));
    elasticStiffness2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 2 row 4"))));
    elasticStiffness2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 2 row 5"))));
    elasticStiffness2.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 2 row 6"))));

    numSlipSystems2=parameter_handler.get_integer("Number of Slip Systems 2");

    latentHardeningRatioFileName2=parameter_handler.get("Latent Hardening Ratio filename 2");

    initialSlipResistance2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance 2")));
    initialHardeningModulus2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus 2")));
    powerLawExponent2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent 2")));
    saturationStress2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress 2")));
    slipDirectionsFile2 = parameter_handler.get("Slip Directions File 2");
    slipNormalsFile2 = parameter_handler.get("Slip Normals File 2");
    enableKinematicHardening2 = parameter_handler.get_bool("Enable Kinematic Hardening 2");

    if(enableKinematicHardening2){
      C_1_slip2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Slip Kinematic Hardening 2")));
      C_2_slip2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Slip Kinematic Hardening 2")));
    }


    enableTwinning2 = parameter_handler.get_bool("Twinning enabled 2");
    if(enableTwinning2){
      numTwinSystems2=parameter_handler.get_integer("Number of Twin Systems 2");
      initialSlipResistanceTwin2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance Twin 2")));
      initialHardeningModulusTwin2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus Twin 2")));
      powerLawExponentTwin2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent Twin 2")));
      saturationStressTwin2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress Twin 2")));
      twinDirectionsFile2 = parameter_handler.get("Twin Directions File 2");
      twinNormalsFile2 = parameter_handler.get("Twin Normals File 2");
      if(enableKinematicHardening2){
        C_1_twin2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Twin Kinematic Hardening 2")));
        C_2_twin2 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Twin Kinematic Hardening 2")));
      }


    }
    else{
      pcout<<"Twinning is not enabled \n";

      numTwinSystems2=1;
      initialSlipResistanceTwin2.push_back(10e14);
      initialHardeningModulusTwin2.push_back(1);
      powerLawExponentTwin2.push_back(0);
      saturationStressTwin2.push_back(10e15);
      C_1_twin2.push_back(0);
      C_2_twin2.push_back(0);

    }

    twinThresholdFraction2=parameter_handler.get_double("Twin Threshold Fraction 2");
    twinSaturationFactor2=parameter_handler.get_double("Twin Saturation Factor 2");
    twinShear2=parameter_handler.get_double("Characteristic Twin Shear 2");

    if (numberofPhases>=3){
      enableUserMaterialModel3 = parameter_handler.get_bool("Enable User Material Model 3");
      numberofUserMatConstants3=parameter_handler.get_integer("Number of User Material Constants 3");
      numberofUserMatStateVar3=parameter_handler.get_integer("Number of User Material State Variables 3");
      UserMatConstants3=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material Constants 3")));
      UserMatStateVar3=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material State Variables Initial Values 3")));
      //elasticStiffness.reinit(6,6);
      elasticStiffness3.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 3 row 1"))));
      elasticStiffness3.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 3 row 2"))));
      elasticStiffness3.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 3 row 3"))));
      elasticStiffness3.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 3 row 4"))));
      elasticStiffness3.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 3 row 5"))));
      elasticStiffness3.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 3 row 6"))));

      numSlipSystems3=parameter_handler.get_integer("Number of Slip Systems 3");

      latentHardeningRatioFileName3=parameter_handler.get("Latent Hardening Ratio filename 3");

      initialSlipResistance3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance 3")));
      initialHardeningModulus3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus 3")));
      powerLawExponent3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent 3")));
      saturationStress3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress 3")));
      slipDirectionsFile3 = parameter_handler.get("Slip Directions File 3");
      slipNormalsFile3 = parameter_handler.get("Slip Normals File 3");
      enableKinematicHardening3 = parameter_handler.get_bool("Enable Kinematic Hardening 3");

      if(enableKinematicHardening3){
        C_1_slip3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Slip Kinematic Hardening 3")));
        C_2_slip3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Slip Kinematic Hardening 3")));
      }


      enableTwinning3 = parameter_handler.get_bool("Twinning enabled 3");
      if(enableTwinning3){
        numTwinSystems3=parameter_handler.get_integer("Number of Twin Systems 3");
        initialSlipResistanceTwin3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance Twin 3")));
        initialHardeningModulusTwin3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus Twin 3")));
        powerLawExponentTwin3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent Twin 3")));
        saturationStressTwin3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress Twin 3")));
        twinDirectionsFile3 = parameter_handler.get("Twin Directions File 3");
        twinNormalsFile3 = parameter_handler.get("Twin Normals File 3");
        if(enableKinematicHardening3){
          C_1_twin3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Twin Kinematic Hardening 3")));
          C_2_twin3 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Twin Kinematic Hardening 3")));
        }


      }
      else{
        pcout<<"Twinning is not enabled \n";

        numTwinSystems3=1;
        initialSlipResistanceTwin3.push_back(10e14);
        initialHardeningModulusTwin3.push_back(1);
        powerLawExponentTwin3.push_back(0);
        saturationStressTwin3.push_back(10e15);
        C_1_twin3.push_back(0);
        C_2_twin3.push_back(0);

      }

      twinThresholdFraction3=parameter_handler.get_double("Twin Threshold Fraction 3");
      twinSaturationFactor3=parameter_handler.get_double("Twin Saturation Factor 3");
      twinShear3=parameter_handler.get_double("Characteristic Twin Shear 3");

      if (numberofPhases>=4){
        enableUserMaterialModel4 = parameter_handler.get_bool("Enable User Material Model 4");
        numberofUserMatConstants4=parameter_handler.get_integer("Number of User Material Constants 4");
        numberofUserMatStateVar4=parameter_handler.get_integer("Number of User Material State Variables 4");
        UserMatConstants4=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material Constants 4")));
        UserMatStateVar4=dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("User Material State Variables Initial Values 4")));

        //elasticStiffness.reinit(6,6);
        elasticStiffness4.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 4 row 1"))));
        elasticStiffness4.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 4 row 2"))));
        elasticStiffness4.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 4 row 3"))));
        elasticStiffness4.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 4 row 4"))));
        elasticStiffness4.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 4 row 5"))));
        elasticStiffness4.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness 4 row 6"))));

        numSlipSystems4=parameter_handler.get_integer("Number of Slip Systems 4");

        latentHardeningRatioFileName4=parameter_handler.get("Latent Hardening Ratio filename 4");

        initialSlipResistance4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance 4")));
        initialHardeningModulus4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus 4")));
        powerLawExponent4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent 4")));
        saturationStress4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress 4")));
        slipDirectionsFile4 = parameter_handler.get("Slip Directions File 4");
        slipNormalsFile4 = parameter_handler.get("Slip Normals File 4");
        enableKinematicHardening4 = parameter_handler.get_bool("Enable Kinematic Hardening 4");

        if(enableKinematicHardening4){
          C_1_slip4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Slip Kinematic Hardening 4")));
          C_2_slip4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Slip Kinematic Hardening 4")));
        }


        enableTwinning4 = parameter_handler.get_bool("Twinning enabled 4");
        if(enableTwinning4){
          numTwinSystems4=parameter_handler.get_integer("Number of Twin Systems 4");
          initialSlipResistanceTwin4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance Twin 4")));
          initialHardeningModulusTwin4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus Twin 4")));
          powerLawExponentTwin4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent Twin 4")));
          saturationStressTwin4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress Twin 4")));
          twinDirectionsFile4 = parameter_handler.get("Twin Directions File 4");
          twinNormalsFile4 = parameter_handler.get("Twin Normals File 4");
          if(enableKinematicHardening4){
            C_1_twin4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_1 Twin Kinematic Hardening 4")));
            C_2_twin4 = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("C_2 Twin Kinematic Hardening 4")));
          }


        }
        else{
          pcout<<"Twinning is not enabled \n";

          numTwinSystems4=1;
          initialSlipResistanceTwin4.push_back(10e14);
          initialHardeningModulusTwin4.push_back(1);
          powerLawExponentTwin4.push_back(0);
          saturationStressTwin4.push_back(10e15);
          C_1_twin4.push_back(0);
          C_2_twin4.push_back(0);

        }

        twinThresholdFraction4=parameter_handler.get_double("Twin Threshold Fraction 4");
        twinSaturationFactor4=parameter_handler.get_double("Twin Saturation Factor 4");
        twinShear4=parameter_handler.get_double("Characteristic Twin Shear 4");
      }
    }

  }

}

void userInputParameters::declare_parameters(dealii::ParameterHandler & parameter_handler){

  parameter_handler.declare_entry("Flag To Use Taylor Model","false",dealii::Patterns::Bool(),"Flag To Use Taylor Model");
  parameter_handler.declare_entry("Number of dimensions","-1",dealii::Patterns::Integer(),"Number of physical dimensions for the simulation");
  parameter_handler.declare_entry("Order of finite elements","-1",dealii::Patterns::Integer(),"Basis function interpolation order (1-linear)");
  parameter_handler.declare_entry("Order of quadrature","-1",dealii::Patterns::Integer(),"Quadrature point order n^3 (2->8 quadrature points)");

  parameter_handler.declare_entry("Domain size X","-1",dealii::Patterns::Double(),"The size of the domain in the x direction.");
  parameter_handler.declare_entry("Domain size Y","-1",dealii::Patterns::Double(),"The size of the domain in the y direction.");
  parameter_handler.declare_entry("Domain size Z","-1",dealii::Patterns::Double(),"The size of the domain in the z direction.");
  parameter_handler.declare_entry("Subdivisions X","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the x direction.");
  parameter_handler.declare_entry("Subdivisions Y","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the y direction.");
  parameter_handler.declare_entry("Subdivisions Z","1",dealii::Patterns::Integer(),"The number of mesh subdivisions in the z direction.");
  parameter_handler.declare_entry("Refine factor","-1",dealii::Patterns::Integer(),"The number of initial refinements of the coarse mesh.");

  parameter_handler.declare_entry("Write Mesh To EPS","false",dealii::Patterns::Bool(),"Only written for serial runs and if number of elements < 10000");

  parameter_handler.declare_entry("Use external mesh","false",dealii::Patterns::Bool(),"Flag to indicate whether to use external mesh");
  parameter_handler.declare_entry("Name of file containing external mesh","",dealii::Patterns::Anything(),"Name of external mesh file");
  parameter_handler.declare_entry("External mesh parameter","0",dealii::Patterns::Double(),"The external mesh parameter: The ratio of defiend region size to the Domain size");

  parameter_handler.declare_entry("Time increments","-1",dealii::Patterns::Double(),"delta T for every increment");
  parameter_handler.declare_entry("critical DeltaF Criteria","10000",dealii::Patterns::Double(),"critical DeltaF Criteria");
  parameter_handler.declare_entry("Total time","-1",dealii::Patterns::Double(),"Total simulation time");

  parameter_handler.declare_entry("Use Simple BCs","true",dealii::Patterns::Bool(),"Flag to indicate whether to use Simple (Basic) BCs");
  parameter_handler.declare_entry("Boundary condition filename","boundaryConditions.txt",dealii::Patterns::Anything(),"File name containing BC information");
  parameter_handler.declare_entry("BC file number of header lines","1",dealii::Patterns::Integer(),"BC file number of header lines");
  parameter_handler.declare_entry("Number of boundary conditions","1",dealii::Patterns::Integer(),"Number of boundary conditions");

  parameter_handler.declare_entry("Use velocity gradient BC","false",dealii::Patterns::Bool(),"Flag to indicate whether to use velocity gradient tensor to apply BCs");
  parameter_handler.declare_entry("Velocity gradient row 1","",dealii::Patterns::List(dealii::Patterns::Double()),"Velocity gradient tensor including the multiplication factor ");
  parameter_handler.declare_entry("Velocity gradient row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Velocity gradient tensor including the multiplication factor ");
  parameter_handler.declare_entry("Velocity gradient row 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Velocity gradient tensor including the multiplication factor ");

  parameter_handler.declare_entry("Use Tabular BCs","false",dealii::Patterns::Bool(),"Flag to indicate whether to use Tabular BCs");
  parameter_handler.declare_entry("Tabular Boundary condition filename","tabularBoundaryConditions.txt",dealii::Patterns::Anything(),"File name containing Tabular BC information");
  parameter_handler.declare_entry("Number of time data for Tabular BCs","0",dealii::Patterns::Integer(),"Number of time data for Tabular BCs (it includes the initial BCs)");
  parameter_handler.declare_entry("Number of tabular boundary conditions","1",dealii::Patterns::Integer(),"Number of tabular boundary conditions");
  parameter_handler.declare_entry("Tabular Time Table","",dealii::Patterns::List(dealii::Patterns::Double()),"Table for Time intervals of Tabular BCs");

  parameter_handler.declare_entry("Use Neumann BCs","false",dealii::Patterns::Bool(),"Flag to indicate whether to Use Neumann BCs");
  parameter_handler.declare_entry("Tabular Neumann Boundary condition filename","tabularBoundaryConditions.txt",dealii::Patterns::Anything(),"File name containing Tabular Neumann Boundary condition");
  parameter_handler.declare_entry("Number of time data for Tabular Neumann BCs","0",dealii::Patterns::Integer(),"Number of time data for Tabular Neumann BCs (it includes the initial BCs)");
  parameter_handler.declare_entry("Number of tabular Neumann boundary conditions","1",dealii::Patterns::Integer(),"Number of tabular Neumann boundary conditions");
  parameter_handler.declare_entry("Boundary IDs of Neumann BCs","",dealii::Patterns::List(dealii::Patterns::Integer()),"Boundary IDs of Neumann BCs");
  parameter_handler.declare_entry("dof of Neumann BCs","",dealii::Patterns::List(dealii::Patterns::Integer()),"dof of Neumann BCs");
  parameter_handler.declare_entry("Tabular Time Neumann BCs","",dealii::Patterns::List(dealii::Patterns::Double()),"Table for Time intervals of Tabular Neumann BCs");


  parameter_handler.declare_entry("Use Torsion BCs","false",dealii::Patterns::Bool(),"Flag to indicate whether to use Torsion BCs");
  parameter_handler.declare_entry("Torsion Axis  x 0  y 1  z 2","2",dealii::Patterns::Integer(),"Torsion Axis  x 0  y 1  z 2");
  parameter_handler.declare_entry("Tabular Time Table for Torsion","",dealii::Patterns::List(dealii::Patterns::Double()),"Table for Time intervals of Tabular Torsion BCs");
  parameter_handler.declare_entry("Tabular Torsion BCs","",dealii::Patterns::List(dealii::Patterns::Double()),"Tabular Torsion BCs (Angular velocity)");
  parameter_handler.declare_entry("Center point for Torsion","",dealii::Patterns::List(dealii::Patterns::Double()),"Center point for Torsion (x,y) for torsion axis=z or (y,z) for torsion axis=x or (z,x) for torsion axis=y");

  parameter_handler.declare_entry("Use DIC pipeline","false",dealii::Patterns::Bool(),"Flag to indicate whether to use DIC experiment pipeline");
  parameter_handler.declare_entry("Number of Input data for DIC experiment","0",dealii::Patterns::Integer(),"Number of Input data for DIC experiment (it includes the initial BCs)");
  parameter_handler.declare_entry("DIC Time Table","",dealii::Patterns::List(dealii::Patterns::Double()),"Table for Time intervals of DIC experiment input");
  parameter_handler.declare_entry("Number of Points in DIC input in X direction","0",dealii::Patterns::Integer(),"Number of Points in DIC input in X direction");
  parameter_handler.declare_entry("Number of Points in DIC input in Y direction","0",dealii::Patterns::Integer(),"Number of Points in DIC input in Y direction");
  parameter_handler.declare_entry("Number of Points in DIC input in Z direction","0",dealii::Patterns::Integer(),"Number of Points in DIC input in Z direction");
  parameter_handler.declare_entry("DIC Boundary condition filename 1","DICboundaryConditions1.txt",dealii::Patterns::Anything(),"DIC Boundary condition filename 1 (x == 0.0)");
  parameter_handler.declare_entry("DIC Boundary condition filename 2","DICboundaryConditions2.txt",dealii::Patterns::Anything(),"DIC Boundary condition filename 2 (x == spanX)");
  parameter_handler.declare_entry("DIC Boundary condition filename 3","DICboundaryConditions3.txt",dealii::Patterns::Anything(),"DIC Boundary condition filename 3 (y == 0.0)");
  parameter_handler.declare_entry("DIC Boundary condition filename 4","DICboundaryConditions4.txt",dealii::Patterns::Anything(),"DIC Boundary condition filename 4 (y == spanY)");


  parameter_handler.declare_entry("Enable cyclic loading","false",dealii::Patterns::Bool(),"Flag to indicate if cyclic loading is enabled");
  parameter_handler.declare_entry("Cyclic loading face","1",dealii::Patterns::Integer(),"Face that is cyclically loaded");
  parameter_handler.declare_entry("Cyclic loading direction","1",dealii::Patterns::Integer(),"Direction that is cyclically loaded");
  parameter_handler.declare_entry("Quarter cycle time","-1",dealii::Patterns::Double(),"Time for finishing quarter of a cyclic loading cycle. One cycle is time taken for starting from 0 displacement and ending at 0 displacement");

  parameter_handler.declare_entry("Enable Nodal Displacement BCs","false",dealii::Patterns::Bool(),"Flag to indicate if Nodal Displacement BCs is enabled");
  parameter_handler.declare_entry("Number of Nodal Displacement BCs","1",dealii::Patterns::Integer(),"Number of Nodal Displacement BCs");
  parameter_handler.declare_entry("Tolerance for Nodal Displacement BCs","0",dealii::Patterns::Double(),"Tolerance for Nodal Displacement BCs: The ratio of defiend tolerance size to the Domain size");
  parameter_handler.declare_entry("Nodal Displacement BCs filename","NodalDisplacementBCs.txt",dealii::Patterns::Anything(),"File name containing Nodal Displacement BCs");

  parameter_handler.declare_entry("Use Periodic BCs","false",dealii::Patterns::Bool(),"Flag to indicate whether to use periodic BCs");
  parameter_handler.declare_entry("Periodic Boundary condition Constraint filename","PeriodicBCsConstraints.txt",dealii::Patterns::Anything(),"File name containing Periodic Boundary condition Constraint");
  parameter_handler.declare_entry("Number of Vertices Constraints","0",dealii::Patterns::Integer(),"Number of Vertices Constraints");
  parameter_handler.declare_entry("Number of Edges Constraints","0",dealii::Patterns::Integer(),"Number of Edges Constraints");
  parameter_handler.declare_entry("Number of Faces Constraints","0",dealii::Patterns::Integer(),"Number of Faces Constraints");
  parameter_handler.declare_entry("Vertices Periodic BCs row 1","",dealii::Patterns::List(dealii::Patterns::Integer()),"Vertices Periodic BCs row 1 ");
  parameter_handler.declare_entry("Vertices Periodic BCs row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Vertices Periodic BCs row 2 ");
  parameter_handler.declare_entry("Use Tabular Periodic BCs","false",dealii::Patterns::Bool(),"Flag to indicate whether to use Tabular Periodic BCs");
  parameter_handler.declare_entry("periodic Tabular time","-1",dealii::Patterns::Double(),"Time for finishing the largest displacement defined in Vertices Periodic BCs row 2");
  parameter_handler.declare_entry("Tabular Periodic Time Table","",dealii::Patterns::List(dealii::Patterns::Double()),"Table for Time intervals of Tabular Periodic BCs");
  parameter_handler.declare_entry("Tabular Periodic Time Table Coefficient","",dealii::Patterns::List(dealii::Patterns::Integer()),"Table for Coefficient of Tabular Periodic BCs ");
    //INDENTATION 2
  parameter_handler.declare_entry("Use Indentation BCs","false",dealii::Patterns::Bool(),"Flag to indicate whether to use indentation BCs");
  parameter_handler.declare_entry("Indentation Active Set Freeze Iterations", "0",dealii::Patterns::Integer(),"How many non linear iterations should active set be frozen for if it is changed");
  parameter_handler.declare_entry("Active Set Lambda Tolerance", "0.0",dealii::Patterns::Double(),"How much to bias against errors in lambda/mass (positive biases into active set)");
  parameter_handler.declare_entry("Indentation Active Set Freeze on Iteration 1","false",dealii::Patterns::Bool(),"Flag to indicate whether to freeze active set for first iteration");
  parameter_handler.declare_entry("Active Set Criterion Coefficient","100",dealii::Patterns::List(dealii::Patterns::Double()),"How much to multiply stiffness by to determine the active set criterion");
  parameter_handler.declare_entry("Indentation Key Frame Number","3",dealii::Patterns::Integer(),"How many key positions does the indenter visit?");
  parameter_handler.declare_entry("Indentation Key Frame Increment Numbers","0, 80, 160",dealii::Patterns::List(dealii::Patterns::Integer()),"Number of times to refine cells in refinement zone");
  parameter_handler.declare_entry("Indentation Boundary condition Constraint filename","IndentationBCsConstraints.txt",dealii::Patterns::Anything(),"File name containing Indentation Boundary condition Constraint");
  parameter_handler.declare_entry("Refinement Zone Center","0.0, 0.0, 0.0",dealii::Patterns::List(dealii::Patterns::Double()),"x, y, z Coordinates of refinement zone center");
  parameter_handler.declare_entry("Refinement Zone Size","0.0",dealii::Patterns::List(dealii::Patterns::Double()),"the radius of the refinement zone for indentation");
  parameter_handler.declare_entry("Refinement Factor","0",dealii::Patterns::List(dealii::Patterns::Integer()),"Number of times to refine cells in refinement zone");
  parameter_handler.declare_entry("Debug Indentation Active Set","false",dealii::Patterns::Bool(),"Flag to indicate whether to debug indentation active set");

  //Isotropic material 2
  parameter_handler.declare_entry("Continuum Isotropic","false",dealii::Patterns::Bool(),"Flag to indicate whether to use isotropic material");
  parameter_handler.declare_entry("lame_lambda","-1",dealii::Patterns::Double(),"Lame' parameter lambda");
  parameter_handler.declare_entry("lame_mu","-1",dealii::Patterns::Double(),"Lame' parameter mu");
  parameter_handler.declare_entry("yield_stress","-1",dealii::Patterns::Double(),"Isotropic yield stress (Kirchhoff)");
  parameter_handler.declare_entry("strain_hardening","-1",dealii::Patterns::Double(),"Linear Isotropic strain hardening coefficient");
  parameter_handler.declare_entry("kinematic_hardening","-1",dealii::Patterns::Double(),"Kinematic strain hardening coefficient");
  parameter_handler.declare_entry("strain_energy_function","quadlog",dealii::Patterns::Anything(),"Strain energy density function (quadlog, stvenkir, or neohook)");
  parameter_handler.declare_entry("yield_function","von_mises",dealii::Patterns::Anything(),"Yield function (currently only von_mises)");
  parameter_handler.declare_entry("iso_hardening_function","linear_hardening",dealii::Patterns::Anything(),"Isotropic hardening function");



  parameter_handler.declare_entry("Write Output","false",dealii::Patterns::Bool(),"Flag to write output vtu and pvtu files");
  parameter_handler.declare_entry("Output Directory",".",dealii::Patterns::Anything(),"Output Directory");

  parameter_handler.declare_entry("Tabular Output","false",dealii::Patterns::Bool(),"Flag to use Tabular Output");
  parameter_handler.declare_entry("Tabular Time Output Table","",dealii::Patterns::List(dealii::Patterns::Double()),"Table for Time Outputs");

  parameter_handler.declare_entry("Skip Output Steps","-1",dealii::Patterns::Integer(),"Skip Output Steps");
  parameter_handler.declare_entry("Write Quadrature Output","false",dealii::Patterns::Bool(),"Flag to write quadrature output");
  parameter_handler.declare_entry("Write Grain Averaged Output","false",dealii::Patterns::Bool(),"Flag to Write Grain Averaged Output");
  parameter_handler.declare_entry("Skip Quadrature Output Steps","-1",dealii::Patterns::Integer(),"Skip Quadrature Output Steps");
  parameter_handler.declare_entry("Skip Grain Averaged Output Steps","-1",dealii::Patterns::Integer(),"Skip Grain Averaged Output Steps");
  parameter_handler.declare_entry("Number of Grain Averaged Output Variables","-1",dealii::Patterns::Integer(),"Number of Grain Averaged Output Variables");
  parameter_handler.declare_entry("Output Equivalent strain","false",dealii::Patterns::Bool(),"Output Equivalent strain");
  parameter_handler.declare_entry("Output Equivalent stress","false",dealii::Patterns::Bool(),"Output Equivalent stress");
  parameter_handler.declare_entry("Output Time","false",dealii::Patterns::Bool(),"Output Time in stressstrain");
  parameter_handler.declare_entry("Output Equivalent plastic strain alpha","false",dealii::Patterns::Bool(),"Output Equivalent plastic strain alpha");

  parameter_handler.declare_entry("Output Indenter Load","false",dealii::Patterns::Bool(),"Output Indenter Load"); //Indentation 3
  parameter_handler.declare_entry("Output Grain ID","false",dealii::Patterns::Bool(),"Output Grain ID");
  parameter_handler.declare_entry("Output Twin fractions","false",dealii::Patterns::Bool(),"Output Twin fractions");

  parameter_handler.declare_entry("Number of Taylor Substeps","1",dealii::Patterns::Integer(),"Number of Taylor Substeps");

  parameter_handler.declare_entry("Output Userdefined Average Variable","false",dealii::Patterns::Bool(),"Flag to enable Output Userdefined Average Variable");
  parameter_handler.declare_entry("Number of Output Userdefined Average Variable","-1",dealii::Patterns::Integer(),"Number of Output Userdefined Average Variable");

  parameter_handler.declare_entry("Output Variable 1","false",dealii::Patterns::Bool(),"Output Variable 1");
  parameter_handler.declare_entry("Output Variable 2","false",dealii::Patterns::Bool(),"Output Variable 2");
  parameter_handler.declare_entry("Output Variable 3","false",dealii::Patterns::Bool(),"Output Variable 3");
  parameter_handler.declare_entry("Output Variable 4","false",dealii::Patterns::Bool(),"Output Variable 4");
  parameter_handler.declare_entry("Output Variable 5","false",dealii::Patterns::Bool(),"Output Variable 5");
  parameter_handler.declare_entry("Output Variable 6","false",dealii::Patterns::Bool(),"Output Variable 6");
  parameter_handler.declare_entry("Output Variable 7","false",dealii::Patterns::Bool(),"Output Variable 7");
  parameter_handler.declare_entry("Output Variable 8","false",dealii::Patterns::Bool(),"Output Variable 8");
  parameter_handler.declare_entry("Output Variable 9","false",dealii::Patterns::Bool(),"Output Variable 9");
  parameter_handler.declare_entry("Output Variable 10","false",dealii::Patterns::Bool(),"Output Variable 10");
  parameter_handler.declare_entry("Output Variable 11","false",dealii::Patterns::Bool(),"Output Variable 11");
  parameter_handler.declare_entry("Output Variable 12","false",dealii::Patterns::Bool(),"Output Variable 12");
  parameter_handler.declare_entry("Output Variable 13","false",dealii::Patterns::Bool(),"Output Variable 13");
  parameter_handler.declare_entry("Output Variable 14","false",dealii::Patterns::Bool(),"Output Variable 14");
  parameter_handler.declare_entry("Output Variable 15","false",dealii::Patterns::Bool(),"Output Variable 15");
  parameter_handler.declare_entry("Output Variable 16","false",dealii::Patterns::Bool(),"Output Variable 16");
  parameter_handler.declare_entry("Output Variable 17","false",dealii::Patterns::Bool(),"Output Variable 17");
  parameter_handler.declare_entry("Output Variable 18","false",dealii::Patterns::Bool(),"Output Variable 18");
  parameter_handler.declare_entry("Output Variable 19","false",dealii::Patterns::Bool(),"Output Variable 19");
  parameter_handler.declare_entry("Output Variable 20","false",dealii::Patterns::Bool(),"Output Variable 20");
  parameter_handler.declare_entry("Output Variable 21","false",dealii::Patterns::Bool(),"Output Variable 21");
  parameter_handler.declare_entry("Output Variable 22","false",dealii::Patterns::Bool(),"Output Variable 22");
  parameter_handler.declare_entry("Output Variable 23","false",dealii::Patterns::Bool(),"Output Variable 23");
  parameter_handler.declare_entry("Output Variable 24","false",dealii::Patterns::Bool(),"Output Variable 24");

  parameter_handler.declare_entry("Output CellCenters Variable 1","false",dealii::Patterns::Bool(),"Output CellCenters Variable 1");
  parameter_handler.declare_entry("Output CellCenters Variable 2","false",dealii::Patterns::Bool(),"Output CellCenters Variable 2");
  parameter_handler.declare_entry("Output CellCenters Variable 3","false",dealii::Patterns::Bool(),"Output CellCenters Variable 3");
  parameter_handler.declare_entry("Output CellCenters Variable 4","false",dealii::Patterns::Bool(),"Output CellCenters Variable 4");
  parameter_handler.declare_entry("Output CellCenters Variable 5","false",dealii::Patterns::Bool(),"Output CellCenters Variable 5");
  parameter_handler.declare_entry("Output CellCenters Variable 6","false",dealii::Patterns::Bool(),"Output CellCenters Variable 6");
  parameter_handler.declare_entry("Output CellCenters Variable 7","false",dealii::Patterns::Bool(),"Output CellCenters Variable 7");
  parameter_handler.declare_entry("Output CellCenters Variable 8","false",dealii::Patterns::Bool(),"Output CellCenters Variable 8");
  parameter_handler.declare_entry("Output CellCenters Variable 9","false",dealii::Patterns::Bool(),"Output CellCenters Variable 9");
  parameter_handler.declare_entry("Output CellCenters Variable 10","false",dealii::Patterns::Bool(),"Output CellCenters Variable 10");
  parameter_handler.declare_entry("Output CellCenters Variable 11","false",dealii::Patterns::Bool(),"Output CellCenters Variable 11");
  parameter_handler.declare_entry("Output CellCenters Variable 12","false",dealii::Patterns::Bool(),"Output CellCenters Variable 12");
  parameter_handler.declare_entry("Output CellCenters Variable 13","false",dealii::Patterns::Bool(),"Output CellCenters Variable 13");
  parameter_handler.declare_entry("Output CellCenters Variable 14","false",dealii::Patterns::Bool(),"Output CellCenters Variable 14");
  parameter_handler.declare_entry("Output CellCenters Variable 15","false",dealii::Patterns::Bool(),"Output CellCenters Variable 15");
  parameter_handler.declare_entry("Output CellCenters Variable 16","false",dealii::Patterns::Bool(),"Output CellCenters Variable 16");
  parameter_handler.declare_entry("Output CellCenters Variable 17","false",dealii::Patterns::Bool(),"Output CellCenters Variable 17");
  parameter_handler.declare_entry("Output CellCenters Variable 18","false",dealii::Patterns::Bool(),"Output CellCenters Variable 18");
  parameter_handler.declare_entry("Output CellCenters Variable 19","false",dealii::Patterns::Bool(),"Output CellCenters Variable 19");
  parameter_handler.declare_entry("Output CellCenters Variable 20","false",dealii::Patterns::Bool(),"Output CellCenters Variable 20");
  parameter_handler.declare_entry("Output CellCenters Variable 21","false",dealii::Patterns::Bool(),"Output CellCenters Variable 21");
  parameter_handler.declare_entry("Output CellCenters Variable 22","false",dealii::Patterns::Bool(),"Output CellCenters Variable 22");
  parameter_handler.declare_entry("Output CellCenters Variable 23","false",dealii::Patterns::Bool(),"Output CellCenters Variable 23");
  parameter_handler.declare_entry("Output CellCenters Variable 24","false",dealii::Patterns::Bool(),"Output CellCenters Variable 24");
  parameter_handler.declare_entry("Output CellCenters Variable 25","false",dealii::Patterns::Bool(),"Output CellCenters Variable 25");
  parameter_handler.declare_entry("Output CellCenters Variable 26","false",dealii::Patterns::Bool(),"Output CellCenters Variable 26");
  parameter_handler.declare_entry("Output CellCenters Variable 27","false",dealii::Patterns::Bool(),"Output CellCenters Variable 27");
  parameter_handler.declare_entry("Output CellCenters Variable 28","false",dealii::Patterns::Bool(),"Output CellCenters Variable 28");

  parameter_handler.declare_entry("Output Buffer Layer Removal Feature","false",dealii::Patterns::Bool(),"Output Buffer Layer Removal Feature");
  parameter_handler.declare_entry("Output Buffer Layer Removal dimension x0 y1 z2","0",dealii::Patterns::Integer(), "MOutput Buffer Layer Removal dimension x0 y1 z2");
  parameter_handler.declare_entry("Output Buffer Layer Removal Lower Bound","0",dealii::Patterns::Double(),"Output Buffer Layer Removal Lower Bound");
  parameter_handler.declare_entry("Output Buffer Layer Removal Upper Bound","1",dealii::Patterns::Double(),"Output Buffer Layer Removal Upper Bound");

  parameter_handler.declare_entry("Maximum linear solver iterations","-1",dealii::Patterns::Integer(), "Maximum iterations for linear solver");
  parameter_handler.declare_entry("Maximum non linear iterations","-1",dealii::Patterns::Integer(),"Maximum no. of non-linear iterations");
  parameter_handler.declare_entry("Relative linear solver tolerance","-1",dealii::Patterns::Double(),"Relative linear solver tolerance");
  parameter_handler.declare_entry("Absolute nonLinear solver tolerance","-1",dealii::Patterns::Double(),"Non-linear solver tolerance");
  parameter_handler.declare_entry("Relative nonLinear solver tolerance","-1",dealii::Patterns::Double(),"Relative non-linear solver tolerance");
  parameter_handler.declare_entry("Stop on convergence failure","false",dealii::Patterns::Bool(),"Flag to stop problem if convergence fails");
  parameter_handler.declare_entry("Enable adaptive Time stepping","false",dealii::Patterns::Bool(),"Flag to enable adaptive time steps");
  parameter_handler.declare_entry("Adaptive load step factor","-1",dealii::Patterns::Double(),"Load step factor");
  parameter_handler.declare_entry("Adaptive load increase Factor","-1",dealii::Patterns::Double(),"adaptive Load Increase Factor");
  parameter_handler.declare_entry("Succesive increment for increasing time step","-1",dealii::Patterns::Double(),"Succesive Inc For Increasing Time Step");


  parameter_handler.declare_entry("Additional Voxel info","0",dealii::Patterns::Integer(),"Number of Additional Voxel info Besides three orientation components and Phase if multiphase is enabled");
  parameter_handler.declare_entry("Enable Multiphase","false",dealii::Patterns::Bool(),"Flag to indicate if Multiphase is enabled");
  parameter_handler.declare_entry("Number of Phases","1",dealii::Patterns::Integer(),"Number of phases in the sample");

  parameter_handler.declare_entry("Enable User Material Model","false",dealii::Patterns::Bool(),"Flag to indicate if User Material Model is enabled");

  parameter_handler.declare_entry("Enable User Material Model 1","false",dealii::Patterns::Bool(),"Flag to indicate if User Material Model is enabled Phase 1");
  parameter_handler.declare_entry("Number of User Material Constants 1","0",dealii::Patterns::Integer(),"Number of User Material Constants in a Material model Phase 1");
  parameter_handler.declare_entry("Number of User Material State Variables 1","0",dealii::Patterns::Integer(),"Number of User Material State Variables in a Material model Phase 1");
  parameter_handler.declare_entry("User Material Constants 1","",dealii::Patterns::List(dealii::Patterns::Double()),"Material Constants in a Material model Phase 1");
  parameter_handler.declare_entry("User Material State Variables Initial Values 1","",dealii::Patterns::List(dealii::Patterns::Double()),"Material State Variables in a Material model Phase 1");

  parameter_handler.declare_entry("Elastic Stiffness row 1","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 3","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 4","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 5","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 6","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");

  parameter_handler.declare_entry("Advanced Rate Dependent Model enabled","false",dealii::Patterns::Bool(),"Flag to indicate if Advanced Rate Dependent Model enabled");
  parameter_handler.declare_entry("Number of Slip Systems","-1",dealii::Patterns::Integer(),"Number of Slip Systems");

  parameter_handler.declare_entry("Latent Hardening Ratio filename","LatentHardeningRatio.txt",dealii::Patterns::Anything(),"Latent Hardening Ratio filename");

  parameter_handler.declare_entry("Initial Slip Resistance","",dealii::Patterns::List(dealii::Patterns::Double()),"RSS of the slip sytems");
  parameter_handler.declare_entry("Initial Hardening Modulus","",dealii::Patterns::List(dealii::Patterns::Double()),"Heardening moduli of slip systems");
  parameter_handler.declare_entry("Power Law Exponent","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law coefficient");
  parameter_handler.declare_entry("Saturation Stress","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress");
  parameter_handler.declare_entry("Slip Directions File","",dealii::Patterns::Anything(),"Slip Directions File");
  parameter_handler.declare_entry("Slip Normals File","",dealii::Patterns::Anything(),"Slip Normals File");
  parameter_handler.declare_entry("Enable Kinematic Hardening","false",dealii::Patterns::Bool(),"Flag to indicate if kinematic hardening is enabled");
  parameter_handler.declare_entry("C_1 Slip Kinematic Hardening","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Slip Kinematic Hardening parameters");
  parameter_handler.declare_entry("C_2 Slip Kinematic Hardening","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Slip Kinematic Hardening parameters");


  parameter_handler.declare_entry("Twinning enabled","false",dealii::Patterns::Bool(),"Flag to indicate if system twins");
  parameter_handler.declare_entry("Advanced Twinning Model enabled","false",dealii::Patterns::Bool(),"Flag to indicate if Advanced Twinning Model enabled");
  parameter_handler.declare_entry("One twin system Reorientation enabled","false",dealii::Patterns::Bool(),"Flag to indicate One twin system Reorientation is allowed");

  parameter_handler.declare_entry("Number of Twin Systems","-1",dealii::Patterns::Integer(),"Number of Twin Systems");
  parameter_handler.declare_entry("Initial Slip Resistance Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Initial CRSS of the twin sytems");
  parameter_handler.declare_entry("Initial Hardening Modulus Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Hardening moduli of twin systems");
  parameter_handler.declare_entry("Power Law Exponent Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law exponents of twin systems");
  parameter_handler.declare_entry("Saturation Stress Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress of twin systems");
  parameter_handler.declare_entry("Twin Directions File","",dealii::Patterns::Anything(),"Twin Directions File");
  parameter_handler.declare_entry("Twin Normals File","",dealii::Patterns::Anything(),"Twin Normals File");

  parameter_handler.declare_entry("C_1 Twin Kinematic Hardening","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Twin Kinematic Hardening parameters");
  parameter_handler.declare_entry("C_2 Twin Kinematic Hardening","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Twin Kinematic Hardening parameters");

  parameter_handler.declare_entry("Twin Threshold Fraction","-1",dealii::Patterns::Double(),"Threshold fraction of characteristic twin shear (<1)");
  parameter_handler.declare_entry("Twin Saturation Factor","-1",dealii::Patterns::Double(),"Twin growth saturation factor  (<(1-twinThresholdFraction))");
  parameter_handler.declare_entry("Characteristic Twin Shear","-1",dealii::Patterns::Double(),"characteristic twin shear");
  parameter_handler.declare_entry("Critical Value for Twin Visualization","1",dealii::Patterns::Double(),"The integration point with Twin volumes larger than this Critical Value is considered twined during visualization");

  parameter_handler.declare_entry("Stress Tolerance","-1",dealii::Patterns::Double(),"Stress tolerance for the yield surface (MPa)");
  parameter_handler.declare_entry("Max Slip Search Iterations","-1",dealii::Patterns::Integer(),"Maximum no. of active slip search iterations");
  parameter_handler.declare_entry("Max Solver Iterations","-1",dealii::Patterns::Integer(),"Maximum no. of iterations to achieve non-linear convergence");
  parameter_handler.declare_entry("Max Plastic Slip L2 Norm","-1",dealii::Patterns::Double(),"L2-Norm of plastic slip strain-used for load-step adaptivity");

  parameter_handler.declare_entry("Grain ID file name","",dealii::Patterns::Anything(),"Grain ID file name");
  parameter_handler.declare_entry("Voxels in X direction","-1",dealii::Patterns::Integer(),"Number of voxels in x direction");
  parameter_handler.declare_entry("Voxels in Y direction","-1",dealii::Patterns::Integer(),"Number of voxels in y direction");
  parameter_handler.declare_entry("Voxels in Z direction","-1",dealii::Patterns::Integer(),"Number of voxels in z direction");

  parameter_handler.declare_entry("Orientations file name","",dealii::Patterns::Anything(),"Grain orientations file name");
  parameter_handler.declare_entry("Header Lines GrainID File","0",dealii::Patterns::Integer(), "Number of header Lines in grain orientations file");

  parameter_handler.declare_entry("Use Element Deletion","false",dealii::Patterns::Bool(),"Flag to indicate if element deletion is enabled");
  parameter_handler.declare_entry("GrainID of Deleted Elements","-1",dealii::Patterns::Integer(),"The grainID for which the element are deleted");


  //if (enableMultiphase&&(numberofPhases>=2)){
  parameter_handler.declare_entry("Enable User Material Model 2","false",dealii::Patterns::Bool(),"Flag to indicate if User Material Model is enabled Phase 2");
  parameter_handler.declare_entry("Number of User Material Constants 2","0",dealii::Patterns::Integer(),"Number of User Material Constants in a Material model Phase 2");
  parameter_handler.declare_entry("Number of User Material State Variables 2","0",dealii::Patterns::Integer(),"Number of User Material State Variables in a Material model Phase 2");
  parameter_handler.declare_entry("User Material Constants 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Material Constants in a Material model Phase 2");
  parameter_handler.declare_entry("User Material State Variables Initial Values 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Material State Variables in a Material model Phase 2");

  parameter_handler.declare_entry("Elastic Stiffness 2 row 1","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 2 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 2 row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 2 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 2 row 3","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 2 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 2 row 4","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 2 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 2 row 5","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 2 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 2 row 6","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 2 Matrix -Voigt Notation (MPa)");

  parameter_handler.declare_entry("Number of Slip Systems 2","-1",dealii::Patterns::Integer(),"Number of Slip Systems Phase 2");

  parameter_handler.declare_entry("Latent Hardening Ratio filename 2","LatentHardeningRatio2.txt",dealii::Patterns::Anything(),"Latent Hardening Ratio filename 2");

  parameter_handler.declare_entry("Initial Slip Resistance 2","",dealii::Patterns::List(dealii::Patterns::Double()),"RSS of the slip sytems Phase 2");
  parameter_handler.declare_entry("Initial Hardening Modulus 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Heardening moduli of slip systems Phase 2");
  parameter_handler.declare_entry("Power Law Exponent 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law coefficient Phase 2");
  parameter_handler.declare_entry("Saturation Stress 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress Phase 2");
  parameter_handler.declare_entry("Slip Directions File 2","",dealii::Patterns::Anything(),"Slip Directions File Phase 2");
  parameter_handler.declare_entry("Slip Normals File 2","",dealii::Patterns::Anything(),"Slip Normals File Phase 2");
  parameter_handler.declare_entry("Enable Kinematic Hardening 2","false",dealii::Patterns::Bool(),"Flag to indicate if kinematic hardening is enabled Phase 2");
  parameter_handler.declare_entry("C_1 Slip Kinematic Hardening 2","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Slip Kinematic Hardening parameters Phase 2");
  parameter_handler.declare_entry("C_2 Slip Kinematic Hardening 2","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Slip Kinematic Hardening parameters Phase 2");


  parameter_handler.declare_entry("Twinning enabled 2","false",dealii::Patterns::Bool(),"Flag to indicate if system twins Phase 2");
  parameter_handler.declare_entry("Number of Twin Systems 2","-1",dealii::Patterns::Integer(),"Number of Twin Systems Phase 2");
  parameter_handler.declare_entry("Initial Slip Resistance Twin 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Initial CRSS of the twin sytems Phase 2");
  parameter_handler.declare_entry("Initial Hardening Modulus Twin 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Hardening moduli of twin systems Phase 2");
  parameter_handler.declare_entry("Power Law Exponent Twin 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law exponents of twin systems Phase 2");
  parameter_handler.declare_entry("Saturation Stress Twin 2","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress of twin systems Phase 2");
  parameter_handler.declare_entry("Twin Directions File 2","",dealii::Patterns::Anything(),"Twin Directions File Phase 2");
  parameter_handler.declare_entry("Twin Normals File 2","",dealii::Patterns::Anything(),"Twin Normals File Phase 2");

  parameter_handler.declare_entry("C_1 Twin Kinematic Hardening 2","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Twin Kinematic Hardening parameters Phase 2");
  parameter_handler.declare_entry("C_2 Twin Kinematic Hardening 2","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Twin Kinematic Hardening parameters Phase 2");


  parameter_handler.declare_entry("Twin Threshold Fraction 2","-1",dealii::Patterns::Double(),"Threshold fraction of characteristic twin shear (<1) Phase 2");
  parameter_handler.declare_entry("Twin Saturation Factor 2","-1",dealii::Patterns::Double(),"Twin growth saturation factor  (<(1-twinThresholdFraction)) Phase 2");
  parameter_handler.declare_entry("Characteristic Twin Shear 2","-1",dealii::Patterns::Double(),"characteristic twin shear Phase 2");

  //  if (numberofPhases>=3){
  parameter_handler.declare_entry("Enable User Material Model 3","false",dealii::Patterns::Bool(),"Flag to indicate if User Material Model is enabled Phase 3");
  parameter_handler.declare_entry("Number of User Material Constants 3","0",dealii::Patterns::Integer(),"Number of User Material Constants in a Material model Phase 3");
  parameter_handler.declare_entry("Number of User Material State Variables 3","0",dealii::Patterns::Integer(),"Number of User Material State Variables in a Material model Phase 3");
  parameter_handler.declare_entry("User Material Constants 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Material Constants in a Material model Phase 3");
  parameter_handler.declare_entry("User Material State Variables Initial Values 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Material State Variables in a Material model Phase 3");

  parameter_handler.declare_entry("Elastic Stiffness 3 row 1","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 3 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 3 row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 3 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 3 row 3","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 3 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 3 row 4","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 3 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 3 row 5","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 3 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 3 row 6","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 3 Matrix -Voigt Notation (MPa)");

  parameter_handler.declare_entry("Number of Slip Systems 3","-1",dealii::Patterns::Integer(),"Number of Slip Systems Phase 3");
  parameter_handler.declare_entry("Latent Hardening Ratio filename 3","LatentHardeningRatio3.txt",dealii::Patterns::Anything(),"Latent Hardening Ratio filename 3");
  parameter_handler.declare_entry("Initial Slip Resistance 3","",dealii::Patterns::List(dealii::Patterns::Double()),"RSS of the slip sytems Phase 3");
  parameter_handler.declare_entry("Initial Hardening Modulus 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Heardening moduli of slip systems Phase 3");
  parameter_handler.declare_entry("Power Law Exponent 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law coefficient Phase 3");
  parameter_handler.declare_entry("Saturation Stress 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress Phase 3");
  parameter_handler.declare_entry("Slip Directions File 3","",dealii::Patterns::Anything(),"Slip Directions File Phase 3");
  parameter_handler.declare_entry("Slip Normals File 3","",dealii::Patterns::Anything(),"Slip Normals File Phase 3");
  parameter_handler.declare_entry("Enable Kinematic Hardening 3","false",dealii::Patterns::Bool(),"Flag to indicate if kinematic hardening is enabled Phase 3");
  parameter_handler.declare_entry("C_1 Slip Kinematic Hardening 3","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Slip Kinematic Hardening parameters Phase 3");
  parameter_handler.declare_entry("C_2 Slip Kinematic Hardening 3","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Slip Kinematic Hardening parameters Phase 3");


  parameter_handler.declare_entry("Twinning enabled 3","false",dealii::Patterns::Bool(),"Flag to indicate if system twins Phase 3");
  parameter_handler.declare_entry("Number of Twin Systems 3","-1",dealii::Patterns::Integer(),"Number of Twin Systems Phase 3");
  parameter_handler.declare_entry("Initial Slip Resistance Twin 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Initial CRSS of the twin sytems Phase 3");
  parameter_handler.declare_entry("Initial Hardening Modulus Twin 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Hardening moduli of twin systems Phase 3");
  parameter_handler.declare_entry("Power Law Exponent Twin 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law exponents of twin systems Phase 3");
  parameter_handler.declare_entry("Saturation Stress Twin 3","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress of twin systems Phase 3");
  parameter_handler.declare_entry("Twin Directions File 3","",dealii::Patterns::Anything(),"Twin Directions File Phase 3");
  parameter_handler.declare_entry("Twin Normals File 3","",dealii::Patterns::Anything(),"Twin Normals File Phase 3");

  parameter_handler.declare_entry("C_1 Twin Kinematic Hardening 3","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Twin Kinematic Hardening parameters Phase 3");
  parameter_handler.declare_entry("C_2 Twin Kinematic Hardening 3","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Twin Kinematic Hardening parameters Phase 3");


  parameter_handler.declare_entry("Twin Threshold Fraction 3","-1",dealii::Patterns::Double(),"Threshold fraction of characteristic twin shear (<1) Phase 3");
  parameter_handler.declare_entry("Twin Saturation Factor 3","-1",dealii::Patterns::Double(),"Twin growth saturation factor  (<(1-twinThresholdFraction)) Phase 3");
  parameter_handler.declare_entry("Characteristic Twin Shear 3","-1",dealii::Patterns::Double(),"characteristic twin shear Phase 3");

  //  if (numberofPhases>=4){
  parameter_handler.declare_entry("Enable User Material Model 4","false",dealii::Patterns::Bool(),"Flag to indicate if User Material Model is enabled Phase 4");
  parameter_handler.declare_entry("Number of User Material Constants 4","0",dealii::Patterns::Integer(),"Number of User Material Constants in a Material model Phase 4");
  parameter_handler.declare_entry("Number of User Material State Variables 4","0",dealii::Patterns::Integer(),"Number of User Material State Variables in a Material model Phase 4");
  parameter_handler.declare_entry("User Material Constants 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Material Constants in a Material model Phase 4");
  parameter_handler.declare_entry("User Material State Variables Initial Values 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Material State Variables in a Material model Phase 4");

  parameter_handler.declare_entry("Elastic Stiffness 4 row 1","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 4 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 4 row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 4 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 4 row 3","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 4 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 4 row 4","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 4 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 4 row 5","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 4 Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness 4 row 6","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Phase 4 Matrix -Voigt Notation (MPa)");

  parameter_handler.declare_entry("Number of Slip Systems 4","-1",dealii::Patterns::Integer(),"Number of Slip Systems Phase 4");
  parameter_handler.declare_entry("Latent Hardening Ratio filename 4","LatentHardeningRatio4.txt",dealii::Patterns::Anything(),"Latent Hardening Ratio filename 4");
  parameter_handler.declare_entry("Initial Slip Resistance 4","",dealii::Patterns::List(dealii::Patterns::Double()),"RSS of the slip sytems Phase 4");
  parameter_handler.declare_entry("Initial Hardening Modulus 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Heardening moduli of slip systems Phase 4");
  parameter_handler.declare_entry("Power Law Exponent 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law coefficient Phase 4");
  parameter_handler.declare_entry("Saturation Stress 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress Phase 4");
  parameter_handler.declare_entry("Slip Directions File 4","",dealii::Patterns::Anything(),"Slip Directions File Phase 4");
  parameter_handler.declare_entry("Slip Normals File 4","",dealii::Patterns::Anything(),"Slip Normals File Phase 4");
  parameter_handler.declare_entry("Enable Kinematic Hardening 4","false",dealii::Patterns::Bool(),"Flag to indicate if kinematic hardening is enabled Phase 4");
  parameter_handler.declare_entry("C_1 Slip Kinematic Hardening 4","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Slip Kinematic Hardening parameters Phase 4");
  parameter_handler.declare_entry("C_2 Slip Kinematic Hardening 4","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Slip Kinematic Hardening parameters Phase 4");


  parameter_handler.declare_entry("Twinning enabled 4","false",dealii::Patterns::Bool(),"Flag to indicate if system twins Phase 4");
  parameter_handler.declare_entry("Number of Twin Systems 4","-1",dealii::Patterns::Integer(),"Number of Twin Systems Phase 4");
  parameter_handler.declare_entry("Initial Slip Resistance Twin 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Initial CRSS of the twin sytems Phase 4");
  parameter_handler.declare_entry("Initial Hardening Modulus Twin 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Hardening moduli of twin systems Phase 4");
  parameter_handler.declare_entry("Power Law Exponent Twin 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law exponents of twin systems Phase 4");
  parameter_handler.declare_entry("Saturation Stress Twin 4","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress of twin systems Phase 4");
  parameter_handler.declare_entry("Twin Directions File 4","",dealii::Patterns::Anything(),"Twin Directions File Phase 4");
  parameter_handler.declare_entry("Twin Normals File 4","",dealii::Patterns::Anything(),"Twin Normals File Phase 4");

  parameter_handler.declare_entry("C_1 Twin Kinematic Hardening 4","",dealii::Patterns::List(dealii::Patterns::Double()),"C_1 Twin Kinematic Hardening parameters Phase 4");
  parameter_handler.declare_entry("C_2 Twin Kinematic Hardening 4","",dealii::Patterns::List(dealii::Patterns::Double()),"C_2 Twin Kinematic Hardening parameters Phase 4");


  parameter_handler.declare_entry("Twin Threshold Fraction 4","-1",dealii::Patterns::Double(),"Threshold fraction of characteristic twin shear (<1) Phase 4");
  parameter_handler.declare_entry("Twin Saturation Factor 4","-1",dealii::Patterns::Double(),"Twin growth saturation factor  (<(1-twinThresholdFraction)) Phase 4");
  parameter_handler.declare_entry("Characteristic Twin Shear 4","-1",dealii::Patterns::Double(),"characteristic twin shear Phase 4");

  //    }

  //  }
  //  }

}
