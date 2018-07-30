// Functions for userInputParameters class

#include "../../include/userInputParameters.h"

userInputParameters::userInputParameters(std::string inputfile, dealii::ParameterHandler & parameter_handler){

  declare_parameters(parameter_handler);

  #if (DEAL_II_VERSION_MAJOR < 9 && DEAL_II_VERSION_MINOR < 5)
  parameter_handler.read_input(inputfile);
  #else
  parameter_handler.parse_input(inputfile);
  #endif

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

  //output parameters
  writeOutput = parameter_handler.get_bool("Write Output");
  outputDirectory = parameter_handler.get("Output Directory");
  skipOutputSteps=parameter_handler.get_integer("Skip Output Steps");

  if(skipOutputSteps<=0)
      skipOutputSteps=1;

  output_Eqv_strain = parameter_handler.get_bool("Output Equivalent strain");
  output_Eqv_stress = parameter_handler.get_bool("Output Equivalent stress");
  output_Grain_ID = parameter_handler.get_bool("Output Grain ID");
  output_Twin = parameter_handler.get_bool("Output Twin fractions");

  totalNumIncrements=parameter_handler.get_integer("Total number of increments");
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

  crystalStructure = parameter_handler.get("Crystal Structure");

  //elasticStiffness.reinit(6,6);
  elasticStiffness.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 1"))));
  elasticStiffness.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 2"))));
  elasticStiffness.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 3"))));
  elasticStiffness.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 4"))));
  elasticStiffness.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 5"))));
  elasticStiffness.push_back(dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Elastic Stiffness row 6"))));

  numSlipSystems=parameter_handler.get_integer("Number of Slip Systems");
  latentHardeningRatio=parameter_handler.get_double("Latent Hardening Ratio");
  initialSlipResistance = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance")));
  initialHardeningModulus = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus")));
  powerLawExponent = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent")));
  saturationStress = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress")));
  slipDirectionsFile = parameter_handler.get("Slip Directions File");
  slipNormalsFile = parameter_handler.get("Slip Normals File");

  enableTwinning = parameter_handler.get_bool("Twinning enabled");
  if(enableTwinning){
    numTwinSystems=parameter_handler.get_integer("Number of Twin Systems");
    initialSlipResistanceTwin = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Slip Resistance Twin")));
    initialHardeningModulusTwin = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Initial Hardening Modulus Twin")));
    powerLawExponentTwin = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Power Law Exponent Twin")));
    saturationStressTwin = dealii::Utilities::string_to_double(dealii::Utilities::split_string_list(parameter_handler.get("Saturation Stress Twin")));
    twinDirectionsFile = parameter_handler.get("Twin Directions File");
    twinNormalsFile = parameter_handler.get("Twin Normals File");}
  else{
    std::cout<<"Twinning is not enabled \n";
    numTwinSystems=1;
    initialSlipResistanceTwin.push_back(10e5);
    initialHardeningModulusTwin.push_back(1);
    powerLawExponentTwin.push_back(0);
    saturationStressTwin.push_back(10e10);
  }

  backstressFactor=parameter_handler.get_double("Backstress Factor");
  twinThresholdFraction=parameter_handler.get_double("Twin Threshold Fraction");
  twinSaturationFactor=parameter_handler.get_double("Twin Saturation Factor");
  twinShear=parameter_handler.get_double("Characteristic Twin Shear");

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

}

void userInputParameters::declare_parameters(dealii::ParameterHandler & parameter_handler){

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

  parameter_handler.declare_entry("Write Output","false",dealii::Patterns::Bool(),"Flag to write output vtu and pvtu files");
  parameter_handler.declare_entry("Output Directory",".",dealii::Patterns::Anything(),"Output Directory");
  parameter_handler.declare_entry("Skip Output Steps","-1",dealii::Patterns::Integer(),"Skip Output Steps");
  parameter_handler.declare_entry("Output Equivalent strain","false",dealii::Patterns::Bool(),"Output Equivalent strain");
  parameter_handler.declare_entry("Output Equivalent stress","false",dealii::Patterns::Bool(),"Output Equivalent stress");
  parameter_handler.declare_entry("Output Grain ID","false",dealii::Patterns::Bool(),"Output Grain ID");
  parameter_handler.declare_entry("Output Twin fractions","false",dealii::Patterns::Bool(),"Output Twin fractions");

  parameter_handler.declare_entry("Total number of increments","-1",dealii::Patterns::Integer(), "No. of increments");

  parameter_handler.declare_entry("Maximum linear solver iterations","-1",dealii::Patterns::Integer(), "Maximum iterations for linear solver");
  parameter_handler.declare_entry("Maximum non linear iterations","-1",dealii::Patterns::Integer(),"Maximum no. of non-linear iterations");
  parameter_handler.declare_entry("Relative linear solver tolerance","-1",dealii::Patterns::Double(),"Relative linear solver tolerance");
  parameter_handler.declare_entry("Absolute nonLinear solver tolerance","-1",dealii::Patterns::Double(),"Non-linear solver tolerance");
  parameter_handler.declare_entry("Relative nonLinear solver tolerance","-1",dealii::Patterns::Double(),"Relative non-linear solver tolerance");
  parameter_handler.declare_entry("Stop on convergence failure","false",dealii::Patterns::Bool(),"Flag to stop problem if convergence fails");
  parameter_handler.declare_entry("Enable adaptive Time stepping","false",dealii::Patterns::Bool(),"lag to enable adaptive time steps");
  parameter_handler.declare_entry("Adaptive load step factor","-1",dealii::Patterns::Double(),"Load step factor");
  parameter_handler.declare_entry("Adaptive load increase Factor","-1",dealii::Patterns::Double(),"adaptive Load Increase Factor");
  parameter_handler.declare_entry("Succesive increment for increasing time step","-1",dealii::Patterns::Double(),"Succesive Inc For Increasing Time Step");

  parameter_handler.declare_entry("Crystal Structure","",dealii::Patterns::Anything(),"Crystal structure of problem");

  parameter_handler.declare_entry("Elastic Stiffness row 1","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 2","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 3","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 4","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 5","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");
  parameter_handler.declare_entry("Elastic Stiffness row 6","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");

  parameter_handler.declare_entry("Number of Slip Systems","-1",dealii::Patterns::Integer(),"Number of Slip Systems");
  parameter_handler.declare_entry("Latent Hardening Ratio","-1",dealii::Patterns::Double(),"Latent Hardening Ratio");
  parameter_handler.declare_entry("Initial Slip Resistance","",dealii::Patterns::List(dealii::Patterns::Double()),"RSS of the slip sytems");
  parameter_handler.declare_entry("Initial Hardening Modulus","",dealii::Patterns::List(dealii::Patterns::Double()),"Heardening moduli of slip systems");
  parameter_handler.declare_entry("Power Law Exponent","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law coefficient");
  parameter_handler.declare_entry("Saturation Stress","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress");
  parameter_handler.declare_entry("Slip Directions File","",dealii::Patterns::Anything(),"Slip Directions File");
  parameter_handler.declare_entry("Slip Normals File","",dealii::Patterns::Anything(),"Slip Normals File");

  parameter_handler.declare_entry("Twinning enabled","false",dealii::Patterns::Bool(),"Flag to indicate if system twins");
  parameter_handler.declare_entry("Number of Twin Systems","-1",dealii::Patterns::Integer(),"Number of Twin Systems");
  parameter_handler.declare_entry("Initial Slip Resistance Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Initial CRSS of the twin sytems");
  parameter_handler.declare_entry("Initial Hardening Modulus Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Hardening moduli of twin systems");
  parameter_handler.declare_entry("Power Law Exponent Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law exponents of twin systems");
  parameter_handler.declare_entry("Saturation Stress Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress of twin systems");
  parameter_handler.declare_entry("Twin Directions File","",dealii::Patterns::Anything(),"Twin Directions File");
  parameter_handler.declare_entry("Twin Normals File","",dealii::Patterns::Anything(),"Twin Normals File");

  parameter_handler.declare_entry("Backstress Factor","-1",dealii::Patterns::Double(),"Ratio between backstress and CRSS during load reversal");
  parameter_handler.declare_entry("Twin Threshold Fraction","-1",dealii::Patterns::Double(),"Threshold fraction of characteristic twin shear (<1)");
  parameter_handler.declare_entry("Twin Saturation Factor","-1",dealii::Patterns::Double(),"Twin growth saturation factor  (<(1-twinThresholdFraction))");
  parameter_handler.declare_entry("Characteristic Twin Shear","-1",dealii::Patterns::Double(),"characteristic twin shear");

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

}
