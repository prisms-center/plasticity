// Functions for userInputParameters class

#include "../../include/userInputParameters.h"

userInputParameters::userInputParameters(std::string inputfile){

  declare_parameters(parameter_handler);

  parameter_handler.read_input(inputfile);

  dim=parameter_handler.get_double("Number of dimensions");

  feOrder=parameter_handler.get_double("Order of finite elements");
  quadOrder=parameter_handler.get_double("Order of quadrature");

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

  refine_factor = parameter_handler.get_integer("Refine factor");

  degree = parameter_handler.get_integer("Element degree");

  // Adaptive meshing parameters
  h_adaptivity = parameter_handler.get_bool("Mesh adaptivity");
  max_refinement_level = parameter_handler.get_integer("Max refinement level");
  min_refinement_level = parameter_handler.get_integer("Min refinement level");

}

userInputParameters::declare_parameters(parameter_handler){

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
  parameter_handler.declare_entry("Write Output","false",dealii::Patterns::Bool(),"Flag to write output vtu and pvtu files");
  parameter_handler.declare_entry("Output Directory","",Patterns::Anything(),"Output Directory");
  parameter_handler.declare_entry("Skip Output Steps","-1",dealii::Patterns::Integer(),"Skip Output Steps");
  parameter_handler.declare_entry("Output Equivalent strain","false",dealii::Patterns::Bool(),"Output Equivalent strain");
  parameter_handler.declare_entry("Output Equivalent stress","false",dealii::Patterns::Bool(),"Output Equivalent stress");
  parameter_handler.declare_entry("Output Grain ID","false",dealii::Patterns::Bool(),"Output Grain ID");
  parameter_handler.declare_entry("Output Twin fractions","false",dealii::Patterns::Bool(),"Output Twin fractions");

  parameter_handler.declare_entry("Total number increments","-1",dealii::Patterns::Integer(), "No. of increments");

  parameter_handler.declare_entry("Maximum linear solver iterations","-1",dealii::Patterns::Integer(), "Maximum iterations for linear solver");
  parameter_handler.declare_entry("Maximum non linear iterations","-1",dealii::Patterns::Integer(),"Maximum no. of non-linear iterations");
  parameter_handler.declare_entry("Relative linear solver tolerance","-1",dealii::Patterns::Double(),"Relative linear solver tolerance");
  parameter_handler.declare_entry("Absolute nonLinear tolerance","-1",dealii::Patterns::Double(),"Non-linear solver tolerance");
  parameter_handler.declare_entry("Relative nonLinear tolerance","-1",dealii::Patterns::Double(),"Relative non-linear solver tolerance");
  parameter_handler.declare_entry("Stop on convergence failure","false",dealii::Patterns::Bool(),"Flag to stop problem if convergence fails");
  parameter_handler.declare_entry("Enable adaptive Time stepping","false",dealii::Patterns::Bool(),"lag to enable adaptive time steps");
  parameter_handler.declare_entry("Adaptive load step factor","-1",dealii::Patterns::Double(),"Load step factor");
  parameter_handler.declare_entry("Adaptive load increase Factor","-1",dealii::Patterns::Double(),"adaptive Load Increase Factor")
  parameter_handler.declare_entry("Succesive increment for increasing time step","-1",dealii::Patterns::Double(),"Succesive Inc For Increasing Time Step")

  parameter_handler.declare_entry("Elastic Stiffness","",dealii::Patterns::List(dealii::Patterns::Double()),"	Elastic Stiffness Matrix -Voigt Notation (MPa)");

  parameter_handler.declare_entry("Number of Slip Systems","-1",dealii::Patterns::Integer(),"Number of Slip Systems")
  parameter_handler.declare_entry("Latent Hardening Ratio","-1",dealii::Patterns::Double(),"Latent Hardening Ratio");
  parameter_handler.declare_entry("Initial Slip Resistance","",dealii::Patterns::List(dealii::Patterns::Double()),"RSS of the slip sytems");
  parameter_handler.declare_entry("Initial Hardening Modulus","",dealii::Patterns::List(dealii::Patterns::Double()),"ardening moduli of slip systems");
  parameter_handler.declare_entry("Power Law Exponent","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law coefficient");
  parameter_handler.declare_entry("Saturation Stress","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress");
  parameter_handler.declare_entry("Slip Directions File","",Patterns::Anything(),"Slip Directions File");
  parameter_handler.declare_entry("Slip Normals File","",Patterns::Anything(),"Slip Normals File");

  parameter_handler.declare_entry("Number of Twin Systems","-1",dealii::Patterns::Integer(),"Number of Twin Systems")
  parameter_handler.declare_entry("Initial Slip Resistance Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Initial CRSS of the twin sytems");
  parameter_handler.declare_entry("Initial Hardening Modulus Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Hardening moduli of twin systems");
  parameter_handler.declare_entry("Power Law Exponent Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Power law exponents of twin systems");
  parameter_handler.declare_entry("Saturation Stress Twin","",dealii::Patterns::List(dealii::Patterns::Double()),"Saturation stress of twin systems");
  parameter_handler.declare_entry("Twin Directions File","",Patterns::Anything(),"Twin Directions File");
  parameter_handler.declare_entry("Twin Normals File","",Patterns::Anything(),"Twin Normals File");

  parameter_handler.declare_entry("Backstress Factor","-1",dealii::Patterns::Double(),"Ratio between backstress and CRSS during load reversal");
  parameter_handler.declare_entry("Twin Threshold Fraction","-1",dealii::Patterns::Double(),"Threshold fraction of characteristic twin shear (<1)");
  parameter_handler.declare_entry("Twin Saturation Factor","-1",dealii::Patterns::Double(),"Twin growth saturation factor  (<(1-twinThresholdFraction))");
  parameter_handler.declare_entry("Characteristic Twin Shear","-1",dealii::Patterns::Double(),"characteristic twin shear");

  parameter_handler.declare_entry("Stress Tolerance","-1",dealii::Patterns::Double(),"Stress tolerance for the yield surface (MPa)");
  parameter_handler.declare_entry("Max Slip Search Iterations","-1",dealii::Patterns::Integer(),"Maximum no. of active slip search iterations");
  parameter_handler.declare_entry("Max Solver Iterations","-1",dealii::Patterns::Integer(),"Maximum no. of iterations to achieve non-linear convergence");
  parameter_handler.declare_entry("Max Plastic Slip L2 Norm","-1",dealii::Patterns::Double(),"L2-Norm of plastic slip strain-used for load-step adaptivity");

  parameter_handler.declare_entry("Grain ID file name","",Patterns::Anything(),"Grain ID file name");
  parameter_handler.declare_entry("Voxels in X","-1",dealii::Patterns::Integer(),"Number of voxels in x direction");
  parameter_handler.declare_entry("Voxels in Y","-1",dealii::Patterns::Integer(),"Number of voxels in y direction");
  parameter_handler.declare_entry("Voxels in Z","-1",dealii::Patterns::Integer(),"Number of voxels in z direction");

  parameter_handler.declare_entry("Orientations file name","",Patterns::Anything(),"Grain orientations file name");
  parameter_handler.declare_entry("Header Lines GrainID File","0",dealii::Patterns::Integer(), "Number of header Lines in grain orientations file");

}
