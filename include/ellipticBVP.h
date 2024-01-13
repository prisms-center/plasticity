//base class for elliptic boundary value problem implementation
#ifndef ELLIPTICBVP_H
#define ELLIPTICBVP_H

//dealii headers
#include "dealIIheaders.h"
#include "userInputParameters.h"
#include "crystalOrientationsIO.h"

using namespace dealii;

//compiler directives to handle warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic pop

//define data types
typedef PETScWrappers::MPI::Vector vectorType;
typedef PETScWrappers::MPI::SparseMatrix matrixType;
//LA::MPI::SparseMatrix
//LA::MPI::Vector

//
//base class for elliptic PDE's
//
template <int dim>
class ellipticBVP : public Subscriptor
{
public:
  ellipticBVP(userInputParameters _userInputs);
  ~ellipticBVP();
  void run   ();
  crystalOrientationsIO<dim> orientations_Mesh;

protected:

  //parallel objects
  MPI_Comm   mpi_communicator;
  IndexSet   locally_owned_dofs;
  IndexSet   locally_owned_dofs_Scalar;
  IndexSet   locally_relevant_dofs,locally_relevant_dofs_Mod, locally_relevant_ghost_dofs;
  IndexSet   locally_relevant_dofs_Scalar;
  IndexSet   dof_FN,dof_FP;
  IndexSet   dof_FXN,dof_FXP;
  IndexSet   dof_FYN,dof_FYP;
  IndexSet   dof_FZN,dof_FZP;
  IndexSet   dof_Boundary_Layer2,vertices_DOFs,global_dof_Boundary_Layer2,Edges_DOFs_1,Edges_DOFs_2,Edges_DOFs_3;
//INDENTATION
  IndexSet      active_set, debug_set, frozen_set;
  unsigned int active_set_size, old_active_set_size, freeze_out_iterations;
  bool  active_set_size_changed;



  unsigned int globalDOF_V000_1,globalDOF_V000_2,globalDOF_V000_3,globalDOF_V001_1,globalDOF_V001_2,globalDOF_V001_3,globalDOF_V010_1,globalDOF_V010_2,globalDOF_V010_3,globalDOF_V100_1,globalDOF_V100_2,globalDOF_V100_3;
  unsigned int globalDOF_V011_1,globalDOF_V011_2,globalDOF_V011_3,globalDOF_V110_1,globalDOF_V110_2,globalDOF_V110_3,globalDOF_V101_1,globalDOF_V101_2,globalDOF_V101_3,globalDOF_V111_1,globalDOF_V111_2,globalDOF_V111_3;

  unsigned int local_globalDOF_V000_1,local_globalDOF_V000_2,local_globalDOF_V000_3,local_globalDOF_V001_1,local_globalDOF_V001_2,local_globalDOF_V001_3,local_globalDOF_V010_1,local_globalDOF_V010_2,local_globalDOF_V010_3,local_globalDOF_V100_1,local_globalDOF_V100_2,local_globalDOF_V100_3;
  unsigned int local_globalDOF_V011_1,local_globalDOF_V011_2,local_globalDOF_V011_3,local_globalDOF_V110_1,local_globalDOF_V110_2,local_globalDOF_V110_3,local_globalDOF_V101_1,local_globalDOF_V101_2,local_globalDOF_V101_3,local_globalDOF_V111_1,local_globalDOF_V111_2,local_globalDOF_V111_3;

  unsigned int global_size_dof_Boundary_Layer2,global_size_dof_Edge;

  std::vector<unsigned int> global_vector_dof_Boundary_Layer2,dofNodalDisplacement;
  std::vector<double> deluNodalDisplacement,initPosIndenter,finalPosIndenter, prevPosIndenter, currentPosIndenter;
  std::vector<Point<dim>> KeyPosIndenter;
  double indenterSize, indenterTolerance, indenterLoad;// depthRefinementMultiple;
  unsigned int indenterShape, indenterFace, loadFace, indentDof, indentationKeyFrames;//, refinementFactor;
  bool roughIndenter;
  std::vector<std::vector<unsigned int> >  vertices_Constraints_Matrix,edges_Constraints_Matrix,faces_Constraints_Matrix, global_Edges_DOFs_Vector_Array;
  std::vector<std::vector<double>>  vertices_Constraints_Coef,edges_Constraints_Coef,faces_Constraints_Coef, global_Edges_DOFs_Coord_Vector_Array,nodalDisplacement;
  std::vector<std::vector<int>> periodicBCsInput; // 	Periodic BCs Input
  std::vector<std::vector<double>> periodicBCsInput2,periodicBCsInput2_Orig; // 	Periodic BCs Input
  std::vector<unsigned int> vertices_Constraint_Known, vertices_DOFs_vector;
  std::vector<std::vector<unsigned int>> edges_DOFs;
  std::vector<IndexSet> faces_dof_Index_vector;

  unsigned int numberVerticesConstraint,numberEdgesConstraint,numberFacesConstraint,totalNumVerticesDOFs,totalNumEdgesDOFs,totalNumEachEdgesNodes,totalNumFacesDOFs;

  Point<dim> vertexNode,node, cellCenter, nodeU, nodedU, nodeU2;
  std::vector<double> displace_local;

  //User input parameters object
  userInputParameters userInputs;

  //FE data structres
  parallel::distributed::Triangulation<dim> triangulation,triangulation2;
  FESystem<dim>      FE;
  FESystem<dim>      FE_Scalar;
  DoFHandler<dim>    dofHandler;
  DoFHandler<dim>    dofHandler_Scalar;
  std::vector<unsigned int> cellOrientationMap_Mesh;

  //methods
  virtual void mesh();
  void init();
  void assemble();
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
  ConstraintMatrix   constraints, constraints_PBCs_Inc0, constraints_PBCs_IncNot0, constraints_PBCs_Inc0Neg;
  ConstraintMatrix   constraintsMassMatrix;
  void solveLinearSystem(ConstraintMatrix& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts);
  void solveLinearSystem2(ConstraintMatrix& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts);
  ///Periodic BCs Implementation
  void setFaceConstraints(ConstraintMatrix& constraintmatrix);
  void setEdgeConstraints(ConstraintMatrix& constraintmatrix);
  void setNodeConstraints(ConstraintMatrix& constraintmatrix);
  #else
  AffineConstraints<double>   constraints, constraints_PBCs_Inc0, constraints_PBCs_IncNot0, constraints_PBCs_Inc0Neg;
  AffineConstraints<double>   constraintsMassMatrix;
  //INDENTATION
  AffineConstraints<double>   indentation_constraints, hanging_constraints;
  void solveLinearSystem(AffineConstraints<double>& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts);
  void solveLinearSystem2(AffineConstraints<double>& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts);
  ///Periodic BCs Implementation
  void setNodeConstraints(AffineConstraints<double>& constraintmatrix);
  void setEdgeConstraints(AffineConstraints<double>& constraintmatrix);
  void setFaceConstraints(AffineConstraints<double>& constraintmatrix);
  #endif

  bool solveNonLinearSystem();
  void solve();
  void output();
  void initProjection();
  void projection();
  void markBoundaries();

  //virtual methods to be implemented in derived class
  //method to calculate elemental Jacobian and Residual,
  //which should be implemented in the derived material model class

  virtual void getElementalValues(FEValues<dim>& fe_values, unsigned int dofs_per_cell, unsigned int num_quad_points, FullMatrix<double>& elementalJacobian, Vector<double>&elementalResidual) = 0;

  //methods to allow for pre/post iteration updates
  virtual void updateBeforeIteration();
  virtual void updateAfterIteration();
  virtual bool testConvergenceAfterIteration();
  //methods to allow for pre/post increment updates
  virtual void updateBeforeIncrement();
  virtual void updateAfterIncrement();
  void updateAfterIncrementBase();
  //methods to apply dirichlet BC's and initial conditions
  void applyDirichletBCs();
  void applyInitialConditions();
  void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  ///////These functions are for Periodic BCs Implementation
  void setPeriodicity();
  void setPeriodicityConstraintsInit();
  void setPeriodicityConstraints();
  void setPeriodicityConstraintsInc0();
  void setPeriodicityConstraintsIncNot0();
  void setPeriodicityConstraintsInc0Neg();
///////These functions are for Indentation BCs Implementation
  void updateIndentPos();
  void setIndentation(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  void setIndentation2(const Point<dim>& node, const unsigned int dof, bool& flag, double& value,
                       double& criterion);
  void setIndentationConstraints();
  void displaceNode(const Point<dim> & p, const Point<dim> & u);
  void meshRefineIndentation();
  bool flagActiveSet(const Point<dim> & p);
  bool flagActiveSetLambda(const Point<dim> & p, double& criterion);
  bool flagActiveSetLambda2(const Point<dim> & p, double& criterion);
  //void updateActiveSet();
  void measureIndentationLoad();
  void setActiveSet();
  void setActiveSet2();
  void setFrozenSet();
  void assemble_mass_matrix_diagonal();
  double Obstacle(const Point<dim> & p, const unsigned int & component, const std::vector<double> & ind);
  ///////These functions are for DIC BCs evaluation
  void bcFunction1(double _yval, double &value_x, double &value_y, double _currentIncr);
  void bcFunction2(double _yval, double &value_x, double &value_y, double _currentIncr);
  void bcFunction3(double _yval, double &value_x, double &value_y, double _currentIncr);
  void bcFunction4(double _yval, double &value_x, double &value_y, double _currentIncr);
  std::map<types::global_dof_index, Point<dim> > supportPoints;
  //parallel data structures
  vectorType solution, oldSolution, residual;
  vectorType solutionWithGhosts, solutionIncWithGhosts;
  //INDENTATION
  vectorType  newton_rhs_uncondensed, newton_rhs_uncondensed_inc, diag_mass_matrix_vector;
  vectorType  lambda;
  matrixType jacobian;
  // Boundary condition variables
  std::vector<std::vector<bool>> faceDOFConstrained;
  std::vector<std::vector<double>> deluConstraint;
  FullMatrix<double> tabularDisplacements,tabularInputNeumannBCs;
  unsigned int timeCounter,timeCounter_Neumann;
  double currentTime;
  /////DIC bc names
  FullMatrix<double> bc_new1,bc_new2,bc_new3,bc_new4;
  FullMatrix<double> Fprev=IdentityMatrix(dim);
  FullMatrix<double> F,deltaF;
  FullMatrix<double> targetVelGrad;
  //misc variables
  double delT,totalT,cycleTime;
  unsigned int currentIteration, currentIncrement;
  unsigned int totalIncrements,periodicTotalIncrements;
  bool resetIncrement;
  double loadFactorSetByModel;
  double totalLoadFactor;
  //INDENTATION
  double currentIndentDisp;
  //parallel message stream
  ConditionalOStream  pcout;
  //compute-time logger
  TimerOutput computing_timer;
  //output variables
  //solution name array
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
  //post processing
  unsigned int numPostProcessedFields;
  unsigned int numPostProcessedFieldsAtCellCenters;
  //postprocessed scalar variable name array (only scalar variables supported currently, will be extended later to vectors and tensors, if required.)
  std::vector<std::string> postprocessed_solution_names, postprocessedFieldsAtCellCenters_solution_names;
  //postprocessing data structures
  std::vector<vectorType*> postFields, postFieldsWithGhosts, postResidual;
  matrixType massMatrix;
  Table<4,double> postprocessValues;
  Table<2,double> postprocessValuesAtCellCenters;
  //user model related variables and methods
  #ifdef enableUserModel
  unsigned int numQuadHistoryVariables;
  Table<3, double> quadHistory;
  virtual void initQuadHistory();
  #endif
};

    #endif
