//base class for elliptic boundary value problem implementation
#ifndef ELLIPTICBVP_H
#define ELLIPTICBVP_H

//dealii headers
#include "dealIIheaders.h"
#include "userInputParameters.h"

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

 protected:

  //parallel objects
  MPI_Comm   mpi_communicator;
  IndexSet   locally_owned_dofs;
  IndexSet   locally_owned_dofs_Scalar;
  IndexSet   locally_relevant_dofs;
  IndexSet   locally_relevant_dofs_Scalar;

  //User input parameters object
  userInputParameters userInputs;

  //FE data structres
  parallel::distributed::Triangulation<dim> triangulation;
  FESystem<dim>      FE;
  FESystem<dim>      FE_Scalar;
  ConstraintMatrix   constraints;
  ConstraintMatrix   constraintsMassMatrix;
  DoFHandler<dim>    dofHandler;
  DoFHandler<dim>    dofHandler_Scalar;

  //methods
  virtual void mesh();
  void init();
  void assemble();
  void solveLinearSystem(ConstraintMatrix& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts);
  void solveLinearSystem2(ConstraintMatrix& constraintmatrix, matrixType& A, vectorType& b, vectorType& x, vectorType& xGhosts, vectorType& dxGhosts);
  bool solveNonLinearSystem();
  void solve();
  void output();
  void initProjection();
  void projection();
  void markBoundaries();

  //virtual methods to be implemented in derived class
  //method to calculate elemental Jacobian and Residual,
  //which should be implemented in the derived material model class
#ifdef enableUserModel
  virtual void getQuadratureValues(unsigned int elementID,
				   unsigned int numElemDofs,
				   unsigned int* componentIndices,
				   double* shapeValues,
				   double* shapeGrads,
				   double* F,
				   double* residual,
				   double* jacobian,
				   double* history) = 0;
#else
  virtual void getElementalValues(FEValues<dim>& fe_values,
				  unsigned int dofs_per_cell,
				  unsigned int num_quad_points,
				  FullMatrix<double>& elementalJacobian,
				  Vector<double>&     elementalResidual) = 0;
#endif
  //methods to allow for pre/post iteration updates
  virtual void updateBeforeIteration();
  virtual void updateAfterIteration();
  virtual bool testConvergenceAfterIteration();
  //methods to allow for pre/post increment updates
  virtual void updateBeforeIncrement();
  virtual void updateAfterIncrement();

  //methods to apply dirichlet BC's and initial conditions
  void applyDirichletBCs();
  void applyInitialConditions();
  virtual void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  std::map<types::global_dof_index,double> boundary_values;
  std::map<types::global_dof_index, Point<dim> > supportPoints;

  //parallel data structures
  vectorType solution, oldSolution, residual;
  vectorType solutionWithGhosts, solutionIncWithGhosts;
  matrixType jacobian;

  //misc variables
  unsigned int currentIteration, currentIncrement;
  unsigned int totalIncrements;
  bool resetIncrement;
  double loadFactorSetByModel;
  double totalLoadFactor;

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
  //postprocessed scalar variable name array (only scalar variables supported currently, will be extended later to vectors and tensors, if required.)
  std::vector<std::string> postprocessed_solution_names;
  //postprocessing data structures
  std::vector<vectorType*> postFields, postFieldsWithGhosts, postResidual;
  matrixType massMatrix;
  Table<4,double> postprocessValues;

  //user model related variables and methods
#ifdef enableUserModel
  unsigned int numQuadHistoryVariables;
  Table<3, double> quadHistory;
  virtual void initQuadHistory();
#endif
};

#endif
