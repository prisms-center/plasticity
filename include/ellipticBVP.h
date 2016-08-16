//base class for elliptic boundary value problem implementation
#ifndef ELLIPTICBVP_H
#define ELLIPTICBVP_H

//general headers
#include <fstream>
#include <sstream>

//dealii headers
#include "dealIIheaders.h"

//compiler directives to handle warnings
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wextra"
#pragma GCC diagnostic pop

using namespace dealii;

//define data types  
typedef PETScWrappers::MPI::Vector vectorType;
typedef PETScWrappers::MPI::SparseMatrix matrixType;
//LA::MPI::SparseMatrix
//LA::MPI::Vector

//
//base class for elliptic PDE's
//
template <int dim>
class ellipticBVP
{
 public:
  ellipticBVP(); 
  ~ellipticBVP(); 
  void run   ();

 protected:
  //parallel objects
  MPI_Comm   mpi_communicator;
  IndexSet   locally_owned_dofs;
  IndexSet   locally_owned_dofs_Scalar;
  IndexSet   locally_relevant_dofs;
  IndexSet   locally_relevant_dofs_Scalar;
  
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
  void initProject();
  void project();

  //virtual methods to be implemented in derived class
  //method to calculate elemental Jacobian and Residual,
  //which should be implemented in the derived material model class
  virtual void getElementalValues(FEValues<dim>& fe_values,
				  unsigned int dofs_per_cell,
				  unsigned int num_quad_points,
				  FullMatrix<double>& elementalJacobian,
				  Vector<double>&     elementalResidual) = 0;
  
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
};

//other ellipticBVP headers 
//(these are source files, which will are temporarily treated as
//header files till library packaging scheme is finalized)
#include "../src/ellipticBVP/ellipticBVP.cc"
#include "../src/ellipticBVP/run.cc"
#include "../src/ellipticBVP/mesh.cc"
#include "../src/ellipticBVP/init.cc"
//#include "../src/ellipticBVP/markBoundaries.cc"
#include "../src/ellipticBVP/initialConditions.cc"
#include "../src/ellipticBVP/boundaryConditions.cc"
#include "../src/ellipticBVP/assemble.cc"
#include "../src/ellipticBVP/solve.cc"
#include "../src/ellipticBVP/solveNonLinearSystem.cc"
#include "../src/ellipticBVP/solveLinearSystem.cc"
#include "../src/ellipticBVP/iterationUpdates.cc"
#include "../src/ellipticBVP/incrementUpdates.cc"
#include "../src/ellipticBVP/output.cc"
#include "../src/ellipticBVP/project.cc"

#endif
