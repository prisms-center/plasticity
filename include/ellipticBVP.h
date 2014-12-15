//base class for elliptic boundary value problem implementation
#ifndef ELLIPTICBVP_H
#define ELLIPTICBVP_H

//general headers
#include <fstream>
#include <sstream>

//dealii headers
#include "dealIIheaders.h"

using namespace dealii;

//define data types  
typedef PETScWrappers::MPI::Vector vectorType;
typedef PETScWrappers::MPI::SparseMatrix matrixType;
//LA::MPI::SparseMatrix
//LA::MPI::Vector

#define feOrder   1
#define quadOrder 2 
#define meshRefineFactor 2
#define writeOutput true
#define solverType PETScWrappers::SolverCG
#define maxSolverIterations 5000
#define relSolverTolerance  1.0e-8

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
  IndexSet   locally_relevant_dofs;
  
  //FE data structres
  parallel::distributed::Triangulation<dim> triangulation;
  FESystem<dim>      FE;
  ConstraintMatrix   constraints;
  DoFHandler<dim>    dofHandler;
    
  //methods
  void mesh();
  void init();
  void assemble();
  void solveIncrement();
  void solve();
  void output();
 
  //virtual methods to be implemented in derived class
  //method to calculate elemental Jacobian and Residual
  virtual void getElementalValues(FullMatrix<double>& elementalJacobian,
				  Vector<double>&     elementalResidual) = 0;

  //methods to apply dirichlet BC's and initial conditions
  virtual void markBoundaries();
  virtual void applyDirichletBCs();
  virtual void applyInitialConditions();
  
  //parallel data structures
  vectorType solution, solutionLocal, residual;
  matrixType jacobian;

  //misc variables
  unsigned int currentIteration, currentIncrement;
  unsigned int totalIncrements;

  //parallel message stream
  ConditionalOStream  pcout;  
  
  //compute-time logger
  TimerOutput computing_timer;

  //output variables
  //solution name array                                                                                                                                                
  std::vector<std::string> nodal_solution_names;
  std::vector<DataComponentInterpretation::DataComponentInterpretation> nodal_data_component_interpretation;
};

//other ellipticBVP headers 
//(these are source files, which will are temporarily treated as
//header files till library packaging scheme is finalized)
#include "../src/ellipticBVP/ellipticBVP.cc"
#include "../src/ellipticBVP/run.cc"
#include "../src/ellipticBVP/mesh.cc"
#include "../src/ellipticBVP/init.cc"
#include "../src/ellipticBVP/markBoundaries.cc"
#include "../src/ellipticBVP/initialConditions.cc"
#include "../src/ellipticBVP/boundaryConditions.cc"
#include "../src/ellipticBVP/assemble.cc"
#include "../src/ellipticBVP/solveIncrement.cc"
#include "../src/ellipticBVP/solve.cc"
#include "../src/ellipticBVP/output.cc"

#endif
