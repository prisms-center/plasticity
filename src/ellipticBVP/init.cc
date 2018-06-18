//initialization method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//initialize all FE objects and data structures
template <int dim>
void ellipticBVP<dim>::init(){
  pcout << "number of MPI processes: "
	<< Utilities::MPI::n_mpi_processes(mpi_communicator)
	<< std::endl;

  //initialize FE objects
  dofHandler.distribute_dofs (FE);
  locally_owned_dofs = dofHandler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dofHandler, locally_relevant_dofs);
  pcout << "number of elements: "
	<< triangulation.n_global_active_cells()
	<< std::endl
	<< "number of degrees of freedom: "
	<< dofHandler.n_dofs()
	<< std::endl;

  //initialize FE objects for scalar field which will be used for post processing
  dofHandler_Scalar.distribute_dofs (FE_Scalar);
  locally_owned_dofs_Scalar = dofHandler_Scalar.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dofHandler_Scalar, locally_relevant_dofs_Scalar);

  //constraints
  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints);
  constraints.close ();

  //constraints for mass matrix used in post-processing for projection
  //(i.e. no constraints, as mass matrix needs no dirichlet BCs)
  constraintsMassMatrix.clear ();
  constraintsMassMatrix.reinit (locally_relevant_dofs_Scalar);
  DoFTools::make_hanging_node_constraints (dofHandler_Scalar, constraintsMassMatrix);
  constraintsMassMatrix.close ();

  //get support points (nodes) for this problem
  DoFTools::map_dofs_to_support_points(MappingQ1<dim, dim>(), dofHandler, supportPoints);

  //initialize global data structures
  solution.reinit (locally_owned_dofs, mpi_communicator); solution=0;
  oldSolution.reinit (locally_owned_dofs, mpi_communicator); oldSolution=0;
  solutionWithGhosts.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  solutionIncWithGhosts.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);
  residual.reinit (locally_owned_dofs, mpi_communicator); residual=0;

  DynamicSparsityPattern dsp (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dofHandler, dsp, constraints, false);
  SparsityTools::distribute_sparsity_pattern (dsp,
					      dofHandler.n_locally_owned_dofs_per_processor(),
					      mpi_communicator,
					      locally_relevant_dofs);
  jacobian.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);

  //apply initial conditions
  applyInitialConditions();
  solutionWithGhosts=solution;
  oldSolution=solution;
}
