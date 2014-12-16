//initialization method for ellipticBVP class

#ifndef INIT_ELLIPTICBVP_H
#define INIT_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

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
  
  //constraints
  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints);
  constraints.close ();

  //initialize global data structures
  solution.reinit (locally_owned_dofs, mpi_communicator); solution=0;
  residual.reinit (locally_owned_dofs, mpi_communicator); residual=0;
  solutionWithGhosts.reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator);

  CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dofHandler, csp, constraints, false);
  SparsityTools::distribute_sparsity_pattern (csp,
					      dofHandler.n_locally_owned_dofs_per_processor(),
					      mpi_communicator,
					      locally_relevant_dofs);
  jacobian.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator);

  //mark boundaries
  markBoundaries();
  //apply initial conditions
  applyInitialConditions();
  solutionWithGhosts=solution;
}

#endif
