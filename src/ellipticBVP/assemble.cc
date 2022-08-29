//assemble method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//FE assemble operation
template <int dim>
void ellipticBVP<dim>::assemble(){

  //local variables
  QGauss<dim>  quadrature(userInputs.quadOrder);
  FEValues<dim> fe_values (FE, quadrature, update_values | update_gradients | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  const unsigned int   num_quad_points = quadrature.size();
  FullMatrix<double>   elementalJacobian (dofs_per_cell, dofs_per_cell);
  Vector<double>       elementalResidual (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //apply Dirichlet BC's
  applyDirichletBCs();

  if(userInputs.enablePeriodicBCs){
    DynamicSparsityPattern dsp (locally_relevant_dofs_Mod);
    DoFTools::make_sparsity_pattern (dofHandler, dsp, constraints, false);
    #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
    SparsityTools::distribute_sparsity_pattern (dsp,
      dofHandler.n_locally_owned_dofs_per_processor(),
      mpi_communicator,
      locally_relevant_dofs_Mod);
      jacobian.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    #else
    SparsityTools::distribute_sparsity_pattern (dsp,
      dofHandler.locally_owned_dofs(),
      mpi_communicator,
      locally_relevant_dofs_Mod);
      jacobian.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    #endif
    }

    if (userInputs.enableNeumannBCs){
      double currentTime;
      currentTime=delT*(currentIncrement+1);
			if (currentIncrement==0){
				timeCounter_Neumann=1;
			}
			if (currentTime>userInputs.tabularTimeNeumannBCs[timeCounter_Neumann]){
				timeCounter_Neumann=timeCounter_Neumann+1;
			}
    }
  //initialize global data structures to zero
  //The additional compress operations are only to flush out data and
  //switch to the correct write state. For  details look at the documentation
  //for PETScWrappers::MPI::Vector()
  jacobian.compress(VectorOperation::add);
  jacobian=0.0;
  residual.compress(VectorOperation::add);
  residual=0.0;
  if (userInputs.enableIndentationBCs){
      newton_rhs_uncondensed_inc.compress(VectorOperation::add);
      newton_rhs_uncondensed_inc = 0.0;
      newton_rhs_uncondensed.compress(VectorOperation::add);
      newton_rhs_uncondensed = 0.0;
  }
  try{
    //parallel loop over all elements
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    unsigned int cellID=0;
    for (; cell!=endc; ++cell) {
      if (cell->is_locally_owned()){
        elementalJacobian = 0;
        elementalResidual = 0;
        cell->set_user_index(cellID);

        //Compute values for the current element
        fe_values.reinit (cell);
        cell->get_dof_indices (local_dof_indices);
        //get elemental jacobian and residual
        getElementalValues(fe_values, dofs_per_cell, num_quad_points, elementalJacobian, elementalResidual);
        //
        constraints.distribute_local_to_global(elementalJacobian,
          elementalResidual,
          local_dof_indices,
          jacobian,
          residual);

        if (userInputs.enableIndentationBCs)
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                newton_rhs_uncondensed_inc(local_dof_indices[i]) += elementalResidual(i);

          cellID++;
        }
      }
    }
    catch (int param){
      std::cout << "skipping assembly and nonlinear solve as resetIncrement==True\n";
    }
    //Check for the state of resetIncrement across all processors and sync the state
    unsigned int resetIncrementFlag= Utilities::MPI::sum((unsigned int) resetIncrement, mpi_communicator);
    if (resetIncrementFlag>0){
      resetIncrement=true;
      loadFactorSetByModel=Utilities::MPI::min(loadFactorSetByModel, mpi_communicator);
    }
    //MPI operation to sync data
    residual.compress(VectorOperation::add);
    jacobian.compress(VectorOperation::add);
  }

  #include "../../include/ellipticBVP_template_instantiations.h"
