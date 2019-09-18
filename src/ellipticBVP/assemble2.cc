//assemble method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//FE assemble operation
template <int dim>
void ellipticBVP<dim>::assemble2(){
  //initialize global data structures to zero
  //The additional compress operations are only to flush out data and
  //switch to the correct write state. For  details look at the documentation
  //for PETScWrappers::MPI::Vector()
  residual.compress(VectorOperation::add); residual=0.0;


  //local variables
  QGauss<dim>  quadrature(userInputs.quadOrder);
  FEValues<dim> fe_values (FE, quadrature, update_values | update_gradients | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  const unsigned int   num_quad_points = quadrature.size();
  Vector<double>       elementalResidual (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //apply Dirichlet BC's
  applyDirichletBCs();

  try{
    //parallel loop over all elements
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    unsigned int cellID=0;
    for (; cell!=endc; ++cell) {
      if (cell->is_locally_owned()){
        elementalResidual = 0;
        cell->set_user_index(cellID);

        //Compute values for the current element
        fe_values.reinit (cell);
        cell->get_dof_indices (local_dof_indices);
        //get elemental jacobian and residual
        getElementalValues2(fe_values, dofs_per_cell, num_quad_points, elementalResidual);
        //
        constraints.distribute_local_to_global(elementalResidual,
          local_dof_indices,
          residual);

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
  }

  #include "../../include/ellipticBVP_template_instantiations.h"
