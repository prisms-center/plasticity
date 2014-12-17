//assemble method for ellipticBVP class

#ifndef ASSEMBLE_H
#define ASSEMBLE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//FE assemble operation
template <int dim>
void ellipticBVP<dim>::assemble(){
  //initialize global data structures to zero
  solution=0.0; 
  residual=0.0; 
  jacobian=0.0;

  //local variables
  QGauss<dim>  quadrature(quadOrder);
  FEValues<dim> fe_values (FE, quadrature, update_values | update_gradients | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  const unsigned int   num_quad_points = quadrature.size();
  FullMatrix<double>   elementalJacobian (dofs_per_cell, dofs_per_cell);
  Vector<double>       elementalResidual (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  //apply Dirichlet BC's
  applyDirichletBCs();  
  
  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
  unsigned int cellID=0;
  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      elementalJacobian = 0;
      elementalResidual = 0;
			cell->set_user_index(cellID++);

      //Compute values for the current element
      fe_values.reinit (cell);
      
      //get elemental jacobian and residual
      getElementalValues(fe_values, dofs_per_cell, num_quad_points, elementalJacobian, elementalResidual);
      
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(elementalJacobian, 
					     elementalResidual,
					     local_dof_indices,
					     jacobian, 
					     residual);

    }
  }
  //MPI operation to sync data 
  residual.compress(VectorOperation::add);
  jacobian.compress(VectorOperation::add);
}

#endif
