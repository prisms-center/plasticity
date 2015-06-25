//project method for ellipticBVP class

#ifndef PROJECT_H
#define PROJECT_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//initialize post processed field projection
template <int dim>
void ellipticBVP<dim>::initProject(){
  //return if no post processing fields
  if (numPostProcessedFields==0) return;

  //create and initialize post processing field vectors
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    //residuals
    vectorType* v=new vectorType;
    v->reinit(locally_owned_dofs, mpi_communicator); (*v)=0;
    postResidual.push_back(v);   
    
    //fields
    v=new vectorType;
    v->reinit(locally_owned_dofs, mpi_communicator); (*v)=0;
    postFields.push_back(v); 

    //fields with ghosts
    v=new vectorType;
    v->reinit (locally_owned_dofs, locally_relevant_dofs, mpi_communicator); (*v)=0;
    postFieldsWithGhosts.push_back(v); 
  }

  //initialize postprocessValues data structure
  QGauss<dim> quadrature(quadOrder);
  const unsigned int num_quad_points = quadrature.size();
  const unsigned int num_local_cells = triangulation.n_locally_owned_active_cells();
  postprocessValues.reinit(TableIndices<4>(num_local_cells, num_quad_points, numPostProcessedFields, dim));

  //create mass matrix
  CompressedSimpleSparsityPattern csp (locally_relevant_dofs);
  DoFTools::make_sparsity_pattern (dofHandler, csp, constraints, false);
  SparsityTools::distribute_sparsity_pattern (csp,
					      dofHandler.n_locally_owned_dofs_per_processor(),
					      mpi_communicator,
					      locally_relevant_dofs);
  massMatrix.reinit (locally_owned_dofs, locally_owned_dofs, csp, mpi_communicator); massMatrix=0.0;

  //local variables
  FEValues<dim> fe_values (FE, quadrature, update_values | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  FullMatrix<double>   elementalMass(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  
  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
  unsigned int cellID=0;
  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      elementalMass = 0;
      cell->set_user_index(cellID++);

      //compute values for the current element
      fe_values.reinit (cell);
      
      //elementalMass=N*N
      for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	for (unsigned int d2=0; d2<dofs_per_cell; ++d2) {
	  for (unsigned int q=0; q<num_quad_points; ++q){
	    elementalMass(d1,d2)+=fe_values.shape_value(d1, q)*fe_values.shape_value(d2, q)*fe_values.JxW(q);
	  }
	}
      }
      
      //assemble
      cell->get_dof_indices (local_dof_indices);
      constraintsMassMatrix.distribute_local_to_global(elementalMass, local_dof_indices, massMatrix);
    }
  }
  //MPI operation to sync data 
  massMatrix.compress(VectorOperation::add);
}

//post processed field projection operation
template <int dim>
void ellipticBVP<dim>::project(){
  //return if no post processing fields
  if (numPostProcessedFields==0) return;
  
  pcout << "projecting post processing fields\n ";
  //initialize global data structures to zero  
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    (*postResidual[field])=0.0;    
  }
  
  //local variables
  QGauss<dim>  quadrature(quadOrder);
  FEValues<dim> fe_values (FE, quadrature, update_values | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  const unsigned int   num_quad_points = quadrature.size();
  Vector<double>       elementalResidual (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    
  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
  unsigned int cellID=0;
  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      cell->set_user_index(cellID);

      //compute values for the current element
      fe_values.reinit (cell);
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int field=0; field<numPostProcessedFields; field++){
	elementalResidual = 0;
	
	//elementalSolution=N*solution
	for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	  unsigned int i = fe_values.get_fe().system_to_component_index(d1).first;
	  for (unsigned int q=0; q<num_quad_points; ++q){
	    elementalResidual(d1)+=fe_values.shape_value(d1, q)*postprocessValues(cellID, q, field, i)*fe_values.JxW(q);
	  }
	}
	//assemble
	constraintsMassMatrix.distribute_local_to_global(elementalResidual, local_dof_indices, *postResidual[field]);
      }
      cellID++;
    }
  }

  //MPI operation to sync data 
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    postResidual[field]->compress(VectorOperation::add);
    
    //L2 projection by solving for Mx=b problem
    solveLinearSystem(constraintsMassMatrix, massMatrix, *postResidual[field], *postFields[field],  *postFieldsWithGhosts[field],  *postFieldsWithGhosts[field]);
  }
}

#endif
