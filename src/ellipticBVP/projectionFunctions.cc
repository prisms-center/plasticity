//projections method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//initialize post processed field projection
template <int dim>
void ellipticBVP<dim>::initProjection(){
  //return if no post processing fields
  if (numPostProcessedFields==0) return;

  //create and initialize post processing field vectors
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    //residuals
    vectorType* v=new vectorType;
    v->reinit(locally_owned_dofs_Scalar, mpi_communicator); (*v)=0;
    postResidual.push_back(v);

    //fields
    v=new vectorType;
    v->reinit(locally_owned_dofs_Scalar, mpi_communicator); (*v)=0;
    postFields.push_back(v);

    //fields with ghosts
    v=new vectorType;
    v->reinit (locally_owned_dofs_Scalar, locally_relevant_dofs_Scalar, mpi_communicator); (*v)=0;
    postFieldsWithGhosts.push_back(v);
  }

  //initialize postprocessValues data structure
  QGauss<dim> quadrature(userInputs.quadOrder);
  const unsigned int num_quad_points = quadrature.size();
  const unsigned int num_local_cells = triangulation.n_locally_owned_active_cells();
  postprocessValues.reinit(TableIndices<4>(num_local_cells, num_quad_points, numPostProcessedFields, dim));
  postprocessValuesAtCellCenters.reinit(TableIndices<2>(num_local_cells, numPostProcessedFieldsAtCellCenters));

  //create mass matrix
  DynamicSparsityPattern dsp (locally_relevant_dofs_Scalar);
  DoFTools::make_sparsity_pattern (dofHandler_Scalar, dsp, constraintsMassMatrix, false);
  SparsityTools::distribute_sparsity_pattern (dsp,
					      dofHandler_Scalar.n_locally_owned_dofs_per_processor(),
					      mpi_communicator,
					      locally_relevant_dofs_Scalar);
  massMatrix.reinit (locally_owned_dofs_Scalar, locally_owned_dofs_Scalar, dsp, mpi_communicator); massMatrix=0.0;

  //local variables
  FEValues<dim> fe_values (FE_Scalar, quadrature, update_values | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE_Scalar.dofs_per_cell;
  FullMatrix<double>   elementalMass(dofs_per_cell, dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler_Scalar.begin_active(), endc = dofHandler_Scalar.end();
  unsigned int cellID=0;
  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      elementalMass = 0;
      cell->set_user_index(cellID++);

      //compute values for the current element
      fe_values.reinit (cell);

      //elementalMass=N*N
      for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	unsigned int i = fe_values.get_fe().system_to_component_index(d1).first;
	for (unsigned int d2=0; d2<dofs_per_cell; ++d2) {
	  unsigned int j = fe_values.get_fe().system_to_component_index(d2).first;
	  for (unsigned int q=0; q<num_quad_points; ++q){
	    if (i==j){
	      if (i==0){
		elementalMass(d1,d2)+=fe_values.shape_value(d1, q)*fe_values.shape_value(d2, q)*fe_values.JxW(q);
	      }
	      else{
		elementalMass(d1,d2)=1.0;
	      }
	    }
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
void ellipticBVP<dim>::projection(){

  //return if no post processing fields
  if (numPostProcessedFields==0) return;

  //check whether to project in current increment
  if ((currentIncrement+1)%userInputs.skipOutputSteps!=0)
    return;

  pcout << "projecting post processing fields\n";

  //initialize global data structures to zero
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    (*postResidual[field])=0.0;
  }

  //local variables
  QGauss<dim>  quadrature(userInputs.quadOrder);
  FEValues<dim> fe_values (FE_Scalar, quadrature, update_values | update_JxW_values);
  const unsigned int   dofs_per_cell   = FE_Scalar.dofs_per_cell;
  const unsigned int   num_quad_points = quadrature.size();
  Vector<double>       elementalResidual (dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);

  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler_Scalar.begin_active(), endc = dofHandler_Scalar.end();
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
    *postFields[field]=0.0;
    solveLinearSystem2(constraintsMassMatrix, massMatrix, *postResidual[field], *postFields[field],  *postFieldsWithGhosts[field],  *postFieldsWithGhosts[field]);
  }
}
#include "../../include/ellipticBVP_template_instantiations.h"
