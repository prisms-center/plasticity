//assemble method for ellipticBVP class
#include "../../include/ellipticBVP.h"

#ifndef ASSEMBLE_H
#define ASSEMBLE_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//FE assemble operation
template <int dim>
void ellipticBVP<dim>::assemble(){
  //initialize global data structures to zero
  //The additional compress operations are only to flush out data and
  //switch to the correct write state. For  details look at the documentation
  //for PETScWrappers::MPI::Vector()
  residual.compress(VectorOperation::add); residual=0.0;
  jacobian.compress(VectorOperation::add); jacobian=0.0;

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

#if enableUserModel == 1
	//fill component indices
	std::vector<unsigned int> componentIndices(dofs_per_cell);
	for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	  componentIndices[d1] = fe_values.get_fe().system_to_component_index(d1).first;
	}
	//quadrature loop
	for (unsigned int q=0; q<num_quad_points; ++q){
	  std::vector<double> quadResidual(dofs_per_cell, 0.0), quadJacobian(dofs_per_cell*dofs_per_cell, 0.0);

	  //load history variables
	  std::vector<double> history(numQuadHistoryVariables);
	  for (unsigned int i=0; i<numQuadHistoryVariables; i++){
	    history[i]= quadHistory(cellID, q, i);
	  }

	  //load shape function values and shape function gradient values
	  std::vector<double> shapeValues(dofs_per_cell), shapeGrads(dofs_per_cell*dim);
	  for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	    shapeValues[d1]=fe_values.shape_value(d1, q);
	    for (unsigned int i=0; i<dim; ++i) {
	      shapeGrads[d1*dim+i]=fe_values.shape_grad(d1, q)[i];
	    }
	  }

	  //compute the deformation gradient at this quad point
	  double gradU[dim*dim], F[dim*dim];
	  for (unsigned int d1=0; d1<dim*dim; ++d1){
	    gradU[d1]=0.0;
	  }
	  for (unsigned int d1=0; d1<dofs_per_cell; ++d1){
	    unsigned int i = fe_values.get_fe().system_to_component_index(d1).first;
	    for (unsigned int j=0; j<dim; ++j){
	      double ULocal=solutionWithGhosts[local_dof_indices[d1]];
	      gradU[i*dim+j]+=ULocal*fe_values.shape_grad(d1, q)[j];
	    }
	  }
	   //F=1+gradU
	  for (unsigned int i=0; i<dim; ++i){
	    for (unsigned int j=0; j<dim; ++j){
	      F[i*dim+j] = (i==j) + gradU[i*dim+j];
	    }
	  }

	  //call the getQuadratureValues method supplied by the user in the userModel
	  getQuadratureValues(cellID,
			      dofs_per_cell,
			      &componentIndices[0],
			      &shapeValues[0],
			      &shapeGrads[0],
			      &F[0],
			      &quadResidual[0],
			      &quadJacobian[0],
			      &history[0]);
	  //update elemental residual and jacobian
	  for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	    elementalResidual[d1]+=quadResidual[d1]*fe_values.JxW(q);
	    for (unsigned int d2=0; d2<dofs_per_cell; ++d2) {
	      elementalJacobian(d1,d2)+=quadJacobian[d1*dofs_per_cell+d2]*fe_values.JxW(q);
	    }
	  }
	  //store history variables
	  for (unsigned int i=0; i<numQuadHistoryVariables; i++){
	   quadHistory(cellID, q, i)= history[i];
	  }
	}
#else
	//get elemental jacobian and residual
	getElementalValues(fe_values, dofs_per_cell, num_quad_points, elementalJacobian, elementalResidual);
#endif
	//
	constraints.distribute_local_to_global(elementalJacobian,
					       elementalResidual,
					       local_dof_indices,
					       jacobian,
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
  jacobian.compress(VectorOperation::add);
  //pcout << "boundary size: " << boundary_values.size() << "\n";
  //MatrixTools::apply_boundary_values (boundary_values, jacobian, solution, residual, false);
  //pcout << "boundary size: " << residual.linfty_norm() << "\n";
}

#endif
