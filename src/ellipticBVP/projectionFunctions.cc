//projections method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//initialize post processed field projection
template <int dim>
void ellipticBVP<dim>::initProjection(){
  //return if no post processing fields
  if (numPostProcessedFields==0) return;
  if (!userInputs.writeOutput) return;

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
  #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
  SparsityTools::distribute_sparsity_pattern (dsp,
					      dofHandler_Scalar.n_locally_owned_dofs_per_processor(),
					      mpi_communicator,
					      locally_relevant_dofs_Scalar);
  #else
  SparsityTools::distribute_sparsity_pattern (dsp, dofHandler_Scalar.locally_owned_dofs(), mpi_communicator, locally_relevant_dofs_Scalar);
  #endif

    if(userInputs.enableIndentationBCs) {
        hanging_constraints.clear ();
        hanging_constraints.reinit (locally_relevant_dofs);
        //pcout<<"lambda -issues: locally_relevant_dofs ascending and 1:1? "<<locally_relevant_dofs.is_ascending_and_one_to_one(mpi_communicator)<<"\n";
        //lambda.reinit(locally_relevant_dofs, mpi_communicator); //THIS ISNT WORKING FOR MPI np>1
        //lambda.reinit(locally_owned_dofs,locally_relevant_ghost_dofs,mpi_communicator);
        //local one to one and ascending is violated
        DoFTools::make_hanging_node_constraints (dofHandler, hanging_constraints);
        hanging_constraints.close ();
        DynamicSparsityPattern dsp2 (locally_relevant_dofs);
        DoFTools::make_sparsity_pattern (dofHandler, dsp2, hanging_constraints, false);
        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        SparsityTools::distribute_sparsity_pattern (dsp2,
                                                    dofHandler.n_locally_owned_dofs_per_processor(),
                                                    mpi_communicator,
                                                    locally_relevant_dofs);
        #else
        SparsityTools::distribute_sparsity_pattern (dsp2, dofHandler.locally_owned_dofs(), mpi_communicator, locally_relevant_dofs);
        #endif
        massMatrix.reinit (locally_owned_dofs, locally_owned_dofs, dsp2, mpi_communicator); massMatrix=0.0;
      //  pcout << "massMatrix.domain " << massMatrix.locally_owned_domain_indices().size() <<std::endl;
        diag_mass_matrix_vector.reinit(locally_owned_dofs, mpi_communicator);
        //diag_mass_matrix_vector = 0;
        //residual.compress();
      //  pcout << "before assemble mass matrix diagonal " << std::endl;
        assemble_mass_matrix_diagonal();
      //  pcout << "after assemble mass matrix diagonal " << std::endl;
        unsigned int start = (residual.local_range().first),end = (residual.local_range().second);
      //  std::cout << "start " << start << "end " << end << std::endl;
        for (unsigned int j = start; j < end; ++j) {
            diag_mass_matrix_vector(j) = massMatrix.diag_element(j);
            //std::cout << "mass matrix vector entry " << j << ": " << diag_mass_matrix_vector(j) << std::endl;
        }


        //pcout << "initializing Indentation arrays\n";
        active_set.clear();
        active_set.set_size(dofHandler.n_dofs());
        indentation_constraints.clear();
        indentation_constraints.reinit(locally_relevant_dofs);
        newton_rhs_uncondensed.reinit(locally_owned_dofs, mpi_communicator);
        newton_rhs_uncondensed_inc.reinit(locally_owned_dofs, mpi_communicator);
      //pcout << "before diagonal mass matrix compress " << std::endl;
        diag_mass_matrix_vector.compress(VectorOperation::insert);
        massMatrix = 0.0;
      //pcout << "after diagonal mass matrix compress " << std::endl;
    }
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
  if (!userInputs.writeOutput) return;
  if (numPostProcessedFields==0) return;



  //////////////////////TabularOutput Start///////////////
  std::vector<unsigned int> tabularTimeInputIncInt;
	std::vector<double> tabularTimeInputInc;
if (userInputs.tabularOutput){

  tabularTimeInputInc=userInputs.tabularTimeOutput;
  for(unsigned int i=0;i<userInputs.tabularTimeOutput.size();i++){
    tabularTimeInputInc[i]=tabularTimeInputInc[i]/delT;
  }
  tabularTimeInputIncInt.resize(userInputs.tabularTimeOutput.size(),0);
  ///Converting to an integer always rounds down, even if the fraction part is 0.99999999.
  //Hence, I add 0.1 to make sure we always get the correct integer.
  for(unsigned int i=0;i<userInputs.tabularTimeOutput.size();i++){
    tabularTimeInputIncInt[i]=int(tabularTimeInputInc[i]+0.1);
  }
}
  //////////////////////TabularOutput Finish///////////////
  //check whether to project in current increment
  if (((!userInputs.tabularOutput)&&((currentIncrement+1)%userInputs.skipOutputSteps!=0))||((userInputs.tabularOutput)&& (std::count(tabularTimeInputIncInt.begin(), tabularTimeInputIncInt.end(), (currentIncrement+1))==0))){
    return;
  }

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
