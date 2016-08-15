//methods to apply Dirichlet boundary conditons 
#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

template <int dim>
void ellipticBVP<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
}

//methods to apply dirichlet BC's
template <int dim>
void ellipticBVP<dim>::applyDirichletBCs(){
  boundary_values.clear();
  constraints.clear();
  
  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  FEValues<dim> fe_values (FE, QGauss<dim>(1), update_values);
  FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(1), update_values);
  
  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
  unsigned int cellID=0;
  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      cell->get_dof_indices (local_dof_indices);
      fe_values.reinit (cell);
      for (unsigned int faceID=0; faceID<2*dim; faceID++){
	if (cell->face(faceID)->at_boundary()){
	  fe_face_values.reinit (cell, faceID);
	  for (unsigned int i=0; i<dofs_per_cell; ++i) {
	    if (fe_face_values.shape_value(i, 0)!=0){
	      const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
	      unsigned int globalDOF=local_dof_indices[i];
	      bool flag=false;
	      double value=0;
	      Point<dim> node=supportPoints[globalDOF];
	      setBoundaryValues(node, dof, flag, value);
	      if (flag){
		//pcout << "node: " << node << " dof:" << dof << " value:" << value << "\n";
		//boundary_values[globalDOF]=value;
		constraints.add_line (globalDOF);
		if (currentIteration==0){
		  value*=loadFactorSetByModel;
		}
		else{
		  value=0.0;
		}
		constraints.set_inhomogeneity(globalDOF, value);
	      }
	    }
	  }
	}
      }
    }
  }
  //
  constraints.close();
}

#endif
