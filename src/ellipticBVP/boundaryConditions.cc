//methods to apply Dirichlet boundary conditons
#include "../../include/ellipticBVP.h"

//Specify Dirichlet boundary conditions
template <int dim>
void ellipticBVP<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
  unsigned int i ;

  //pcout<<node[0]<<" "<<node[1]<<" "<<node[2]<<" "<<dof<<std::endl;

  if(userInputs.enableCyclicLoading){
    if(dof==(userInputs.cyclicLoadingDOF-1)){
      switch (userInputs.cyclicLoadingFace){
        case 1:
        if (node[0] == 0.0)
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 2:
        if (node[0] == userInputs.span[0])
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            //pcout<<"Positive"<<std::endl;
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            //pcout<<"negative"<<std::endl;
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{//pcout<<"Positive"<<std::endl;flag=true;
          value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 3:
        if (node[1] == 0.0)
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 4:
        if (node[1] == userInputs.span[1])
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 5:
        if (node[2] == 0.0)
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 6:
        if (node[2] == userInputs.span[2])
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
      }
    }
  }

  for (i=0;i<2*dim;i++){
    if(faceDOFConstrained[i][dof])
      switch (i+1){
        case 1:
        if (node[0] == 0.0)
            {//pcout<<i<<" "<<dof<<std::endl;
              flag=true; value=deluConstraint[i][dof];return;}
        break;
        case 2:
        if (node[0] == userInputs.span[0])
            {//pcout<<i<<" "<<dof<<std::endl;
              flag=true; value=deluConstraint[i][dof];return;}
        break;
        case 3:
        if (node[1] == 0.0)
            {//pcout<<i<<" "<<dof<<std::endl;
              flag=true; value=deluConstraint[i][dof];return;}
        break;
        case 4:
        if (node[1] == userInputs.span[1])
            {//pcout<<i<<" "<<dof<<std::endl;
              flag=true; value=deluConstraint[i][dof];return;}
        break;
        case 5:
        if (node[2] == 0.0)
            {//pcout<<i<<" "<<dof<<std::endl;
              flag=true; value=deluConstraint[i][dof];return;}
        break;
        case 6:
        if (node[2] == userInputs.span[2])
            {//pcout<<i<<" "<<dof<<std::endl;
              flag=true; value=deluConstraint[i][dof];return;}
      }
  }
}

//methods to apply dirichlet BC's
template <int dim>
void ellipticBVP<dim>::applyDirichletBCs(){
  constraints.clear();

  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  FEValues<dim> fe_values (FE, QGauss<dim>(1), update_values);
  FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(1), update_values);

  //parallel loop over all elements
  typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
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
#include "../../include/ellipticBVP_template_instantiations.h"
