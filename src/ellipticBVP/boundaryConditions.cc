//methods to apply Dirichlet boundary conditons
#include "../../include/ellipticBVP.h"

//Specify Dirichlet boundary conditions
template <int dim>
void ellipticBVP<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
  unsigned int i,dof_1, dof_2 ;

  Vector<double> externalMeshParameterBCs(dim),nodalDisplacementBCsToleranceVector(dim);
  double alpha_Torsion;
  for (unsigned int i=0; i<dim; ++i) {
    externalMeshParameterBCs(i)=userInputs.externalMeshParameter*userInputs.span[i];
    nodalDisplacementBCsToleranceVector(i)=userInputs.nodalDisplacementBCsTolerance*userInputs.span[i];
  }

  if(userInputs.enableCyclicLoading){
    if(dof==(userInputs.cyclicLoadingDOF-1)){
      switch (userInputs.cyclicLoadingFace){
        case 1:
        if (node[0] <= externalMeshParameterBCs(0))
        {
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
	}
      	break;

        case 2:
        if (node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))
        {
          //pcout<<"Positive12"<<std::endl;
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            //pcout<<"Positive"<<std::endl;
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            //pcout<<"negative"<<std::endl;
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{//pcout<<"Positive12"<<std::endl;
          flag=true;value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        }
        break;
        case 3:
        if (node[1] <= externalMeshParameterBCs(1))
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 4:
        if (node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 5:
        if (node[2] <= externalMeshParameterBCs(2))
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
        break;
        case 6:
        if (node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))
          if(fmod((currentIncrement*delT),cycleTime)<userInputs.quarterCycleTime){
            flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else if(fmod((currentIncrement*delT),cycleTime)<3*userInputs.quarterCycleTime){
            flag=true; value=-deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
          else{flag=true; value=deluConstraint[userInputs.cyclicLoadingFace-1][dof];return;}
      }
    }
  }

  if((userInputs.enableSimpleBCs)||(userInputs.enableCyclicLoading)){
    for (i=0;i<2*dim;i++){
      if(faceDOFConstrained[i][dof])
        switch (i+1){
          case 1:
          if (node[0] <= externalMeshParameterBCs(0))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 2:
          if (node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 3:
          if (node[1] <= externalMeshParameterBCs(1))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 4:
          if (node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 5:
          if (node[2] <= externalMeshParameterBCs(2))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
          break;
          case 6:
          if (node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true; value=deluConstraint[i][dof];return;}
        }
    }
  }

  if(userInputs.enableTorsionBCs){
    currentTime=delT*(currentIncrement+1);

    if (currentIncrement==0){
      timeCounter=1;
    }
    if (currentTime>userInputs.tabularTimeInputTorsion[timeCounter]){
      timeCounter=timeCounter+1;
    }
    alpha_Torsion=userInputs.tabularTorsionBCsInput[timeCounter]*delT;

    if (userInputs.torsionAxis==2){ //z-axis is torsion axis
      dof_1=0; dof_2=1;
    }
    else if (userInputs.torsionAxis==0){ //x-axis is torsion axis
      dof_1=1; dof_2=2;
    }
    else{ //y-axis is torsion axis
      dof_1=2; dof_2=0;
    }

    if ((node[userInputs.torsionAxis] >= (userInputs.span[userInputs.torsionAxis]-externalMeshParameterBCs(userInputs.torsionAxis)))&&(dof==dof_1)){
                flag=true; value=-(node[dof_2]-userInputs.centerTorsion[1])*alpha_Torsion;}
    if ((node[userInputs.torsionAxis] >= (userInputs.span[userInputs.torsionAxis]-externalMeshParameterBCs(userInputs.torsionAxis)))&&(dof==dof_2)){
                flag=true; value=(node[dof_1]-userInputs.centerTorsion[0])*alpha_Torsion;}
  }

  if(userInputs.enableNodalDisplacementBCs){
    for (i=0;i<userInputs.numberOfNodalBCs;i++){
      if ((fabs(node[0]-nodalDisplacement[i][0])<= nodalDisplacementBCsToleranceVector(0))&&(fabs(node[1]-nodalDisplacement[i][1])<= nodalDisplacementBCsToleranceVector(1))&&(fabs(node[2]-nodalDisplacement[i][2])<= nodalDisplacementBCsToleranceVector(2))){
        if (dof==dofNodalDisplacement[i]){
          flag=true; value=deluNodalDisplacement[i];return;
        }
      }
    }
  }

  if(userInputs.enableTabularBCs){

    for (i=0;i<2*dim;i++){
      if(faceDOFConstrained[i][dof])
        switch (i+1){
          case 1:
          if (node[0] <= externalMeshParameterBCs(0))
              {//pcout<<i<<" "<<dof<<std::endl;
                flag=true;

                currentTime=delT*(currentIncrement+1);

                if (currentIncrement==0){
                  timeCounter=1;
                }
                if (currentTime>userInputs.tabularTimeInput[timeCounter]){
                  timeCounter=timeCounter+1;
                }

                value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs.tabularTimeInput[timeCounter]-userInputs.tabularTimeInput[timeCounter-1])*delT ;

                return;}
          break;
          case 2:
          if (node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime=delT*(currentIncrement+1);

              if (currentIncrement==0){
                timeCounter=1;
              }
              if (currentTime>userInputs.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs.tabularTimeInput[timeCounter]-userInputs.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 3:
          if (node[1] <= externalMeshParameterBCs(1))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime=delT*(currentIncrement+1);

              if (currentIncrement==0){
                timeCounter=1;
              }
              if (currentTime>userInputs.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs.tabularTimeInput[timeCounter]-userInputs.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 4:
          if (node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime=delT*(currentIncrement+1);

              if (currentIncrement==0){
                timeCounter=1;
              }
              if (currentTime>userInputs.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs.tabularTimeInput[timeCounter]-userInputs.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 5:
          if (node[2] <= externalMeshParameterBCs(2))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime=delT*(currentIncrement+1);

              if (currentIncrement==0){
                timeCounter=1;
              }
              if (currentTime>userInputs.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs.tabularTimeInput[timeCounter]-userInputs.tabularTimeInput[timeCounter-1])*delT ;

              return;}
          break;
          case 6:
          if (node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))
              {//pcout<<i<<" "<<dof<<std::endl;
              flag=true;

              currentTime=delT*(currentIncrement+1);

              if (currentIncrement==0){
                timeCounter=1;
              }
              if (currentTime>userInputs.tabularTimeInput[timeCounter]){
                timeCounter=timeCounter+1;
              }

              value=(-tabularDisplacements[3*i+dof][timeCounter-1]+tabularDisplacements[3*i+dof][timeCounter])/(userInputs.tabularTimeInput[timeCounter]-userInputs.tabularTimeInput[timeCounter-1])*delT ;

              return;}
        }
    }
  }

  if(userInputs.useVelocityGrad){
    value = 0;
    for(i=0;i<dim;i++)
      value+=deltaF[dof][i]*node[i];
    flag=true;
    return;
  }

  if(userInputs.enableDICpipeline){

  //back boundary
  if (node[0] <= externalMeshParameterBCs(0)){
    double value_x, value_y;
    bcFunction1(node[1],value_x,value_y, currentIncrement);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}
  }


  //front boundary
  if (node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0))){
    double value_x, value_y;
    bcFunction2(node[1],value_x,value_y, currentIncrement);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}
  }


  //left boundary
  if (node[1] <= externalMeshParameterBCs(1)){
    double value_x, value_y;
    bcFunction3(node[0],value_x,value_y, currentIncrement);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}

  }

  //left boundary
  if (node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1))){
    double value_x, value_y;
    bcFunction4(node[0],value_x,value_y, currentIncrement);
    if (dof==0) {flag=true; value=value_x;}
    if (dof==1) {flag=true; value=value_y;}

  }


  //bottom boundary: u_z=0
  if (node[2] <= externalMeshParameterBCs(2)){
    if(node[0]<= externalMeshParameterBCs(0) && node[1] <= externalMeshParameterBCs(1)){
      if (dof==2) {flag=true; value=0.0;}
    }
  }
  }

}


template <int dim>
void ellipticBVP<dim>::bcFunction1(double yval, double &value_x, double &value_y, double currentIncr){

  unsigned int n_div=userInputs.Y_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(yval-bc_new1[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);

  currentTime=delT*(currentIncrement+1);

  if (currentIncrement==0){
    timeCounter=1;
  }
  if (currentTime>userInputs.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new1[coord_pos][1+(timeCounter-1)*2]+bc_new1[coord_pos][1+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new1[coord_pos][2+(timeCounter-1)*2]+bc_new1[coord_pos][2+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;
  //return value1 ; // displacement along X-Direction

}

template <int dim>
void ellipticBVP<dim>::bcFunction2(double yval, double &value_x, double &value_y,double currentIncr){

  unsigned int n_div=userInputs.Y_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(yval-bc_new2[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);

  currentTime=delT*(currentIncrement+1);

  if (currentIncrement==0){
    timeCounter=1;
  }
  if (currentTime>userInputs.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new2[coord_pos][1+(timeCounter-1)*2]+bc_new2[coord_pos][1+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new2[coord_pos][2+(timeCounter-1)*2]+bc_new2[coord_pos][2+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;



  //return value1 ; // displacement along X-Direction

}

template <int dim>
void ellipticBVP<dim>::bcFunction3(double xval, double &value_x, double &value_y,double currentIncr){

  unsigned int n_div=userInputs.X_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(xval-bc_new3[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);

  currentTime=delT*(currentIncrement+1);

  if (currentIncrement==0){
    timeCounter=1;
  }
  if (currentTime>userInputs.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new3[coord_pos][1+(timeCounter-1)*2]+bc_new3[coord_pos][1+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new3[coord_pos][2+(timeCounter-1)*2]+bc_new3[coord_pos][2+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;


}

template <int dim>
void ellipticBVP<dim>::bcFunction4(double xval, double &value_x, double &value_y,double currentIncr){

  unsigned int n_div=userInputs.X_dic;
  std::vector<double> coord(n_div);

  for (unsigned int i=0; i<(n_div); i++){
    coord[i]=fabs(xval-bc_new4[i][0]);
  }

  std::vector<double>::iterator result;
  result = std::min_element(coord.begin(), coord.end());
  int coord_pos;
  coord_pos= std::distance(coord.begin(), result);
  currentTime=delT*(currentIncrement+1);

  if (currentIncrement==0){
    timeCounter=1;
  }
  if (currentTime>userInputs.timeInputDIC[timeCounter]){
    timeCounter=timeCounter+1;
  }

  value_x=(-bc_new4[coord_pos][1+(timeCounter-1)*2]+bc_new4[coord_pos][1+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;
  value_y=(-bc_new4[coord_pos][2+(timeCounter-1)*2]+bc_new4[coord_pos][2+timeCounter*2])/(userInputs.timeInputDIC[timeCounter]-userInputs.timeInputDIC[timeCounter-1])*delT ;


}



//methods to apply dirichlet BC's
template <int dim>
void ellipticBVP<dim>::applyDirichletBCs(){
    //pcout<<"debug setDirichlet 0\n";
  if(!userInputs.enablePeriodicBCs){
     // pcout<<"debug setDirichlet b1\n";
    constraints.clear();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dofHandler, constraints);
    const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
     // pcout<<"debug setDirichlet b2\n";
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim> fe_values (FE, QGauss<dim>(userInputs.quadOrder), update_values);
    FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs.quadOrder), update_values);
     // pcout<<"debug setDirichlet b3\n";
    //parallel loop over all elements
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    for (; cell!=endc; ++cell) {
      if (cell->is_locally_owned()){
        cell->get_dof_indices (local_dof_indices);
        fe_values.reinit (cell);
        for (unsigned int faceID=0; faceID<GeometryInfo<dim>::faces_per_cell; faceID++){
          if (cell->face(faceID)->at_boundary()){
            fe_face_values.reinit (cell, faceID);
            for (unsigned int i=0; i<dofs_per_cell; ++i) {
              if (fe_face_values.shape_value(i, 0)!=0){
                const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
                unsigned int globalDOF=local_dof_indices[i];
                bool flag=false;
                double value=0;
                node=supportPoints[globalDOF];
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
  }
  else{
    setPeriodicityConstraints();
  }
  if (userInputs.enableIndentationBCs){
      //pcout<<"debug setDirichlet c1\n";
      setIndentationConstraints();
      //pcout<<"debug setDirichlet c2\n";
      constraints.merge(indentation_constraints,dealii::AffineConstraints<>::right_object_wins);
      //pcout<<"debug setDirichlet c3\n";
      //indentation conditions will overwrite dirichlet conditions from other sources (as we prefer)
  }
  constraints.close ();
  //pcout<<"debug setDirichlet 4\n";

}
#include "../../include/ellipticBVP_template_instantiations.h"
