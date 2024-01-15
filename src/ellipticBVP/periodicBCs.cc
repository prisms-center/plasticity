//solve method for ellipticBVP class
#include "../../include/ellipticBVP.h"
#include <fstream>
#include <iostream>

//loop over increments and solve each increment
template <int dim>
void ellipticBVP<dim>::setPeriodicity(){
  //This functions connects the periodic faces to each other, and add those dofs as ghost cells of the other face.
  std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > periodicity_vector;

///Here, we assumed that in the face constraints lines, we always start with FXP, FYP, or FZP.
  GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 1, /*b_id2*/ 0, /*direction*/ 0, periodicity_vector);
  GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 3, /*b_id2*/ 2, /*direction*/ 1, periodicity_vector);
  GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 5, /*b_id2*/ 4, /*direction*/ 2, periodicity_vector);


  triangulation.add_periodicity(periodicity_vector);
  pcout << "periodic facepairs: " << periodicity_vector.size() << std::endl;
}

template <int dim>
void ellipticBVP<dim>::setPeriodicityConstraintsInit(){
  // This function generates the required inputs for applying Periodic BCs constraints.

  // #set Vertices Periodic BCs row order:
  // # 0=V000_1;1=V000_2;2=V000_3; 3=V100_1;4=V100_2;5=V100_3; 6=V010_1;7=V010_2;8=V010_3; 9=V001_1;10=V001_2;11=V001_3;
  // # 12=V110_1;13=V110_2;14=V110_3; 15=V101_1;16=V101_2;17=V101_3; 18=V011_1;19=V011_2;20=V011_3; 21=V111_1;22=V111_2;23=V111_3;
  // #set Edges Periodic BCs row order:
  // # 0=EX000_1;1=EX000_2;2=EX000_3; 3=EX001_1;4=EX001_2;5=EX001_3; 6=EX010_1;7=EX010_2;8=EX010_3; 9=EX011_1;10=EX011_2;11=EX011_3;
  // # 12=EY000_1;13=EY000_2;14=EY000_3; 15=EY001_1;16=EY001_2;17=EY001_3; 18=EY100_1;19=EY100_2;20=EY100_3; 21=EY101_1;22=EY101_2;23=EY101_3;
  // # 24=EZ000_1;25=EZ000_2;26=EZ000_3; 27=EZ010_1;28=EZ010_2;29=EZ010_3; 30=EZ100_1;31=EZ100_2;32=EZ100_3; 33=EZ110_1;34=EZ110_2;35=EZ110_3;

  // #set Faces Periodic BCs row order:
  // # 0=FXN_1;1=FXN_2;2=FXN_3; 3=FXP_1;4=FXP_2;5=FXP_3; 6=FYN_1;7=FYN_2;8=FYN_3; 9=FYP_1;10=FYP_2;11=FYP_3;
  // # 12=FZN_1;13=FZN_2;14=FZN_3; 15=FZP_1;16=FZP_2;17=FZP_3;

//////////////////////////Initialization Start////////////////////////////////
//vertices_DOFs: This index set has the global dofs for vertices.
  vertices_DOFs.set_size(1000000000);
  dof_Boundary_Layer2.set_size(1000000000);
  Edges_DOFs_1.set_size(1000000000);
  Edges_DOFs_2.set_size(1000000000);
  Edges_DOFs_3.set_size(1000000000);

  totalNumEdgesDOFs=36;
  periodicBCsInput=userInputs.periodicBCsInput;
  periodicBCsInput2_Orig=userInputs.periodicBCsInput2;
  totalNumVerticesDOFs=24;
  totalNumEachEdgesNodes=userInputs.subdivisions[1]-1;
  totalNumFacesDOFs=18;

  IndexSet dof_FXN_1,dof_FXN_2,dof_FXN_3, dof_FXP_1, dof_FXP_2, dof_FXP_3;
  IndexSet dof_FYN_1,dof_FYN_2,dof_FYN_3, dof_FYP_1, dof_FYP_2, dof_FYP_3;
  IndexSet dof_FZN_1,dof_FZN_2,dof_FZN_3, dof_FZP_1, dof_FZP_2, dof_FZP_3;

  unsigned int flagVertices;
  unsigned int check1,check2;
  unsigned int kFlag=0;
  unsigned int local_num_Edges_nodes;
  unsigned int size_dof_Edge_processor;
  unsigned int checkBoundary_id;
  unsigned int n_processes;
  unsigned int size_dof_Boundary_Layer2;
  unsigned int total_num_edges_node_processor;

  double node0,node1,node2;

  std::vector<std::vector<unsigned int> >  global_Edges_DOFs_Vector,Edges_DOFs_Vector,local_Edges_DOFs_Vector;
  std::vector<std::vector<double>>  global_Edges_DOFs_Coord_Vector,Edges_DOFs_Coord_Vector,local_Edges_DOFs_Coord_Vector;


  const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
  std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
  FEValues<dim> fe_values (FE, QGauss<dim>(userInputs.quadOrder), update_values);
  FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs.quadOrder), update_values);

//vertices_DOFs_vector: This vector has the the global dofs for vertices.
  vertices_DOFs_vector.resize(totalNumVerticesDOFs,0);
//vertices_Constraint_Known: This vector will give 0 or 1 for each vertices constraint. 1: means
//that all dofs invlovled in the that constraint line are defined. 0: means not 1!
  vertices_Constraint_Known.resize(numberVerticesConstraint,0);

//edges_DOFs: In each row of this vector, there is a vector including the dofs
//of an edge (for example the first row is the all dofs for EX000_1);
  edges_DOFs.resize(totalNumEdgesDOFs,std::vector<unsigned int>(totalNumEachEdgesNodes,0));

//faces_dof_Index_vector:In each row of this vector, an indexset of dof for each face is saved.
  faces_dof_Index_vector.resize(totalNumFacesDOFs);

  globalDOF_V000_1=0;globalDOF_V000_2=0;globalDOF_V000_3=0;globalDOF_V001_1=0;
  globalDOF_V001_2=0;globalDOF_V001_3=0;globalDOF_V010_1=0;globalDOF_V010_2=0;
  globalDOF_V010_3=0;globalDOF_V100_1=0;globalDOF_V100_2=0;globalDOF_V100_3=0;

  globalDOF_V011_1=0;globalDOF_V011_2=0;globalDOF_V011_3=0;globalDOF_V110_1=0;
  globalDOF_V110_2=0;globalDOF_V110_3=0;globalDOF_V101_1=0;globalDOF_V101_2=0;
  globalDOF_V101_3=0;globalDOF_V111_1=0;globalDOF_V111_2=0;globalDOF_V111_3=0;

  local_globalDOF_V000_1=0;local_globalDOF_V000_2=0;local_globalDOF_V000_3=0;local_globalDOF_V001_1=0;
  local_globalDOF_V001_2=0;local_globalDOF_V001_3=0;local_globalDOF_V010_1=0;local_globalDOF_V010_2=0;
  local_globalDOF_V010_3=0;local_globalDOF_V100_1=0;local_globalDOF_V100_2=0;local_globalDOF_V100_3=0;

  local_globalDOF_V011_1=0;local_globalDOF_V011_2=0;local_globalDOF_V011_3=0;local_globalDOF_V110_1=0;
  local_globalDOF_V110_2=0;local_globalDOF_V110_3=0;local_globalDOF_V101_1=0;local_globalDOF_V101_2=0;
  local_globalDOF_V101_3=0;local_globalDOF_V111_1=0;local_globalDOF_V111_2=0;local_globalDOF_V111_3=0;

  if(userInputs.enableTabularPeriodicBCs){
    for(unsigned int i=0;i<totalNumVerticesDOFs;i++){
      periodicBCsInput2_Orig[0][i]=periodicBCsInput2_Orig[0][i]/periodicTotalIncrements;
    }
  }
  else {
    for(unsigned int i=0;i<totalNumVerticesDOFs;i++){
      periodicBCsInput2_Orig[0][i]=periodicBCsInput2_Orig[0][i]/totalIncrements;
    }
  }


//////////////////////////Initialization Finish///////////////////////////////


////////////////////////Redefining boundary id Start////////////////////////
//Here, The external boundary conditions of the faces of the elements at vertices and edges
//will be changed from its initial to (100).
//The reason to do this, later on, when we want to apply the face constraint using
//collect_periodic_faces and make_periodicity_constraints in setFaceConstraintsInc0 and setFaceConstraintsIncNot0,
//we want to separatethe Vertices and Edges from the faces.
  Vector<double> externalMeshParameterBCs(dim);
  for (unsigned int i=0; i<dim; ++i) {
    externalMeshParameterBCs(i)=userInputs.externalMeshParameter*userInputs.span[i];
  }

  for (const auto &cell : triangulation.cell_iterators()){

    for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
    {
      vertexNode=cell->vertex(v);
      if ((((vertexNode[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))||(vertexNode[0] <= externalMeshParameterBCs(0)))&&
      ((vertexNode[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))||(vertexNode[1] <= externalMeshParameterBCs(1))))||
      (((vertexNode[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))||(vertexNode[0] <= externalMeshParameterBCs(0)))&&
      ((vertexNode[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))||(vertexNode[2] <= externalMeshParameterBCs(2))))||
      (((vertexNode[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))||(vertexNode[2] <= externalMeshParameterBCs(2)))&&
      ((vertexNode[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))||(vertexNode[1] <= externalMeshParameterBCs(1))))){
        for (unsigned int face_number = 0; face_number < GeometryInfo<dim>::faces_per_cell;++face_number){
          if (cell->face(face_number)->at_boundary()){
            cell->face(face_number)->set_boundary_id(100);
          }
        }
      }
    }
  }

////////////////////////Redefining boundary id End////////////////////////

///////////////////////Defining different dofs on the faces Start//////////////
//After the boundary id of vertices and edges are separated from the face, we group all nodes on different faces
//based on the direction of the normal to the face. Next, for each node, we have
// 3 dofs in 3d, and the number define which of these dofs we mention. If it doesn't have number, it means it includes all three dofs.
// P and N  in the name defines the positive face or negative face.
//Number of boundary_ids are based on the deal.ii notation for, which can be fined in:
// https://www.dealii.org/current/doxygen/deal.II/namespaceGridGenerator.html under the command of subdivided_hyper_rectangle.

  const FEValuesExtractors::Scalar  displacementX(0);
  const FEValuesExtractors::Scalar  displacementY(1);
  const FEValuesExtractors::Scalar  displacementZ(2);

  std::set< types::boundary_id > boundary_ids1;
  boundary_ids1.insert(1);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementX),dof_FXP_1,boundary_ids1);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementY),dof_FXP_2,boundary_ids1);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementZ),dof_FXP_3,boundary_ids1);
  DoFTools::extract_boundary_dofs(dofHandler, ComponentMask(),dof_FXP,boundary_ids1);

  std::set< types::boundary_id > boundary_ids2;
  boundary_ids2.insert(3);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementX),dof_FYP_1,boundary_ids2);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementY),dof_FYP_2,boundary_ids2);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementZ),dof_FYP_3,boundary_ids2);
  DoFTools::extract_boundary_dofs(dofHandler, ComponentMask(),dof_FYP,boundary_ids2);

  std::set< types::boundary_id > boundary_ids3;
  boundary_ids3.insert(5);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementX),dof_FZP_1,boundary_ids3);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementY),dof_FZP_2,boundary_ids3);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementZ),dof_FZP_3,boundary_ids3);
  DoFTools::extract_boundary_dofs(dofHandler, ComponentMask(),dof_FZP,boundary_ids3);

//Here, We merge dof_FXP,dof_FYP, and dof_FZP as dof_FP to include all dofs of positive faces.
  dof_FP=dof_FXP;
  if (dof_FYP.n_elements()>0){
    IndexSet::ElementIterator dof_FYP_counter = dof_FYP.begin();
    if (!dof_FP.is_element(*dof_FYP_counter)){
      dof_FP.add_index(*dof_FYP_counter);
    }
    for(unsigned int j = 1; j < dof_FYP.n_elements(); j++){
      dof_FYP_counter++;
      if (!dof_FP.is_element(*dof_FYP_counter)){
        dof_FP.add_index(*dof_FYP_counter);
      }
    }
  }


  if (dof_FZP.n_elements()>0){
    IndexSet::ElementIterator dof_FZP_counter = dof_FZP.begin();
    if (!dof_FP.is_element(*dof_FZP_counter)){
      dof_FP.add_index(*dof_FZP_counter);
    }
    for(unsigned int j = 1; j < dof_FZP.n_elements(); j++){
      dof_FZP_counter++;
      if (!dof_FP.is_element(*dof_FZP_counter)){
        dof_FP.add_index(*dof_FZP_counter);
      }
    }
  }


//Here, we are doing the same for negative faces.
  std::set< types::boundary_id > boundary_ids4;
  boundary_ids4.insert(0);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementX),dof_FXN_1,boundary_ids4);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementY),dof_FXN_2,boundary_ids4);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementZ),dof_FXN_3,boundary_ids4);
  DoFTools::extract_boundary_dofs(dofHandler, ComponentMask(),dof_FXN,boundary_ids4);

  std::set< types::boundary_id > boundary_ids5;
  boundary_ids5.insert(2);
  // const FEValuesExtractors::Scalar  displacementY(1);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementX),dof_FYN_1,boundary_ids5);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementY),dof_FYN_2,boundary_ids5);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementZ),dof_FYN_3,boundary_ids5);
  DoFTools::extract_boundary_dofs(dofHandler, ComponentMask(),dof_FYN,boundary_ids5);

  std::set< types::boundary_id > boundary_ids6;
  boundary_ids6.insert(4);
  // const FEValuesExtractors::Scalar  displacementZ(2);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementX),dof_FZN_1,boundary_ids6);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementY),dof_FZN_2,boundary_ids6);
  DoFTools::extract_boundary_dofs(dofHandler, FE.component_mask(displacementZ),dof_FZN_3,boundary_ids6);
  DoFTools::extract_boundary_dofs(dofHandler, ComponentMask(),dof_FZN,boundary_ids6);

//Here, We merge dof_FXN,dof_FYN, and dof_FZN as dof_FN to include all dofs of negative faces.
  dof_FN=dof_FXN;
  if (dof_FYN.n_elements()>0){
    IndexSet::ElementIterator dof_FYN_counter = dof_FYN.begin();
    if (!dof_FN.is_element(*dof_FYN_counter)){
      dof_FN.add_index(*dof_FYN_counter);
    }
    for(unsigned int j = 1; j < dof_FYN.n_elements(); j++){
      dof_FYN_counter++;
      if (!dof_FN.is_element(*dof_FYN_counter)){
        dof_FN.add_index(*dof_FYN_counter);
      }
    }
  }

  if (dof_FZN.n_elements()>0){
    IndexSet::ElementIterator dof_FZN_counter = dof_FZN.begin();
    if (!dof_FN.is_element(*dof_FZN_counter)){
      dof_FN.add_index(*dof_FZN_counter);
    }
    for(unsigned int j = 1; j < dof_FZN.n_elements(); j++){
      dof_FZN_counter++;
      if (!dof_FN.is_element(*dof_FZN_counter)){
        dof_FN.add_index(*dof_FZN_counter);
      }
    }
  }


/////Assigning the dof for each face into a vector faces_dof_Index_vector.
  faces_dof_Index_vector[0]=dof_FXN_1; faces_dof_Index_vector[1]=dof_FXN_2; faces_dof_Index_vector[2]=dof_FXN_3;
  faces_dof_Index_vector[3]=dof_FXP_1; faces_dof_Index_vector[4]=dof_FXP_2; faces_dof_Index_vector[5]=dof_FXP_3;

  faces_dof_Index_vector[6]=dof_FYN_1; faces_dof_Index_vector[7]=dof_FYN_2; faces_dof_Index_vector[8]=dof_FYN_3;
  faces_dof_Index_vector[9]=dof_FYP_1; faces_dof_Index_vector[10]=dof_FYP_2; faces_dof_Index_vector[11]=dof_FYP_3;

  faces_dof_Index_vector[12]=dof_FZN_1; faces_dof_Index_vector[13]=dof_FZN_2; faces_dof_Index_vector[14]=dof_FZN_3;
  faces_dof_Index_vector[15]=dof_FZP_1; faces_dof_Index_vector[16]=dof_FZP_2; faces_dof_Index_vector[17]=dof_FZP_3;

///////////////////////Defining different dofs on the faces End//////////////

///////////////////////Defining global dofs for Vertices Start//////////////
//(DOFs on Edeges are also defined here which will be used later)
  //Here, for each processor, it goes over all the locally_owned cells and try to find the global dofs for Vertices
  //and Edges. Since it iis for each processor, we call it local_globalDOF.
  //Later one, we merged all the numbers for all processors and get the global and transfer it to all processors to have it.
  //It is important to know that it is assumed here that the sample is a cuboid with x->[0,userInputs.span[0]]
  //y->[0,userInputs.span[1]], and z->[0,userInputs.span[2]].
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
              node=supportPoints[globalDOF];
              flagVertices=0;

              if ((node[0] <= externalMeshParameterBCs(0))&&(node[1] <= externalMeshParameterBCs(1))&&(node[2] <= externalMeshParameterBCs(2))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V000_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V000_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V000_3=globalDOF;
                }
              }

              if ((node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node[1] <= externalMeshParameterBCs(1))&&(node[2] <= externalMeshParameterBCs(2))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V100_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V100_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V100_3=globalDOF;
                }
              }

              if ((node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))&&(node[0] <= externalMeshParameterBCs(0))&&(node[2] <= externalMeshParameterBCs(2))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V010_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V010_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V010_3=globalDOF;
                }
              }

              if ((node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))&&(node[0] <= externalMeshParameterBCs(0))&&(node[1] <= externalMeshParameterBCs(1))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V001_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V001_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V001_3=globalDOF;
                }
              }

              if ((node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))&&(node[2] <= externalMeshParameterBCs(2))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V110_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V110_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V110_3=globalDOF;
                }
              }

              if ((node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))&&(node[1] <= externalMeshParameterBCs(1))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V101_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V101_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V101_3=globalDOF;
                }
              }

              if ((node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))&&(node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))&&(node[0] <= externalMeshParameterBCs(0))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V011_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V011_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V011_3=globalDOF;
                }
              }

              if ((node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))&&(node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))){
                flagVertices=1;
                if (dof==0){
                  local_globalDOF_V111_1=globalDOF;
                }
                if (dof==1){
                  local_globalDOF_V111_2=globalDOF;
                }
                if (dof==2){
                  local_globalDOF_V111_3=globalDOF;
                }
              }
              if (flagVertices==0){
                if ((((node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))||(node[0] <= externalMeshParameterBCs(0)))&&
                ((node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))||(node[1] <= externalMeshParameterBCs(1))))||
                (((node[0] >= (userInputs.span[0]-externalMeshParameterBCs(0)))||(node[0] <= externalMeshParameterBCs(0)))&&
                ((node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))||(node[2] <= externalMeshParameterBCs(2))))||
                (((node[2] >= (userInputs.span[2]-externalMeshParameterBCs(2)))||(node[2] <= externalMeshParameterBCs(2)))&&
                ((node[1] >= (userInputs.span[1]-externalMeshParameterBCs(1)))||(node[1] <= externalMeshParameterBCs(1))))){
                  if (dof==0){
                    Edges_DOFs_1.add_index(globalDOF);
                  }
                  if (dof==1){
                    Edges_DOFs_2.add_index(globalDOF);
                  }
                  if (dof==2){
                    Edges_DOFs_3.add_index(globalDOF);
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  //Here, it goes over all processors. The processors which does not include
  //any of the globalDOF, they initialized zero. accordingly, the max among all processors
  //will give the real global_DOF.
  globalDOF_V000_1=Utilities::MPI::max(local_globalDOF_V000_1,mpi_communicator);
  globalDOF_V000_2=Utilities::MPI::max(local_globalDOF_V000_2,mpi_communicator);
  globalDOF_V000_3=Utilities::MPI::max(local_globalDOF_V000_3,mpi_communicator);
  globalDOF_V001_1=Utilities::MPI::max(local_globalDOF_V001_1,mpi_communicator);
  globalDOF_V001_2=Utilities::MPI::max(local_globalDOF_V001_2,mpi_communicator);
  globalDOF_V001_3=Utilities::MPI::max(local_globalDOF_V001_3,mpi_communicator);
  globalDOF_V010_1=Utilities::MPI::max(local_globalDOF_V010_1,mpi_communicator);
  globalDOF_V010_2=Utilities::MPI::max(local_globalDOF_V010_2,mpi_communicator);
  globalDOF_V010_3=Utilities::MPI::max(local_globalDOF_V010_3,mpi_communicator);
  globalDOF_V100_1=Utilities::MPI::max(local_globalDOF_V100_1,mpi_communicator);
  globalDOF_V100_2=Utilities::MPI::max(local_globalDOF_V100_2,mpi_communicator);
  globalDOF_V100_3=Utilities::MPI::max(local_globalDOF_V100_3,mpi_communicator);

  globalDOF_V011_1=Utilities::MPI::max(local_globalDOF_V011_1,mpi_communicator);
  globalDOF_V011_2=Utilities::MPI::max(local_globalDOF_V011_2,mpi_communicator);
  globalDOF_V011_3=Utilities::MPI::max(local_globalDOF_V011_3,mpi_communicator);
  globalDOF_V110_1=Utilities::MPI::max(local_globalDOF_V110_1,mpi_communicator);
  globalDOF_V110_2=Utilities::MPI::max(local_globalDOF_V110_2,mpi_communicator);
  globalDOF_V110_3=Utilities::MPI::max(local_globalDOF_V110_3,mpi_communicator);
  globalDOF_V101_1=Utilities::MPI::max(local_globalDOF_V101_1,mpi_communicator);
  globalDOF_V101_2=Utilities::MPI::max(local_globalDOF_V101_2,mpi_communicator);
  globalDOF_V101_3=Utilities::MPI::max(local_globalDOF_V101_3,mpi_communicator);
  globalDOF_V111_1=Utilities::MPI::max(local_globalDOF_V111_1,mpi_communicator);
  globalDOF_V111_2=Utilities::MPI::max(local_globalDOF_V111_2,mpi_communicator);
  globalDOF_V111_3=Utilities::MPI::max(local_globalDOF_V111_3,mpi_communicator);

//vertices_DOFs: This index set has the global dofs for vertices.
  vertices_DOFs.add_index(globalDOF_V000_1);vertices_DOFs.add_index(globalDOF_V000_2);vertices_DOFs.add_index(globalDOF_V000_3);
  vertices_DOFs.add_index(globalDOF_V100_1);vertices_DOFs.add_index(globalDOF_V100_2);vertices_DOFs.add_index(globalDOF_V100_3);
  vertices_DOFs.add_index(globalDOF_V010_1);vertices_DOFs.add_index(globalDOF_V010_2);vertices_DOFs.add_index(globalDOF_V010_3);
  vertices_DOFs.add_index(globalDOF_V001_1);vertices_DOFs.add_index(globalDOF_V001_2);vertices_DOFs.add_index(globalDOF_V001_3);
  vertices_DOFs.add_index(globalDOF_V110_1);vertices_DOFs.add_index(globalDOF_V110_2);vertices_DOFs.add_index(globalDOF_V110_3);
  vertices_DOFs.add_index(globalDOF_V101_1);vertices_DOFs.add_index(globalDOF_V101_2);vertices_DOFs.add_index(globalDOF_V101_3);
  vertices_DOFs.add_index(globalDOF_V011_1);vertices_DOFs.add_index(globalDOF_V011_2);vertices_DOFs.add_index(globalDOF_V011_3);
  vertices_DOFs.add_index(globalDOF_V111_1);vertices_DOFs.add_index(globalDOF_V111_2);vertices_DOFs.add_index(globalDOF_V111_3);

//vertices_DOFs_vector: This vector has the the global dofs for vertices.
  vertices_DOFs_vector[0]=globalDOF_V000_1;vertices_DOFs_vector[1]=globalDOF_V000_2;vertices_DOFs_vector[2]=globalDOF_V000_3;
  vertices_DOFs_vector[3]=globalDOF_V100_1;vertices_DOFs_vector[4]=globalDOF_V100_2;vertices_DOFs_vector[5]=globalDOF_V100_3;
  vertices_DOFs_vector[6]=globalDOF_V010_1;vertices_DOFs_vector[7]=globalDOF_V010_2;vertices_DOFs_vector[8]=globalDOF_V010_3;
  vertices_DOFs_vector[9]=globalDOF_V001_1;vertices_DOFs_vector[10]=globalDOF_V001_2;vertices_DOFs_vector[11]=globalDOF_V001_3;
  vertices_DOFs_vector[12]=globalDOF_V110_1;vertices_DOFs_vector[13]=globalDOF_V110_2;vertices_DOFs_vector[14]=globalDOF_V110_3;
  vertices_DOFs_vector[15]=globalDOF_V101_1;vertices_DOFs_vector[16]=globalDOF_V101_2;vertices_DOFs_vector[17]=globalDOF_V101_3;
  vertices_DOFs_vector[18]=globalDOF_V011_1;vertices_DOFs_vector[19]=globalDOF_V011_2;vertices_DOFs_vector[20]=globalDOF_V011_3;
  vertices_DOFs_vector[21]=globalDOF_V111_1;vertices_DOFs_vector[22]=globalDOF_V111_2;vertices_DOFs_vector[23]=globalDOF_V111_3;

///////////////////////Defining global dofs for Vertices End//////////////

///////////////////////Checking the Vertices constraint lines Start//////////////
//Here, we checked the constraint lines, if all except one dof is not defined,
//we obtain that dof and update the vertices_Constraint_Known vector to 1,
//which means all the dofs for this line is defined.
//We'll go over all constraint again and again (numberVerticesConstraint times to be precise)
//to make sure that we doesn't miss anything, because it is possible that one
//dofs which is not defined in line 1,
//is later on defined in line 3 from other dofs,while we counted it as undefined in line 1.
// So we want to make sure that we do not miss anything.

   for(unsigned int i=0;i<numberVerticesConstraint;i++){
    for(unsigned int j=0;j<numberVerticesConstraint;j++){
      check1=0;
      for(unsigned int k=0;k<vertices_Constraints_Matrix[j][0];k++){
        if (periodicBCsInput[0][vertices_Constraints_Matrix[j][k+1]]==1){
          check1=check1+1;
        }
      }
      if (check1==(vertices_Constraints_Matrix[j][0])){
        vertices_Constraint_Known[j]=1;
      }
      if (check1==(vertices_Constraints_Matrix[j][0]-1)){
        check2=0;
        for(unsigned int k=0;k<vertices_Constraints_Matrix[j][0];k++){
          if (periodicBCsInput[0][vertices_Constraints_Matrix[j][k+1]]==1){
            check2=check2+periodicBCsInput2_Orig[0][vertices_Constraints_Matrix[j][k+1]]*vertices_Constraints_Coef[j][k+1];
          }
          if (periodicBCsInput[0][vertices_Constraints_Matrix[j][k+1]]==0){
            kFlag=k;
          }
        }
        periodicBCsInput[0][vertices_Constraints_Matrix[j][kFlag+1]]=1;
        periodicBCsInput2_Orig[0][vertices_Constraints_Matrix[j][kFlag+1]]=check2*(-1/vertices_Constraints_Coef[j][kFlag+1]);
        vertices_Constraint_Known[j]=1;
      }
    }
  }

  periodicBCsInput2=periodicBCsInput2_Orig;

///////////////////////Checking the Vertices constraint lines End//////////////

//////////////////////Defining the global boundaryLayer2 DOFs Start/////////////////
//The object of this function is to find   global_vector_dof_Boundary_Layer2=(all dofs for elements with a boundary face).
//The first step is to find the local dof_Boundary_Layer2, which is done below.
  cell = dofHandler.begin_active(), endc = dofHandler.end();
  unsigned int CellId=0;
  for (; cell!=endc; ++cell) {
    if ((cell->is_locally_owned())&&(cell->has_boundary_lines())){
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i) {
        unsigned int globalDOF=local_dof_indices[i];
        if ((!dof_Boundary_Layer2.is_element(globalDOF))&&(!vertices_DOFs.is_element(globalDOF))){
          dof_Boundary_Layer2.add_index(globalDOF);
        }
        CellId=CellId+1;
      }
    }
  }

//Here, we search inside the dof of negative faces. If they have not been assigned to dof_Boundary_Layer2
//before, we'll add them.
  if (dof_FN.n_elements()>0){
    IndexSet::ElementIterator dof_FN_counter = dof_FN.begin();
    if (!dof_Boundary_Layer2.is_element(*dof_FN_counter)){
      dof_Boundary_Layer2.add_index(*dof_FN_counter);
    }
    for(unsigned int j = 1; j < dof_FN.n_elements(); j++){
      dof_FN_counter++;
      if (!dof_Boundary_Layer2.is_element(*dof_FN_counter)){
        dof_Boundary_Layer2.add_index(*dof_FN_counter);
      }
    }
  }

  //Here, we gather the local dof_Boundary_Layer2 for all processors as a vector dof_Boundary_Layer2_vector which
  //can be accessed by the processor number 0.
  std::vector<IndexSet> dof_Boundary_Layer2_vector=Utilities::MPI::all_gather(mpi_communicator,dof_Boundary_Layer2);
  n_processes=Utilities::MPI::n_mpi_processes(mpi_communicator);
  global_dof_Boundary_Layer2=dof_Boundary_Layer2_vector[0];

  //Here, we go over all components of dof_Boundary_Layer2_vector, which is local dof_Boundary_Layer2 for all processors
  //and add them to a global_dof_Boundary_Layer2.
  for (unsigned int i=1; i<n_processes; ++i) {
    if (dof_Boundary_Layer2_vector[i].n_elements()>0){
      IndexSet::ElementIterator dofs_boundary = dof_Boundary_Layer2_vector[i].begin();
      if (!global_dof_Boundary_Layer2.is_element(*dofs_boundary)){
        global_dof_Boundary_Layer2.add_index(*dofs_boundary);
      }
      for(unsigned int j = 1; j < dof_Boundary_Layer2_vector[i].n_elements(); j++){
        dofs_boundary++;
        if (!global_dof_Boundary_Layer2.is_element(*dofs_boundary)){
          global_dof_Boundary_Layer2.add_index(*dofs_boundary);
        }
      }
    }
  }

  //Here, we find out the number of dofs for global_dof_Boundary_Layer2 by getting the maximum of the length
  //of global_dof_Boundary_Layer2 over all processors. Only processor rank 0 has a number and all the other ones will end up o.
  //Accordingly, the max becomes equla to the dofs for global_dof_Boundary_Layer2 of processor 0.
  size_dof_Boundary_Layer2=global_dof_Boundary_Layer2.n_elements();
  global_size_dof_Boundary_Layer2=Utilities::MPI::max(size_dof_Boundary_Layer2,mpi_communicator);
  std::vector<unsigned int> vector_dof_Boundary_Layer2;
  vector_dof_Boundary_Layer2.resize(global_size_dof_Boundary_Layer2,0);
  global_vector_dof_Boundary_Layer2.resize(global_size_dof_Boundary_Layer2,0);

  if (size_dof_Boundary_Layer2>0){
    IndexSet::ElementIterator global_dofs_boundary = global_dof_Boundary_Layer2.begin();
    vector_dof_Boundary_Layer2[0]=*global_dofs_boundary;
    for(unsigned int j = 1; j < size_dof_Boundary_Layer2; j++){
      global_dofs_boundary++;
      vector_dof_Boundary_Layer2[j]=*global_dofs_boundary;
    }
  }
//Finally, we copy all components of global_size_dof_Boundary_Layer2 from processor rank 0 to all other processors by applying
//the max operation. One should not for all other processors, the value is 0, after applying
//max operator, the value of global_vector_dof_Boundary_Layer2 will be the same for all processors which is equal to the
// dof_Boundary_Layer2 gathered from all processors.
  for(unsigned int j = 0; j < global_size_dof_Boundary_Layer2; j++){
    global_vector_dof_Boundary_Layer2[j]=Utilities::MPI::max(vector_dof_Boundary_Layer2[j],mpi_communicator);
  }

//////////////////////Defining the global boundaryLayer2 DOFs End/////////////////


//////////////////////Defining the global dofs for Edges Start/////////////////
//The object of this part is to find  global_Edges_DOFs_Vector_Array=(all dofs for on the edges).
//The first step is to find the local dof_Boundary_Layer2, which is done below.

//Edges_DOFs_1 is a indexSet we saved befor, which includes allthe dof=0 on the all edges.
//Edges_DOFs_2 is the same for dof=1. Edges_DOFs_3 is the same for dof=2.
  local_num_Edges_nodes=Edges_DOFs_1.n_elements();
  if (local_num_Edges_nodes>0){
//Edges_DOFs_Vector: is a vector which gather the dofs on the edge for the local processor for three dof of 0, 1, 2 (in 3D).
//Edges_DOFs_Coord_Vector: is a vector that include the x,y, and z of the edge nodes for the local processor.
    Edges_DOFs_Vector.resize(3,std::vector<unsigned int>(local_num_Edges_nodes,0));
    Edges_DOFs_Coord_Vector.resize(3,std::vector<double>(local_num_Edges_nodes,0));

    IndexSet::ElementIterator Edges_DOFs_1_counter = Edges_DOFs_1.begin();
    IndexSet::ElementIterator Edges_DOFs_2_counter = Edges_DOFs_2.begin();
    IndexSet::ElementIterator Edges_DOFs_3_counter = Edges_DOFs_3.begin();

    Edges_DOFs_Vector[0][0]=*Edges_DOFs_1_counter;
    Edges_DOFs_Vector[1][0]=*Edges_DOFs_2_counter;
    Edges_DOFs_Vector[2][0]=*Edges_DOFs_3_counter;
    node=supportPoints[*Edges_DOFs_1_counter];
    Edges_DOFs_Coord_Vector[0][0]=node[0];
    Edges_DOFs_Coord_Vector[1][0]=node[1];
    Edges_DOFs_Coord_Vector[2][0]=node[2];

    for(unsigned int i=1;i<local_num_Edges_nodes;i++){
      Edges_DOFs_1_counter++;Edges_DOFs_2_counter++;Edges_DOFs_3_counter++;
      Edges_DOFs_Vector[0][i]=*Edges_DOFs_1_counter;
      Edges_DOFs_Vector[1][i]=*Edges_DOFs_2_counter;
      Edges_DOFs_Vector[2][i]=*Edges_DOFs_3_counter;
      node=supportPoints[*Edges_DOFs_1_counter];
      Edges_DOFs_Coord_Vector[0][i]=node[0];
      Edges_DOFs_Coord_Vector[1][i]=node[1];
      Edges_DOFs_Coord_Vector[2][i]=node[2];
    }
  }
  //Here, we gather the local Edges_DOFs_Vector and Edges_DOFs_Coord_Vector for all processors as a vector Edges_DOFs_Vector_Array and
  // Edges_DOFs_Coord_Vector_Array which can be accessed by the processor number 0.
  std::vector<std::vector<std::vector<unsigned int>>> Edges_DOFs_Vector_Array=Utilities::MPI::all_gather(mpi_communicator,Edges_DOFs_Vector);
  std::vector<std::vector<std::vector<double>>> Edges_DOFs_Coord_Vector_Array=Utilities::MPI::all_gather(mpi_communicator,Edges_DOFs_Coord_Vector);

  //Here, it is assumed that we are dealing with a cube with equal descritization in X, Y, and Z (totalNumEachEdgesNodes_X=totalNumEachEdgesNodes_Y=totalNumEachEdgesNodes_Z)
  //Each cube has 12 edges->TotalNumEdgesNodes=12*totalNumEachEdgesNodes;
  global_Edges_DOFs_Vector.resize(3,std::vector<unsigned int>(totalNumEachEdgesNodes*12,0));
  global_Edges_DOFs_Coord_Vector.resize(3,std::vector<double>(totalNumEachEdgesNodes*12,0));
  local_Edges_DOFs_Vector=Edges_DOFs_Vector_Array[0];
  local_Edges_DOFs_Coord_Vector=Edges_DOFs_Coord_Vector_Array[0];
  total_num_edges_node_processor=local_Edges_DOFs_Vector[0].size();
  if (total_num_edges_node_processor>0){
    for (unsigned int i=0; i<total_num_edges_node_processor; ++i) {
      global_Edges_DOFs_Vector[0][i]=local_Edges_DOFs_Vector[0][i];
      global_Edges_DOFs_Vector[1][i]=local_Edges_DOFs_Vector[1][i];
      global_Edges_DOFs_Vector[2][i]=local_Edges_DOFs_Vector[2][i];
      global_Edges_DOFs_Coord_Vector[0][i]=local_Edges_DOFs_Coord_Vector[0][i];
      global_Edges_DOFs_Coord_Vector[1][i]=local_Edges_DOFs_Coord_Vector[1][i];
      global_Edges_DOFs_Coord_Vector[2][i]=local_Edges_DOFs_Coord_Vector[2][i];
    }
  }

  local_Edges_DOFs_Vector.clear();
  local_Edges_DOFs_Coord_Vector.clear();

  unsigned int num_edges_node_processor,local_edge_dof;
  for (unsigned int i=1; i<n_processes; ++i) {
    if (!Edges_DOFs_Vector_Array[i].empty()){
      local_Edges_DOFs_Vector=Edges_DOFs_Vector_Array[i];
      local_Edges_DOFs_Coord_Vector=Edges_DOFs_Coord_Vector_Array[i];
      num_edges_node_processor=local_Edges_DOFs_Vector[0].size();
      if (num_edges_node_processor>0){
        for(unsigned int j = 0; j < num_edges_node_processor; j++){
          local_edge_dof=local_Edges_DOFs_Vector[0][j];
          if (std::count(global_Edges_DOFs_Vector[0].begin(), global_Edges_DOFs_Vector[0].end(), local_edge_dof)==0){
            global_Edges_DOFs_Vector[0][total_num_edges_node_processor]=local_Edges_DOFs_Vector[0][j];
            global_Edges_DOFs_Vector[1][total_num_edges_node_processor]=local_Edges_DOFs_Vector[1][j];
            global_Edges_DOFs_Vector[2][total_num_edges_node_processor]=local_Edges_DOFs_Vector[2][j];
            global_Edges_DOFs_Coord_Vector[0][total_num_edges_node_processor]=local_Edges_DOFs_Coord_Vector[0][j];
            global_Edges_DOFs_Coord_Vector[1][total_num_edges_node_processor]=local_Edges_DOFs_Coord_Vector[1][j];
            global_Edges_DOFs_Coord_Vector[2][total_num_edges_node_processor]=local_Edges_DOFs_Coord_Vector[2][j];
            total_num_edges_node_processor=total_num_edges_node_processor+1;
          }
        }
      }
      local_Edges_DOFs_Vector.clear();
      local_Edges_DOFs_Coord_Vector.clear();
    }
  }

  size_dof_Edge_processor=global_Edges_DOFs_Vector[0].size();

  local_Edges_DOFs_Vector.resize(3,std::vector<unsigned int>(totalNumEachEdgesNodes*12,0));
  local_Edges_DOFs_Coord_Vector.resize(3,std::vector<double>(totalNumEachEdgesNodes*12,0));

  global_Edges_DOFs_Vector_Array.resize(3,std::vector<unsigned int>(totalNumEachEdgesNodes*12,0));
  global_Edges_DOFs_Coord_Vector_Array.resize(3,std::vector<double>(totalNumEachEdgesNodes*12,0));

  if (size_dof_Edge_processor>0){
    for(unsigned int i = 0; i < size_dof_Edge_processor; i++){
      local_Edges_DOFs_Vector[0][i]=global_Edges_DOFs_Vector[0][i];
      local_Edges_DOFs_Vector[1][i]=global_Edges_DOFs_Vector[1][i];
      local_Edges_DOFs_Vector[2][i]=global_Edges_DOFs_Vector[2][i];
      local_Edges_DOFs_Coord_Vector[0][i]=global_Edges_DOFs_Coord_Vector[0][i];
      local_Edges_DOFs_Coord_Vector[1][i]=global_Edges_DOFs_Coord_Vector[1][i];
      local_Edges_DOFs_Coord_Vector[2][i]=global_Edges_DOFs_Coord_Vector[2][i];
    }
  }

  for(unsigned int i = 0; i < (totalNumEachEdgesNodes*12); i++){

    //Finally, we copy all components of global_Edges_DOFs_Vector_Array and global_Edges_DOFs_Coord_Vector_Array
    // from processor rank 0 to all other processors by applying
    //the max operation. One should not for all other processors, the value is 0, after applying
    //max operator, the value of global_Edges_DOFs_Vector_Array and global_Edges_DOFs_Coord_Vector_Array
    //will be the same for all processors which is equal to the local_Edges_DOFs_Vector and local_Edges_DOFs_Coord_Vector
    // gathered from all processors. Here, it is assumed that all the coordinates are bigger than zero, i.e., x,y,z>=0.
    global_Edges_DOFs_Vector_Array[0][i]=Utilities::MPI::max(local_Edges_DOFs_Vector[0][i],mpi_communicator);
    global_Edges_DOFs_Vector_Array[1][i]=Utilities::MPI::max(local_Edges_DOFs_Vector[1][i],mpi_communicator);
    global_Edges_DOFs_Vector_Array[2][i]=Utilities::MPI::max(local_Edges_DOFs_Vector[2][i],mpi_communicator);

    global_Edges_DOFs_Coord_Vector_Array[0][i]=Utilities::MPI::max(local_Edges_DOFs_Coord_Vector[0][i],mpi_communicator);
    global_Edges_DOFs_Coord_Vector_Array[1][i]=Utilities::MPI::max(local_Edges_DOFs_Coord_Vector[1][i],mpi_communicator);
    global_Edges_DOFs_Coord_Vector_Array[2][i]=Utilities::MPI::max(local_Edges_DOFs_Coord_Vector[2][i],mpi_communicator);
  }
  //////////////////////Defining the global dofs for Edges End/////////////////

  ////////////Defining the global EX000,EX001,EX010,.......... Start///////////
  //The object of this part is to find   the global dofs for EX000,EX001,EX010,
  //EX011, EY000,EY001,EY100,EY101, EZ000,EZ010,EZ100,EZ110 are gathered.
  //Here, we have global_Edges_DOFs_Vector_Array from previous step as the global dofs on edges gathered from all processors.
  //Two steps are required here: First, to define out of all global edges nodes, which ones belong to EX000,EX001,...
  //Second: In each of the EX000, EX001,.... arrangment of nodes based on the corresponding coordinate.
  //For example in EX000 based on the x coordinate. For EY000 based on y coordinate,
  //and for EZ000 based on z coordinate.

// EX000,EX001,.... are defined as a pair. The reason is we want to be able to arrange the nodes in them based on the coordinates.
//Here, the second component of pair is dof0 of node, and the first component is the required coordinate for arranging the nodes.
//We only saved the first dof (dof-0) from each node, because we know that dof=1 and dof=2 are (dof=1)=(dof=0)+1 and (dof=2)=(dof=0)+2.
  std::vector< std::pair <double,int> > EX000,EX001,EX010,EX011;
  std::vector< std::pair <double,int> > EY000,EY001,EY100,EY101;
  std::vector< std::pair <double,int> > EZ000,EZ010,EZ100,EZ110;

//Here, we find the dof0 belong to the nodes of each EX000,EX001,.....
  for(unsigned int i = 0; i < (totalNumEachEdgesNodes*12); i++){

    node0=global_Edges_DOFs_Coord_Vector_Array[0][i];
    node1=global_Edges_DOFs_Coord_Vector_Array[1][i];
    node2=global_Edges_DOFs_Coord_Vector_Array[2][i];

    if ((node1 <= externalMeshParameterBCs(1))&&(node2 <= externalMeshParameterBCs(2))){
      EX000.push_back(std::make_pair(node0,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node1 <= externalMeshParameterBCs(1))&&(node2 >= (userInputs.span[2]-externalMeshParameterBCs(2)))){
      EX001.push_back(std::make_pair(node0,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node1 >= (userInputs.span[1]-externalMeshParameterBCs(1)))&&(node2 <= externalMeshParameterBCs(2))){
      EX010.push_back(std::make_pair(node0,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node1 >= (userInputs.span[1]-externalMeshParameterBCs(1)))&&(node2 >= (userInputs.span[2]-externalMeshParameterBCs(2)))){
      EX011.push_back(std::make_pair(node0,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 <= externalMeshParameterBCs(0))&&(node2 <= externalMeshParameterBCs(2))){
      EY000.push_back(std::make_pair(node1,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 <= externalMeshParameterBCs(0))&&(node2 >= (userInputs.span[2]-externalMeshParameterBCs(2)))){
      EY001.push_back(std::make_pair(node1,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node2 <= externalMeshParameterBCs(2))){
      EY100.push_back(std::make_pair(node1,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node2 >= (userInputs.span[2]-externalMeshParameterBCs(2)))){
      EY101.push_back(std::make_pair(node1,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 <= externalMeshParameterBCs(0))&&(node1 <= externalMeshParameterBCs(1))){
      EZ000.push_back(std::make_pair(node2,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 <= externalMeshParameterBCs(0))&&(node1 >= (userInputs.span[1]-externalMeshParameterBCs(1)))){
      EZ010.push_back(std::make_pair(node2,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node1 <= externalMeshParameterBCs(1))){
      EZ100.push_back(std::make_pair(node2,global_Edges_DOFs_Vector_Array[0][i]));
    }
    if ((node0 >= (userInputs.span[0]-externalMeshParameterBCs(0)))&&(node1 >= (userInputs.span[1]-externalMeshParameterBCs(1)))){
      EZ110.push_back(std::make_pair(node2,global_Edges_DOFs_Vector_Array[0][i]));
    }
  }

  //Here, we sort the dof0 (second compondent of the pair) based on the cooresponding coordinates (first component of the pair)
  std::sort(EX000.begin(), EX000.end());
  std::sort(EX001.begin(), EX001.end());
  std::sort(EX010.begin(), EX010.end());
  std::sort(EX011.begin(), EX011.end());
  std::sort(EY000.begin(), EY000.end());
  std::sort(EY001.begin(), EY001.end());
  std::sort(EY100.begin(), EY100.end());
  std::sort(EY101.begin(), EY101.end());
  std::sort(EZ000.begin(), EZ000.end());
  std::sort(EZ010.begin(), EZ010.end());
  std::sort(EZ100.begin(), EZ100.end());
  std::sort(EZ110.begin(), EZ110.end());

//Here, we assume that for each node, dof=1 and dof=2 are (dof=1)=(dof=0)+1 and (dof=2)=(dof=0)+2.
  for(unsigned int i = 0; i < totalNumEachEdgesNodes; i++){
    edges_DOFs[0][i]=EX000[i].second;
    edges_DOFs[1][i]=EX000[i].second+1;
    edges_DOFs[2][i]=EX000[i].second+2;

    edges_DOFs[3][i]=EX001[i].second;
    edges_DOFs[4][i]=EX001[i].second+1;
    edges_DOFs[5][i]=EX001[i].second+2;

    edges_DOFs[6][i]=EX010[i].second;
    edges_DOFs[7][i]=EX010[i].second+1;
    edges_DOFs[8][i]=EX010[i].second+2;

    edges_DOFs[9][i]=EX011[i].second;
    edges_DOFs[10][i]=EX011[i].second+1;
    edges_DOFs[11][i]=EX011[i].second+2;

    edges_DOFs[12][i]=EY000[i].second;
    edges_DOFs[13][i]=EY000[i].second+1;
    edges_DOFs[14][i]=EY000[i].second+2;

    edges_DOFs[15][i]=EY001[i].second;
    edges_DOFs[16][i]=EY001[i].second+1;
    edges_DOFs[17][i]=EY001[i].second+2;

    edges_DOFs[18][i]=EY100[i].second;
    edges_DOFs[19][i]=EY100[i].second+1;
    edges_DOFs[20][i]=EY100[i].second+2;

    edges_DOFs[21][i]=EY101[i].second;
    edges_DOFs[22][i]=EY101[i].second+1;
    edges_DOFs[23][i]=EY101[i].second+2;

    edges_DOFs[24][i]=EZ000[i].second;
    edges_DOFs[25][i]=EZ000[i].second+1;
    edges_DOFs[26][i]=EZ000[i].second+2;

    edges_DOFs[27][i]=EZ010[i].second;
    edges_DOFs[28][i]=EZ010[i].second+1;
    edges_DOFs[29][i]=EZ010[i].second+2;

    edges_DOFs[30][i]=EZ100[i].second;
    edges_DOFs[31][i]=EZ100[i].second+1;
    edges_DOFs[32][i]=EZ100[i].second+2;

    edges_DOFs[33][i]=EZ110[i].second;
    edges_DOFs[34][i]=EZ110[i].second+1;
    edges_DOFs[35][i]=EZ110[i].second+2;
  }

  ////////////Defining the global EX000,EX001,EX010,.......... End///////////
}

template <int dim>
#if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
void ellipticBVP<dim>::setNodeConstraints(ConstraintMatrix& constraintmatrix){
#else
void ellipticBVP<dim>::setNodeConstraints(AffineConstraints<double>& constraintmatrix){
#endif
//This function define the required Vertices constraints for increment 0 of each timestep.

//The reason that increment 0 (Inc0) and the rest of increaments (IncNot0) in the case of constraints is
//at Inc0, the nonhomegenous constraints nonequal to zero, if it
//is available based on the constraints, should be applied as the external load.
//However, in the following increments (IncNot0), they are just for equilibrating the residuals
//and no nonhomogenous constraints should be added.

  unsigned int globalDOF1,globalDOF2;
  double globalDOF1_Coef;
  double inhomogeneity_Total;
  for(unsigned int j = 0; j < numberVerticesConstraint; j++){
    if (vertices_Constraint_Known[j]==0){
      globalDOF1=vertices_DOFs_vector[vertices_Constraints_Matrix[j][1]];
      // pcout << "globalDOF1 " << globalDOF1 << " j number " << j << std::endl;
      if (locally_relevant_dofs_Mod.is_element(globalDOF1)){

        globalDOF1_Coef=vertices_Constraints_Coef[j][1];
        constraintmatrix.add_line (globalDOF1);
        inhomogeneity_Total=0;
        for(unsigned int k=1;k<vertices_Constraints_Matrix[j][0];k++){
          globalDOF2=vertices_DOFs_vector[vertices_Constraints_Matrix[j][k+1]];
          // pcout << " j " << j << " vertices_Constraints_Matrix " << vertices_Constraints_Matrix[j][k+1] << " check2 " << periodicBCsInput[0][vertices_Constraints_Matrix[j][k+1]] << " globalDOF2 " << globalDOF2 << std::endl;
          if (periodicBCsInput[0][vertices_Constraints_Matrix[j][k+1]]==1){
            inhomogeneity_Total=inhomogeneity_Total+vertices_Constraints_Coef[j][k+1]*periodicBCsInput2[0][vertices_Constraints_Matrix[j][k+1]];
          }
          else {
            constraintmatrix.add_entry(globalDOF1,globalDOF2, (-1/globalDOF1_Coef)*vertices_Constraints_Coef[j][k+1]);
          }

        }
        inhomogeneity_Total=(-1/globalDOF1_Coef)*inhomogeneity_Total;

        constraintmatrix.set_inhomogeneity(globalDOF1, inhomogeneity_Total);
      }
    }
  }

  for(unsigned int i=0;i<totalNumVerticesDOFs;i++){
    if (periodicBCsInput[0][i]==1){
      globalDOF1=vertices_DOFs_vector[i];
      if (locally_relevant_dofs_Mod.is_element(globalDOF1)){
        constraintmatrix.add_line (globalDOF1);
        constraintmatrix.set_inhomogeneity(globalDOF1, periodicBCsInput2[0][i]);
      }
    }
  }

}

template <int dim>
#if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
void ellipticBVP<dim>::setEdgeConstraints(ConstraintMatrix& constraintmatrix){
#else
void ellipticBVP<dim>::setEdgeConstraints(AffineConstraints<double>& constraintmatrix){
#endif
  //This function define the required Edges constraints for increment 0 of each timestep.

  //The reason that increment 0 (Inc0) and the rest of increaments (IncNot0) in the case of constraints is
  //at Inc0, the nonhomegenous constraints nonequal to zero, if it
  //is available based on the constraints, should be applied as the external load.
  //However, in the following increments (IncNot0), they are just for equilibrating the residuals
  //and no nonhomogenous constraints should be added.

  unsigned int globalVerticesDOF1,globalVerticesDOF2;
  unsigned int globalEdgesDOF1,globalEdgesDOF2;
  double globalVerticesDOF1_Coef,globalVerticesDOF2_Coef,globalEdgesDOF1_Coef,globalEdgesDOF2_Coef;
  double inhomogeneity_Total=0;
  unsigned int FlagDOF1,FlagDOF2,case_DOF;

  for(unsigned int j = 0; j < numberEdgesConstraint; j++){
    FlagDOF1=0;
    FlagDOF2=0;
    case_DOF=0;
    globalVerticesDOF1=vertices_DOFs_vector[edges_Constraints_Matrix[j][2]];
    globalVerticesDOF2=vertices_DOFs_vector[edges_Constraints_Matrix[j][3]];
    globalVerticesDOF1_Coef=edges_Constraints_Coef[j][2];
    globalVerticesDOF2_Coef=edges_Constraints_Coef[j][3];
    globalEdgesDOF1_Coef=edges_Constraints_Coef[j][0];
    globalEdgesDOF2_Coef=edges_Constraints_Coef[j][1];

    if (periodicBCsInput[0][edges_Constraints_Matrix[j][2]]==1){
      FlagDOF1=1;
    }

    if (periodicBCsInput[0][edges_Constraints_Matrix[j][3]]==1){
      FlagDOF2=1;
    }
    inhomogeneity_Total=0;
    if ((FlagDOF1==1)&&(FlagDOF2==1)){
      inhomogeneity_Total=(-1/globalEdgesDOF1_Coef)*globalVerticesDOF1_Coef*periodicBCsInput2[0][edges_Constraints_Matrix[j][2]]+
      (-1/globalEdgesDOF1_Coef)*globalVerticesDOF2_Coef*periodicBCsInput2[0][edges_Constraints_Matrix[j][3]];
      case_DOF=1;
    }
    else if ((FlagDOF1==1)&&(FlagDOF2==0)) {
      inhomogeneity_Total=(-1/globalEdgesDOF1_Coef)*globalVerticesDOF1_Coef*periodicBCsInput2[0][edges_Constraints_Matrix[j][2]];
      case_DOF=2;
    }
    else if ((FlagDOF1==0)&&(FlagDOF2==1)) {
      inhomogeneity_Total=(-1/globalEdgesDOF1_Coef)*globalVerticesDOF2_Coef*periodicBCsInput2[0][edges_Constraints_Matrix[j][3]];
      case_DOF=3;
    }
    else {
      case_DOF=4;
    }

    for(unsigned int i = 0; i < totalNumEachEdgesNodes; i++){
      globalEdgesDOF1=edges_DOFs[edges_Constraints_Matrix[j][0]][i];
      globalEdgesDOF2=edges_DOFs[edges_Constraints_Matrix[j][1]][i];
      // pcout << "globalDOF1 " << globalDOF1 << " j number " << j << std::endl;
      if (locally_relevant_dofs_Mod.is_element(globalEdgesDOF1)){
        constraintmatrix.add_line (globalEdgesDOF1);
        constraintmatrix.add_entry(globalEdgesDOF1,globalEdgesDOF2, (-1/globalEdgesDOF1_Coef)*globalEdgesDOF2_Coef);
        if (case_DOF==1){
          constraintmatrix.set_inhomogeneity(globalEdgesDOF1, inhomogeneity_Total);
        }
        else if (case_DOF==2){
          constraintmatrix.set_inhomogeneity(globalEdgesDOF1, inhomogeneity_Total);
          constraintmatrix.add_entry(globalEdgesDOF1,globalVerticesDOF2, (-1/globalEdgesDOF1_Coef)*globalVerticesDOF2_Coef);
        }
        else if (case_DOF==3){
          constraintmatrix.set_inhomogeneity(globalEdgesDOF1, inhomogeneity_Total);
          constraintmatrix.add_entry(globalEdgesDOF1,globalVerticesDOF1, (-1/globalEdgesDOF1_Coef)*globalVerticesDOF1_Coef);
        }
        else if (case_DOF==4){
          constraintmatrix.add_entry(globalEdgesDOF1,globalVerticesDOF2, (-1/globalEdgesDOF1_Coef)*globalVerticesDOF2_Coef);
          constraintmatrix.add_entry(globalEdgesDOF1,globalVerticesDOF1, (-1/globalEdgesDOF1_Coef)*globalVerticesDOF1_Coef);
        }
      }
    }
  }

}

template <int dim>
#if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 1)&&(DEAL_II_VERSION_MAJOR==9)))
void ellipticBVP<dim>::setFaceConstraints(ConstraintMatrix& constraintmatrix){
#else
void ellipticBVP<dim>::setFaceConstraints(AffineConstraints<double>& constraintmatrix){
#endif

  //This function define the required Faces constraints for increment 0 of each timestep.
  ///Here, we assumed that in the face constraints lines, we always start with FXP, FYP, or FZP.

  //The reason that increment 0 (Inc0) and the rest of increaments (IncNot0) in the case of constraints is
  //at Inc0, the nonhomegenous constraints nonequal to zero, if it
  //is available based on the constraints, should be applied as the external load.
  //However, in the following increments (IncNot0), they are just for equilibrating the residuals
  //and no nonhomogenous constraints should be added.

  IndexSet currentIndexSet;
  unsigned int nb_dofs_CurrentFace;
  unsigned int globalVerticesDOF1,globalVerticesDOF2;
  unsigned int globalFacesDOF1;
  double globalVerticesDOF1_Coef,globalVerticesDOF2_Coef,globalFacesDOF1_Coef,globalFacesDOF2_Coef;
  double inhomogeneity_Total=0;
  unsigned int FlagDOF1,FlagDOF2,case_DOF;

  unsigned int nb_dofs_face_FXN = dof_FXN.n_elements();
  unsigned int nb_dofs_face_FYN = dof_FYN.n_elements();
  unsigned int nb_dofs_face_FZP = dof_FZP.n_elements();


  const FEValuesExtractors::Scalar  displacementX(0);
  const FEValuesExtractors::Scalar  displacementY(1);
  const FEValuesExtractors::Scalar  displacementZ(2);

  std::vector<GridTools::PeriodicFacePair<typename DoFHandler<dim>::cell_iterator> > periodicity_vectorX,periodicity_vectorY,periodicity_vectorZ;

///Here, we assumed that in the face constraints lines, we always start with FXP, FYP, or FZP.
  if (nb_dofs_face_FXN>0){
    GridTools::collect_periodic_faces(dofHandler, /*b_id1*/ 1, /*b_id2*/ 0, /*direction*/ 0, periodicity_vectorX);
  }

  if (nb_dofs_face_FYN>0){
    GridTools::collect_periodic_faces(dofHandler, /*b_id1*/ 3, /*b_id2*/ 2, /*direction*/ 1, periodicity_vectorY);
  }

  if (nb_dofs_face_FZP>0){
    GridTools::collect_periodic_faces(dofHandler, /*b_id1*/ 5, /*b_id2*/ 4, /*direction*/ 2, periodicity_vectorZ);
  }

  // #set Faces Periodic BCs row order:
  // # 0=FXN_1;1=FXN_2;2=FXN_3; 3=FXP_1;4=FXP_2;5=FXP_3; 6=FYN_1;7=FYN_2;8=FYN_3; 9=FYP_1;10=FYP_2;11=FYP_3;
  // # 12=FZN_1;13=FZN_2;14=FZN_3; 15=FZP_1;16=FZP_2;17=FZP_3;

  for (unsigned int i=0; i <numberFacesConstraint;i++){
    currentIndexSet=faces_dof_Index_vector[faces_Constraints_Matrix[i][0]];
    nb_dofs_CurrentFace=currentIndexSet.n_elements();
    if (nb_dofs_CurrentFace>0){
      #if (DEAL_II_VERSION_MAJOR == 9 && DEAL_II_VERSION_MINOR >= 4)
        if (faces_Constraints_Matrix[i][0]==3){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorX, constraintmatrix, FE.component_mask(displacementX));
        }
        if (faces_Constraints_Matrix[i][0]==4){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorX, constraintmatrix, FE.component_mask(displacementY));
        }
        if (faces_Constraints_Matrix[i][0]==5){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorX, constraintmatrix, FE.component_mask(displacementZ));
        }
        if (faces_Constraints_Matrix[i][0]==9){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorY, constraintmatrix, FE.component_mask(displacementX));
        }
        if (faces_Constraints_Matrix[i][0]==10){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorY, constraintmatrix, FE.component_mask(displacementY));
        }
        if (faces_Constraints_Matrix[i][0]==11){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorY, constraintmatrix, FE.component_mask(displacementZ));
        }
        if (faces_Constraints_Matrix[i][0]==15){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorZ, constraintmatrix, FE.component_mask(displacementX));
        }
        if (faces_Constraints_Matrix[i][0]==16){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorZ, constraintmatrix, FE.component_mask(displacementY));
        }
        if (faces_Constraints_Matrix[i][0]==17){
          DoFTools::make_periodicity_constraints<dim, dim>(periodicity_vectorZ, constraintmatrix, FE.component_mask(displacementZ));
        }
      #else
        if (faces_Constraints_Matrix[i][0]==3){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorX, constraintmatrix, FE.component_mask(displacementX));
        }
        if (faces_Constraints_Matrix[i][0]==4){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorX, constraintmatrix, FE.component_mask(displacementY));
        }
        if (faces_Constraints_Matrix[i][0]==5){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorX, constraintmatrix, FE.component_mask(displacementZ));
        }
        if (faces_Constraints_Matrix[i][0]==9){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorY, constraintmatrix, FE.component_mask(displacementX));
        }
        if (faces_Constraints_Matrix[i][0]==10){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorY, constraintmatrix, FE.component_mask(displacementY));
        }
        if (faces_Constraints_Matrix[i][0]==11){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorY, constraintmatrix, FE.component_mask(displacementZ));
        }
        if (faces_Constraints_Matrix[i][0]==15){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorZ, constraintmatrix, FE.component_mask(displacementX));
        }
        if (faces_Constraints_Matrix[i][0]==16){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorZ, constraintmatrix, FE.component_mask(displacementY));
        }
        if (faces_Constraints_Matrix[i][0]==17){
          DoFTools::make_periodicity_constraints<DoFHandler<dim> >(periodicity_vectorZ, constraintmatrix, FE.component_mask(displacementZ));
        }
      #endif
    }
    currentIndexSet.clear();
  }

  for (unsigned int i=0; i <numberFacesConstraint;i++){
    FlagDOF1=0;
    FlagDOF2=0;
    case_DOF=0;
    globalVerticesDOF1=vertices_DOFs_vector[faces_Constraints_Matrix[i][2]];
    globalVerticesDOF2=vertices_DOFs_vector[faces_Constraints_Matrix[i][3]];
    globalVerticesDOF1_Coef=faces_Constraints_Coef[i][2];
    globalVerticesDOF2_Coef=faces_Constraints_Coef[i][3];
    globalFacesDOF1_Coef=faces_Constraints_Coef[i][0];

    if (periodicBCsInput[0][faces_Constraints_Matrix[i][2]]==1){
      FlagDOF1=1;
    }

    if (periodicBCsInput[0][faces_Constraints_Matrix[i][3]]==1){
      FlagDOF2=1;
    }
    inhomogeneity_Total=0;
    if ((FlagDOF1==1)&&(FlagDOF2==1)){
      inhomogeneity_Total=(-1/globalFacesDOF1_Coef)*globalVerticesDOF1_Coef*periodicBCsInput2[0][faces_Constraints_Matrix[i][2]]+
      (-1/globalFacesDOF1_Coef)*globalVerticesDOF2_Coef*periodicBCsInput2[0][faces_Constraints_Matrix[i][3]];
      case_DOF=1;
    }
    else if ((FlagDOF1==1)&&(FlagDOF2==0)) {
      inhomogeneity_Total=(-1/globalFacesDOF1_Coef)*globalVerticesDOF1_Coef*periodicBCsInput2[0][faces_Constraints_Matrix[i][2]];
      case_DOF=2;
    }
    else if ((FlagDOF1==0)&&(FlagDOF2==1)) {
      inhomogeneity_Total=(-1/globalFacesDOF1_Coef)*globalVerticesDOF2_Coef*periodicBCsInput2[0][faces_Constraints_Matrix[i][3]];
      case_DOF=3;
    }
    else {
      case_DOF=4;
    }


    currentIndexSet=faces_dof_Index_vector[faces_Constraints_Matrix[i][0]];
    nb_dofs_CurrentFace=currentIndexSet.n_elements();
    if (nb_dofs_CurrentFace>0){
      IndexSet::ElementIterator dofs_currentFace = currentIndexSet.begin();
      if (constraintmatrix.is_constrained(*dofs_currentFace)){
        if (case_DOF==1){
          constraintmatrix.set_inhomogeneity(*dofs_currentFace, inhomogeneity_Total);
        }
        else if (case_DOF==2){
          constraintmatrix.set_inhomogeneity(*dofs_currentFace, inhomogeneity_Total);
          constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF2, (-1/globalFacesDOF1_Coef)*globalVerticesDOF2_Coef);
        }
        else if (case_DOF==3){
          constraintmatrix.set_inhomogeneity(*dofs_currentFace, inhomogeneity_Total);
          constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF1, (-1/globalFacesDOF1_Coef)*globalVerticesDOF1_Coef);
        }
        else if (case_DOF==4){
          constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF2, (-1/globalFacesDOF1_Coef)*globalVerticesDOF2_Coef);
          constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF1, (-1/globalFacesDOF1_Coef)*globalVerticesDOF1_Coef);
        }
      }

      for(unsigned int j = 1; j < nb_dofs_CurrentFace; j++){
        dofs_currentFace++;
        if (constraintmatrix.is_constrained(*dofs_currentFace)){
          if (case_DOF==1){
            constraintmatrix.set_inhomogeneity(*dofs_currentFace, inhomogeneity_Total);
          }
          else if (case_DOF==2){
            constraintmatrix.set_inhomogeneity(*dofs_currentFace, inhomogeneity_Total);
            constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF2, (-1/globalFacesDOF1_Coef)*globalVerticesDOF2_Coef);
          }
          else if (case_DOF==3){
            constraintmatrix.set_inhomogeneity(*dofs_currentFace, inhomogeneity_Total);
            constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF1, (-1/globalFacesDOF1_Coef)*globalVerticesDOF1_Coef);
          }
          else if (case_DOF==4){
            constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF2, (-1/globalFacesDOF1_Coef)*globalVerticesDOF2_Coef);
            constraintmatrix.add_entry(*dofs_currentFace,globalVerticesDOF1, (-1/globalFacesDOF1_Coef)*globalVerticesDOF1_Coef);
          }
        }
      }
    }
    currentIndexSet.clear();
  }
}

template <int dim>
void ellipticBVP<dim>::setPeriodicityConstraintsInc0(){
//This functions apply all Vertices, Edges, and Faces constraint for increment 0 of each timestep.
//This function also generate the positive loading for increment 0
//of each timestep in the case of loading when TabularPeriodicBCs option is used.

  periodicBCsInput2=periodicBCsInput2_Orig;


//The reason that increment 0 (Inc0) and the rest of increaments (IncNot0) in the case of constraints is
//at Inc0, the nonhomegenous constraints nonequal to zero, if it
//is available based on the constraints, should be applied as the external load.
//However, in the following increments (IncNot0), they are just for equilibrating the residuals
//and no nonhomogenous constraints should be added.
  constraints_PBCs_Inc0.clear();
  constraints_PBCs_Inc0.reinit (locally_relevant_dofs_Mod);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints_PBCs_Inc0);
  setFaceConstraints(constraints_PBCs_Inc0);
  setEdgeConstraints(constraints_PBCs_Inc0);
  setNodeConstraints(constraints_PBCs_Inc0);
  constraints_PBCs_Inc0. close ();
}

template <int dim>
void ellipticBVP<dim>::setPeriodicityConstraintsInc0Neg(){
//This functions apply all Vertices, Edges, and Faces constraint for increment 0
//of each timestep in the case of unloading when TabularPeriodicBCs option is used.

//Here, the BCs increment is multiplied by negative sign
for(unsigned int i=0;i<totalNumVerticesDOFs;i++){
  periodicBCsInput2[0][i]=-periodicBCsInput2_Orig[0][i];
}

//The reason that increment 0 (Inc0) and the rest of increaments (IncNot0) in the case of constraints is
//at Inc0, the nonhomegenous constraints nonequal to zero, if it
//is available based on the constraints, should be applied as the external load.
//However, in the following increments (IncNot0), they are just for equilibrating the residuals
//and no nonhomogenous constraints should be added.
  constraints_PBCs_Inc0Neg.clear();
  constraints_PBCs_Inc0Neg.reinit (locally_relevant_dofs_Mod);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints_PBCs_Inc0Neg);
  setFaceConstraints(constraints_PBCs_Inc0Neg);
  setEdgeConstraints(constraints_PBCs_Inc0Neg);
  setNodeConstraints(constraints_PBCs_Inc0Neg);
  constraints_PBCs_Inc0Neg. close ();
}

// Set constraints to enforce periodic boundary conditions
template <int dim>
void ellipticBVP<dim>::setPeriodicityConstraintsIncNot0(){
  //This functions apply all Vertices, Edges, and Faces constraint for increment after 0 (IncNot0) of each timestep.

  for(unsigned int i=0;i<totalNumVerticesDOFs;i++){
    periodicBCsInput2[0][i]=0;
  }

  //The reason that increment 0 (Inc0) and the rest of increaments (IncNot0) in the case of constraints is
  //at Inc0, the nonhomegenous constraints nonequal to zero, if it
  //is available based on the constraints, should be applied as the external load.
  //However, in the following increments (IncNot0), they are just for equilibrating the residuals
  //and no nonhomogenous constraints should be added.
  constraints_PBCs_IncNot0.clear();
  constraints_PBCs_IncNot0.reinit (locally_relevant_dofs_Mod);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints_PBCs_IncNot0);
  setFaceConstraints(constraints_PBCs_IncNot0);
  setEdgeConstraints(constraints_PBCs_IncNot0);
  setNodeConstraints(constraints_PBCs_IncNot0);
  constraints_PBCs_IncNot0. close ();
}

template <int dim>
void ellipticBVP<dim>::setPeriodicityConstraints(){
  constraints.clear();
  constraints.reinit (locally_relevant_dofs_Mod);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints);
  if (currentIteration==0){
    if(userInputs.enableTabularPeriodicBCs){
      //I added delT/1000 as some small value to make sure we have BCs applied correct.
      currentTime=delT*(currentIncrement+1)-delT/1000;
      if (currentIncrement==0){
        timeCounter=1;
      }
      if (currentTime>userInputs.tabularPeriodicTimeInput[0][timeCounter]){
        timeCounter=timeCounter+1;
      }
      ////////Loading case
      if (userInputs.tabularPeriodicCoef[0][timeCounter-1]==1){
        constraints.copy_from(constraints_PBCs_Inc0);
      }
      ////////Unloading case
      else if (userInputs.tabularPeriodicCoef[0][timeCounter-1]==-1) {
        constraints.copy_from(constraints_PBCs_Inc0Neg);
      }
      ////////Neutral= No Loading or Unloading
      else if (userInputs.tabularPeriodicCoef[0][timeCounter-1]==0){
        constraints.copy_from(constraints_PBCs_IncNot0);
      }
    }
    else {
      constraints.copy_from(constraints_PBCs_Inc0);
    }
  }
  else{
    constraints.copy_from(constraints_PBCs_IncNot0);
  }
}

        #include "../../include/ellipticBVP_template_instantiations.h"
