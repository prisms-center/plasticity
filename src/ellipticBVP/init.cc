//initialization method for ellipticBVP class
#include "../../include/ellipticBVP.h"
#include <fstream>


//initialize all FE objects and data structures
template <int dim>
void ellipticBVP<dim>::init(){
  std::string line;
  std::ifstream bcDataFile;
  unsigned int id;
  double totalU;
  unsigned int faceID,dof;

  for(unsigned int i=0;i<2*dim;i++){
    faceDOFConstrained.push_back({false,false,false});
    deluConstraint.push_back({0.0,0.0,0.0});
  }


  pcout << "number of MPI processes: "
  << Utilities::MPI::n_mpi_processes(mpi_communicator)
  << std::endl;

  //initialize FE objects
  dofHandler.distribute_dofs (FE);
  locally_owned_dofs = dofHandler.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dofHandler, locally_relevant_dofs);

  //In the case of Periodic BCs, locally_relevant_dofs_Mod will be updated to
  //include the required DOFs to implement the periodicity over different processors.
  //If the BCs is not Periodic, locally_relevant_dofs_Mod will remian similar to locally_relevant_dofs.
  locally_relevant_dofs_Mod=locally_relevant_dofs;

  pcout << "number of elements: "
  << triangulation.n_global_active_cells()
  << std::endl
  << "number of degrees of freedom: "
  << dofHandler.n_dofs()
  << std::endl;

  //initialize FE objects for scalar field which will be used for post processing
  dofHandler_Scalar.distribute_dofs (FE_Scalar);
  locally_owned_dofs_Scalar = dofHandler_Scalar.locally_owned_dofs ();
  DoFTools::extract_locally_relevant_dofs (dofHandler_Scalar, locally_relevant_dofs_Scalar);

  if(!userInputs.enablePeriodicBCs){
    //constraints
    constraints.clear ();
    constraints.reinit (locally_relevant_dofs);
    DoFTools::make_hanging_node_constraints (dofHandler, constraints);
    constraints.close ();
  }

  //constraints for mass matrix used in post-processing for projection
  //(i.e. no constraints, as mass matrix needs no dirichlet BCs)
  constraintsMassMatrix.clear ();
  constraintsMassMatrix.reinit (locally_relevant_dofs_Scalar);
  DoFTools::make_hanging_node_constraints (dofHandler_Scalar, constraintsMassMatrix);
  constraintsMassMatrix.close ();

  //get support points (nodes) for this problem
  DoFTools::map_dofs_to_support_points(MappingQ1<dim, dim>(), dofHandler, supportPoints);


  if(userInputs.enablePeriodicBCs){
    numberVerticesConstraint=userInputs.numberVerticesConstraint;
    numberEdgesConstraint=userInputs.numberEdgesConstraint;
    numberFacesConstraint=userInputs.numberFacesConstraint;
    vertices_Constraints_Matrix.resize(numberVerticesConstraint,std::vector<unsigned int>(5,0));
    vertices_Constraints_Coef.resize(numberVerticesConstraint,std::vector<double>(5,0));
    edges_Constraints_Matrix.resize(numberEdgesConstraint,std::vector<unsigned int>(4,0));
    edges_Constraints_Coef.resize(numberEdgesConstraint,std::vector<double>(4,0));
    faces_Constraints_Matrix.resize(numberFacesConstraint,std::vector<unsigned int>(4,0));
    faces_Constraints_Coef.resize(numberFacesConstraint,std::vector<double>(4,0));

    //open data file to boundary displacements
    bcDataFile.open(userInputs.Periodic_BCfilename);
    id=0;
    //read data
    if (bcDataFile.is_open()){
      //It has 4 lines of headers as explanation
      for (unsigned int i=0; i<4; i++) std::getline (bcDataFile,line);

      while (id<(numberVerticesConstraint)){
        getline (bcDataFile,line);
        std::stringstream ss(line);
        ss >>vertices_Constraints_Matrix[id][0] >> vertices_Constraints_Matrix[id][1] >> vertices_Constraints_Matrix[id][2] >> vertices_Constraints_Matrix[id][3] >> vertices_Constraints_Matrix[id][4];
        std::getline (bcDataFile,line);
        std::stringstream ss2(line);
        ss2 >>vertices_Constraints_Coef[id][0] >> vertices_Constraints_Coef[id][1] >> vertices_Constraints_Coef[id][2] >> vertices_Constraints_Coef[id][3] >> vertices_Constraints_Coef[id][4];

        id=id+1;
      }

      id=0;
      //It has 4 lines of headers as explanation
      for (unsigned int i=0; i<5; i++) std::getline (bcDataFile,line);

      while (id<(numberEdgesConstraint)){
        getline (bcDataFile,line);
        std::stringstream ss(line);
        ss >>edges_Constraints_Matrix[id][0] >> edges_Constraints_Matrix[id][1] >> edges_Constraints_Matrix[id][2] >> edges_Constraints_Matrix[id][3] ;
        std::getline (bcDataFile,line);
        std::stringstream ss2(line);
        ss2 >>edges_Constraints_Coef[id][0] >> edges_Constraints_Coef[id][1] >> edges_Constraints_Coef[id][2] >> edges_Constraints_Coef[id][3] ;

        id=id+1;
      }

      id=0;
      //It has 4 lines of headers as explanation
      for (unsigned int i=0; i<5; i++) std::getline (bcDataFile,line);

      while (id<(numberFacesConstraint)){
        getline (bcDataFile,line);
        std::stringstream ss(line);
        ss >>faces_Constraints_Matrix[id][0] >> faces_Constraints_Matrix[id][1] >> faces_Constraints_Matrix[id][2] >> faces_Constraints_Matrix[id][3] ;
        std::getline (bcDataFile,line);
        std::stringstream ss2(line);
        ss2 >>faces_Constraints_Coef[id][0] >> faces_Constraints_Coef[id][1] >> faces_Constraints_Coef[id][2] >> faces_Constraints_Coef[id][3] ;

        id=id+1;
      }

    }
    else{
      pcout << "Unable to open Periodic BCfilename \n";
      exit(1);
    }
    bcDataFile.close();

    setPeriodicityConstraintsInit();
    //In the case of Periodic BCs, locally_relevant_dofs_Mod will be updated to
    //include the required DOFs, including global_vector_dof_Boundary_Layer2=(all dofs for elements with a boundary face),
    // edges, and vertices to implement the periodicity over different processors.
    unsigned int mpi_related_dof,globalDOF1;
    for(unsigned int i=0;i<global_size_dof_Boundary_Layer2;i++){
      mpi_related_dof=global_vector_dof_Boundary_Layer2[i];
      for(unsigned int j=0;j<totalNumVerticesDOFs;j++){
        globalDOF1=vertices_DOFs_vector[j];
        if ((locally_relevant_dofs_Mod.is_element(globalDOF1))&&(!locally_relevant_dofs_Mod.is_element(mpi_related_dof))){
          locally_relevant_dofs_Mod.add_index(mpi_related_dof);
          break;
        }
      }
    }

    for(unsigned int i=0;i<totalNumEachEdgesNodes*12;i++){
      for (unsigned int j=0;j<dim;j++){
        mpi_related_dof=global_Edges_DOFs_Vector_Array[j][i];
        if (!locally_relevant_dofs_Mod.is_element(mpi_related_dof)){
          locally_relevant_dofs_Mod.add_index(mpi_related_dof);
        }
      }
    }

    IndexSet::ElementIterator dofs_vertices = vertices_DOFs.begin();
    if (!locally_relevant_dofs_Mod.is_element(*dofs_vertices)){
      locally_relevant_dofs_Mod.add_index(*dofs_vertices);
    }
    for(unsigned int j = 1; j < vertices_DOFs.n_elements(); j++){
      dofs_vertices++;
      if (!locally_relevant_dofs_Mod.is_element(*dofs_vertices)){
        locally_relevant_dofs_Mod.add_index(*dofs_vertices);
      }
    }

    setPeriodicityConstraintsInc0();
    if(userInputs.enableTabularPeriodicBCs){
      setPeriodicityConstraintsInc0Neg();
    }
    setPeriodicityConstraintsIncNot0();
  }

  //initialize global data structures
  solution.reinit (locally_owned_dofs, mpi_communicator); solution=0;
  oldSolution.reinit (locally_owned_dofs, mpi_communicator); oldSolution=0;
  solutionWithGhosts.reinit (locally_owned_dofs, locally_relevant_dofs_Mod, mpi_communicator);solutionWithGhosts=0;
  solutionIncWithGhosts.reinit (locally_owned_dofs, locally_relevant_dofs_Mod, mpi_communicator);solutionIncWithGhosts=0;
  residual.reinit (locally_owned_dofs, mpi_communicator); residual=0;

  if(!userInputs.enablePeriodicBCs){
    DynamicSparsityPattern dsp (locally_relevant_dofs);
    DoFTools::make_sparsity_pattern (dofHandler, dsp, constraints, false);
    #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
    SparsityTools::distribute_sparsity_pattern (dsp,
      dofHandler.n_locally_owned_dofs_per_processor(),
      mpi_communicator,
      locally_relevant_dofs);
      jacobian.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    #else
    SparsityTools::distribute_sparsity_pattern (dsp,
      dofHandler.locally_owned_dofs(),
      mpi_communicator,
      locally_relevant_dofs);
      jacobian.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
    #endif
    }
//INDENTATION
    if(userInputs.enableIndentationBCs){

      locally_relevant_ghost_dofs = locally_relevant_dofs;
      locally_relevant_ghost_dofs.subtract_set(locally_owned_dofs);

      std::ifstream IndentationBCfile(userInputs.Indentation_BCfilename);

      //read 5 header lines for documentation for now!!!
      // read 4 data lines for now!!

      if (IndentationBCfile.is_open()){
            KeyPosIndenter.resize(userInputs.indentationKeyFrames);
            initPosIndenter.resize(dim);
            finalPosIndenter.resize(dim);
            currentPosIndenter.resize(dim);
            prevPosIndenter.resize(dim);
            pcout << "Reading Indentation boundary conditions\n";
            for (unsigned int i=0; i<4; i++) std::getline (IndentationBCfile,line);
            //pcout << "After Header lines\n";

            unsigned int temp;
            unsigned int roughTemp;
            unsigned int lines =userInputs.indentationKeyFrames + 1;
            for (unsigned int i=0; i < lines; i++){
                std::getline (IndentationBCfile,line);
                std::stringstream ss(line);
                if (i == 0) {
                        ss>>indenterShape>>indenterSize>>temp>>indenterTolerance>>roughTemp;
                        indenterFace = temp - 1;
                        indentDof = (int)(temp - 1) / (int)2;
                        if (roughTemp == 0) {
                            roughIndenter = false;
                        }
                        else {
                            roughIndenter = true;
                        }
                    }
                else {
                    if (dim==3){
                        double tem1=0;
                        double tem2=0;
                        double tem3=0;
                        std::string debugstring;
                        ss>>tem1>>tem2>>tem3;
                        KeyPosIndenter[i-1](0)=tem1;
                        KeyPosIndenter[i-1](1)=tem2;
                        KeyPosIndenter[i-1](2)=tem3;
                    }
                    else if (dim==2){
                        double tem1=0, tem2=0;
                        ss>>tem1>>tem2;
                        KeyPosIndenter[i-1](0)=tem1;
                        KeyPosIndenter[i-1](1)=tem2;
                    }
                    pcout<<KeyPosIndenter[i-1](0)<<" "<<KeyPosIndenter[i-1](1)<<" "<<KeyPosIndenter[i-1](2)<<"\n";
                }
            }
            for (unsigned int i=0; i<3; i++){
//                currentPosIndenter[i] = initPosIndenter[i];
//                prevPosIndenter[i] = initPosIndenter[i];
                currentPosIndenter[i] = KeyPosIndenter[0](i);
                prevPosIndenter[i] = KeyPosIndenter[0](i);
            }

        }
      if (indenterFace == 5) loadFace = 4;
      if (indenterFace == 4) loadFace = 5;
      if (indenterFace == 3) loadFace = 2;
      if (indenterFace == 2) loadFace = 3;
      if (indenterFace == 1) loadFace = 0;
      if (indenterFace == 0) loadFace = 1;

      if(userInputs.readExternalMesh && false){ // DISABLE FOR IMPORTING BOUNDARIES AS PHYSICAL SURFACES
          // Set which (if any) faces of the triangulation are indentation
          QGaussLobatto<dim - 1> face_quadrature_formula(FE.degree + 1);
          FEFaceValues<dim> fe_face_values(FE,face_quadrature_formula,update_values | update_JxW_values);
          const unsigned int dofs_per_cell   = FE.n_dofs_per_cell();
          std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
          Vector<double> externalMeshParameterBCs(dim);
          Point<dim> node_BoundaryID;
          unsigned int globalDOF;
          for (unsigned int i=0; i<dim; ++i) {
            externalMeshParameterBCs(i)=this->userInputs.externalMeshParameter*this->userInputs.span[i];
          }
          for (const auto &cell : dofHandler.active_cell_iterators()){
            if (cell->is_locally_owned()) {
              //std::cout << "cell locally owned! "<< std::endl;
              for (unsigned int faceID = 0; faceID < GeometryInfo<dim>::faces_per_cell; faceID++) { //(const auto &face: cell->face_iterators()){
                if (cell->face(faceID)->at_boundary()) { //(face->at_boundary() && face->boundary_id()==indenterFace) {
                    fe_face_values.reinit(cell, faceID);
                  cell->get_dof_indices (local_dof_indices);
                  for (unsigned int i=0; i<dofs_per_cell; ++i) {
                    if (fe_face_values.shape_value(i, 0)!=0){ //skip cell support points not on face
                      globalDOF=local_dof_indices[i];
                      node_BoundaryID=this->supportPoints[globalDOF];
                      unsigned int boundary= dim * 2 + 1;
                      for (unsigned int i2=0; i2<dim; ++i2)
                      {
                          if (node_BoundaryID[i2] <= externalMeshParameterBCs(0)) {
                              if (boundary == dim * 2 + 1) boundary = 2 * i2;
                              else boundary = 2 * dim;
                              //cell->face(faceID)->set_boundary_id(2 * i2);
                              //break;
                          }

                          if (node_BoundaryID[i2] >= (this->userInputs.span[i2]-externalMeshParameterBCs(0))) {
                              if (boundary == dim * 2 + 1) boundary = 2 * i2 + 1;
                              else boundary = 2 * dim;
                              //break;
                          }
                      }
                      if (boundary < dim * 2){
                          cell->face(faceID)->set_boundary_id(boundary);
                          break;
                      }

                    }

                  }
                }
              }
            }
          }
        }
    }
    // Read boundary conditions
    if((userInputs.enableSimpleBCs)||(userInputs.enableCyclicLoading)){
      std::ifstream BCfile(userInputs.BCfilename);

      //read data
      if (BCfile.is_open()){
        pcout << "Reading boundary conditions\n";
        //skip header lines
        for (unsigned int i=0; i<userInputs.BCheaderLines; i++) std::getline (BCfile,line);
        for (unsigned int i=0; i<userInputs.NumberofBCs; i++){
          std::getline (BCfile,line);
          std::stringstream ss(line);
          ss>>faceID>>dof;
          faceDOFConstrained[faceID-1][dof-1]=true;
          ss>>totalU;
          deluConstraint[faceID-1][dof-1]=totalU/totalIncrements;
        }

        if(userInputs.enableCyclicLoading){
          deluConstraint[userInputs.cyclicLoadingFace-1][userInputs.cyclicLoadingDOF-1]=deluConstraint[userInputs.cyclicLoadingFace-1][userInputs.cyclicLoadingDOF-1]*totalIncrements*userInputs.delT/userInputs.quarterCycleTime;
        }
      }
    }


    if(userInputs.enableNodalDisplacementBCs){


      std::ifstream BCfileNodal(userInputs.nodalDisplacement_BCfilename);
      nodalDisplacement.resize(userInputs.numberOfNodalBCs,std::vector<double>(3,0));
      dofNodalDisplacement.resize(userInputs.numberOfNodalBCs);
      deluNodalDisplacement.resize(userInputs.numberOfNodalBCs);

      //read data
      if (BCfileNodal.is_open()){
        pcout << "Reading Nodal boundary conditions\n";
        for (unsigned int i=0; i<userInputs.numberOfNodalBCs; i++){
          std::getline (BCfileNodal,line);
          std::stringstream ss(line);
          ss>>nodalDisplacement[i][0]>>nodalDisplacement[i][1]>>nodalDisplacement[i][2]>>dofNodalDisplacement[i];
          ss>>totalU;
          deluNodalDisplacement[i]=totalU/totalIncrements;
        }
      }
    }

// Isotropic Continuum
//    if(userInputs.continuum_Isotropic){
//        numPostProcessedFields = 2;
//        numPostProcessedFieldsAtCellCenters = 1;
//    }

    if(userInputs.useVelocityGrad){
      targetVelGrad.reinit(3,3); targetVelGrad=0.0;

      for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
          targetVelGrad[i][j] = userInputs.targetVelGrad[i][j];
        }
      }
    }


    if(userInputs.enableTabularBCs){
      tabularDisplacements.reinit(2*dim*dim,userInputs.tabularBCs_InputStepNumber);
      tabularDisplacements=0;
      //open data file to boundary displacements
      bcDataFile.open(userInputs.Tabular_BCfilename);
      //read data
      if (bcDataFile.is_open()){
        //read data
        id=0;
        while (getline (bcDataFile,line) && id<(userInputs.tabularNumberofBCs)){
          std::stringstream ss(line);
          ss>>faceID>>dof;
          faceDOFConstrained[faceID-1][dof-1]=true;

          for (unsigned int i=0; i<userInputs.tabularBCs_InputStepNumber; i++){
            ss >> tabularDisplacements[3*(faceID-1)+(dof-1)][i];
          }
          //cout<<id<<'\t'<<bc_new[id][0]<<'\t'<<bc_new[id][1]<<'\t'<<bc_new[id][2]<<'\n';
          id=id+1;
        }
      }
      else{
        pcout << "Unable to open Tabular_BCfilename \n";
        exit(1);
      }

      bcDataFile.close();
    }

    if(userInputs.enableNeumannBCs){
      tabularInputNeumannBCs.reinit(userInputs.neumannBCsNumber,userInputs.tabularNeumannBCs_InputStepNumber);
      tabularInputNeumannBCs=0;
      //open data file to boundary displacements
      bcDataFile.open(userInputs.Tabular_NeumannBCfilename);
      //read data
      if (bcDataFile.is_open()){
        //read data
        id=0;
        while (getline (bcDataFile,line) && id<(userInputs.neumannBCsNumber)){
          std::stringstream ss(line);
          for (unsigned int i=0; i<userInputs.tabularNeumannBCs_InputStepNumber; i++){
            ss >> tabularInputNeumannBCs[id][i];
          }
          //cout<<id<<'\t'<<bc_new[id][0]<<'\t'<<bc_new[id][1]<<'\t'<<bc_new[id][2]<<'\n';
          id=id+1;
        }
      }
      else{
        pcout << "Unable to open Tabular_NeumannBCfilename \n";
        exit(1);
      }

      bcDataFile.close();
    }

    if(userInputs.enableDICpipeline){
      bc_new1.reinit(userInputs.Y_dic,1+2*userInputs.DIC_InputStepNumber);
      bc_new2.reinit(userInputs.Y_dic,1+2*userInputs.DIC_InputStepNumber);
      bc_new3.reinit(userInputs.X_dic,1+2*userInputs.DIC_InputStepNumber);
      bc_new4.reinit(userInputs.X_dic,1+2*userInputs.DIC_InputStepNumber);
      //open data file to boundary displacements
      bcDataFile.open(userInputs.DIC_BCfilename1);
      //read data
      id=0;
      if (bcDataFile.is_open()){
        //read data
        while (getline (bcDataFile,line) && id<(userInputs.Y_dic)){
          std::stringstream ss(line);
          ss >> bc_new1[id][0];
          for (unsigned int i=0; i<userInputs.DIC_InputStepNumber; i++){
            ss >> bc_new1[id][1+i*2];
            ss >> bc_new1[id][2+i*2];
          }
          //cout<<id<<'\t'<<bc_new[id][0]<<'\t'<<bc_new[id][1]<<'\t'<<bc_new[id][2]<<'\n';
          id=id+1;
        }
      }
      else{
        pcout << "Unable to open DIC_BCfilename1 \n";
        exit(1);
      }

      bcDataFile.close();



      //open data file to boundary displacements
      bcDataFile.open(userInputs.DIC_BCfilename2);
      //read data
      id=0;
      if (bcDataFile.is_open()){
        //read data
        while (getline (bcDataFile,line) && id<(userInputs.Y_dic)){
          std::stringstream ss(line);
          ss >> bc_new2[id][0];
          for (unsigned int i=0; i<userInputs.DIC_InputStepNumber; i++){
            ss >> bc_new2[id][1+i*2];
            ss >> bc_new2[id][2+i*2];
          }
          //cout<<id<<'\t'<<bc_new[id][0]<<'\t'<<bc_new[id][1]<<'\t'<<bc_new[id][2]<<'\n';
          id=id+1;
        }
      }
      else{
        pcout << "Unable to open DIC_BCfilename2 \n";
        exit(1);
      }

      bcDataFile.close();


      //open data file to boundary displacements
      bcDataFile.open(userInputs.DIC_BCfilename3);
      //read data
      id=0;
      if (bcDataFile.is_open()){
        //read data
        while (getline (bcDataFile,line) && id<(userInputs.X_dic)){
          std::stringstream ss(line);
          ss >> bc_new3[id][0];
          for (unsigned int i=0; i<userInputs.DIC_InputStepNumber; i++){
            ss >> bc_new3[id][1+i*2];
            ss >> bc_new3[id][2+i*2];
          }
          //cout<<id<<'\t'<<bc_new[id][0]<<'\t'<<bc_new[id][1]<<'\t'<<bc_new[id][2]<<'\n';
          id=id+1;
        }
      }
      else{
        pcout << "Unable to open DIC_BCfilename3 \n";
        exit(1);
      }

      bcDataFile.close();


      //open data file to boundary displacements
      bcDataFile.open(userInputs.DIC_BCfilename4);
      //read data
      id=0;
      if (bcDataFile.is_open()){
        //read data
        while (getline (bcDataFile,line) && id<(userInputs.X_dic)){
          std::stringstream ss(line);
          ss >> bc_new4[id][0];
          for (unsigned int i=0; i<userInputs.DIC_InputStepNumber; i++){
            ss >> bc_new4[id][1+i*2];
            ss >> bc_new4[id][2+i*2];
          }
          //cout<<id<<'\t'<<bc_new[id][0]<<'\t'<<bc_new[id][1]<<'\t'<<bc_new[id][2]<<'\n';
          id=id+1;
        }
      }
      else{
        pcout << "Unable to open DIC_BCfilename4 \n";
        exit(1);
      }

      bcDataFile.close();

    }

    Fprev=IdentityMatrix(dim);

    //apply initial conditions
    applyInitialConditions();
    solutionWithGhosts=solution;
    oldSolution=solution;

  }
  #include "../../include/ellipticBVP_template_instantiations.h"
