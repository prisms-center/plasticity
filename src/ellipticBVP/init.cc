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
    SparsityTools::distribute_sparsity_pattern (dsp,
      dofHandler.n_locally_owned_dofs_per_processor(),
      mpi_communicator,
      locally_relevant_dofs);
      jacobian.reinit (locally_owned_dofs, locally_owned_dofs, dsp, mpi_communicator);
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
