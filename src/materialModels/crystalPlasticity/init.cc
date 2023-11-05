#include "../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

template <int dim>
void crystalPlasticity<dim>::init(unsigned int num_quad_points)
{

  //call loadOrientations to load material orientations
  loadOrientations();

  local_strain.reinit(dim,dim);
  local_stress.reinit(dim,dim);
  global_strain.reinit(dim,dim);
  global_stress.reinit(dim,dim);
  local_strain=0.0;
  local_stress=0.0;
  global_strain=0.0;
  global_stress=0.0;
  FullMatrix<double> CauchyStress_init(dim,dim), TinterStress_init(dim,dim), TinterStress_diff_init(dim,dim);
  CauchyStress_init=0;
  TinterStress_init = 0 ;
  TinterStress_diff_init = 0 ;
  unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();
  F.reinit(dim, dim);

  double m_norm , n_norm ;
  unsigned int n_Tslip_systems_Real_SinglePhase,n_Tslip_systems_Real;
  n_slip_systems_SinglePhase=this->userInputs.numSlipSystems1;
  n_Tslip_systems_SinglePhase=n_slip_systems_SinglePhase;
  if(this->userInputs.enableTwinning1){
    n_Tslip_systems_SinglePhase+=this->userInputs.numTwinSystems1;
    n_twin_systems_SinglePhase=this->userInputs.numTwinSystems1;
    n_Tslip_systems_Real_SinglePhase=n_Tslip_systems_SinglePhase;

  }
  else{
    n_Tslip_systems_SinglePhase+=1;
    n_twin_systems_SinglePhase=1;
    n_Tslip_systems_Real_SinglePhase=n_Tslip_systems_SinglePhase-1;
  }


  n_alpha_SinglePhase.reinit(n_Tslip_systems_SinglePhase,3);
  m_alpha_SinglePhase.reinit(n_Tslip_systems_SinglePhase,3);



  std::string line;

  q_phase1.reinit(n_Tslip_systems_SinglePhase,n_Tslip_systems_SinglePhase);
  q_phase1=0;
  //open data file to read latent hardening ratios
  std::ifstream latentHardeningratioFile(this->userInputs.latentHardeningRatioFileName1);
  //read data
  unsigned int id=0;
  if (latentHardeningratioFile.is_open()){
    while (getline (latentHardeningratioFile,line) && id<n_Tslip_systems_Real_SinglePhase){
      std::stringstream ss(line);
      for (unsigned int i=0; i<n_Tslip_systems_Real_SinglePhase; i++){
        ss >> q_phase1[id][i];
      }
      id=id+1;
    }
  }
  else{
    std::cout << "Unable to open latent hardening ratio file \n";
    exit(1);
  }


  //open data file to read slip normals
  std::ifstream slipNormalsDataFile(this->userInputs.slipNormalsFile1);
  //read data
  id=0;
  if (slipNormalsDataFile.is_open()){
    while (getline (slipNormalsDataFile,line) && id<n_slip_systems_SinglePhase){
      std::stringstream ss(line);
      ss >> n_alpha_SinglePhase[id][0];
      ss >> n_alpha_SinglePhase[id][1];
      ss >> n_alpha_SinglePhase[id][2];
      n_norm = 0 ;
      n_norm = n_norm + n_alpha_SinglePhase[id][0]*n_alpha_SinglePhase[id][0] ;
      n_norm = n_norm + n_alpha_SinglePhase[id][1]*n_alpha_SinglePhase[id][1] ;
      n_norm = n_norm + n_alpha_SinglePhase[id][2]*n_alpha_SinglePhase[id][2] ;
      n_norm = sqrt(n_norm) ;
      n_alpha_SinglePhase[id][0] = n_alpha_SinglePhase[id][0]/n_norm ;
      n_alpha_SinglePhase[id][1] = n_alpha_SinglePhase[id][1]/n_norm ;
      n_alpha_SinglePhase[id][2] = n_alpha_SinglePhase[id][2]/n_norm ;

      id=id+1;
    }
  }
  else{
    std::cout << "Unable to open slip normals file \n";
    exit(1);
  }

  //open data file to read slip directions
  std::ifstream slipDirectionsDataFile(this->userInputs.slipDirectionsFile1);
  //read data
  id=0;
  if (slipDirectionsDataFile.is_open()){
    while (getline (slipDirectionsDataFile,line)&& id<n_slip_systems_SinglePhase){
      std::stringstream ss(line);
      ss >> m_alpha_SinglePhase[id][0];
      ss >> m_alpha_SinglePhase[id][1];
      ss >> m_alpha_SinglePhase[id][2];
      m_norm = 0 ;
      m_norm = m_norm + m_alpha_SinglePhase[id][0]*m_alpha_SinglePhase[id][0] ;
      m_norm = m_norm + m_alpha_SinglePhase[id][1]*m_alpha_SinglePhase[id][1] ;
      m_norm = m_norm + m_alpha_SinglePhase[id][2]*m_alpha_SinglePhase[id][2] ;
      m_norm = sqrt(m_norm) ;
      m_alpha_SinglePhase[id][0] = m_alpha_SinglePhase[id][0]/m_norm ;
      m_alpha_SinglePhase[id][1] = m_alpha_SinglePhase[id][1]/m_norm ;
      m_alpha_SinglePhase[id][2] = m_alpha_SinglePhase[id][2]/m_norm ;
      id=id+1;
    }
  }
  else{
    std::cout << "Unable to open slip directions file \n";
    exit(1);
  }

  if(this->userInputs.enableTwinning1){
    //open data file to read twin normals
    std::ifstream twinNormalsDataFile(this->userInputs.twinNormalsFile1);
    //read data
    id=n_slip_systems_SinglePhase;
    if (twinNormalsDataFile.is_open())
    while (getline (twinNormalsDataFile,line) && id<n_Tslip_systems_SinglePhase){
      std::stringstream ss(line);
      ss >> n_alpha_SinglePhase[id][0];
      ss >> n_alpha_SinglePhase[id][1];
      ss >> n_alpha_SinglePhase[id][2];
      n_norm = 0 ;
      n_norm = n_norm + n_alpha_SinglePhase[id][0]*n_alpha_SinglePhase[id][0] ;
      n_norm = n_norm + n_alpha_SinglePhase[id][1]*n_alpha_SinglePhase[id][1] ;
      n_norm = n_norm + n_alpha_SinglePhase[id][2]*n_alpha_SinglePhase[id][2] ;
      n_norm = sqrt(n_norm) ;
      n_alpha_SinglePhase[id][0] = n_alpha_SinglePhase[id][0]/n_norm ;
      n_alpha_SinglePhase[id][1] = n_alpha_SinglePhase[id][1]/n_norm ;
      n_alpha_SinglePhase[id][2] = n_alpha_SinglePhase[id][2]/n_norm ;

      id=id+1;
    }
    else{
      std::cout << "Unable to open twin normals file\n";
      exit(1);
    }

    //open data file to read twin directions
    std::ifstream twinDirectionsDataFile(this->userInputs.twinDirectionsFile1);
    //read data
    id=n_slip_systems_SinglePhase;
    if (twinDirectionsDataFile.is_open())
    //read data
    while (getline (twinDirectionsDataFile,line)&& id<n_Tslip_systems_SinglePhase){
      std::stringstream ss(line);
      ss >> m_alpha_SinglePhase[id][0];
      ss >> m_alpha_SinglePhase[id][1];
      ss >> m_alpha_SinglePhase[id][2];
      m_norm = 0 ;
      m_norm = m_norm + m_alpha_SinglePhase[id][0]*m_alpha_SinglePhase[id][0] ;
      m_norm = m_norm + m_alpha_SinglePhase[id][1]*m_alpha_SinglePhase[id][1] ;
      m_norm = m_norm + m_alpha_SinglePhase[id][2]*m_alpha_SinglePhase[id][2] ;
      m_norm = sqrt(m_norm) ;
      m_alpha_SinglePhase[id][0] = m_alpha_SinglePhase[id][0]/m_norm ;
      m_alpha_SinglePhase[id][1] = m_alpha_SinglePhase[id][1]/m_norm ;
      m_alpha_SinglePhase[id][2] = m_alpha_SinglePhase[id][2]/m_norm ;
      id=id+1;
    }
    else{
      std::cout << "Unable to open twin directions file \n";
      exit(1);
    }
  }
  else{
    n_alpha_SinglePhase[id][0]=1;
    n_alpha_SinglePhase[id][1]=0;
    n_alpha_SinglePhase[id][2]=0;
    m_alpha_SinglePhase[id][0]=0;
    m_alpha_SinglePhase[id][1]=1;
    m_alpha_SinglePhase[id][2]=0;
  }



  //Elastic Stiffness Matrix Dmat
  Dmat_SinglePhase.reinit(6,6); Dmat=0.0;

  for(unsigned int i=0;i<6;i++){
    for(unsigned int j=0;j<6;j++){
      Dmat_SinglePhase[i][j] = this->userInputs.elasticStiffness1[i][j];
    }
  }

  Vector<double> rot_init(dim),rotnew_init(dim),VoxelData_init;
  unsigned int numberOfAdditionalVoxelInfo=this->userInputs.additionalVoxelInfo;

  for (unsigned int i=0;i<dim;i++){
    rot_init(i)=0.0;
    rotnew_init(i)=0.0;
  }

  rot_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rot_init));
  rotnew_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rotnew_init));
  rot_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rot_init));
  rotnew_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rotnew_init));
  phase.resize(num_local_cells,std::vector<unsigned int>(num_quad_points,1));
  if (this->userInputs.enableMultiphase){
    numberofPhases=this->userInputs.numberofPhases;
  }

  if (numberOfAdditionalVoxelInfo>0){
    VoxelData_init.reinit(numberOfAdditionalVoxelInfo);
    for (unsigned int i=0;i<this->userInputs.additionalVoxelInfo;i++){
      VoxelData_init(i)=0.0;
    }
    VoxelData.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,VoxelData_init));
  }

  unsigned int counter=0;
  //load rot, rotnew and VoxelData
  for (unsigned int cell=0; cell<num_local_cells; cell++){
    unsigned int materialID=cellOrientationMap[cell];
    for (unsigned int q=0; q<num_quad_points; q++){
      for (unsigned int i=0; i<dim; i++){
        rot_iter[cell][q][i]=orientations.eulerAngles[materialID][i];
        rotnew_iter[cell][q][i]=orientations.eulerAngles[materialID][i];
      }
      if (this->userInputs.enableMultiphase){
        phase[cell][q]=orientations.eulerAngles[materialID][3];
        counter=1;
      }
      if (numberOfAdditionalVoxelInfo>0){
        for (unsigned int i=dim+counter; i<dim+counter+this->userInputs.additionalVoxelInfo; i++){
          VoxelData[cell][q][i-dim-counter]=orientations.eulerAngles[materialID][i];
        }
      }
    }
  }

  rot_conv=rot_iter;
  rotnew_conv=rotnew_iter;



  Vector<double> s0_init (n_Tslip_systems_SinglePhase);
  std::vector<double> twin_init(n_twin_systems_SinglePhase),slip_init(n_slip_systems_SinglePhase);
  Vector<double> s0_init1,s0_init2,s0_init3,s0_init4;
  Vector<double> stateVar_init,stateVar_init1, stateVar_init2, stateVar_init3, stateVar_init4;
  std::vector<double> twin_init1,slip_init1;
  Vector<double> W_kh_init1;

  for (unsigned int i=0;i<n_slip_systems_SinglePhase;i++){
    s0_init(i)=this->userInputs.initialSlipResistance1[i];
  }

  for (unsigned int i=0;i<n_twin_systems_SinglePhase;i++){
    s0_init(i+n_slip_systems_SinglePhase)=this->userInputs.initialSlipResistanceTwin1[i];
  }


  for (unsigned int i=0;i<n_slip_systems_SinglePhase;i++){
    slip_init[i]=0.0;
  }

  Vector<double> W_kh_init(n_Tslip_systems_SinglePhase);
  for (unsigned int i = 0;i<n_Tslip_systems_SinglePhase;i++) {
    W_kh_init[i] = 0.0;
  }

  for (unsigned int i=0;i<n_twin_systems_SinglePhase;i++){
    twin_init[i]=0.0;
  }

  if (this->userInputs.enableUserMaterialModel){
    if (this->userInputs.enableUserMaterialModel1){
      if (this->userInputs.numberofUserMatStateVar1==0){
        n_UserMatStateVar_SinglePhase=1;
        stateVar_init.reinit(n_UserMatStateVar_SinglePhase);
        stateVar_init=0;
      }
      else{
        n_UserMatStateVar_SinglePhase=this->userInputs.numberofUserMatStateVar1;
        stateVar_init.reinit(n_UserMatStateVar_SinglePhase);
        for (unsigned int i=0;i<n_UserMatStateVar_SinglePhase;i++){
          stateVar_init(i)=this->userInputs.UserMatStateVar1[i];
        }
      }
    }
    else{
      n_UserMatStateVar_SinglePhase=1;
      stateVar_init.reinit(n_UserMatStateVar_SinglePhase);
      stateVar_init=0;
    }
  }



  if (!this->userInputs.enableMultiphase){
    //Resize the vectors of history variables
    Fp_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init));
    W_kh_conv.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, W_kh_init));
    W_kh_iter.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, W_kh_init));
    Fp_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    CauchyStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));
    F_lastIter_Global.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    
    FirstPiolaStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));
    SecondPiolaStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));
    workDensity1.reinit(num_local_cells); workDensity1 = 0.0;
    workDensity2.reinit(num_local_cells); workDensity2 = 0.0;
    workDensityTotal1.reinit(num_local_cells); workDensityTotal1 = 0.0;
    workDensityTotal2.reinit(num_local_cells); workDensityTotal2 = 0.0;
    workDensity1_Tr.reinit(num_local_cells); workDensity1_Tr = 0.0;
    workDensity2_Tr.reinit(num_local_cells); workDensity2_Tr = 0.0;
    workDensityTotal1_Tr.reinit(num_local_cells); workDensityTotal1_Tr = 0.0;
    workDensityTotal2_Tr.reinit(num_local_cells); workDensityTotal2_Tr = 0.0;
	 TinterStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,TinterStress_init));
	 TinterStress_diff.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,TinterStress_diff_init));
    s_alpha_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init));
    twinfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
    slipfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,slip_init));
    twinfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
    slipfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,slip_init));
    twin_ouput.resize(num_local_cells, std::vector<double>(num_quad_points,0.0));
    twin_conv.resize(num_local_cells,std::vector<unsigned int>(num_quad_points,0));
    twin_iter.resize(num_local_cells,std::vector<unsigned int>(num_quad_points,0));

    if (this->userInputs.enableUserMaterialModel){
      stateVar_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,stateVar_init));
      stateVar_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,stateVar_init));
    }
  }


  unsigned int Max_n_slip_systems_MultiPhase,Max_n_twin_systems_MultiPhase,Max_n_Tslip_systems_MultiPhase, Max_n_UserMatStateVar_MultiPhase;
  if(this->userInputs.enableMultiphase){
    n_slip_systems_MultiPhase.resize(this->userInputs.numberofPhases);
    n_twin_systems_MultiPhase.resize(this->userInputs.numberofPhases);
    n_Tslip_systems_MultiPhase.resize(this->userInputs.numberofPhases);
    n_Tslip_systems_Real_MultiPhase.resize(this->userInputs.numberofPhases);



    n_slip_systems_MultiPhase[0]=n_slip_systems_SinglePhase;
    n_twin_systems_MultiPhase[0]=n_twin_systems_SinglePhase;
    n_Tslip_systems_MultiPhase[0]=n_Tslip_systems_SinglePhase;
    n_Tslip_systems_Real_MultiPhase[0]=n_Tslip_systems_Real_SinglePhase;

    n_UserMatStateVar_MultiPhase.resize(this->userInputs.numberofPhases);
    n_UserMatStateVar_MultiPhase[0]=n_UserMatStateVar_SinglePhase;


    if (numberofPhases>=2){
      n_slip_systems=this->userInputs.numSlipSystems2;
      n_Tslip_systems=n_slip_systems;
      if(this->userInputs.enableTwinning2){
        n_Tslip_systems+=this->userInputs.numTwinSystems2;
        n_twin_systems=this->userInputs.numTwinSystems2;
        n_Tslip_systems_Real=n_Tslip_systems;
      }
      else{
        n_Tslip_systems+=1;
        n_twin_systems=1;
        n_Tslip_systems_Real=n_Tslip_systems-1;
      }
      n_slip_systems_MultiPhase[1]=n_slip_systems;
      n_twin_systems_MultiPhase[1]=n_twin_systems;
      n_Tslip_systems_MultiPhase[1]=n_Tslip_systems;
      n_Tslip_systems_Real_MultiPhase[1]=n_Tslip_systems_Real;

      if (this->userInputs.enableUserMaterialModel2){
        if (this->userInputs.numberofUserMatStateVar2==0){
          n_UserMatStateVar=1;
        }
        else{
          n_UserMatStateVar=this->userInputs.numberofUserMatStateVar2;
        }
      }
      else{
        n_UserMatStateVar=1;
      }
      n_UserMatStateVar_MultiPhase[1]=n_UserMatStateVar;

      if (numberofPhases>=3){
        n_slip_systems=this->userInputs.numSlipSystems3;
        n_Tslip_systems=n_slip_systems;
        if(this->userInputs.enableTwinning3){
          n_Tslip_systems+=this->userInputs.numTwinSystems3;
          n_twin_systems=this->userInputs.numTwinSystems3;
          n_Tslip_systems_Real=n_Tslip_systems;
        }
        else{
          n_Tslip_systems+=1;
          n_twin_systems=1;
          n_Tslip_systems_Real=n_Tslip_systems-1;
        }
        n_slip_systems_MultiPhase[2]=n_slip_systems;
        n_twin_systems_MultiPhase[2]=n_twin_systems;
        n_Tslip_systems_MultiPhase[2]=n_Tslip_systems;
        n_Tslip_systems_Real_MultiPhase[2]=n_Tslip_systems_Real;

        if (this->userInputs.enableUserMaterialModel3){
          if (this->userInputs.numberofUserMatStateVar3==0){
            n_UserMatStateVar=1;
          }
          else{
            n_UserMatStateVar=this->userInputs.numberofUserMatStateVar3;
          }
        }
        else{
          n_UserMatStateVar=1;
        }
        n_UserMatStateVar_MultiPhase[2]=n_UserMatStateVar;

        if (numberofPhases>=4){
          n_slip_systems=this->userInputs.numSlipSystems4;
          n_Tslip_systems=n_slip_systems;
          if(this->userInputs.enableTwinning4){
            n_Tslip_systems+=this->userInputs.numTwinSystems4;
            n_twin_systems=this->userInputs.numTwinSystems4;
            n_Tslip_systems_Real=n_Tslip_systems;
          }
          else{
            n_Tslip_systems+=1;
            n_twin_systems=1;
            n_Tslip_systems_Real=n_Tslip_systems-1;
          }
          n_slip_systems_MultiPhase[3]=n_slip_systems;
          n_twin_systems_MultiPhase[3]=n_twin_systems;
          n_Tslip_systems_MultiPhase[3]=n_Tslip_systems;
          n_Tslip_systems_Real_MultiPhase[3]=n_Tslip_systems_Real;

          if (this->userInputs.enableUserMaterialModel4){
            if (this->userInputs.numberofUserMatStateVar4==0){
              n_UserMatStateVar=1;
            }
            else{
              n_UserMatStateVar=this->userInputs.numberofUserMatStateVar4;
            }
          }
          else{
            n_UserMatStateVar=1;
          }
          n_UserMatStateVar_MultiPhase[3]=n_UserMatStateVar;

        }

      }

    }
    unsigned int n_Tslip_systems_Total=0;
    for(unsigned int i=0;i<this->userInputs.numberofPhases;i++){
      n_Tslip_systems_Total=n_Tslip_systems_Total+n_Tslip_systems_MultiPhase[i];
    }

    Max_n_slip_systems_MultiPhase=std::max(n_slip_systems_MultiPhase[0],n_slip_systems_MultiPhase[1]);
    Max_n_slip_systems_MultiPhase=std::max(Max_n_slip_systems_MultiPhase,n_slip_systems_MultiPhase[2]);
    Max_n_slip_systems_MultiPhase=std::max(Max_n_slip_systems_MultiPhase,n_slip_systems_MultiPhase[3]);

    Max_n_Tslip_systems_MultiPhase=std::max(n_Tslip_systems_MultiPhase[0],n_Tslip_systems_MultiPhase[1]);
    Max_n_Tslip_systems_MultiPhase=std::max(Max_n_Tslip_systems_MultiPhase,n_Tslip_systems_MultiPhase[2]);
    Max_n_Tslip_systems_MultiPhase=std::max(Max_n_Tslip_systems_MultiPhase,n_Tslip_systems_MultiPhase[3]);

    Max_n_twin_systems_MultiPhase=std::max(n_twin_systems_MultiPhase[0],n_twin_systems_MultiPhase[1]);
    Max_n_twin_systems_MultiPhase=std::max(Max_n_twin_systems_MultiPhase,n_twin_systems_MultiPhase[2]);
    Max_n_twin_systems_MultiPhase=std::max(Max_n_twin_systems_MultiPhase,n_twin_systems_MultiPhase[3]);

    if (this->userInputs.enableUserMaterialModel){
      Max_n_UserMatStateVar_MultiPhase=std::max(n_UserMatStateVar_MultiPhase[0],n_UserMatStateVar_MultiPhase[1]);
      Max_n_UserMatStateVar_MultiPhase=std::max(Max_n_UserMatStateVar_MultiPhase,n_UserMatStateVar_MultiPhase[2]);
      Max_n_UserMatStateVar_MultiPhase=std::max(Max_n_UserMatStateVar_MultiPhase,n_UserMatStateVar_MultiPhase[3]);
    }

    n_alpha_MultiPhase.reinit(n_Tslip_systems_Total,3);
    m_alpha_MultiPhase.reinit(n_Tslip_systems_Total,3);

    for(unsigned int i=0;i<n_Tslip_systems_MultiPhase[0];i++){
      m_norm = 0 ; n_norm = 0 ;
      for(unsigned int j=0;j<dim;j++){
        m_alpha_MultiPhase[i][j]=m_alpha_SinglePhase[i][j];
        n_alpha_MultiPhase[i][j]=n_alpha_SinglePhase[i][j];
        m_norm += m_alpha_MultiPhase[i][j]*m_alpha_MultiPhase[i][j] ;
        n_norm += n_alpha_MultiPhase[i][j]*n_alpha_MultiPhase[i][j] ;
      }
     m_norm = sqrt(m_norm) ; n_norm = sqrt(n_norm) ;

      for(unsigned int j=0;j<dim;j++){
        m_alpha_MultiPhase[i][j]=m_alpha_MultiPhase[i][j]/m_norm;
        n_alpha_MultiPhase[i][j]=n_alpha_MultiPhase[i][j]/n_norm;
      }

    }

    if (numberofPhases>=2){

      q_phase2.reinit(n_Tslip_systems_MultiPhase[1],n_Tslip_systems_MultiPhase[1]);
      q_phase2=0;
      //open data file to read latent hardening ratios
      std::ifstream latentHardeningratioFile2(this->userInputs.latentHardeningRatioFileName2);
      //read data
      id=0;
      if (latentHardeningratioFile2.is_open()){
        while (getline (latentHardeningratioFile2,line) && id<n_Tslip_systems_Real_MultiPhase[1]){
          std::stringstream ss(line);
          for (unsigned int i=0; i<n_Tslip_systems_Real_MultiPhase[1]; i++){
            ss >> q_phase2[id][i];
          }
          id=id+1;
        }
      }
      else{
        std::cout << "Unable to open latent hardening ratio file 2 \n";
        exit(1);
      }


      //open data file to read slip normals
      std::ifstream slipNormalsDataFile2(this->userInputs.slipNormalsFile2);
      //read data
      id=n_Tslip_systems_MultiPhase[0];
      if (slipNormalsDataFile2.is_open()){
        while (getline (slipNormalsDataFile2,line) && id<n_Tslip_systems_MultiPhase[0]+n_slip_systems_MultiPhase[1]){
          std::stringstream ss(line);
          ss >> n_alpha_MultiPhase[id][0];
          ss >> n_alpha_MultiPhase[id][1];
          ss >> n_alpha_MultiPhase[id][2];
          n_norm = 0 ;
          n_norm = n_norm + n_alpha_MultiPhase[id][0]*n_alpha_MultiPhase[id][0] ;
          n_norm = n_norm + n_alpha_MultiPhase[id][1]*n_alpha_MultiPhase[id][1] ;
          n_norm = n_norm + n_alpha_MultiPhase[id][2]*n_alpha_MultiPhase[id][2] ;
          n_norm = sqrt(n_norm) ;
          n_alpha_MultiPhase[id][0] = n_alpha_MultiPhase[id][0]/n_norm ;
          n_alpha_MultiPhase[id][1] = n_alpha_MultiPhase[id][1]/n_norm ;
          n_alpha_MultiPhase[id][2] = n_alpha_MultiPhase[id][2]/n_norm ;
          id=id+1;
        }
      }
      else{
        std::cout << "Unable to open slip normals file \n";
        exit(1);
      }

      //open data file to read slip directions
      std::ifstream slipDirectionsDataFile2(this->userInputs.slipDirectionsFile2);
      //read data
      id=n_Tslip_systems_MultiPhase[0];
      if (slipDirectionsDataFile2.is_open()){
        while (getline (slipDirectionsDataFile2,line)&& id<n_Tslip_systems_MultiPhase[0]+n_slip_systems_MultiPhase[1]){
          std::stringstream ss(line);
          ss >> m_alpha_MultiPhase[id][0];
          ss >> m_alpha_MultiPhase[id][1];
          ss >> m_alpha_MultiPhase[id][2];
          m_norm = 0 ;
          m_norm = m_norm + m_alpha_MultiPhase[id][0]*m_alpha_MultiPhase[id][0] ;
          m_norm = m_norm + m_alpha_MultiPhase[id][1]*m_alpha_MultiPhase[id][1] ;
          m_norm = m_norm + m_alpha_MultiPhase[id][2]*m_alpha_MultiPhase[id][2] ;
          m_norm = sqrt(m_norm) ;
          m_alpha_MultiPhase[id][0] = m_alpha_MultiPhase[id][0]/m_norm ;
          m_alpha_MultiPhase[id][1] = m_alpha_MultiPhase[id][1]/m_norm ;
          m_alpha_MultiPhase[id][2] = m_alpha_MultiPhase[id][2]/m_norm ;
          id=id+1;
        }
      }
      else{
        std::cout << "Unable to open slip directions file \n";
        exit(1);
      }

      if(this->userInputs.enableTwinning2){
        //open data file to read twin normals
        std::ifstream twinNormalsDataFile2(this->userInputs.twinNormalsFile2);
        //read data
        id=n_slip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
        if (twinNormalsDataFile2.is_open())
        while (getline (twinNormalsDataFile2,line) && id<n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]){
          std::stringstream ss(line);
          ss >> n_alpha_MultiPhase[id][0];
          ss >> n_alpha_MultiPhase[id][1];
          ss >> n_alpha_MultiPhase[id][2];
          n_norm = 0 ;
          n_norm = n_norm + n_alpha_MultiPhase[id][0]*n_alpha_MultiPhase[id][0] ;
          n_norm = n_norm + n_alpha_MultiPhase[id][1]*n_alpha_MultiPhase[id][1] ;
          n_norm = n_norm + n_alpha_MultiPhase[id][2]*n_alpha_MultiPhase[id][2] ;
          n_norm = sqrt(n_norm) ;
          n_alpha_MultiPhase[id][0] = n_alpha_MultiPhase[id][0]/n_norm ;
          n_alpha_MultiPhase[id][1] = n_alpha_MultiPhase[id][1]/n_norm ;
          n_alpha_MultiPhase[id][2] = n_alpha_MultiPhase[id][2]/n_norm ;
          id=id+1;
        }
        else{
          std::cout << "Unable to open twin normals file\n";
          exit(1);
        }

        //open data file to read twin directions
        std::ifstream twinDirectionsDataFile2(this->userInputs.twinDirectionsFile2);
        //read data
        id=n_slip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
        if (twinDirectionsDataFile2.is_open())
        //read data
        while (getline (twinDirectionsDataFile2,line)&& id<n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]){
          std::stringstream ss(line);
          ss >> m_alpha_MultiPhase[id][0];
          ss >> m_alpha_MultiPhase[id][1];
          ss >> m_alpha_MultiPhase[id][2];
          m_norm = 0 ;
          m_norm = m_norm + m_alpha_MultiPhase[id][0]*m_alpha_MultiPhase[id][0] ;
          m_norm = m_norm + m_alpha_MultiPhase[id][1]*m_alpha_MultiPhase[id][1] ;
          m_norm = m_norm + m_alpha_MultiPhase[id][2]*m_alpha_MultiPhase[id][2] ;
          m_norm = sqrt(m_norm) ;
          m_alpha_MultiPhase[id][0] = m_alpha_MultiPhase[id][0]/m_norm ;
          m_alpha_MultiPhase[id][1] = m_alpha_MultiPhase[id][1]/m_norm ;
          m_alpha_MultiPhase[id][2] = m_alpha_MultiPhase[id][2]/m_norm ;

          id=id+1;
        }
        else{
          std::cout << "Unable to open twin directions file \n";
          exit(1);
        }
      }
      else{
        n_alpha_MultiPhase[id][0]=1;
        n_alpha_MultiPhase[id][1]=0;
        n_alpha_MultiPhase[id][2]=0;
        m_alpha_MultiPhase[id][0]=0;
        m_alpha_MultiPhase[id][1]=1;
        m_alpha_MultiPhase[id][2]=0;
      }



      if (numberofPhases>=3){

        q_phase3.reinit(n_Tslip_systems_MultiPhase[2],n_Tslip_systems_MultiPhase[2]);
        q_phase3=0;
        //open data file to read latent hardening ratios
        std::ifstream latentHardeningratioFile3(this->userInputs.latentHardeningRatioFileName3);
        //read data
        id=0;
        if (latentHardeningratioFile3.is_open()){
          while (getline (latentHardeningratioFile3,line) && id<n_Tslip_systems_Real_MultiPhase[2]){
            std::stringstream ss(line);
            for (unsigned int i=0; i<n_Tslip_systems_Real_MultiPhase[2]; i++){
              ss >> q_phase3[id][i];
            }
            id=id+1;
          }
        }
        else{
          std::cout << "Unable to open latent hardening ratio file 3 \n";
          exit(1);
        }


        //open data file to read slip normals
        std::ifstream slipNormalsDataFile3(this->userInputs.slipNormalsFile3);
        //read data
        id=n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1];
        if (slipNormalsDataFile3.is_open()){
          while (getline (slipNormalsDataFile3,line) && id<n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]+n_slip_systems_MultiPhase[2]){
            std::stringstream ss(line);
            ss >> n_alpha_MultiPhase[id][0];
            ss >> n_alpha_MultiPhase[id][1];
            ss >> n_alpha_MultiPhase[id][2];
            n_norm = 0 ;
            n_norm = n_norm + n_alpha_MultiPhase[id][0]*n_alpha_MultiPhase[id][0] ;
            n_norm = n_norm + n_alpha_MultiPhase[id][1]*n_alpha_MultiPhase[id][1] ;
            n_norm = n_norm + n_alpha_MultiPhase[id][2]*n_alpha_MultiPhase[id][2] ;
            n_norm = sqrt(n_norm) ;
            n_alpha_MultiPhase[id][0] = n_alpha_MultiPhase[id][0]/n_norm ;
            n_alpha_MultiPhase[id][1] = n_alpha_MultiPhase[id][1]/n_norm ;
            n_alpha_MultiPhase[id][2] = n_alpha_MultiPhase[id][2]/n_norm ;
            id=id+1;
          }
        }
        else{
          std::cout << "Unable to open slip normals file \n";
          exit(1);
        }

        //open data file to read slip directions
        std::ifstream slipDirectionsDataFile3(this->userInputs.slipDirectionsFile3);
        //read data
        id=n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1];
        if (slipDirectionsDataFile3.is_open()){
          while (getline (slipDirectionsDataFile3,line)&& id<n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]+n_slip_systems_MultiPhase[2]){
            std::stringstream ss(line);
            ss >> m_alpha_MultiPhase[id][0];
            ss >> m_alpha_MultiPhase[id][1];
            ss >> m_alpha_MultiPhase[id][2];
            m_norm = 0 ;
            m_norm = m_norm + m_alpha_MultiPhase[id][0]*m_alpha_MultiPhase[id][0] ;
            m_norm = m_norm + m_alpha_MultiPhase[id][1]*m_alpha_MultiPhase[id][1] ;
            m_norm = m_norm + m_alpha_MultiPhase[id][2]*m_alpha_MultiPhase[id][2] ;
            m_norm = sqrt(m_norm) ;
            m_alpha_MultiPhase[id][0] = m_alpha_MultiPhase[id][0]/m_norm ;
            m_alpha_MultiPhase[id][1] = m_alpha_MultiPhase[id][1]/m_norm ;
            m_alpha_MultiPhase[id][2] = m_alpha_MultiPhase[id][2]/m_norm ;

            id=id+1;
          }
        }
        else{
          std::cout << "Unable to open slip directions file \n";
          exit(1);
        }

        if(this->userInputs.enableTwinning3){
          //open data file to read twin normals
          std::ifstream twinNormalsDataFile3(this->userInputs.twinNormalsFile3);
          //read data
          id=n_slip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
          if (twinNormalsDataFile3.is_open())
          while (getline (twinNormalsDataFile3,line) && id<n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0]){
            std::stringstream ss(line);
            ss >> n_alpha_MultiPhase[id][0];
            ss >> n_alpha_MultiPhase[id][1];
            ss >> n_alpha_MultiPhase[id][2];
            n_norm = 0 ;
            n_norm = n_norm + n_alpha_MultiPhase[id][0]*n_alpha_MultiPhase[id][0] ;
            n_norm = n_norm + n_alpha_MultiPhase[id][1]*n_alpha_MultiPhase[id][1] ;
            n_norm = n_norm + n_alpha_MultiPhase[id][2]*n_alpha_MultiPhase[id][2] ;
            n_norm = sqrt(n_norm) ;
            n_alpha_MultiPhase[id][0] = n_alpha_MultiPhase[id][0]/n_norm ;
            n_alpha_MultiPhase[id][1] = n_alpha_MultiPhase[id][1]/n_norm ;
            n_alpha_MultiPhase[id][2] = n_alpha_MultiPhase[id][2]/n_norm ;
            id=id+1;
          }
          else{
            std::cout << "Unable to open twin normals file\n";
            exit(1);
          }

          //open data file to read twin directions
          std::ifstream twinDirectionsDataFile3(this->userInputs.twinDirectionsFile3);
          //read data
          id=n_slip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
          if (twinDirectionsDataFile3.is_open())
          //read data
          while (getline (twinDirectionsDataFile3,line)&& id<n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0]){
            std::stringstream ss(line);
            ss >> m_alpha_MultiPhase[id][0];
            ss >> m_alpha_MultiPhase[id][1];
            ss >> m_alpha_MultiPhase[id][2];
            m_norm = 0 ;
            m_norm = m_norm + m_alpha_MultiPhase[id][0]*m_alpha_MultiPhase[id][0] ;
            m_norm = m_norm + m_alpha_MultiPhase[id][1]*m_alpha_MultiPhase[id][1] ;
            m_norm = m_norm + m_alpha_MultiPhase[id][2]*m_alpha_MultiPhase[id][2] ;
            m_norm = sqrt(m_norm) ;
            m_alpha_MultiPhase[id][0] = m_alpha_MultiPhase[id][0]/m_norm ;
            m_alpha_MultiPhase[id][1] = m_alpha_MultiPhase[id][1]/m_norm ;
            m_alpha_MultiPhase[id][2] = m_alpha_MultiPhase[id][2]/m_norm ;

            id=id+1;
          }
          else{
            std::cout << "Unable to open twin directions file \n";
            exit(1);
          }
        }
        else{
          n_alpha_MultiPhase[id][0]=1;
          n_alpha_MultiPhase[id][1]=0;
          n_alpha_MultiPhase[id][2]=0;
          m_alpha_MultiPhase[id][0]=0;
          m_alpha_MultiPhase[id][1]=1;
          m_alpha_MultiPhase[id][2]=0;
        }



        if (numberofPhases>=4){

          q_phase4.reinit(n_Tslip_systems_MultiPhase[3],n_Tslip_systems_MultiPhase[3]);
          q_phase4=0;
          //open data file to read latent hardening ratios
          std::ifstream latentHardeningratioFile4(this->userInputs.latentHardeningRatioFileName4);
          //read data
          id=0;
          if (latentHardeningratioFile4.is_open()){
            while (getline (latentHardeningratioFile4,line) && id<n_Tslip_systems_Real_MultiPhase[3]){
              std::stringstream ss(line);
              for (unsigned int i=0; i<n_Tslip_systems_Real_MultiPhase[3]; i++){
                ss >> q_phase4[id][i];
              }
              id=id+1;
            }
          }
          else{
            std::cout << "Unable to open latent hardening ratio file 4\n";
            exit(1);
          }


          //open data file to read slip normals
          std::ifstream slipNormalsDataFile4(this->userInputs.slipNormalsFile4);
          //read data
          id=n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
          if (slipNormalsDataFile4.is_open()){
            while (getline (slipNormalsDataFile4,line) && id<n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0]+n_slip_systems_MultiPhase[3]){
              std::stringstream ss(line);
              ss >> n_alpha_MultiPhase[id][0];
              ss >> n_alpha_MultiPhase[id][1];
              ss >> n_alpha_MultiPhase[id][2];
              n_norm = 0 ;
              n_norm = n_norm + n_alpha_MultiPhase[id][0]*n_alpha_MultiPhase[id][0] ;
              n_norm = n_norm + n_alpha_MultiPhase[id][1]*n_alpha_MultiPhase[id][1] ;
              n_norm = n_norm + n_alpha_MultiPhase[id][2]*n_alpha_MultiPhase[id][2] ;
              n_norm = sqrt(n_norm) ;
              n_alpha_MultiPhase[id][0] = n_alpha_MultiPhase[id][0]/n_norm ;
              n_alpha_MultiPhase[id][1] = n_alpha_MultiPhase[id][1]/n_norm ;
              n_alpha_MultiPhase[id][2] = n_alpha_MultiPhase[id][2]/n_norm ;
              id=id+1;
            }
          }
          else{
            std::cout << "Unable to open slip normals file \n";
            exit(1);
          }

          //open data file to read slip directions
          std::ifstream slipDirectionsDataFile4(this->userInputs.slipDirectionsFile4);
          //read data
          id=n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
          if (slipDirectionsDataFile4.is_open()){
            while (getline (slipDirectionsDataFile4,line)&& id<n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0]+n_slip_systems_MultiPhase[3]){
              std::stringstream ss(line);
              ss >> m_alpha_MultiPhase[id][0];
              ss >> m_alpha_MultiPhase[id][1];
              ss >> m_alpha_MultiPhase[id][2];
              m_norm = 0 ;
              m_norm = m_norm + m_alpha_MultiPhase[id][0]*m_alpha_MultiPhase[id][0] ;
              m_norm = m_norm + m_alpha_MultiPhase[id][1]*m_alpha_MultiPhase[id][1] ;
              m_norm = m_norm + m_alpha_MultiPhase[id][2]*m_alpha_MultiPhase[id][2] ;
              m_norm = sqrt(m_norm) ;
              m_alpha_MultiPhase[id][0] = m_alpha_MultiPhase[id][0]/m_norm ;
              m_alpha_MultiPhase[id][1] = m_alpha_MultiPhase[id][1]/m_norm ;
              m_alpha_MultiPhase[id][2] = m_alpha_MultiPhase[id][2]/m_norm ;

              id=id+1;
            }
          }
          else{
            std::cout << "Unable to open slip directions file \n";
            exit(1);
          }

          if(this->userInputs.enableTwinning4){
            //open data file to read twin normals
            std::ifstream twinNormalsDataFile4(this->userInputs.twinNormalsFile4);
            //read data
            id=n_slip_systems_MultiPhase[3]+n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
            if (twinNormalsDataFile4.is_open())
            while (getline (twinNormalsDataFile4,line) && id<n_Tslip_systems_MultiPhase[3]+n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0]){
              std::stringstream ss(line);
              ss >> n_alpha_MultiPhase[id][0];
              ss >> n_alpha_MultiPhase[id][1];
              ss >> n_alpha_MultiPhase[id][2];
              n_norm = 0 ;
              n_norm = n_norm + n_alpha_MultiPhase[id][0]*n_alpha_MultiPhase[id][0] ;
              n_norm = n_norm + n_alpha_MultiPhase[id][1]*n_alpha_MultiPhase[id][1] ;
              n_norm = n_norm + n_alpha_MultiPhase[id][2]*n_alpha_MultiPhase[id][2] ;
              n_norm = sqrt(n_norm) ;
              n_alpha_MultiPhase[id][0] = n_alpha_MultiPhase[id][0]/n_norm ;
              n_alpha_MultiPhase[id][1] = n_alpha_MultiPhase[id][1]/n_norm ;
              n_alpha_MultiPhase[id][2] = n_alpha_MultiPhase[id][2]/n_norm ;
              id=id+1;
            }
            else{
              std::cout << "Unable to open twin normals file\n";
              exit(1);
            }

            //open data file to read twin directions
            std::ifstream twinDirectionsDataFile4(this->userInputs.twinDirectionsFile4);
            //read data
            id=n_slip_systems_MultiPhase[3]+n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0];
            if (twinDirectionsDataFile4.is_open())
            //read data
            while (getline (twinDirectionsDataFile4,line)&& id<n_Tslip_systems_MultiPhase[3]+n_Tslip_systems_MultiPhase[2]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[0]){
              std::stringstream ss(line);
              ss >> m_alpha_MultiPhase[id][0];
              ss >> m_alpha_MultiPhase[id][1];
              ss >> m_alpha_MultiPhase[id][2];
              m_norm = 0 ;
              m_norm = m_norm + m_alpha_MultiPhase[id][0]*m_alpha_MultiPhase[id][0] ;
              m_norm = m_norm + m_alpha_MultiPhase[id][1]*m_alpha_MultiPhase[id][1] ;
              m_norm = m_norm + m_alpha_MultiPhase[id][2]*m_alpha_MultiPhase[id][2] ;
              m_norm = sqrt(m_norm) ;
              m_alpha_MultiPhase[id][0] = m_alpha_MultiPhase[id][0]/m_norm ;
              m_alpha_MultiPhase[id][1] = m_alpha_MultiPhase[id][1]/m_norm ;
              m_alpha_MultiPhase[id][2] = m_alpha_MultiPhase[id][2]/m_norm ;

              id=id+1;
            }
            else{
              std::cout << "Unable to open twin directions file \n";
              exit(1);
            }
          }
          else{
            n_alpha_MultiPhase[id][0]=1;
            n_alpha_MultiPhase[id][1]=0;
            n_alpha_MultiPhase[id][2]=0;
            m_alpha_MultiPhase[id][0]=0;
            m_alpha_MultiPhase[id][1]=1;
            m_alpha_MultiPhase[id][2]=0;
          }



        }



      }



    }

    Dmat_MultiPhase.reinit(6*this->userInputs.numberofPhases,6*this->userInputs.numberofPhases); Dmat=0.0;

    for(unsigned int i=0;i<6;i++){
      for(unsigned int j=0;j<6;j++){
        Dmat_MultiPhase[i][j] = Dmat_SinglePhase[i][j];
      }
    }

    if (numberofPhases>=2){

      for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
          Dmat_MultiPhase[i+6][j] = this->userInputs.elasticStiffness2[i][j];
        }
      }

      if (numberofPhases>=3){

        for(unsigned int i=0;i<6;i++){
          for(unsigned int j=0;j<6;j++){
            Dmat_MultiPhase[i+6*2][j] = this->userInputs.elasticStiffness3[i][j];
          }
        }

        if (numberofPhases>=4){

          for(unsigned int i=0;i<6;i++){
            for(unsigned int j=0;j<6;j++){
              Dmat_MultiPhase[i+6*3][j] = this->userInputs.elasticStiffness4[i][j];
            }
          }

        }
      }


    }


    s0_init1.reinit(Max_n_Tslip_systems_MultiPhase);
    W_kh_init1.reinit(Max_n_Tslip_systems_MultiPhase);
    twin_init1.resize(Max_n_twin_systems_MultiPhase);
    slip_init1.resize(Max_n_slip_systems_MultiPhase);


    for (unsigned int i=0;i<Max_n_Tslip_systems_MultiPhase;i++){
      s0_init1(i)=0;
    }

    for (unsigned int i=0;i<Max_n_slip_systems_MultiPhase;i++){
      slip_init1[i]=0.0;
    }


    for (unsigned int i = 0;i<Max_n_Tslip_systems_MultiPhase;i++) {
      W_kh_init1[i] = 0.0;
    }

    for (unsigned int i=0;i<Max_n_twin_systems_MultiPhase;i++){
      twin_init1[i]=0.0;
    }

    if (this->userInputs.enableUserMaterialModel){
      stateVar_init1.reinit(Max_n_UserMatStateVar_MultiPhase);
      for (unsigned int i=0;i<Max_n_UserMatStateVar_MultiPhase;i++){
        stateVar_init1(i)=0.0;
      }
      stateVar_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,stateVar_init1));
      stateVar_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,stateVar_init1));
      for (unsigned int i=0;i<n_UserMatStateVar_MultiPhase[0];i++){
        stateVar_init1(i)=stateVar_init(i);
      }
    }




    Fp_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fp_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    CauchyStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));
    F_lastIter_Global.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));

    FirstPiolaStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));
    SecondPiolaStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));
    workDensity1.reinit(num_local_cells); workDensity1 = 0.0;
    workDensity2.reinit(num_local_cells); workDensity2 = 0.0;
    workDensityTotal1.reinit(num_local_cells); workDensityTotal1 = 0.0;
    workDensityTotal2.reinit(num_local_cells); workDensityTotal2 = 0.0;
    workDensity1_Tr.reinit(num_local_cells); workDensity1_Tr = 0.0;
    workDensity2_Tr.reinit(num_local_cells); workDensity2_Tr = 0.0;
    workDensityTotal1_Tr.reinit(num_local_cells); workDensityTotal1_Tr = 0.0;
    workDensityTotal2_Tr.reinit(num_local_cells); workDensityTotal2_Tr = 0.0;
	TinterStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,TinterStress_init));
	TinterStress_diff.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,TinterStress_diff_init));
    s_alpha_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init1));
    W_kh_conv.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, W_kh_init1));
    W_kh_iter.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, W_kh_init1));
    s_alpha_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init1));
    twinfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init1));
    slipfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,slip_init1));
    twinfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init1));
    slipfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,slip_init1));
    twin_ouput.resize(num_local_cells, std::vector<double>(num_quad_points,0.0));
    twin_conv.resize(num_local_cells,std::vector<unsigned int>(num_quad_points,0));
    twin_iter.resize(num_local_cells,std::vector<unsigned int>(num_quad_points,0));


    for (unsigned int i=0;i<n_Tslip_systems_MultiPhase[0];i++){
      s0_init1(i)=s0_init(i);
    }



    if (numberofPhases>=2){
      s0_init2.reinit(Max_n_Tslip_systems_MultiPhase);

      for (unsigned int i=0;i<Max_n_Tslip_systems_MultiPhase;i++){
        s0_init2(i)=0;
      }

      for (unsigned int i=0;i<n_slip_systems_MultiPhase[1];i++){
        s0_init2(i)=this->userInputs.initialSlipResistance2[i];
      }

      for (unsigned int i=0;i<n_twin_systems_MultiPhase[1];i++){
        s0_init2(i+n_slip_systems_MultiPhase[1])=this->userInputs.initialSlipResistanceTwin2[i];
      }

      if (this->userInputs.enableUserMaterialModel){
        stateVar_init2.reinit(Max_n_UserMatStateVar_MultiPhase);
        for (unsigned int i=0;i<Max_n_UserMatStateVar_MultiPhase;i++){
          stateVar_init2(i)=0.0;
        }
        if (this->userInputs.enableUserMaterialModel2){
          if (this->userInputs.numberofUserMatStateVar2==0){
            stateVar_init2(1)=0;
          }
          else{
            for (unsigned int i=0;i<n_UserMatStateVar_MultiPhase[1];i++){
              stateVar_init2(i)=this->userInputs.UserMatStateVar2[i];
            }
          }
        }
      }

      if (numberofPhases>=3){
        s0_init3.reinit(Max_n_Tslip_systems_MultiPhase);

        for (unsigned int i=0;i<Max_n_Tslip_systems_MultiPhase;i++){
          s0_init3(i)=0;
        }

        for (unsigned int i=0;i<n_slip_systems_MultiPhase[2];i++){
          s0_init3(i)=this->userInputs.initialSlipResistance3[i];
        }

        for (unsigned int i=0;i<n_twin_systems_MultiPhase[2];i++){
          s0_init3(i+n_slip_systems_MultiPhase[2])=this->userInputs.initialSlipResistanceTwin3[i];
        }

        if (this->userInputs.enableUserMaterialModel){
          stateVar_init3.reinit(Max_n_UserMatStateVar_MultiPhase);
          for (unsigned int i=0;i<Max_n_UserMatStateVar_MultiPhase;i++){
            stateVar_init3(i)=0.0;
          }
          if (this->userInputs.enableUserMaterialModel3){
            if (this->userInputs.numberofUserMatStateVar3==0){
              stateVar_init3(1)=0;
            }
            else{
              for (unsigned int i=0;i<n_UserMatStateVar_MultiPhase[2];i++){
                stateVar_init3(i)=this->userInputs.UserMatStateVar3[i];
              }
            }
          }
        }

        if (numberofPhases>=4){
          s0_init4.reinit(Max_n_Tslip_systems_MultiPhase);

          for (unsigned int i=0;i<Max_n_Tslip_systems_MultiPhase;i++){
            s0_init4(i)=0;
          }

          for (unsigned int i=0;i<n_slip_systems_MultiPhase[3];i++){
            s0_init4(i)=this->userInputs.initialSlipResistance4[i];
          }

          for (unsigned int i=0;i<n_twin_systems_MultiPhase[3];i++){
            s0_init4(i+n_slip_systems_MultiPhase[3])=this->userInputs.initialSlipResistanceTwin4[i];
          }

          if (this->userInputs.enableUserMaterialModel){
            stateVar_init4.reinit(Max_n_UserMatStateVar_MultiPhase);
            for (unsigned int i=0;i<Max_n_UserMatStateVar_MultiPhase;i++){
              stateVar_init4(i)=0.0;
            }
            if (this->userInputs.enableUserMaterialModel4){
              if (this->userInputs.numberofUserMatStateVar4==0){
                stateVar_init4(1)=0;
              }
              else{
                for (unsigned int i=0;i<n_UserMatStateVar_MultiPhase[3];i++){
                  stateVar_init4(i)=this->userInputs.UserMatStateVar4[i];
                }
              }
            }
          }
        }
      }
    }




    for (unsigned int cell=0; cell<num_local_cells; cell++){
      for (unsigned int q=0; q<num_quad_points; q++){
        phaseMaterial=phase[cell][q];
        if (phaseMaterial==1){
          for (unsigned int i=0; i<Max_n_Tslip_systems_MultiPhase; i++){
            s_alpha_conv[cell][q][i]=s0_init1(i);
            s_alpha_iter[cell][q][i]=s0_init1(i);
          }
          if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel1)){
            for (unsigned int i=0; i<Max_n_UserMatStateVar_MultiPhase; i++){
              stateVar_conv[cell][q][i]=stateVar_init1(i);
              stateVar_iter[cell][q][i]=stateVar_init1(i);
            }
          }

        }
        else if (phaseMaterial==2){
          for (unsigned int i=0; i<Max_n_Tslip_systems_MultiPhase; i++){
            s_alpha_conv[cell][q][i]=s0_init2(i);
            s_alpha_iter[cell][q][i]=s0_init2(i);
          }
          if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel2)){
            for (unsigned int i=0; i<Max_n_UserMatStateVar_MultiPhase; i++){
              stateVar_conv[cell][q][i]=stateVar_init2(i);
              stateVar_iter[cell][q][i]=stateVar_init2(i);
            }
          }
        }
        else if (phaseMaterial==3){
          for (unsigned int i=0; i<Max_n_Tslip_systems_MultiPhase; i++){
            s_alpha_conv[cell][q][i]=s0_init3(i);
            s_alpha_iter[cell][q][i]=s0_init3(i);
          }
          if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel3)){
            for (unsigned int i=0; i<Max_n_UserMatStateVar_MultiPhase; i++){
              stateVar_conv[cell][q][i]=stateVar_init3(i);
              stateVar_iter[cell][q][i]=stateVar_init3(i);
            }
          }
        }
        else if (phaseMaterial==4){
          for (unsigned int i=0; i<Max_n_Tslip_systems_MultiPhase; i++){
            s_alpha_conv[cell][q][i]=s0_init4(i);
            s_alpha_iter[cell][q][i]=s0_init4(i);
          }
          if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel4)){
            for (unsigned int i=0; i<Max_n_UserMatStateVar_MultiPhase; i++){
              stateVar_conv[cell][q][i]=stateVar_init4(i);
              stateVar_iter[cell][q][i]=stateVar_init4(i);
            }
          }
        }
      }
    }


  }

  N_qpts=num_quad_points;
  initCalled=true;

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
