#include "../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

template <int dim>
void crystalPlasticity<dim>::init2(unsigned int num_quad_points)
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

  FullMatrix<double> CauchyStress_init(dim,dim);
  CauchyStress_init=0;

  unsigned int n_twin_systems = 0;
  double m_norm , n_norm ;
  unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();
  F.reinit(dim, dim);

  unsigned int n_slip_systemsWOtwin = this->userInputs.numSlipSystems1;
  unsigned int n_slip_systems= n_slip_systemsWOtwin;
  if(this->userInputs.enableTwinning1){
    n_slip_systems+=this->userInputs.numTwinSystems1*2;
    n_twin_systems =this->userInputs.numTwinSystems1*2;
  }
  else{
    n_slip_systems+=1;
    n_twin_systems =1;
  }

  n_alpha.reinit(n_slip_systems,3);
  m_alpha.reinit(n_slip_systems,3);

  std::string line;
  q_phase1.reinit(n_slip_systems,n_slip_systems);
  q_phase1=0;
  //open data file to read latent hardening ratios
  std::ifstream latentHardeningratioFile(this->userInputs.latentHardeningRatioFileName1);
  //read data
  unsigned int id=0;
  if (latentHardeningratioFile.is_open()){
    while (getline (latentHardeningratioFile,line) && id<n_slip_systems){
      std::stringstream ss(line);
      for (unsigned int i=0; i<n_slip_systems; i++){
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
    while (getline (slipNormalsDataFile,line) && id<n_slip_systemsWOtwin){
      std::stringstream ss(line);
      ss >> n_alpha[id][0];
      ss >> n_alpha[id][1];
      ss >> n_alpha[id][2];
      n_norm = 0 ;
      n_norm = n_norm + n_alpha[id][0]*n_alpha[id][0] ;
      n_norm = n_norm + n_alpha[id][1]*n_alpha[id][1] ;
      n_norm = n_norm + n_alpha[id][2]*n_alpha[id][2] ;
      n_norm = sqrt(n_norm) ;
      n_alpha[id][0] = n_alpha[id][0]/n_norm ;
      n_alpha[id][1] = n_alpha[id][1]/n_norm ;
      n_alpha[id][2] = n_alpha[id][2]/n_norm ;

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
    while (getline (slipDirectionsDataFile,line)&& id<n_slip_systemsWOtwin){
      std::stringstream ss(line);
      ss >> m_alpha[id][0];
      ss >> m_alpha[id][1];
      ss >> m_alpha[id][2];
      m_norm = 0 ;
      m_norm = m_norm + m_alpha[id][0]*m_alpha[id][0] ;
      m_norm = m_norm + m_alpha[id][1]*m_alpha[id][1] ;
      m_norm = m_norm + m_alpha[id][2]*m_alpha[id][2] ;
      m_norm = sqrt(m_norm) ;
      m_alpha[id][0] = m_alpha[id][0]/m_norm ;
      m_alpha[id][1] = m_alpha[id][1]/m_norm ;
      m_alpha[id][2] = m_alpha[id][2]/m_norm ;

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
    id= n_slip_systemsWOtwin;
    if (twinNormalsDataFile.is_open()){
      while (getline (twinNormalsDataFile,line) && id<n_slip_systems){
        std::stringstream ss(line);
        ss >> n_alpha[id][0];
        ss >> n_alpha[id][1];
        ss >> n_alpha[id][2];
         n_norm = 0 ;
         n_norm = n_norm + n_alpha[id][0]*n_alpha[id][0] ;
         n_norm = n_norm + n_alpha[id][1]*n_alpha[id][1] ;
         n_norm = n_norm + n_alpha[id][2]*n_alpha[id][2] ;
         n_norm = sqrt(n_norm) ;
         n_alpha[id][0] = n_alpha[id][0]/n_norm ;
         n_alpha[id][1] = n_alpha[id][1]/n_norm ;
         n_alpha[id][2] = n_alpha[id][2]/n_norm ;

        id=id+1;
      }
      for(unsigned int i=n_slip_systemsWOtwin;i<this->userInputs.numTwinSystems1+n_slip_systemsWOtwin;i++){
        for(unsigned int j=0;j<dim;j++){
          n_alpha[i+this->userInputs.numTwinSystems1][j]=n_alpha[i][j];
        }
      }
    }
    else{
      std::cout << "Unable to open twin normals file\n";
      exit(1);
    }

    //open data file to read twin directions
    std::ifstream twinDirectionsDataFile(this->userInputs.twinDirectionsFile1);
    //read data
    id= n_slip_systemsWOtwin;
    if (twinDirectionsDataFile.is_open()){
      //read data
      while (getline (twinDirectionsDataFile,line)&& id<n_slip_systems){
        std::stringstream ss(line);
        ss >> m_alpha[id][0];
        ss >> m_alpha[id][1];
        ss >> m_alpha[id][2];
        m_norm = 0 ;
        m_norm = m_norm + m_alpha[id][0]*m_alpha[id][0] ;
        m_norm = m_norm + m_alpha[id][1]*m_alpha[id][1] ;
        m_norm = m_norm + m_alpha[id][2]*m_alpha[id][2] ;
        m_norm = sqrt(m_norm) ;
        m_alpha[id][0] = m_alpha[id][0]/m_norm ;
        m_alpha[id][1] = m_alpha[id][1]/m_norm ;
        m_alpha[id][2] = m_alpha[id][2]/m_norm ;

        id=id+1;
      }
      for(unsigned int i=n_slip_systemsWOtwin;i<this->userInputs.numTwinSystems1+n_slip_systemsWOtwin;i++){
        for(unsigned int j=0;j<dim;j++){
          m_alpha[i+this->userInputs.numTwinSystems1][j]=m_alpha[i][j];
        }
      }
    }
    else{
      std::cout << "Unable to open twin directions file \n";
      exit(1);
    }
  }
  else{
    n_alpha[id][0]=1;
    n_alpha[id][1]=0;
    n_alpha[id][2]=0;
    m_alpha[id][0]=0;
    m_alpha[id][1]=1;
    m_alpha[id][2]=0;
  }

  //q is a parameter in the hardening model
  //q definition is transfered in calculatePlasticity.cc

  //Elastic Stiffness Matrix Dmat
  Dmat.reinit(6,6); Dmat=0.0;

  for(unsigned int i=0;i<6;i++){
    for(unsigned int j=0;j<6;j++){
      Dmat[i][j] = this->userInputs.elasticStiffness1[i][j];
    }
  }


  for(unsigned int i=0;i<6;i++){
    for(unsigned int j=3;j<6;j++){
      Dmat[i][j] = 2*Dmat[i][j];
    }
  }

  Vector<double> s0_init (n_slip_systems*((n_twin_systems / 2) + 1)),rot_init(dim*((n_twin_systems / 2) + 1)),rotnew_init(dim*((n_twin_systems / 2) + 1));
  Vector<double> stateVar_init;
  std::vector<double> twin_init(n_twin_systems/2),TwinOutput_init(n_twin_systems*2),slip_init(n_slip_systemsWOtwin*((n_twin_systems / 2) + 1));
  std::vector<unsigned int> twin_init2(n_twin_systems / 2);
  for (unsigned int i=0;i<n_slip_systemsWOtwin;i++){
    for (unsigned int j = 0;j < (n_twin_systems / 2) + 1;j++) {
      s0_init(i+ n_slip_systems*j) = this->userInputs.initialSlipResistance1[i];
    }
  }

  for (unsigned int i=0;i<n_twin_systems;i++){
    for (unsigned int j = 0;j < (n_twin_systems / 2) + 1;j++) {
      s0_init(i + n_slip_systemsWOtwin + n_slip_systems*j) = this->userInputs.initialSlipResistanceTwin1[i];
    }
  }


  for (unsigned int i=0;i<n_slip_systemsWOtwin;i++){
    for (unsigned int j = 0;j < (n_twin_systems / 2) + 1;j++) {
      slip_init[i+j*n_slip_systemsWOtwin] = 0.0;
    }
  }

  for (unsigned int i=0;i<n_twin_systems/2;i++){
    twin_init[i]=0.0;
    twin_init2[i] = 0;
  }
  for (unsigned int i=0;i<n_twin_systems*2;i++){
    TwinOutput_init[i] = 0;
  }


  for (unsigned int i=0;i<dim*((n_twin_systems / 2) + 1);i++){
    rot_init(i)=0.0;
    rotnew_init(i)=0.0;
  }
  FullMatrix<double> Fp_conv_init(dim, (dim)*((n_twin_systems / 2) + 1));

  for (unsigned int i = 0;i < (n_twin_systems / 2) + 1;i++) {
    unsigned int alpha = (i)*(dim);
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        if (j == k) {
          Fp_conv_init[j][k+alpha] = 1;
        }
        else {
          Fp_conv_init[j][k+alpha] = 0;
        }
      }
    }
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


  //Resize the vectors of history variables
  Fp_conv.resize(num_local_cells, std::vector<FullMatrix<double> >(num_quad_points, Fp_conv_init));
  Fe_conv.resize(num_local_cells, std::vector<FullMatrix<double> >(num_quad_points, Fp_conv_init));
  Fp_iter.resize(num_local_cells, std::vector<FullMatrix<double> >(num_quad_points, Fp_conv_init));
  Fe_iter.resize(num_local_cells, std::vector<FullMatrix<double> >(num_quad_points, Fp_conv_init));
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
  s_alpha_conv.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, s0_init));
  s_alpha_iter.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, s0_init));
  slipfraction_iter.resize(num_local_cells, std::vector<std::vector<double> >(num_quad_points, slip_init));
  slipfraction_conv.resize(num_local_cells, std::vector<std::vector<double> >(num_quad_points, slip_init));
  TwinOutputfraction_iter.resize(num_local_cells, std::vector<std::vector<double> >(num_quad_points, TwinOutput_init));
  TwinOutputfraction_conv.resize(num_local_cells, std::vector<std::vector<double> >(num_quad_points, TwinOutput_init));
  rot.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, rot_init));
  rotnew_conv.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, rotnew_init));
  rotnew_iter.resize(num_local_cells, std::vector<Vector<double> >(num_quad_points, rotnew_init));

  twinfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
  twinfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
  twin_ouput.resize(num_local_cells, std::vector<double>(num_quad_points,0.0));
  TwinMaxFlag_conv.resize(num_local_cells, std::vector<unsigned int>(num_quad_points, 1));
  TwinFlag_conv.resize(num_local_cells, std::vector<std::vector<unsigned int> >(num_quad_points, twin_init2));
  ActiveTwinSystems_conv.resize(num_local_cells, std::vector<std::vector<unsigned int> >(num_quad_points, twin_init2));
  NumberOfTwinnedRegion_conv.resize(num_local_cells, std::vector<unsigned int>(num_quad_points, 0));
  TotaltwinvfK.resize(num_local_cells, std::vector<double>(num_quad_points, 0));
  TwinMaxFlag_iter.resize(num_local_cells, std::vector<unsigned int>(num_quad_points, 1));
  TwinFlag_iter.resize(num_local_cells, std::vector<std::vector<unsigned int> >(num_quad_points, twin_init2));
  ActiveTwinSystems_iter.resize(num_local_cells, std::vector<std::vector<unsigned int> >(num_quad_points, twin_init2));
  NumberOfTwinnedRegion_iter.resize(num_local_cells, std::vector<unsigned int>(num_quad_points, 0));

  if (this->userInputs.enableUserMaterialModel){
      stateVar_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,stateVar_init));
      stateVar_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,stateVar_init));
  }

  double s0_twin=this->userInputs.initialSlipResistanceTwin1[n_twin_systems-1];
  for (unsigned int cell = 0; cell<num_local_cells; cell++) {
    for (unsigned int region = 1; region<(n_twin_systems / 2) + 1; region++) {
      for (unsigned int q = 0; q<num_quad_points; q++) {
        s_alpha_conv[cell][q][region*n_slip_systems+n_slip_systemsWOtwin] = s0_twin;
        s_alpha_iter[cell][q][region*n_slip_systems+n_slip_systemsWOtwin] = s0_twin;
        s_alpha_conv[cell][q][region*n_slip_systems+n_slip_systemsWOtwin+1] = s0_twin;
        s_alpha_iter[cell][q][region*n_slip_systems+n_slip_systemsWOtwin+1] = s0_twin;
        for (unsigned int i = n_slip_systemsWOtwin + 2; i < n_slip_systems; i++) {
          s_alpha_conv[cell][q][region*n_slip_systems+i] = 0;
          s_alpha_iter[cell][q][region*n_slip_systems+i] = 0;
        }
      }
    }
  }




  //load rot and rotnew
  for (unsigned int cell=0; cell<num_local_cells; cell++){
    unsigned int materialID=cellOrientationMap[cell];
    for (unsigned int q=0; q<num_quad_points; q++){
      for (unsigned int i = 0; i<dim; i++){
        rot[cell][q][i]=orientations.eulerAngles[materialID][i];
        rotnew_iter[cell][q][i]=orientations.eulerAngles[materialID][i];
        rotnew_conv[cell][q][i] = orientations.eulerAngles[materialID][i];
      }
      for (unsigned int Region = 1; Region<(n_twin_systems / 2) + 1; Region++) {

        Vector<double> quat1(4), rod(3), quat2(4), quatprod(4);
        rod(0) = rot[cell][q][0];rod(1) = rot[cell][q][1];rod(2) = rot[cell][q][2];
        rod2quat(quat2, rod);
        quat1(0) = 0;quat1(1) = n_alpha[n_slip_systemsWOtwin + Region - 1][0];
        quat1(2) = n_alpha[n_slip_systemsWOtwin + Region - 1][1];quat1(3) = n_alpha[n_slip_systemsWOtwin + Region - 1][2];
        quatproduct(quatprod, quat2, quat1);
        quat2rod(quatprod, rod);
        rot[cell][q][Region*dim] = rod(0);rot[cell][q][Region*dim+1] = rod(1);rot[cell][q][Region*dim+2] = rod(2);
        rotnew_iter[cell][q][Region*dim] = rod(0);rotnew_iter[cell][q][Region*dim+1] = rod(1);rotnew_iter[cell][q][Region*dim+2] = rod(2);
        rotnew_conv[cell][q][Region*dim] = rod(0);rotnew_conv[cell][q][Region*dim+1] = rod(1);rotnew_conv[cell][q][Region*dim+2] = rod(2);
      }

    }
  }

  N_qpts=num_quad_points;
  initCalled=true;


}

#include "../../../include/crystalPlasticity_template_instantiations.h"
