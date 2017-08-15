
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
    
    unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();
    F.reinit(dim, dim);
    
    n_slip_systems1=numSlipSystems1;
    n_slip_systems2=numSlipSystems2;
    n_slip_systems2+=numTwinSystems;
    n_twin_systems=numTwinSystems;
    
    
    
    n_alpha1.reinit(n_slip_systems1,3);
    m_alpha1.reinit(n_slip_systems1,3);
    string line;
    
    //open data file to read slip normals
    ifstream slipNormalsDataFile(slipDirectionsFile1);
    //read data
    unsigned int id=0;
    if (slipNormalsDataFile.is_open()){
        cout << "reading slip Normals file\n";
        //read data
        while (getline (slipNormalsDataFile,line) && id<n_slip_systems1){
            stringstream ss(line);
            ss >> n_alpha1[id][0];
            ss >> n_alpha1[id][1];
            ss >> n_alpha1[id][2];
            //cout<<id<<'\t'<<n_alpha[id][0]<<'\t'<<n_alpha[id][1]<<'\t'<<n_alpha[id][2]<<'\n';
            id=id+1;
        }
    }
    else{
        cout << "Unable to open slipNormals.txt \n";
        exit(1);
    }
    
    //open data file to read slip directions
    ifstream slipDirectionsDataFile(slipNormalsFile1);
    //read data
    id=0;
    if (slipDirectionsDataFile.is_open()){
        cout << "reading slip Directions file\n";
        //read data
        while (getline (slipDirectionsDataFile,line)&& id<n_slip_systems1){
            stringstream ss(line);
            ss >> m_alpha1[id][0];
            ss >> m_alpha1[id][1];
            ss >> m_alpha1[id][2];
            //cout<<id<<'\t'<<m_alpha[id][0]<<'\t'<<m_alpha[id][1]<<'\t'<<m_alpha[id][2]<<'\n';
            id=id+1;
        }
    }
    else{
        cout << "Unable to open slipDirections.txt \n";
        exit(1);
    }
    
    
    
    
    //q is a parameter in the hardening model
    q1.reinit(n_slip_systems1,n_slip_systems1);
    for(unsigned int i=0;i<n_slip_systems1;i++){
        for(unsigned int j=0;j<n_slip_systems1;j++){
            q1[i][j] = latentHardeningRatio1;
        }
    }
    
    for(unsigned int i=0;i<n_slip_systems1;i++){
        q1[i][i] = 1.0;
    }
    
    //Elastic Stiffness Matrix Dmat
    Dmat11.reinit(6,6); Dmat11=0.0;
    
    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
            Dmat11[i][j] = elasticStiffness1[i][j];
        }
    }
    
    
    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=3;j<6;j++){
            Dmat11[i][j] = 2*Dmat11[i][j];
        }
    }
    
    Vector<double> s0_init1 (n_slip_systems1);
    
    for (unsigned int i=0;i<n_slip_systems1;i++){
        s0_init1(i)=initialSlipResistance1[i];
    }
    
    
    std::vector<double> slip_init1(numSlipSystems1);
    
    for (unsigned int i=0;i<numSlipSystems1;i++){
        slip_init1[i]=0.0;
    }
    
    
    
    
    
    
   
    n_alpha2.reinit(n_slip_systems2,3);
      m_alpha2.reinit(n_slip_systems2,3);

      
      //open data file to read slip normals
      slipNormalsDataFile.open(slipDirectionsFile2);
      //read data
       id=0;
      if (slipNormalsDataFile.is_open()){
	//cout << "reading slip Normals file\n";
	//read data
	while (getline (slipNormalsDataFile,line) && id<numSlipSystems2){
	  stringstream ss(line);
	  ss >> n_alpha2[id][0];
	  ss >> n_alpha2[id][1];
	  ss >> n_alpha2[id][2];
	  //cout<<id<<'\t'<<n_alpha[id][0]<<'\t'<<n_alpha[id][1]<<'\t'<<n_alpha[id][2]<<'\n';
	  id=id+1;
	}
      }
      else{
	cout << "Unable to open slipNormals.txt \n";
	exit(1);
      }
      
      //open data file to read slip directions
      slipDirectionsDataFile.open(slipNormalsFile2);
      //read data
      id=0;
      if (slipDirectionsDataFile.is_open()){
	//cout << "reading slip Directions file\n";
	//read data
	while (getline (slipDirectionsDataFile,line)&& id<numSlipSystems2){
	  stringstream ss(line);
	  ss >> m_alpha2[id][0];
	  ss >> m_alpha2[id][1];
	  ss >> m_alpha2[id][2];
	  //cout<<id<<'\t'<<m_alpha[id][0]<<'\t'<<m_alpha[id][1]<<'\t'<<m_alpha[id][2]<<'\n';
	  id=id+1;
	}
      }
      else{
	cout << "Unable to open slipDirections.txt \n";
	exit(1);
      }


	//open data file to read twin normals
      ifstream twinNormalsDataFile(twinDirectionsFile);
      //read data
      id=numSlipSystems2;
      if (twinNormalsDataFile.is_open()){
	cout << "reading slip Normals file\n";
	//read data
	while (getline (twinNormalsDataFile,line) && id<n_slip_systems2){
	  stringstream ss(line);
	  ss >> n_alpha2[id][0];
	  ss >> n_alpha2[id][1];
	  ss >> n_alpha2[id][2];
	  //cout<<id<<'\t'<<n_alpha[id][0]<<'\t'<<n_alpha[id][1]<<'\t'<<n_alpha[id][2]<<'\n';
	  id=id+1;
	}
      }
      else{
	cout << "Unable to open slipNormals.txt \n";
	exit(1);
      }
      
      //open data file to read twin directions
      ifstream twinDirectionsDataFile(twinNormalsFile);
      //read data
      id=numSlipSystems2;
      if (twinDirectionsDataFile.is_open()){
	cout << "reading slip Directions file\n";
	//read data
	while (getline (twinDirectionsDataFile,line)&& id<n_slip_systems2){
	  stringstream ss(line);
	  ss >> m_alpha2[id][0];
	  ss >> m_alpha2[id][1];
	  ss >> m_alpha2[id][2];
	  //cout<<id<<'\t'<<m_alpha[id][0]<<'\t'<<m_alpha[id][1]<<'\t'<<m_alpha[id][2]<<'\n';
	  id=id+1;
	}
      }
      else{
	cout << "Unable to open slipDirections.txt \n";
	exit(1);
      }
      
        
    
    //q is a parameter in the hardening model
    q2.reinit(n_slip_systems2,n_slip_systems2);
    for(unsigned int i=0;i<n_slip_systems2;i++){
        for(unsigned int j=0;j<n_slip_systems2;j++){
            q2[i][j] = latentHardeningRatio2;
        }
    }
    
    for(unsigned int i=0;i<n_slip_systems2;i++){
        q2[i][i] = 1.0;
    }
    
    //Basal Slip
    for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            q2[i][j] = 1.0;
        }
    }
    
    //Elastic Stiffness Matrix Dmat
    Dmat12.reinit(6,6); Dmat12=0.0;
    
    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
            Dmat12[i][j] = elasticStiffness2[i][j];
        }
    }
    
    
    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=3;j<6;j++){
            Dmat12[i][j] = 2*Dmat12[i][j];
        }
    }
    
    Vector<double> s0_init2 (n_slip_systems2),rot_init(dim),rotnew_init(dim);
    std::vector<double> twin_init(numTwinSystems),slip_init2(numSlipSystems2);
    
    for (unsigned int i=0;i<numSlipSystems2;i++){
        s0_init2(i)=initialSlipResistance2[i];
    }
    
    for (unsigned int i=0;i<numTwinSystems;i++){
        s0_init2(i+numSlipSystems2)=initialSlipResistanceTwin[i];
    }
    
    
    for (unsigned int i=0;i<numSlipSystems2;i++){
        slip_init2[i]=0.0;
    }
    
    for (unsigned int i=0;i<numTwinSystems;i++){
        twin_init[i]=0.0;
    }
    
    for (unsigned int i=0;i<dim;i++){
        rot_init(i)=0.0;
        rotnew_init(i)=0.0;
    }
    

    //Resize the vectors of history variables
    Fp_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_conv2.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init2));
    Fp_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_iter2.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init2));
    twinfraction_iter.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,twin_init));
    slipfraction_iter2.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,slip_init2));
    twinfraction_conv.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,twin_init));
    slipfraction_conv2.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,slip_init2));
    rot.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rot_init));
    rotnew.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rotnew_init));
    twin.resize(num_local_cells,std::vector<double>(num_quad_points,0.0));
    phaseID.resize(num_local_cells,std::vector<double>(num_quad_points,1.0));
    
    s_alpha_conv1.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init1));
    s_alpha_iter1.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init1));
    slipfraction_iter1.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,slip_init1));
    slipfraction_conv1.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,slip_init1));
   
    
    
    //load rot and rotnew
    for (unsigned int cell=0; cell<num_local_cells; cell++){
        for (unsigned int q=0; q<num_quad_points; q++){
            unsigned int materialID=quadratureOrientationsMap[cell][q];
            for (unsigned int i=0; i<dim; i++){
                rot[cell][q][i]=orientations.eulerAngles[materialID][i];
                rotnew[cell][q][i]=orientations.eulerAngles[materialID][i];
            }
            phaseID[cell][q]=orientations.eulerAngles[materialID][dim];
        }
    }
    N_qpts=num_quad_points;
    initCalled=true;
    

}
