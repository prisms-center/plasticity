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
    
    n_slip_systems=properties.n_slip_systems;
    n_twin_systems=properties.n_twin_systems;
    m_alpha.reinit(properties.n_slip_systems,dim);
    n_alpha.reinit(properties.n_slip_systems,dim);
    
    m_alpha.fill(properties.m_alpha);
    n_alpha.fill(properties.n_alpha);
    
    
    //q is a parameter in the hardening model
    q.reinit(n_slip_systems,n_slip_systems);
    for(unsigned int i=0;i<n_slip_systems;i++){
        for(unsigned int j=0;j<n_slip_systems;j++){
            q[i][j] = properties.q1;
        }
    }
    
    for(unsigned int i=0;i<n_slip_systems;i++){
        q[i][i] = properties.q2;
    }
    
    //Basal Slip
    for(unsigned int i=0;i<3;i++){
        for(unsigned int j=0;j<3;j++){
            q[i][j] = properties.q2;
        }
    }
    
    //Elastic Stiffness Matrix Dmat
    Dmat.reinit(6,6); Dmat=0.0;
    Dmat[0][0]=properties.C11; Dmat[0][1]=properties.C12; Dmat[0][2]=properties.C13; Dmat[1][0]=properties.C12; Dmat[1][1]=properties.C11; Dmat[1][2]=properties.C13; Dmat[2][0]=properties.C13; Dmat[2][1]=properties.C13; Dmat[2][2]=properties.C33;
    Dmat[3][3]=2*properties.C44; Dmat[4][4]=2*properties.C44; Dmat[5][5]=(properties.C11-properties.C12);
    
    Vector<double> s0_init (n_slip_systems),rot_init(dim),rotnew_init(dim);
    std::vector<double> twin_init(6),slip_init(18);
    
    for (unsigned int i=0;i<n_slip_systems;i++){
        if(i<3)
            s0_init(i)=properties.s01;
        else if(i<6)
            s0_init(i)=properties.s02;
        else if(i<12)
            s0_init(i)=properties.s03;
        else if(i<18)
            s0_init(i)=properties.s04;
        else if(i<24)
            s0_init(i)=properties.s05;
    }
    
    for (unsigned int i=0;i<6;i++){
        twin_init[i]=0.0;
    }
    
    
    for (unsigned int i=0;i<18;i++){
        slip_init[i]=0.0;
    }
    
    for (unsigned int i=0;i<dim;i++){
        rot_init(i)=0.0;
        rotnew_init(i)=0.0;
    }
    
    //Resize the vectors of history variables
    Fp_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init));
    Fp_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init));
    twinfraction_iter.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,twin_init));
    slipfraction_iter.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,slip_init));
    twinfraction_conv.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,twin_init));
    slipfraction_conv.resize(num_local_cells,std::vector<vector<double> >(num_quad_points,slip_init));
    rot.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rot_init));
    rotnew.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rotnew_init));
    twin.resize(num_local_cells,std::vector<double>(num_quad_points,0.0));
    
    //load rot and rotnew
    for (unsigned int cell=0; cell<num_local_cells; cell++){
        for (unsigned int q=0; q<num_quad_points; q++){
            unsigned int materialID=quadratureOrientationsMap[cell][q];
            for (unsigned int i=0; i<dim; i++){
                rot[cell][q][i]=orientations.eulerAngles[materialID][i];
                rotnew[cell][q][i]=orientations.eulerAngles[materialID][i];
            }
        }
    }
    N_qpts=num_quad_points;
    initCalled=true;
}
