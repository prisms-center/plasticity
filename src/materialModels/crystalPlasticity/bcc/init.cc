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
    
    //Elastic Stiffness Matrix Dmat
    Dmat.reinit(6,6); Dmat=0.0;
    Dmat[0][0]=properties.C11; Dmat[0][1]=properties.C12; Dmat[0][2]=properties.C12; Dmat[1][0]=properties.C12; Dmat[1][1]=properties.C11; Dmat[1][2]=properties.C12; Dmat[2][0]=properties.C12; Dmat[2][1]=properties.C12; Dmat[2][2]=properties.C11;
    Dmat[3][3]=2*properties.C44; Dmat[4][4]=2*properties.C44; Dmat[5][5]=2*properties.C44;
    
    Vector<double> s0_init (n_slip_systems),rot_init(dim),rotnew_init(dim);
    
    for (unsigned int i=0;i<n_slip_systems;i++){
        s0_init(i)=properties.s0;
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
    rot.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rot_init));
    rotnew.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rotnew_init));
    
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