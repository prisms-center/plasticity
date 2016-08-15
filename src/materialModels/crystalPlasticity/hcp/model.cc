//constructor
template <int dim>
crystalPlasticity<dim>::crystalPlasticity() :
ellipticBVP<dim>(),
F(dim,dim),
F_tau(dim,dim),
FP_tau(dim,dim),
FE_tau(dim,dim),
T(dim,dim),
P(dim,dim)
{
    initCalled = false;
    
    //post processing
    ellipticBVP<dim>::numPostProcessedFields=4;
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Grain_ID");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("twin");
    
    
}
        
//implementation of the getElementalValues method
template <int dim>
void crystalPlasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
                                                unsigned int dofs_per_cell,
                                                unsigned int num_quad_points,
                                                FullMatrix<double>& elementalJacobian,
                                                Vector<double>&     elementalResidual)
{
    
    //Initialized history variables and pfunction variables if unititialized
    if(initCalled == false){
        init(num_quad_points);
    }
    
    unsigned int cellID = fe_values.get_cell()->user_index();
    std::vector<unsigned int> local_dof_indices(dofs_per_cell);
    Vector<double> Ulocal(dofs_per_cell);
    
    typename DoFHandler<dim>::active_cell_iterator cell(& this->triangulation,
                                                        fe_values.get_cell()->level(),
                                                        fe_values.get_cell()->index(),
                                                        & this->dofHandler);
    cell->set_user_index(fe_values.get_cell()->user_index());
    cell->get_dof_indices (local_dof_indices);
    for(unsigned int i=0; i<dofs_per_cell; i++){
        Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
    }
    
    //local data structures
    FullMatrix<double> K_local(dofs_per_cell,dofs_per_cell),CE_tau(dim,dim),E_tau(dim,dim),temp,temp2,temp3;
    Vector<double> Rlocal (dofs_per_cell);
    K_local = 0.0; Rlocal = 0.0;
    
    
    //loop over quadrature points
    for (unsigned int q=0; q<num_quad_points; ++q){
        //Get deformation gradient
        F=0.0;
        for (unsigned int d=0; d<dofs_per_cell; ++d){
            unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
            for (unsigned int j=0; j<dim; ++j){
                F[i][j]+=Ulocal(d)*fe_values.shape_grad(d, q)[j]; // u_{i,j}= U(d)*N(d)_{,j}, where d is the DOF correonding to the i'th dimension
            }
        }
        for (unsigned int i=0; i<dim; ++i){
            F[i][i]+=1;
        }
        
        
        //Update strain, stress, and tangent for current time step/quadrature point
        calculatePlasticity(cellID, q);
        
        //Fill local residual
        for (unsigned int d=0; d<dofs_per_cell; ++d) {
            unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
            for (unsigned int j = 0; j < dim; j++){
                Rlocal(d) -=  fe_values.shape_grad(d, q)[j]*P[i][j]*fe_values.JxW(q);
                
                
            }
            //if(q==7)
            // this->pcout<<Rlocal(d)<<'\n';
        }
        
        temp.reinit(dim,dim); temp=0.0;
        temp2.reinit(dim,dim); temp2=0.0;
        temp3.reinit(dim,dim); temp3=0.0;
        CE_tau=0.0;
        temp=F;
        F.Tmmult(CE_tau,temp);
        E_tau=CE_tau;
        temp=IdentityMatrix(dim);
        for(unsigned int i=0;i<dim;i++){
            for(unsigned int j=0;j<dim;j++){
                E_tau[i][j] = 0.5*(E_tau[i][j]-temp[i][j]);
                temp2[i][j]=E_tau[i][j]*fe_values.JxW(q);
                temp3[i][j]=T[i][j]*fe_values.JxW(q);
            }
        }
        //cout<<E_tau[0][0]<<"\t"<<T[0][0]<<"\t"<<fe_values.JxW(q)<<"\n";
        local_strain.add(1.0,temp2);
        local_stress.add(1.0,temp3);
        local_microvol=local_microvol+fe_values.JxW(q);
        
        double traceE, traceT,vonmises,eqvstrain;
        FullMatrix<double> deve(dim,dim),devt(dim,dim);
        
        
        traceE=E_tau.trace();
        traceT=T.trace();
        temp=IdentityMatrix(3);
        temp.equ(traceE/3,temp);
        
        deve=E_tau;
        deve.add(-1.0,temp);
        
        temp=IdentityMatrix(3);
        temp.equ(traceT/3,temp);
        
        devt=T;
        devt.add(-1.0,temp);
        
        vonmises= devt.frobenius_norm();
        vonmises=sqrt(3.0/2.0)*vonmises;
        eqvstrain=deve.frobenius_norm();
        eqvstrain=sqrt(2.0/3.0)*eqvstrain;
        
        
        //fill in post processing field values
        
        this->postprocessValues(cellID, q, 0, 0)=vonmises;
        this->postprocessValues(cellID, q, 1, 0)=eqvstrain;
        this->postprocessValues(cellID, q, 2, 0)=quadratureOrientationsMap[cellID][q];
        this->postprocessValues(cellID, q, 3, 0)=twin[cellID][q];
        
        
        
        
        std::cout.precision(3);
        
        //evaluate elemental stiffness matrix, K_{ij} = N_{i,k}*C_{mknl}*F_{im}*F{jn}*N_{j,l} + N_{i,k}*F_{kl}*N_{j,l}*del{ij} dV 
        for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
            unsigned int i = fe_values.get_fe().system_to_component_index(d1).first;
            for (unsigned int d2=0; d2<dofs_per_cell; ++d2) {
                unsigned int j = fe_values.get_fe().system_to_component_index(d2).first;
                for (unsigned int k = 0; k < dim; k++){
                    for (unsigned int l= 0; l< dim; l++){
                        K_local(d1,d2) +=  fe_values.shape_grad(d1, q)[k]*dP_dF[i][k][j][l]*fe_values.shape_grad(d2, q)[l]*fe_values.JxW(q);
                        
                    }
                }
                //if(q==7)
                //this->pcout<<K_local(d1,d2)<<'\t';
            }
            //if(q==7)
            //this->pcout<<'\n';
        }
    }
    elementalJacobian = K_local;
    elementalResidual = Rlocal;
}



 template <int dim>
 void crystalPlasticity<dim>::updateBeforeIteration()
 {
     local_strain=0.0;
     local_stress=0.0;
     local_microvol=0.0;

     //call base class project() function to project post processed fields
     //ellipticBVP<dim>::project();
 }

 template <int dim>
 void crystalPlasticity<dim>::updateBeforeIncrement()
 {
     microvol=0.0;
     //call base class project() function to project post processed fields
     //ellipticBVP<dim>::project();
 }



//implementation of the getElementalValues method
template <int dim>
void crystalPlasticity<dim>::updateAfterIncrement()
{
    reorient();
    
    twinfraction_conv=twinfraction_iter;
    slipfraction_conv=slipfraction_iter;
    
    
    //copy rotnew to output
    orientations.outputOrientations.clear();
    QGauss<dim>  quadrature(quadOrder);
    const unsigned int num_quad_points = quadrature.size();
    FEValues<dim> fe_values (this->FE, quadrature, update_quadrature_points | update_JxW_values);
    //loop over elements
    unsigned int cellID=0;
    typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){
            fe_values.reinit(cell);
            //loop over quadrature points
            for (unsigned int q=0; q<num_quad_points; ++q){
                std::vector<double> temp;
                temp.push_back(fe_values.get_quadrature_points()[q][0]);
                temp.push_back(fe_values.get_quadrature_points()[q][1]);
                temp.push_back(fe_values.get_quadrature_points()[q][2]);
                temp.push_back(rotnew[cellID][q][0]);
                temp.push_back(rotnew[cellID][q][1]);
                temp.push_back(rotnew[cellID][q][2]);
                temp.push_back(fe_values.JxW(q));
                temp.push_back(quadratureOrientationsMap[cellID][q]);
                temp.push_back(slipfraction_conv[cellID][q][0]);
                temp.push_back(slipfraction_conv[cellID][q][1]);
                temp.push_back(slipfraction_conv[cellID][q][2]);
                temp.push_back(slipfraction_conv[cellID][q][3]);
                temp.push_back(slipfraction_conv[cellID][q][4]);
                temp.push_back(slipfraction_conv[cellID][q][5]);
                
                temp.push_back(twinfraction_conv[cellID][q][0]);
                temp.push_back(twinfraction_conv[cellID][q][1]);
                temp.push_back(twinfraction_conv[cellID][q][2]);
                temp.push_back(twinfraction_conv[cellID][q][3]);
                temp.push_back(twinfraction_conv[cellID][q][4]);
                temp.push_back(twinfraction_conv[cellID][q][5]);
                temp.push_back(twin[cellID][q]);
                orientations.addToOutputOrientations(temp);
                local_F_e=local_F_e+twin[cellID][q]*fe_values.JxW(q);
                for(unsigned int i=0;i<6;i++){
                    local_F_r=local_F_r+twinfraction_conv[cellID][q][i]*fe_values.JxW(q);
                }
                
            }
            cellID++;
        }
    }
    orientations.writeOutputOrientations();
    
    //Update the history variables when convergence is reached for the current increment
    Fe_conv=Fe_iter;
    Fp_conv=Fp_iter;
    s_alpha_conv=s_alpha_iter;
    
    //double temp4,temp5;
    //temp4=Lambda[0][0];
    //temp5=Utilities::MPI::sum(temp4,this->mpi_communicator);
    //cout << temp5<<"\n";
    microvol=Utilities::MPI::sum(local_microvol,this->mpi_communicator);
    
    for(unsigned int i=0;i<dim;i++){
        for(unsigned int j=0;j<dim;j++){
            global_strain[i][j]=Utilities::MPI::sum(local_strain[i][j]/microvol,this->mpi_communicator);
            global_stress[i][j]=Utilities::MPI::sum(local_stress[i][j]/microvol,this->mpi_communicator);
        }
        
    }
    F_e=Utilities::MPI::sum(local_F_e/microvol,this->mpi_communicator);
    F_r=Utilities::MPI::sum(local_F_r/microvol,this->mpi_communicator);
    
    
    
    
    ofstream outputFile;
    
    if(this->currentIncrement==0){
        outputFile.open("stressstrain.txt");
        outputFile.close();
    }
    outputFile.open("stressstrain.txt",ios::app);
    if(Utilities::MPI::this_mpi_process(this->mpi_communicator)==0){
        outputFile << global_strain[0][0]<<'\t'<<global_strain[1][1]<<'\t'<<global_strain[2][2]<<'\t'<<global_strain[1][2]<<'\t'<<global_strain[0][2]<<'\t'<<global_strain[0][1]<<'\t'<<global_stress[0][0]<<'\t'<<global_stress[1][1]<<'\t'<<global_stress[2][2]<<'\t'<<global_stress[1][2]<<'\t'<<global_stress[0][2]<<'\t'<<global_stress[0][1]<<'\n';
    }
    outputFile.close();
    global_strain=0.0;
    global_stress=0.0;
    
    
    cellID=0;
    cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){
            fe_values.reinit(cell);
            //loop over quadrature points
            for (unsigned int q=0; q<num_quad_points; ++q){
                std::vector<double> local_twin;
                local_twin.resize(6,0.0);
                local_twin=twinfraction_conv[cellID][q];
                std::vector<double>::iterator result;
                result = std::max_element(local_twin.begin(), local_twin.end());
                double twin_pos, twin_max;
                twin_pos= std::distance(local_twin.begin(), result);
                twin_max=local_twin[twin_pos];
                
                if(twin_max>=(twinThresholdFraction+twinSaturationFactor*F_e/F_r)){
                    
                    FullMatrix<double> FE_t(dim,dim), FP_t(dim,dim),Twin_T(dim,dim),temp(dim,dim);
                    
                    FE_t=Fe_conv[cellID][q];
                    FP_t=Fp_conv[cellID][q];
                    
                    
                    Twin_image(twin_pos,cellID,q);
                    double s_alpha_twin=s_alpha_conv[cellID][q][numSlipSystems+twin_pos];
                    for(unsigned int i=0;i<numTwinSystems;i++){
                        twinfraction_conv[cellID][q][i]=0;
                        s_alpha_conv[cellID][q][numSlipSystems+twin_pos]=s_alpha_twin;
                        
                    }
                    
                    Vector<double> n(dim);
                    n(0)=n_alpha[numSlipSystems+twin_pos][0];
                    n(1)=n_alpha[numSlipSystems+twin_pos][1];
                    n(2)=n_alpha[numSlipSystems+twin_pos][2];
                    
                    for(unsigned int i=0;i<dim;i++){
                        for(unsigned int j=0;j<dim;j++){
                            Twin_T(i,j)=2.0*n(i)*n(j)-1.0*(i==j);
                            
                        }
                        
                    }
                    
                    //this->pcout<<Twin_T(0,0)<<'\t'<<Twin_T(1,1)<<'\t'<<Twin_T(2,2)<<'\t'<<Twin_T(0,1)<<'\n';
                    
                    FE_t.mmult(temp,Twin_T);
                    Twin_T.mmult(FE_t,temp);
                    
                    FP_t.mmult(temp,Twin_T);
                    Twin_T.mmult(FP_t,temp);
                    
                    Fe_conv[cellID][q]=FE_t;
                    Fp_conv[cellID][q]=FP_t;
                    
                    
                    twin[cellID][q]=1.0;
                    
                }
                
                
            }
            cellID++;
        }
    }
    
    
    
    //call base class project() function to project post processed fields
    ellipticBVP<dim>::project();
}



template <int dim>
void crystalPlasticity<dim>::Twin_image(double twin_pos,unsigned int cellID,
                                        unsigned int quadPtID)
{
    Vector<double> quat1(4),rod(3),quat2(4),quatprod(4);
    rod(0) = rot[cellID][quadPtID][0];rod(1) = rot[cellID][quadPtID][1];rod(2) = rot[cellID][quadPtID][2];
    
    rod2quat(quat2,rod);
    
    quat1(0) = 0;quat1(1) = n_alpha[18+twin_pos][0];
    quat1(2) = n_alpha[18+twin_pos][1];quat1(3) = n_alpha[18+twin_pos][2];
    
    
    quatproduct(quatprod,quat2,quat1);
    
    quat2rod(quatprod,rod);
    
    //this->pcout<<rod(0)<<'\t'<<rod(1)<<'\t'<<rod(2)<<'\n';
    
    rot[cellID][quadPtID][0]=rod(0);rot[cellID][quadPtID][1]=rod(1);rot[cellID][quadPtID][2]=rod(2);
    rotnew[cellID][quadPtID][0]=rod(0);rotnew[cellID][quadPtID][1]=rod(1);rotnew[cellID][quadPtID][2]=rod(2);
    
    
}


//------------------------------------------------------------------------
template <int dim>
void crystalPlasticity<dim>::rod2quat(Vector<double> &quat,Vector<double> &rod)
{
    double dotrod = rod(0)*rod(0) + rod(1)*rod(1) + rod(2)*rod(2);
    double cphiby2   = cos(atan(sqrt(dotrod)));
    quat(0) = cphiby2;
    quat(1) = cphiby2* rod(0);
    quat(2) = cphiby2* rod(1);
    quat(3) = cphiby2* rod(2);
}
//------------------------------------------------------------------------
template <int dim>
void crystalPlasticity<dim>::quatproduct(Vector<double> &quatp,Vector<double> &quat2,Vector<double> &quat1)
//R(qp) = R(q2)R(q1)
{
    double a = quat2(0);
    double b = quat1(0);
    double dot1 = quat1(1)*quat2(1) + quat1(2)*quat2(2) + quat1(3)*quat2(3);
    quatp(0) = (a*b) - dot1;
    quatp(1) = a*quat1(1) + b*quat2(1)+ quat2(2)*quat1(3) - quat1(2)*quat2(3);
    quatp(2) = a*quat1(2) + b*quat2(2)- quat2(1)*quat1(3) + quat1(1)*quat2(3);
    quatp(3) = a*quat1(3) + b*quat2(3)+ quat2(1)*quat1(2) - quat1(1)*quat2(2);
    //if (quatp(0) < 0)
    //quatp.mult(-1);
}
//------------------------------------------------------------------------
template <int dim>
void crystalPlasticity<dim>::quat2rod(Vector<double> &quat,Vector<double> &rod)
{
    double invquat1 = 1/quat(0);
    
    for (int i = 0;i <= 2;i++)
        rod(i) = quat(i+1)*invquat1;
    
}