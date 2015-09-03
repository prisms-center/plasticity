template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
                                                 unsigned int quadPtID)
{
    
    F_tau=F; // Deformation Gradient
    FullMatrix<double> FE_t(dim,dim),FP_t(dim,dim);  //Elastic and Plastic deformation gradient
    Vector<double> s_alpha_t(n_slip_systems); // Slip resistance
    Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)
    
    int old_precision = std::cout.precision();
    
    
    std::cout.precision(16);
    
    FE_t=Fe_conv[cellID][quadPtID];
    FP_t=Fp_conv[cellID][quadPtID];
    s_alpha_t=s_alpha_conv[cellID][quadPtID];
    rot1=rot[cellID][quadPtID];
    
    // Rotation matrix of the crystal orientation
    FullMatrix<double> rotmat(dim,dim);
    rotmat=0.0;
    odfpoint(rotmat,rot1);
    
    FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim); // Temporary matrices
    
    //convert to crystal coordinates F_tau=R'*F_tau*R
    temp=0.0;
    rotmat.Tmmult(temp,F_tau);
    temp.mmult(F_tau,rotmat);
    
    //FE_tau_trial=F_tau*inv(FP_t)
    FullMatrix<double> Fpn_inv(dim,dim),FE_tau_trial(dim,dim),F_trial(dim,dim),CE_tau_trial(dim,dim),Ee_tau_trial(dim,dim);
    Fpn_inv=0.0; Fpn_inv.invert(FP_t);
    FE_tau_trial=0.0;
    F_trial=0.0;
    F_tau.mmult(FE_tau_trial,Fpn_inv);F_trial = FE_tau_trial;
    
    //% % % % % STEP 1 % % % % %
    //Calculate Trial Elastic Strain Ee_tau_trial
    //Ee_tau_trial=0.5(FE_tau_trial'*FE_tau_trial-I)
    temp.reinit(dim,dim); temp=0.0;
    temp=FE_tau_trial;
    FE_tau_trial.Tmmult(CE_tau_trial,temp);
    Ee_tau_trial=CE_tau_trial;
    temp=IdentityMatrix(dim);
    for(unsigned int i=0;i<dim;i++){
        for(unsigned int j=0;j<dim;j++){
            Ee_tau_trial[i][j] = 0.5*(Ee_tau_trial[i][j]-temp[i][j]);
        }
    }
    
    // Calculation of Schmid Tensors  and B= symm(FE_tau_trial'*FE_tau_trial*S_alpha)
    FullMatrix<double> SCHMID_TENSOR1(n_slip_systems*dim,dim),B(n_slip_systems*dim,dim);
    Vector<double> m1(dim),n1(dim);
    
    for(unsigned int i=0;i<n_slip_systems;i++){
        for (unsigned int j=0;j<dim;j++){
            m1(j)=m_alpha[i][j];
            n1(j)=n_alpha[i][j];
        }
        
        for (unsigned int j=0;j<dim;j++){
            for (unsigned int k=0;k<dim;k++){
                temp[j][k]=m1(j)*n1(k);
                SCHMID_TENSOR1[dim*i+j][k]=m1(j)*n1(k);
            }
        }
        
        CE_tau_trial.mmult(temp2,temp);
        temp2.symmetrize();
        for (unsigned int j=0;j<dim;j++){
            for (unsigned int k=0;k<dim;k++){
                B[dim*i+j][k]=2*temp2[j][k];
            }
        }
    }
    
    //% % % % % STEP 2 % % % % %
    // Calculate the trial stress T_star_tau_trial
    Vector<double> tempv1(6),tempv2(6);
    FullMatrix<double> T_star_tau_trial(dim,dim);
    tempv1=0.0;
    Dmat.vmult(tempv1, vecform(Ee_tau_trial));
    matform(T_star_tau_trial,tempv1);
    
    //% % % % % STEP 3 % % % % %
    // Calculate the trial resolved shear stress resolved_shear_tau_trial for each slip system
    int n_PA=0;	// Number of active slip systems
    Vector<double> PA, PA_temp(1);
    Vector<double> resolved_shear_tau_trial(n_slip_systems),b(n_slip_systems),resolved_shear_tau(n_slip_systems);
    resolved_shear_tau_trial=0.0;
    
    
    for(unsigned int i=0;i<n_slip_systems;i++){
        
        for (unsigned int j=0;j<dim;j++){
            for (unsigned int k=0;k<dim;k++){
                resolved_shear_tau_trial(i)+=T_star_tau_trial[j][k]*SCHMID_TENSOR1[dim*i+j][k];
            }
        }
        //% % % % % STEP 4 % % % % %
        //Determine the set set of the n potentially active slip systems
        b(i)=fabs(resolved_shear_tau_trial(i))-s_alpha_t(i);
        if( b(i)>=0){
            if(n_PA==0){
                n_PA=n_PA+1;
                PA.reinit(n_PA);
                PA(0)=i;
            }
            else{
                PA_temp=PA;
                n_PA=n_PA+1;
                PA.reinit(n_PA);
                for (unsigned int j=0;j<(n_PA-1);j++){
                    PA(j)=PA_temp(j);
                }
                PA(n_PA-1)=i;
                PA_temp.reinit(n_PA);     //%%%%% Potentially active slip systems
            }
        }
        resolved_shear_tau(i)=fabs(resolved_shear_tau_trial(i));
    }
    
    
    //% % % % % STEP 5 % % % % %
    //Calculate the shear increments from the consistency condition
    Vector<double> s_beta(n_slip_systems),h_beta(n_slip_systems),h0(n_slip_systems),a_pow(n_slip_systems),s_s(n_slip_systems);
    FullMatrix<double> h_alpha_beta_t(n_slip_systems,n_slip_systems),A(n_slip_systems,n_slip_systems);
    s_beta=s_alpha_t;
    
    
    for (unsigned int i=0;i<n_slip_systems;i++){
        if(i<3){
            h0(i)=properties.h01;
            a_pow(i)=properties.a1;
            s_s(i)=properties.s_s1;
        }
        else if(i<6){
            h0(i)=properties.h02;
            a_pow(i)=properties.a2;
            s_s(i)=properties.s_s2;
        }
        else if(i<12){
            h0(i)=properties.h03;
            a_pow(i)=properties.a3;
            s_s(i)=properties.s_s3;
        }
        else if(i<18){
            h0(i)=properties.h04;
            a_pow(i)=properties.a4;
            s_s(i)=properties.s_s4;
        }
        else if(i<24){
            h0(i)=properties.h05;
            a_pow(i)=properties.a5;
            s_s(i)=properties.s_s5;
        }
        
    }
    
    
    // Single slip hardening rate
    for(unsigned int i=0;i<n_slip_systems;i++){
        if(s_beta(i)>s_s(i))
            s_beta(i)=0.98*s_s(i);
        h_beta(i)=h0(i)*pow((1-s_beta(i)/s_s(i)),a_pow(i));
    }
    
    
    for(unsigned int i=0;i<n_slip_systems;i++){
        for(unsigned int j=0;j<n_slip_systems;j++){
            h_alpha_beta_t[i][j] = q[i][j]*h_beta(j);
            A[i][j]=h_alpha_beta_t[i][j];
        }
    }
    
    // Calculate the Stiffness Matrix A
    for(unsigned int i=0;i<n_slip_systems;i++){
        for(unsigned int j=0;j<n_slip_systems;j++){
            temp1.reinit(dim,dim); temp1=0.0;
            
            for(unsigned int k=0;k<dim;k++){
                for(unsigned int l=0;l<dim;l++){
                    temp[k][l]=SCHMID_TENSOR1(dim*j+k,l);
                }
            }
            temp2.reinit(dim,dim); CE_tau_trial.mmult(temp2,temp);
            temp2.symmetrize();
            tempv1=0.0; Dmat.vmult(tempv1, vecform(temp2));
            temp3=0.0; matform(temp3,tempv1);
            
            for(unsigned int k=0;k<dim;k++){
                for(unsigned int l=0;l<dim;l++){
                    if((resolved_shear_tau_trial(i)<0.0)^(resolved_shear_tau_trial(j)<0.0))
                        A[i][j]-=SCHMID_TENSOR1(dim*i+k,l)*temp3[k][l];
                    else
                        A[i][j]+=SCHMID_TENSOR1(dim*i+k,l)*temp3[k][l];
                    
                }
            }
        }
    }
    
    
    // Caluclate the trial Cauchy Stress T_tau and trial PK1 Stress P_tau
    FullMatrix<double> T_tau(dim,dim),P_tau(dim,dim);
    Vector<double> s_alpha_tau;
    FP_tau=FP_t;
    FE_tau.reinit(dim,dim);
    F_tau.mmult(FE_tau,Fpn_inv);
    
    double det_FE_tau,det_F_tau, det_FP_tau;
    det_FE_tau=FE_tau.determinant();
    temp.reinit(dim,dim); FE_tau.mmult(temp,T_star_tau_trial);
    temp.equ(1.0/det_FE_tau,temp); temp.mTmult(T_tau,FE_tau);
    det_F_tau=F_tau.determinant();
    temp=0.0;
    temp.invert(F_tau); T_tau.mTmult(P_tau,temp);
    P_tau.equ(det_F_tau,P_tau);
    s_alpha_tau=s_alpha_t;
    
    double gamma = 0;
    int iter=0,iter1=0,iter2=0,iter3=0;
    Vector<double> active;
    Vector<double> x_beta;
    FullMatrix<double> A_PA;
    
    
    // Determination of active slip systems and shear increments
    if (n_PA > 0){
        
        Vector<double> inactive(n_slip_systems-n_PA); // Set of inactive slip systems
        
        
        inactive_slip_removal(inactive,active,x_beta,n_PA, PA, b, A, A_PA);
        
        // % % % % % STEP 6 % % % % %
        temp.reinit(dim,dim);
        for (unsigned int i=0;i<n_slip_systems;i++){
            for (unsigned int j=0;j<dim;j++){
                for (unsigned int k=0;k<dim;k++){
                    temp[j][k]=SCHMID_TENSOR1[dim*i+j][k];
                }
            }
            
            temp.mmult(temp2,FP_t);
            for (unsigned int j=0;j<dim;j++){
                for (unsigned int k=0;k<dim;k++){
                    if(resolved_shear_tau_trial(i)>0)
                        FP_tau[j][k]=FP_tau[j][k]+x_beta(i)*temp2[j][k];
                    else
                        FP_tau[j][k]=FP_tau[j][k]-x_beta(i)*temp2[j][k];
                }
            }
            
        }
        
        //% % % % % STEP 7 % % % % %
        //%     FP_tau = FP_tau/(det(FP_tau)^(1/3));
        
        det_FP_tau=FP_tau.determinant();
        FP_tau.equ(1.0/pow(det_FP_tau,1.0/3.0),FP_tau);
        
        
        // % % % % % STEP 8 % % % % %
        temp.invert(FP_tau);
        F_tau.mmult(FE_tau,temp);
        FullMatrix<double> Ce_tau(dim,dim),T_star_tau(dim,dim);
        FE_tau.Tmmult(Ce_tau,FE_tau);
        T_star_tau=0.0;
        
        for(unsigned int i=0;i<n_slip_systems;i++){
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    temp[j][k]=SCHMID_TENSOR1(dim*i+j,k);
                }
            }
            
            CE_tau_trial.mmult(temp2,temp);
            temp2.symmetrize();
            tempv1=0.0; Dmat.vmult(tempv1, vecform(temp2));
            matform(temp3,tempv1);
            
            
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    if(resolved_shear_tau_trial(i)>0.0)
                        T_star_tau[j][k]=  T_star_tau[j][k] - x_beta(i)*temp3[j][k];
                    else
                        T_star_tau[j][k]=  T_star_tau[j][k] + x_beta(i)*temp3[j][k];
                }
            }
        }
        
        for (unsigned int i=0;i<6;i++){
            twinfraction_iter[cellID][quadPtID][i]=twinfraction_conv[cellID][quadPtID][i]+x_beta[i+18]/0.129;
        }
        
        for (unsigned int i=0;i<18;i++){
            slipfraction_iter[cellID][quadPtID][i]=slipfraction_conv[cellID][quadPtID][i]+x_beta[i]/0.129;
        }
        
        
        
        T_star_tau.add(1.0,T_star_tau_trial);
        
        // % % % % % STEP 9 % % % % %
        
        temp.reinit(dim,dim);
        det_FE_tau=FE_tau.determinant();
        FE_tau.mmult(temp,T_star_tau); temp.equ(1.0/det_FE_tau,temp);
        temp.mTmult(T_tau,FE_tau);
        
        det_F_tau=F_tau.determinant();
        temp=0.0;
        temp.invert(F_tau); T_tau.mTmult(P_tau,temp);
        P_tau.equ(det_F_tau,P_tau);
        
        double h1=0;
        for(unsigned int i=0;i<n_slip_systems;i++){
            h1=0;
            for(unsigned int j=0;j<n_slip_systems;j++){
                h1=h1+h_alpha_beta_t(i,j)*x_beta(j);
            }
            s_alpha_tau(i)=s_alpha_t(i)+h1;
        }
        
        //% see whether shear resistance is passed
        for(unsigned int i=0;i<n_slip_systems;i++){
            for(unsigned int j=0;j<n_slip_systems;j++){
                temp.reinit(dim,dim);
                for(unsigned int k=0;k<dim;k++){
                    for(unsigned int l=0;l<dim;l++){
                        temp[k][l]=SCHMID_TENSOR1(dim*j+k,l);
                    }
                }
                temp2.reinit(dim,dim); CE_tau_trial.mmult(temp2,temp);
                temp2.symmetrize();
                tempv1.reinit(2*dim); tempv1=0.0;
                Dmat.vmult(tempv1, vecform(temp2));
                temp3.reinit(dim,dim); 	matform(temp3,tempv1);
                temp.reinit(dim,dim);
                
                for(unsigned int k=0;k<dim;k++){
                    
                    for(unsigned int l=0;l<dim;l++){
                        
                        temp[k][l]=SCHMID_TENSOR1(dim*i+k,l);
                    }
                }
                
                //Update the resolved shear stress
                for(unsigned int k=0;k<dim;k++){
                    for(unsigned int l=0;l<dim;l++){
                        if((resolved_shear_tau_trial(i)<0.0)^(resolved_shear_tau_trial(j)<0.0))
                            resolved_shear_tau(i) = resolved_shear_tau(i)+temp[k][l]*temp3[k][l]*x_beta(j);
                        else
                            resolved_shear_tau(i) = resolved_shear_tau(i)-temp[k][l]*temp3[k][l]*x_beta(j);
                    }
                }
                
            }
            
            gamma=gamma+resolved_shear_tau(i)*x_beta(i);
        }
        
        
    }
    
    FullMatrix<double> PK1_Stiff(dim*dim,dim*dim);
    
    tangent_modulus(F_trial, Fpn_inv, SCHMID_TENSOR1,A,A_PA,B,T_tau, PK1_Stiff, active, resolved_shear_tau_trial, x_beta, PA, n_PA,det_F_tau,det_FE_tau );
    
    
    
    temp.reinit(dim,dim);
    temp=0.0;
    T_tau.mTmult(temp,rotmat);
    T_tau=0.0;
    rotmat.mmult(T_tau,temp);
    temp=0.0;
    temp.reinit(dim,dim); P_tau.mTmult(temp,rotmat);
    P_tau=0.0;
    rotmat.mmult(P_tau,temp);
    
    dP_dF=0.0;
    FullMatrix<double> L(dim,dim),mn(dim,dim);
    L=0.0;
    temp1.reinit(dim,dim); temp1=IdentityMatrix(dim);
    rotmat.Tmmult(L,temp1);
    
    // Transform the tangent modulus back to crystal frame
    
    
    for(unsigned int m=0;m<dim;m++){
        for(unsigned int n=0;n<dim;n++){
            for(unsigned int o=0;o<dim;o++){
                for(unsigned int p=0;p<dim;p++){
                    for(unsigned int i=0;i<dim;i++){
                        for(unsigned int j=0;j<dim;j++){
                            for(unsigned int k=0;k<dim;k++){
                                for(unsigned int l=0;l<dim;l++){
                                    dP_dF[m][n][o][p]=dP_dF[m][n][o][p]+PK1_Stiff(dim*i+j,dim*k+l)*L(i,m)*L(j,n)*L(k,o)*L(l,p);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    
    
    
    P.reinit(dim,dim);
    P=P_tau;
    T=T_tau;
    
    
    sres_tau.reinit(n_slip_systems);
    sres_tau = s_alpha_tau;
    
    // Update the history variables
    Fe_iter[cellID][quadPtID]=FE_tau;
    Fp_iter[cellID][quadPtID]=FP_tau;
    s_alpha_iter[cellID][quadPtID]=sres_tau;
    
}

