#include "../../../include/crystalPlasticity.h"


template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
                                                 unsigned int quadPtID)
{

    F_tau=F; // Deformation Gradient
    FullMatrix<double> FE_t(dim,dim),FP_t(dim,dim);  //Elastic and Plastic deformation gradient
    Vector<double> s_alpha_t(n_slip_systems); // Slip resistance
    Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)

    // Tolerance

    double tol1=this->userInputs.modelStressTolerance;



    std::cout.precision(16);

    FE_t=Fe_conv[cellID][quadPtID];
    FP_t=Fp_conv[cellID][quadPtID];
    s_alpha_t=s_alpha_conv[cellID][quadPtID];
    rot1=rot[cellID][quadPtID];



    // Rotation matrix of the crystal orientation
    FullMatrix<double> rotmat(dim,dim);
    rotmat=0.0;
    odfpoint(rotmat,rot1);


    FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim),temp4(dim,dim),temp5(dim,dim),temp6(dim,dim); // Temporary matrices
    FullMatrix<double> T_tau(dim,dim),P_tau(dim,dim);
    FullMatrix<double> Fpn_inv(dim,dim),FE_tau_trial(dim,dim),F_trial(dim,dim),CE_tau_trial(dim,dim),FP_t2(dim,dim),Ee_tau_trial(dim,dim);

    // Calculation of Schmid Tensors  and B= symm(FE_tau_trial'*FE_tau_trial*S_alpha)
    FullMatrix<double> SCHMID_TENSOR1(n_slip_systems*dim,dim),B(n_slip_systems*dim,dim);
    Vector<double> m1(dim),n1(dim);


    // Elastic Modulus

    FullMatrix<double> Dmat2(2*dim,2*dim),TM(dim*dim,dim*dim), ElasticityTensor(2*dim,2*dim), elasticStiffnessMatrix(2 * dim, 2 * dim);
    Vector<double> vec1(2*dim),vec2(dim*dim);
	for (unsigned int i = 0;i<6;i++) {
		for (unsigned int j = 0;j<6;j++) {
			elasticStiffnessMatrix[i][j] = this->userInputs.elasticStiffness[i][j];
		}
	}

	elasticmoduli(Dmat2, rotmat, elasticStiffnessMatrix);

	//Elastic Stiffness Matrix Dmat
		Dmat.reinit(6, 6); Dmat = 0.0;

	for (unsigned int i = 0;i<6;i++) {
		for (unsigned int j = 0;j<6;j++) {
			Dmat[i][j] = Dmat2[i][j];
		}
	}


	for (unsigned int i = 0;i<6;i++) {
		for (unsigned int j = 3;j<6;j++) {
			Dmat[i][j] = 2 * Dmat[i][j];
		}
	}

    vec1(0)=0;vec1(1)=5;vec1(2)=4;vec1(3)=1;vec1(4)=3;vec1(5)=2;
    vec2(0)=0;vec2(1)=5;vec2(2)=4;vec2(3)=5;vec2(4)=1;vec2(5)=3;vec2(6)=4;vec2(7)=3;vec2(8)=2;



    for(unsigned int i=0;i<9;i++){
        for(unsigned int j=0;j<9;j++){
            TM[i][j]=Dmat2(vec2(i),vec2(j));
        }
    }

	for (unsigned int i = 0;i<6;i++) {
		for (unsigned int j = 0;j<6;j++) {
			ElasticityTensor[i][j] = Dmat(vec1(i), vec1(j));
		}
	}


    Vector<double> s_alpha_tau;
    Vector<double> s_beta(n_slip_systems),h_beta(n_slip_systems),delh_beta_dels(n_slip_systems),h0(n_slip_systems),a_pow(n_slip_systems),s_s(n_slip_systems);
    FullMatrix<double> h_alpha_beta_t(n_slip_systems,n_slip_systems),A(n_slip_systems,n_slip_systems);
    FullMatrix<double> del_FP(dim,dim);
    FullMatrix<double> A_PA;
    Vector<double> active, active_First;
    Vector<double> PA, PA_temp(1);
    Vector<double> resolved_shear_tau_trial(n_slip_systems),b(n_slip_systems),resolved_shear_tau(n_slip_systems);
    Vector<double> x_beta_old(n_slip_systems);

    Vector<double> x_beta(n_slip_systems);

    FullMatrix<double> PK1_Stiff(dim*dim,dim*dim);
    FullMatrix<double> CE_tau(dim,dim),T_star_tau(dim,dim);
    FullMatrix<double> T_star_tau_trial(dim,dim);



	double det_FE_tau, det_F_tau, det_FP_tau;
    unsigned int n_PA=0;	// Number of active slip systems


	    FP_tau = FP_t;
	    Fpn_inv = 0.0; Fpn_inv.invert(FP_t);
	    s_alpha_tau = s_alpha_t;
        FE_tau_trial=0.0;
        F_trial=0.0;
        F_tau.mmult(FE_tau_trial,Fpn_inv);F_trial = FE_tau_trial;

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

		for (unsigned int i = 0;i<n_slip_systems;i++) {
			for (unsigned int j = 0;j<dim;j++) {
				m1(j) = m_alpha[i][j];
				n1(j) = n_alpha[i][j];
			}

			temp = 0.0;
			temp2 = 0.0;
			for (unsigned int j = 0;j<dim;j++) {
				for (unsigned int k = 0;k<dim;k++) {
					temp[j][k] = m1(j)*n1(k);
				}
			}
			//convert the Shmitch tensor to Sample coordinates
			rotmat.mmult(temp2, temp);
			temp2.mTmult(temp, rotmat);
			for (unsigned int j = 0;j<dim;j++) {
				for (unsigned int k = 0;k<dim;k++) {

					SCHMID_TENSOR1[dim*i + j][k] = temp[j][k];
				}
			}
			CE_tau_trial.mmult(temp2, temp);
			temp2.symmetrize();
			for (unsigned int j = 0;j<dim;j++) {
				for (unsigned int k = 0;k<dim;k++) {
					B[dim*i + j][k] = 2 * temp2[j][k];
				}
			}
		}


        //% % % % % STEP 2 % % % % %
        // Calculate the trial stress T_star_tau_trial
        Vector<double> tempv1(6),tempv2(6);
        tempv1=0.0;
        Dmat.vmult(tempv1, vecform(Ee_tau_trial));
        matform(T_star_tau_trial,tempv1);
        T_star_tau.equ(1.0,T_star_tau_trial);




        //% % % % % STEP 3 % % % % %
        // Calculate the trial resolved shear stress resolved_shear_tau_trial for each slip system

        resolved_shear_tau_trial=0.0;

		for (unsigned int i = 0;i < n_slip_systems;i++) {

			for (unsigned int j = 0;j < dim;j++) {
				for (unsigned int k = 0;k < dim;k++) {
						resolved_shear_tau_trial(i) += T_star_tau_trial[j][k] * SCHMID_TENSOR1[dim*i + j][k];
				}
			}

			if (i > this->userInputs.numSlipSystems - 1) {
				if (resolved_shear_tau_trial(i) < 0)
					resolved_shear_tau_trial(i) = 0;
			}
		}


			x_beta_old = 0.0;
			unsigned int iter1 = 1;
			unsigned int flag2 = 0;
			resolved_shear_tau = resolved_shear_tau_trial;

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////Start Nonlinear iteration for Slip increments////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


			while (iter1) {
				x_beta = 0.0;

				if (iter1>this->userInputs.modelMaxSlipSearchIterations) {
					flag2 = 1;
					break;
				}


            //% % % % % STEP 4 % % % % %
            n_PA = 0;	// Number of active slip systems
			b.reinit(n_slip_systems);
            //Determine the set set of the n potentially active slip systems
			for (unsigned int i = 0;i<n_slip_systems;i++) {
				b(i)=fabs(resolved_shear_tau(i))-s_alpha_tau(i);
				if (b(i) >= tol1) {
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
        }



        if(n_PA==0)
            break;

		Vector<double> b_PA(n_PA);
		b_PA.reinit(n_PA);

		for (unsigned int i = 0;i<n_PA;i++) {
			b_PA(i) = b(PA(i));
		}

		if ((b_PA.linfty_norm())<tol1)
			break;

		if (iter1 > 1) {
			temp.reinit(dim, dim);
		}

	//	resolved_shear_tau_trial = resolved_shear_tau;
		FP_t2 = FP_tau;
		Fpn_inv = 0.0; Fpn_inv.invert(FP_t2);
		FE_tau = 0.0;
		F_tau.mmult(FE_tau, Fpn_inv);
		temp.reinit(dim, dim); temp = 0.0;
		temp = FE_tau;
		FE_tau.Tmmult(CE_tau, temp);

        //% % % % % STEP 5 % % % % %
        //Calculate the shear increments from the consistency condition
        s_beta=s_alpha_tau;

        // Single slip hardening rate
        for(unsigned int i=0;i<this->userInputs.numSlipSystems;i++){
            h_beta(i)=this->userInputs.initialHardeningModulus[i]*pow((1-s_beta(i)/this->userInputs.saturationStress[i]),this->userInputs.powerLawExponent[i]);
        }


        for(unsigned int i=0;i<this->userInputs.numTwinSystems;i++){
            h_beta(this->userInputs.numSlipSystems+i)=this->userInputs.initialHardeningModulusTwin[i]*pow((1-s_beta(this->userInputs.numSlipSystems+i)/this->userInputs.saturationStressTwin[i]),this->userInputs.powerLawExponentTwin[i]);
        }


        for(unsigned int i=0;i<n_slip_systems;i++){
            for(unsigned int j=0;j<n_slip_systems;j++){
                h_alpha_beta_t[i][j] = q[i][j]*h_beta(j);
                A[i][j]=h_alpha_beta_t[i][j];
            }
        }

        for(unsigned int i=0;i<n_slip_systems;i++){
           temp1.reinit(dim,dim);temp1=0.0;
                for(unsigned int k=0;k<dim;k++){
                    for(unsigned int l=0;l<dim;l++){
                        temp1[k][l]=SCHMID_TENSOR1(dim*i+k,l);
                    }
                }
             for(unsigned int j=0;j<n_slip_systems;j++){
                temp.reinit(dim,dim);temp=0.0;
                for(unsigned int k=0;k<dim;k++){
                    for(unsigned int l=0;l<dim;l++){
                        temp[k][l]=SCHMID_TENSOR1(dim*j+k,l);
                    }
                }
                temp2.reinit(dim,dim); CE_tau.mmult(temp2,temp);
                temp2.symmetrize();
                tempv1=0.0; Dmat.vmult(tempv1, vecform(temp2));
                temp3=0.0; matform(temp3,tempv1);
                //temp3 is symm in Matlab line 94 of constitutive.m


                for(unsigned int k=0;k<dim;k++){
                    for(unsigned int l=0;l<dim;l++){
                        if((resolved_shear_tau(i)*resolved_shear_tau(j))<0.0)
                            A[i][j]-=temp1[k][l]*temp3[k][l];
                        else
                            A[i][j]+=temp1[k][l]*temp3[k][l];

                    }
                }
            }
        }




            //Modified slip system search for adding corrective term
            inactive_slip_removal(active,x_beta_old,x_beta,n_PA,PA,b,A,A_PA);

            temp.reinit(dim,dim);
            del_FP.reinit(dim,dim);
            del_FP=0.0;
            for (unsigned int i=0;i<n_slip_systems;i++){
                for (unsigned int j=0;j<dim;j++){
                    for (unsigned int k=0;k<dim;k++){
                        temp[j][k]=SCHMID_TENSOR1[dim*i+j][k];
                    }
                }

                temp2.reinit(dim,dim);
                temp.mmult(temp2,FP_t2);
                    for (unsigned int j=0;j<dim;j++){
                    for (unsigned int k=0;k<dim;k++){
                       if(resolved_shear_tau(i)>0)
                           FP_tau[j][k]=FP_tau[j][k]+x_beta(i)*temp2[j][k];
                        else
                           FP_tau[j][k]=FP_tau[j][k]-x_beta(i)*temp2[j][k];
           }
          }

          }

                det_FP_tau=FP_tau.determinant();
                for (unsigned int j=0;j<dim;j++){
                    for (unsigned int k=0;k<dim;k++){
                        FP_tau[j][k]=FP_tau[j][k]*pow(det_FP_tau,-1/3);
                       }
                }


            // % % % % % STEP 8 % % % % %

				for (unsigned int i = 0;i < n_slip_systems;i++) {
					for (unsigned int j = 0;j < dim;j++) {
						for (unsigned int k = 0;k < dim;k++) {
							temp[j][k] = SCHMID_TENSOR1[dim*i + j][k];
						}
					}

					temp2.reinit(dim, dim);
					CE_tau.mmult(temp2, temp);
					temp3.reinit(dim, dim);
					for (unsigned int j = 0;j < dim;j++) {
						for (unsigned int k = 0;k < dim;k++) {
							temp3[j][k] = temp2[j][k] + temp2[k][j];
						}
					}
					temp4.reinit(dim, dim);
					tempv1 = 0.0;
					Dmat.vmult(tempv1, vecform(temp3));
					matform(temp4, tempv1);
					for (unsigned int j = 0;j < dim;j++) {
						for (unsigned int k = 0;k < dim;k++) {
							if (resolved_shear_tau(i) > 0)
								T_star_tau[j][k] = T_star_tau[j][k] - 0.5*x_beta(i)*temp4[j][k];
							else
								T_star_tau[j][k] = T_star_tau[j][k] + 0.5*x_beta(i)*temp4[j][k];
						}
					}
				}


            resolved_shear_tau=0.0;
            for(unsigned int i=0;i<n_slip_systems;i++){

                for (unsigned int j=0;j<dim;j++){
                    for (unsigned int k=0;k<dim;k++){
            resolved_shear_tau(i)+=T_star_tau[j][k]*SCHMID_TENSOR1[dim*i+j][k];
                    }
                }
				if (i > this->userInputs.numSlipSystems - 1) {
					if (resolved_shear_tau(i) < 0)
						resolved_shear_tau(i) = 0;
				}
            }


            // % % % % % STEP 9 % % % % %



            double h1=0;
            for(unsigned int i=0;i<n_slip_systems;i++){
                h1=0;
                for(unsigned int j=0;j<n_slip_systems;j++){
                    h1=h1+h_alpha_beta_t(i,j)*x_beta(j);
                }
                s_alpha_tau(i)=s_alpha_tau(i)+h1;
            }


//           for (unsigned int i=0;i<this->userInputs.numSlipSystems;i++){//

//                    //if(s_alpha_tau(i)>(this->userInputs.saturationStress[i])){
//                      // s_alpha_tau(i)=this->userInputs.saturationStress[i];
////abort();//

//              }//

//            }//

//     for (unsigned int i=0;i<this->userInputs.numTwinSystems;i++){//

//                if(s_alpha_tau(this->userInputs.numSlipSystems+i)>(this->userInputs.saturationStressTwin[i])){
//                s_alpha_tau(this->userInputs.numSlipSystems+i)=this->userInputs.saturationStressTwin[i];

//abort();
//}


//    }


		iter1 = iter1 + 1;

  }

  ////////////////////////////////////End Nonlinear iteration for Slip increments////////////////////////////////////

  Fpn_inv = 0.0; Fpn_inv.invert(FP_tau);
  FE_tau = 0.0;
  F_tau.mmult(FE_tau, Fpn_inv);
  temp.reinit(dim, dim);
  det_FE_tau = FE_tau.determinant();
  FE_tau.mmult(temp, T_star_tau); temp.equ(1.0 / det_FE_tau, temp);
  temp.mTmult(T_tau, FE_tau);

  det_F_tau = F_tau.determinant();
  temp.invert(F_tau); T_tau.mTmult(P_tau, temp);
  P_tau.equ(det_F_tau, P_tau);

         for (unsigned int i=0;i<this->userInputs.numTwinSystems;i++){
            twinfraction_iter[cellID][quadPtID][i]=twinfraction_conv[cellID][quadPtID][i]+x_beta_old[i+this->userInputs.numSlipSystems]/this->userInputs.twinShear;
            twinfraction_iter_Twin[cellID][quadPtID][i] = twinfraction_conv_Twin[cellID][quadPtID][i] + x_beta_old[i + this->userInputs.numSlipSystems] / this->userInputs.twinShear;
}

        for (unsigned int i=0;i<this->userInputs.numSlipSystems;i++){
            slipfraction_iter[cellID][quadPtID][i]=slipfraction_conv[cellID][quadPtID][i]+x_beta_old[i];
        }

		n_PA = 0;
		for (unsigned int i = 0;i<n_slip_systems;i++) {
			if (x_beta_old(i) > 0) {
				if (n_PA == 0) {
					n_PA = n_PA + 1;
					PA.reinit(n_PA);
					PA(0) = i;
				}
				else {
					PA_temp = PA;
					n_PA = n_PA + 1;
					PA.reinit(n_PA);
					for (unsigned int j = 0;j<(n_PA - 1);j++) {
						PA(j) = PA_temp(j);
					}
					PA(n_PA - 1) = i;
					PA_temp.reinit(n_PA);     //%%%%% Potentially active slip systems
				}
			}
		}
		active.reinit(n_PA); active = PA;


			FullMatrix<double> F_temp(dim, dim);
			int iter = 0, iter2 = 0, iter3 = 0;
			Fpn_inv = 0.0; Fpn_inv.invert(FP_t);
			F_temp = Fpn_inv;
			FullMatrix<double> scratch_1, scratch_2, EtF;
			right(scratch_1, F_temp);
			for (unsigned int i = 0;i<dim;i++) {
				for (unsigned int j = 0;j<dim;j++) {
					F_temp[i][j] = F_trial[j][i];
				}
			}

			left(scratch_2, F_temp);
			F_temp.reinit(9, 9);
			scratch_1.mmult(F_temp, scratch_2);
			symmf(EtF, F_temp);




			FullMatrix<double> s(dim, dim), p(n_PA, dim*dim), dgammadEmat(dim, dim), dgammadEtrial(n_PA, dim*dim), dgammadF(n_PA, dim*dim);
			Vector<double> P_as_vec3(dim*dim), s1(dim*dim), s2(dim*dim), temps2(dim*dim);

			if (n_PA > 0) {

				s2 = 0.0;
				for (unsigned int alpha = 0;alpha<n_PA;alpha++) {
					iter = active(alpha);
					iter2 = 0;
					for (unsigned int j = 0;j<dim;j++) {
						for (unsigned int k = 0;k<dim;k++) {
							s(j, k) = SCHMID_TENSOR1(dim*iter + j, k);
							P_as_vec3(iter2) = s(j, k);
							iter2++;
						}
					}

					TM.Tvmult(s1, P_as_vec3);
					if (resolved_shear_tau_trial(iter)<0) {
						s1.equ(-1.0, s1);
					}
					s2 = 0.0;

					for (unsigned int beta = 0;beta<n_PA;beta++) {
						iter3 = active(beta);

						temp3.reinit(dim, dim);
						temp.reinit(dim, dim);
						for (unsigned int k = 0;k<dim;k++) {
							for (unsigned int l = 0;l<dim;l++) {
								temp3(k, l) = SCHMID_TENSOR1(dim*iter3 + k, l);
								temp(l, k) = temp3(k, l);
							}
						}

						temp1.reinit(dim*dim, dim*dim); temp2.reinit(dim*dim, dim*dim);
						right(temp1, temp3); left(temp2, temp);
						temp1.add(1.0, temp2);
						TM.mmult(temp2, temp1); temp2.Tvmult(temps2, P_as_vec3);

						if ((resolved_shear_tau_trial(iter)<0.0) ^ (resolved_shear_tau_trial(iter3)<0.0))
							temps2.equ(-1.0*x_beta_old(iter3), temps2);
						else
							temps2.equ(x_beta_old(iter3), temps2);

						s2+=temps2;

					}

					for (unsigned int index1 = 0;index1<dim*dim;index1++) {
						p(alpha, index1) = s1(index1) - s2(index1);
					}

				}

				A_PA.reinit(n_PA, n_PA);
				for (unsigned int i = 0;i<(n_PA);i++) {

					for (unsigned int j = 0;j<n_PA;j++) {
						A_PA[i][j] = A[PA(i)][PA(j)];
					}
				}

				temp.reinit(n_PA, n_PA); temp.invert(A_PA);
				temp.mmult(dgammadEtrial, p); dgammadEtrial.mmult(dgammadF, EtF);

			}

			s.reinit(dim*dim, dim*dim); s = 0.0;

			int alpha = 0;
			if (n_PA>0) {
				for (unsigned int i = 0;i<n_PA;i++) {
					alpha = active(i);
					iter = 0;
					temp3.reinit(dim, dim);temp.reinit(dim, dim);
					for (unsigned int j = 0;j<dim;j++) {
						for (unsigned int k = 0;k<dim;k++) {
							dgammadEmat[k][j] = dgammadEtrial[i][dim*j + k];
							temp[j][k] = B[dim*alpha + j][k];

						}
					}

					temp1.reinit(dim*dim, dim*dim); right(temp1, dgammadEmat);
					ElasticProd(temp3, temp, ElasticityTensor);
					temp2.reinit(dim*dim, dim*dim); tracev(temp2, temp1, temp3);

					if (resolved_shear_tau_trial(alpha)>0) {
						temp2.equ(-0.5, temp2);
					}
					else {
						temp2.equ(0.5, temp2);
					}
					s.add(1.0, temp2);
				}
			}

			FullMatrix<double> smat1(dim*dim, dim*dim), smat2(dim*dim, dim*dim), smat3(dim*dim, dim*dim);

			smat1 = 0.0; smat1.add(1.0, s); smat1.add(1.0, TM);
			smat3 = 0.0;

			if (n_PA>0) {

				for (unsigned int i = 0;i<n_PA;i++) {
					alpha = active(i);
					F_temp.reinit(dim, dim);
					for (unsigned int j = 0;j<dim;j++) {
						for (unsigned int k = 0;k<dim;k++) {
							F_temp[j][k] = SCHMID_TENSOR1[dim*alpha + j][k];
							temp3[k][j] = F_temp[j][k];
						}
					}

					right(temp1, F_temp);
					left(temp2, temp3);
					temp1.add(1.0, temp2);
					TM.mmult(temp2, temp1);
					if (resolved_shear_tau_trial(alpha)>0) {
						temp2.equ(-1.0*x_beta_old(alpha), temp2);
					}
					else {
						temp2.equ(1.0*x_beta_old(alpha), temp2);
					}
					smat3.add(1.0, temp2);
				}
			}


			FullMatrix<double> tangent_moduli(dim*dim, dim*dim), dgammadFmat(dim, dim);

			tangent_moduli = 0.0; tangent_moduli.add(1.0, smat1); tangent_moduli.add(1.0, smat3);
			smat2 = 0.0;

			if (n_PA>0) {

				for (unsigned int i = 0;i<n_PA;i++) {
					alpha = active(i);
					temp2.reinit(dim, dim);
					for (unsigned int j = 0;j<dim;j++) {
						for (unsigned int k = 0;k<dim;k++) {
							dgammadFmat[k][j] = dgammadF[i][dim*j + k];
							temp2[j][k] = SCHMID_TENSOR1[dim*alpha + j][k];
						}
					}

					right(temp1, dgammadFmat);
					F_trial.mmult(temp3, temp2);
					tracev(temp2, temp1, temp3);
					if (resolved_shear_tau_trial(alpha)>0) {
						temp2.equ(-1.0, temp2);
					}
					smat2.add(1.0, temp2);
				}
			}

			temp2.reinit(dim, dim);
			smat1.reinit(dim*dim, dim*dim);
			temp2.invert(FP_tau);
			right(smat1, temp2);


			FullMatrix<double> FeF(dim*dim, dim*dim);
			FeF = 0.0; FeF.add(1.0, smat1); FeF.add(1.0, smat2);
			temp1.reinit(dim, dim); temp1.invert(FE_tau);
			F_temp.reinit(dim, dim);


			//term 1; -trace(dFe Fe_inv) T * (1/det(Fe)) = [Term1]{dF}
			F_temp = temp1;
			FullMatrix<double> C4(dim*dim, dim*dim);

			right(C4, F_temp);

			FullMatrix<double> Term1(dim*dim, dim*dim), Term2(dim*dim, dim*dim), Term3(dim*dim, dim*dim), Term24(dim*dim, dim*dim);
			tracev(Term3, C4, T_tau);
			Term3.mmult(Term1, FeF); Term1.equ(-1.0, Term1);

			//term 3; Fe * dT_bar * FeT * (1/det(Fe)) = [Term3]{dF}

			F_temp.reinit(dim, dim);
			F_temp = FE_tau;

			FullMatrix<double> cauchy_Stiff(dim*dim, dim*dim);
			scratch_1.reinit(dim*dim, dim*dim); left(scratch_1, F_temp);
			temp1.reinit(dim, dim); temp1 = IdentityMatrix(dim);

			temp2.reinit(dim, dim); temp2 = F_temp;
			temp2.Tmmult(F_temp, temp1);

			right(C4, F_temp); C4.mmult(scratch_2, scratch_1);
			scratch_1.reinit(dim*dim, dim*dim); scratch_2.mmult(scratch_1, tangent_moduli);
			Term3.reinit(dim*dim, dim*dim); scratch_1.mmult(Term3, EtF);
			Term3.equ(1.0 / det_FE_tau, Term3);




			//term 2&4: (dFe * Fe_inv) * T + T * (dFe * Fe_inv)^T = 2Symm(dFe * Fe_inv * T) = [Term24]{dF}

			temp2.invert(FE_tau); temp2.mmult(F_temp, T_tau);
			right(C4, F_temp);
			temp3.reinit(dim*dim, dim*dim); temp3 = C4;
			symmf(C4, temp3);
			Term24 = 0.0; C4.mmult(Term24, FeF); Term24.equ(2.0, Term24);
			cauchy_Stiff = 0.0; cauchy_Stiff.add(1.0, Term1); cauchy_Stiff.add(1.0, Term24); cauchy_Stiff.add(1.0, Term3);
			temp1.reinit(dim, dim); temp1.invert(F_tau);

			C4.reinit(dim*dim, dim*dim); right(C4, temp1);
			Term1.reinit(dim*dim, dim*dim); tracev(Term1, C4, T_tau);

			temp2.reinit(dim, dim);temp2.invert(F_tau); temp2.mmult(F_temp, T_tau);
			scratch_1.reinit(dim*dim, dim*dim); right(scratch_1, F_temp);
			Term2.reinit(dim*dim, dim*dim); symmf(Term2, scratch_1);
			Term2.equ(2.0, Term2); Term2.add(-1.0, scratch_1); Term2.equ(-1.0, Term2);




			Term3.reinit(dim*dim, dim*dim); Term3 = cauchy_Stiff;

			PK1_Stiff = 0.0; PK1_Stiff.add(1.0, Term1); PK1_Stiff.add(1.0, Term2); PK1_Stiff.add(1.0, Term3);
			temp2.reinit(dim, dim); temp2 = IdentityMatrix(dim);
			temp3.reinit(dim, dim); temp3 = temp1; temp3.Tmmult(temp1, temp2);
			temp2.reinit(dim*dim, dim*dim); right(temp2, temp1);
			temp3.reinit(dim*dim, dim*dim); temp3 = PK1_Stiff;
			temp3.mmult(PK1_Stiff, temp2);
			PK1_Stiff.equ(det_F_tau, PK1_Stiff);



			dP_dF = 0.0;
			FullMatrix<double> L(dim, dim);
			L = IdentityMatrix(dim);
			for (unsigned int m = 0;m<dim;m++) {
				for (unsigned int n = 0;n<dim;n++) {
					for (unsigned int o = 0;o<dim;o++) {
						for (unsigned int p = 0;p<dim;p++) {
							for (unsigned int i = 0;i<dim;i++) {
								for (unsigned int j = 0;j<dim;j++) {
									for (unsigned int k = 0;k<dim;k++) {
										for (unsigned int l = 0;l<dim;l++) {
											dP_dF[m][n][o][p] = dP_dF[m][n][o][p] + PK1_Stiff(dim*i + j, dim*k + l)*L(i, m)*L(j, n)*L(k, o)*L(l, p);;
										}
									}
								}
							}
						}
					}
				}
			}



			P.reinit(dim, dim);
			P = P_tau;
			T = T_tau;


    sres_tau.reinit(n_slip_systems);
    sres_tau = s_alpha_tau;

    // Update the history variables
    Fe_iter[cellID][quadPtID]=FE_tau;
    Fp_iter[cellID][quadPtID]=FP_tau;
    s_alpha_iter[cellID][quadPtID]=sres_tau;


}

#include "../../../include/crystalPlasticity_template_instantiations.h"
