#include "../../../include/crystalPlasticity.h"


//////////////////////////////////////////////////////////////////////////
//calculatePlasticity.cc numerically integrates the constitive model.
//This calculatePlasticity.cc is the modified version of the following rate-independent crystal plasticity model:
//Mohammadreza Yaghoobi, John E. Allison, Veera Sundararaghavan,
//Multiscale modeling of twinning and detwinning behavior of HCP polycrystals,
// International Journal of Plasticity, December 2019, 102653.
//
//To use this file, one should copy it into the following folder (replacing the original calculatePlasticity.cc inside the following folder with this new one):
//    plasticity/src/materialModels/crystalPlasticity/
// Finaly, one should recompile PRISMS-Plasticity.
//
//This model includes a multiscale scheme to capture the twinning and detwinning mechanisms during
//cyclic loading of HCP polycrystals.
//////////////////////////////////////////////////////////////////////////

template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
  unsigned int quadPtID, unsigned int StiffnessCalFlag)
  {

    unsigned int n_slip_systemsWOtwin = this->userInputs.numSlipSystems1;
    n_Tslip_systems = n_slip_systemsWOtwin;
    n_twin_systems = 0;
    if (this->userInputs.enableTwinning1) {
      n_Tslip_systems += this->userInputs.numTwinSystems1 * 2;
      n_twin_systems = this->userInputs.numTwinSystems1 * 2;
    }
    else {
      n_Tslip_systems += 1;
      n_twin_systems = 1;
    }
    std::vector<double> ttwinvf(this->userInputs.numTwinSystems1);
    ttwinvf = twinfraction_conv[cellID][quadPtID];

    unsigned int tTwinMaxFlag= TwinMaxFlag_conv[cellID][quadPtID];
    unsigned int n_twin_systems_Size = this->userInputs.numTwinSystems1;
    unsigned int alpha=0;
    n_slip_systems = n_Tslip_systems;
    double det_F_tau;
    Vector<double> initialHardeningM(n_slip_systemsWOtwin), powerLawExp(n_slip_systemsWOtwin), saturationStressV(n_slip_systemsWOtwin);
    Vector<double> initialHardeningMTwin(n_twin_systems), powerLawExpTwin(n_twin_systems), saturationStressVTwin(n_twin_systems);
    FullMatrix<double> M_alpha(n_Tslip_systems, n_Tslip_systems), N_alpha(n_Tslip_systems, n_Tslip_systems);

    F_tau=F; // Deformation Gradient
    std::vector<unsigned int> activeTwinSystems(n_twin_systems_Size), tTwinFlag(n_twin_systems_Size);

    unsigned int NumberOfTwinnedRegionK = NumberOfTwinnedRegion_conv[cellID][quadPtID];
    activeTwinSystems =ActiveTwinSystems_conv[cellID][quadPtID];
    tTwinFlag = TwinFlag_conv[cellID][quadPtID];
    unsigned int RegionSize = NumberOfTwinnedRegionK + 1;
    FullMatrix<double> PK1_Stiff_RVE(dim*dim, dim*dim);
    FullMatrix<double> temp(dim, dim), temp1(dim, dim), temp2(dim, dim), temp3(dim, dim), temp4(dim, dim), temp5(dim, dim), temp6(dim, dim), temp7(dim*dim,dim*dim); // Temporary matrices
    FullMatrix<double> T_tau(dim, dim);
    FullMatrix<double> Fpn_inv(dim, dim), FE_tau_trial(dim, dim), F_trial(dim, dim), CE_tau_trial(dim, dim), FP_t2(dim, dim), Ee_tau_trial(dim, dim);
    temp = 0;temp7 = 0;

    std::vector<std::vector<FullMatrix<double> > >	dgammadEmat1;
    std::vector<FullMatrix<double> >	PK1_Stiff_Region, P_tau, T_star_tau_Region;
    dgammadEmat1.resize(RegionSize, std::vector<FullMatrix<double> >(n_twin_systems_Size, temp));
    PK1_Stiff_Region.resize(RegionSize, FullMatrix<double>(temp7));
    P_tau.resize(RegionSize, FullMatrix<double>(temp));
    T_star_tau_Region.resize(RegionSize, FullMatrix<double>(temp));

    FullMatrix<double> tslipvfsys(n_twin_systems_Size + 1, n_Tslip_systems);
    tslipvfsys = 0;





    for (unsigned int Region = 0; Region < NumberOfTwinnedRegionK + 1; Region++) {
      if (Region == 0) {
        n_twin_systems = n_twin_systems_Size * 2;
        M_alpha.reinit(n_Tslip_systems, n_Tslip_systems); N_alpha.reinit(n_Tslip_systems, n_Tslip_systems);
        M_alpha = m_alpha;
        N_alpha = n_alpha;
        alpha = 0;
        n_slip_systems = n_Tslip_systems;

        q.reinit(n_slip_systems,n_slip_systems);
        q=q_phase1;

        initialHardeningMTwin.reinit(n_twin_systems); powerLawExpTwin.reinit(n_twin_systems); saturationStressVTwin.reinit(n_twin_systems);
        for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {
          initialHardeningM[i] = this->userInputs.initialHardeningModulus1[i];
          saturationStressV[i] = this->userInputs.saturationStress1[i];
          powerLawExp[i] = this->userInputs.powerLawExponent1[i];
        }
        for (unsigned int i = 0;i < n_twin_systems;i++) {
          initialHardeningMTwin[i] = this->userInputs.initialHardeningModulusTwin1[i];
          saturationStressVTwin[i] = this->userInputs.saturationStressTwin1[i];
          powerLawExpTwin[i] = this->userInputs.powerLawExponentTwin1[i];
        }
      }
      else {
        n_twin_systems = 2;
        M_alpha.reinit(n_slip_systemsWOtwin + 2, n_slip_systemsWOtwin + 2); N_alpha.reinit(n_slip_systemsWOtwin + 2, n_slip_systemsWOtwin + 2);
        alpha = activeTwinSystems[Region - 1];
        for (unsigned int j = 0;j < n_slip_systemsWOtwin;j++) {
          for (unsigned int i = 0;i < dim;i++) {
            M_alpha[j][i] = m_alpha[j][i];
            N_alpha[j][i] = n_alpha[j][i];
          }
        }
        for (unsigned int i = 0;i < dim;i++) {
          M_alpha[n_slip_systemsWOtwin ][i] = m_alpha[n_slip_systemsWOtwin + alpha-1][i];
          N_alpha[n_slip_systemsWOtwin ][i] = n_alpha[n_slip_systemsWOtwin + alpha-1][i];
          M_alpha[n_slip_systemsWOtwin + 1][i] = m_alpha[n_slip_systemsWOtwin + alpha-1][i];
          N_alpha[n_slip_systemsWOtwin + 1][i] = n_alpha[n_slip_systemsWOtwin + alpha-1][i];
        }
        n_slip_systems = n_slip_systemsWOtwin + 2;

        q.reinit(n_slip_systems,n_slip_systems);
        //Here, we assumed that twinning and detwinning all have equal latent hardening ratios.
        for(unsigned int i=0;i<n_slip_systems;i++){
          for(unsigned int j=0;j<n_slip_systems;j++){
            q[i][j] = q_phase1[i][j];
          }
        }

        initialHardeningMTwin.reinit(n_twin_systems); powerLawExpTwin.reinit(n_twin_systems); saturationStressVTwin.reinit(n_twin_systems);
        for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {
          initialHardeningM[i] = this->userInputs.initialHardeningModulus1[i];
          saturationStressV[i] = this->userInputs.saturationStress1[i];
          powerLawExp[i] = this->userInputs.powerLawExponent1[i];
        }
        for (unsigned int i = 0;i < n_twin_systems;i++) {
          initialHardeningMTwin[i] = this->userInputs.initialHardeningModulusTwin1[i+ (n_twin_systems_Size )];
          saturationStressVTwin[i] = this->userInputs.saturationStressTwin1[i + (n_twin_systems_Size)];
          powerLawExpTwin[i] = this->userInputs.powerLawExponentTwin1[i + (n_twin_systems_Size)];
        }
      }

      FullMatrix<double> FE_t(dim, dim), FP_t(dim, dim);  //Elastic and Plastic deformation gradient
      Vector<double> s_alpha_t(n_slip_systems); // Slip resistance
      Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)

      // Tolerance

      double tol1 = this->userInputs.modelStressTolerance;


      std::cout.precision(16);
      for (unsigned int i = 0;i < dim;i++) {
        for (unsigned int j = 0;j < dim;j++) {
          FE_t[i][j] = Fe_conv[cellID][quadPtID][i][j + alpha*dim];
          FP_t[i][j] = Fp_conv[cellID][quadPtID][i][j + alpha*dim];
        }
      }
      for (unsigned int i = 0;i < n_slip_systems;i++) {
        s_alpha_t[i] = s_alpha_conv[cellID][quadPtID][alpha*n_Tslip_systems + i];
      }

      for (unsigned int i = 0;i < dim;i++) {
        rot1[i] = rot[cellID][quadPtID][alpha*dim + i];
      }



      // Rotation matrix of the crystal orientation
      FullMatrix<double> rotmat(dim, dim);
      rotmat = 0.0;
      odfpoint(rotmat, rot1);

      // Calculation of Schmid Tensors  and B= symm(FE_tau_trial'*FE_tau_trial*S_alpha)
      FullMatrix<double> SCHMID_TENSOR1(n_slip_systems*dim, dim), B(n_slip_systems*dim, dim);
      Vector<double> m1(dim), n1(dim);


      // Elastic Modulus

      FullMatrix<double> Dmat2(2 * dim, 2 * dim), TM(dim*dim, dim*dim), ElasticityTensor(2 * dim, 2 * dim), elasticStiffnessMatrix(2 * dim, 2 * dim);
      Vector<double> vec1(2 * dim), vec2(dim*dim);
      for (unsigned int i = 0;i < 6;i++) {
        for (unsigned int j = 0;j < 6;j++) {
          elasticStiffnessMatrix[i][j] = this->userInputs.elasticStiffness1[i][j];
        }
      }

      elasticmoduli(Dmat2, rotmat, elasticStiffnessMatrix);

      //Elastic Stiffness Matrix Dmat
      Dmat.reinit(6, 6); Dmat = 0.0;

      for (unsigned int i = 0;i < 6;i++) {
        for (unsigned int j = 0;j < 6;j++) {
          Dmat[i][j] = Dmat2[i][j];
        }
      }


      for (unsigned int i = 0;i < 6;i++) {
        for (unsigned int j = 3;j < 6;j++) {
          Dmat[i][j] = 2 * Dmat[i][j];
        }
      }

      vec1(0) = 0;vec1(1) = 5;vec1(2) = 4;vec1(3) = 1;vec1(4) = 3;vec1(5) = 2;
      vec2(0) = 0;vec2(1) = 5;vec2(2) = 4;vec2(3) = 5;vec2(4) = 1;vec2(5) = 3;vec2(6) = 4;vec2(7) = 3;vec2(8) = 2;



      for (unsigned int i = 0;i < 9;i++) {
        for (unsigned int j = 0;j < 9;j++) {
          TM[i][j] = Dmat2(vec2(i), vec2(j));
        }
      }

      for (unsigned int i = 0;i < 6;i++) {
        for (unsigned int j = 0;j < 6;j++) {
          ElasticityTensor[i][j] = Dmat(vec1(i), vec1(j));
        }
      }


      Vector<double> s_alpha_tau;
      Vector<double> s_beta(n_slip_systems), h_beta(n_slip_systems), delh_beta_dels(n_slip_systems), h0(n_slip_systems), a_pow(n_slip_systems), s_s(n_slip_systems);
      FullMatrix<double> h_alpha_beta_t(n_slip_systems, n_slip_systems), A(n_slip_systems, n_slip_systems);
      FullMatrix<double> del_FP(dim, dim);
      FullMatrix<double> A_PA;
      Vector<double> active, active_First;
      Vector<double> PA, PA_temp(1);
      Vector<double> resolved_shear_tau_trial(n_slip_systems), b(n_slip_systems), resolved_shear_tau(n_slip_systems);
      Vector<double> x_beta_old(n_slip_systems);

      Vector<double> x_beta(n_slip_systems);

      FullMatrix<double> PK1_Stiff(dim*dim, dim*dim);
      FullMatrix<double> CE_tau(dim, dim), T_star_tau(dim, dim);
      FullMatrix<double> T_star_tau_trial(dim, dim);



      double det_FE_tau,det_FP_tau;
      unsigned int n_PA = 0;	// Number of active slip systems


      FP_tau = FP_t;
      Fpn_inv = 0.0; Fpn_inv.invert(FP_t);
      s_alpha_tau = s_alpha_t;

      if (this->userInputs.enableOneTwinSys_Reorien){
        for (unsigned int i = 0;i < n_twin_systems;i++) {//
          if (s_alpha_tau(n_slip_systemsWOtwin + i) > 10000*saturationStressVTwin[i]) {
            saturationStressVTwin[i]=s_alpha_tau(n_slip_systemsWOtwin + i)*10;
          }
        }
      }

      FE_tau_trial = 0.0;
      F_trial = 0.0;
      F_tau.mmult(FE_tau_trial, Fpn_inv);F_trial = FE_tau_trial;
      temp2.reinit(dim, dim);
      temp.reinit(dim, dim); temp = 0.0;
      temp = FE_tau_trial;
      FE_tau_trial.Tmmult(CE_tau_trial, temp);
      Ee_tau_trial = CE_tau_trial;
      temp = IdentityMatrix(dim);
      for (unsigned int i = 0;i < dim;i++) {
        for (unsigned int j = 0;j < dim;j++) {
          Ee_tau_trial[i][j] = 0.5*(Ee_tau_trial[i][j] - temp[i][j]);
        }
      }

      for (unsigned int i = 0;i < n_slip_systems;i++) {
        for (unsigned int j = 0;j < dim;j++) {
          m1(j) = M_alpha[i][j];
          n1(j) = N_alpha[i][j];
        }

        temp = 0.0;
        temp2 = 0.0;
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            temp[j][k] = m1(j)*n1(k);
          }
        }
        //convert the Shmitch tensor to Sample coordinates
        rotmat.mmult(temp2, temp);
        temp2.mTmult(temp, rotmat);
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {

            SCHMID_TENSOR1[dim*i + j][k] = temp[j][k];
          }
        }
        CE_tau_trial.mmult(temp2, temp);
        temp2.symmetrize();
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            B[dim*i + j][k] = 2 * temp2[j][k];
          }
        }
      }


      //% % % % % STEP 2 % % % % %
      // Calculate the trial stress T_star_tau_trial
      Vector<double> tempv1(6), tempv2(6);
      tempv1 = 0.0;
      Dmat.vmult(tempv1, vecform(Ee_tau_trial));
      matform(T_star_tau_trial, tempv1);
      T_star_tau.equ(1.0, T_star_tau_trial);



      //% % % % % STEP 3 % % % % %
      // Calculate the trial resolved shear stress resolved_shear_tau_trial for each slip system

      resolved_shear_tau_trial = 0.0;

      for (unsigned int i = 0;i < n_slip_systems;i++) {

        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            resolved_shear_tau_trial[i] += T_star_tau_trial[j][k] * SCHMID_TENSOR1[dim*i + j][k];
          }
        }

    }

    if (Region == 0) {
      for (unsigned int i = 0;i < (n_twin_systems / 2);i++) {
        unsigned int j = i + n_slip_systemsWOtwin;
        if ((resolved_shear_tau_trial[j] < 0) || (tTwinMaxFlag == 0)) {
          resolved_shear_tau_trial[j] = 0;
        }
      }

      /////////////The modification compared to the original model is applied to the following lines//////////////////
      for (unsigned int i = (n_twin_systems / 2);i < n_twin_systems;i++) {
        unsigned int j = i + n_slip_systemsWOtwin;
        if ((resolved_shear_tau_trial[j] > 0) || (ttwinvf[i - (n_twin_systems / 2)] <= 0)) {
          resolved_shear_tau_trial[j] = 0;
        }
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    }
    else {

      if (resolved_shear_tau_trial[n_slip_systemsWOtwin] < 0) {
        resolved_shear_tau_trial[n_slip_systemsWOtwin] = 0;
      }
      if ((resolved_shear_tau_trial[1+n_slip_systemsWOtwin] > 0)||(tTwinMaxFlag == 0)) {
        resolved_shear_tau_trial[1+n_slip_systemsWOtwin] = 0;
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

      if (iter1 > this->userInputs.modelMaxSlipSearchIterations) {
        flag2 = 1;
        break;
      }


      //% % % % % STEP 4 % % % % %
      n_PA = 0;	// Number of active slip systems
      b.reinit(n_slip_systems);
      //Determine the set set of the n potentially active slip systems
      for (unsigned int i = 0;i < n_slip_systems;i++) {
        b(i) = fabs(resolved_shear_tau(i)) - s_alpha_tau(i);
        if (b(i) >= tol1) {
          if (n_PA == 0) {
            n_PA = n_PA + 1;
            PA.reinit(n_PA);
            PA(0) = i;
          }
          else {
            PA_temp = PA;
            n_PA = n_PA + 1;
            PA.reinit(n_PA);
            for (unsigned int j = 0;j < (n_PA - 1);j++) {
              PA(j) = PA_temp(j);
            }
            PA(n_PA - 1) = i;
            PA_temp.reinit(n_PA);     //%%%%% Potentially active slip systems
          }
        }
      }



      if (n_PA == 0)
      break;

      Vector<double> b_PA(n_PA);
      b_PA.reinit(n_PA);

      for (unsigned int i = 0;i < n_PA;i++) {
        b_PA(i) = b(PA(i));
      }

      if ((b_PA.linfty_norm()) < tol1)
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
      s_beta = s_alpha_tau;

      // Single slip hardening rate
      for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {
        h_beta(i) = initialHardeningM[i] * pow((1 - s_beta(i) / saturationStressV[i]), powerLawExp[i]);
      }

      for (unsigned int i = 0;i < n_twin_systems;i++) {
        h_beta(n_slip_systemsWOtwin + i) = initialHardeningMTwin[i] * pow((1 - s_beta(n_slip_systemsWOtwin + i) / saturationStressVTwin[i]), powerLawExpTwin[i]);
      }


      for (unsigned int i = 0;i < n_slip_systems;i++) {
        for (unsigned int j = 0;j < n_slip_systems;j++) {
          h_alpha_beta_t[i][j] = q[i][j] * h_beta(j);
          A[i][j] = h_alpha_beta_t[i][j];
        }
      }




      for (unsigned int j = 0;j < n_slip_systems;j++) {
        temp2.reinit(dim, dim);temp = 0.0;
        for (unsigned int k = 0;k < dim;k++) {
          for (unsigned int l = 0;l < dim;l++) {
            temp2[k][l] = 0.5* B(dim*j + k, l);
          }
        }
        tempv1 = 0.0; Dmat.vmult(tempv1, vecform(temp2));
        temp3 = 0.0; matform(temp3, tempv1);
        //temp3 is symm in Matlab line 94 of constitutive.m

        for (unsigned int i = 0;i < n_slip_systems;i++) {
          temp1.reinit(dim, dim);temp1 = 0.0;
          for (unsigned int k = 0;k < dim;k++) {
            for (unsigned int l = 0;l < dim;l++) {
              temp1[k][l] = SCHMID_TENSOR1(dim*i + k, l);
            }
          }
          for (unsigned int k = 0;k < dim;k++) {
            for (unsigned int l = 0;l < dim;l++) {
              if ((resolved_shear_tau(i)*resolved_shear_tau(j)) < 0.0)
              A[i][j] -= temp1[k][l] * temp3[k][l];
              else
              A[i][j] += temp1[k][l] * temp3[k][l];

            }
          }
        }
      }

      //Modified slip system search for adding corrective term
      inactive_slip_removal(active, x_beta_old, x_beta, n_PA, n_slip_systems, PA, b, A, A_PA);

      temp.reinit(dim, dim);
      del_FP.reinit(dim, dim);
      del_FP = 0.0;
      for (unsigned int i = 0;i < n_slip_systems;i++) {
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            temp[j][k] = SCHMID_TENSOR1[dim*i + j][k];
          }
        }

        temp2.reinit(dim, dim);
        temp.mmult(temp2, FP_t2);
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            if (resolved_shear_tau(i) > 0)
            FP_tau[j][k] = FP_tau[j][k] + x_beta(i)*temp2[j][k];
            else
            FP_tau[j][k] = FP_tau[j][k] - x_beta(i)*temp2[j][k];
          }
        }

      }

      det_FP_tau = FP_tau.determinant();
      for (unsigned int j = 0;j < dim;j++) {
        for (unsigned int k = 0;k < dim;k++) {
          FP_tau[j][k] = FP_tau[j][k] * pow(det_FP_tau, -1 / 3);
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


      resolved_shear_tau = 0.0;
      for (unsigned int i = 0;i < n_slip_systems;i++) {

        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            resolved_shear_tau(i) += T_star_tau[j][k] * SCHMID_TENSOR1[dim*i + j][k];
          }
        }
        /*if (i > this->userInputs.numSlipSystems - 1) {
        if (resolved_shear_tau(i) < 0)
        resolved_shear_tau(i) = 0;
      }*/
    }

    if (Region == 0) {
      for (unsigned int i = 0;i < (n_twin_systems / 2);i++) {
        unsigned int j = i + n_slip_systemsWOtwin;
        if ((resolved_shear_tau_trial[j] < 0) || (tTwinMaxFlag == 0)) {
          resolved_shear_tau_trial[j] = 0;
        }
      }

      /////////////The modification compared to the original model is applied to the following lines//////////////////
      for (unsigned int i = (n_twin_systems / 2);i < n_twin_systems;i++) {
        unsigned int j = i + n_slip_systemsWOtwin;
        if ((resolved_shear_tau_trial[j] > 0) || (ttwinvf[i - (n_twin_systems / 2)] <= 0)) {
          resolved_shear_tau_trial[j] = 0;
        }
      }
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    }
    else {

      if (resolved_shear_tau_trial[n_slip_systemsWOtwin] < 0) {
        resolved_shear_tau_trial[n_slip_systemsWOtwin] = 0;
      }
      if ((resolved_shear_tau_trial[1 + n_slip_systemsWOtwin] > 0) || (tTwinMaxFlag == 0)) {
        resolved_shear_tau_trial[1 + n_slip_systemsWOtwin] = 0;
      }

    }


    // % % % % % STEP 9 % % % % %



    double h1 = 0;
    for (unsigned int i = 0;i < n_slip_systems;i++) {
      h1 = 0;
      for (unsigned int j = 0;j < n_slip_systems;j++) {
        h1 = h1 + h_alpha_beta_t(i, j)*x_beta(j);
      }
      s_alpha_tau(i) = s_alpha_tau(i) + h1;
    }

    for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {//

      if (s_alpha_tau(i) > saturationStressV[i]) {
        s_alpha_tau(i) = saturationStressV[i];

      }//

    }//

    for (unsigned int i = 0;i < n_twin_systems;i++) {//

      if (s_alpha_tau(n_slip_systemsWOtwin + i) > saturationStressVTwin[i]) {
        s_alpha_tau(n_slip_systemsWOtwin + i) = saturationStressVTwin[i];

        //abort();
      }

    }


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
  temp2.reinit(dim, dim);
  det_F_tau = F_tau.determinant();
  temp.invert(F_tau); T_tau.mTmult(temp2, temp);
  temp2.equ(det_F_tau, temp2);

  for (unsigned int i = 0;i < dim;i++) {
    for (unsigned int j = 0;j < dim;j++) {
      T_star_tau_Region[Region][i][j] = T_star_tau[i][j];
    }
  }

  for (unsigned int i = 0;i < dim;i++) {
    for (unsigned int j = 0;j < dim;j++) {
      P_tau[Region][i][j] = temp2[i][j];
    }
  }

FullMatrix<double> F_temp(dim, dim);
int iter = 0, iter2 = 0, iter3 = 0;
Fpn_inv = 0.0; Fpn_inv.invert(FP_t);
F_temp = Fpn_inv;
FullMatrix<double> scratch_1, scratch_2, EtF;
right(scratch_1, F_temp);
for (unsigned int i = 0;i < dim;i++) {
  for (unsigned int j = 0;j < dim;j++) {
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
  for (unsigned int alpha1 = 0;alpha1 < n_PA;alpha1++) {
    iter = active(alpha1);
    iter2 = 0;
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        s(j, k) = SCHMID_TENSOR1(dim*iter + j, k);
        P_as_vec3(iter2) = s(j, k);
        iter2++;
      }
    }

    TM.Tvmult(s1, P_as_vec3);
    if (resolved_shear_tau_trial(iter) < 0) {
      s1.equ(-1.0, s1);
    }
    s2 = 0.0;

    for (unsigned int beta = 0;beta < n_PA;beta++) {
      iter3 = active(beta);

      temp3.reinit(dim, dim);
      temp.reinit(dim, dim);
      for (unsigned int k = 0;k < dim;k++) {
        for (unsigned int l = 0;l < dim;l++) {
          temp3(k, l) = SCHMID_TENSOR1(dim*iter3 + k, l);
          temp(l, k) = temp3(k, l);
        }
      }

      temp1.reinit(dim*dim, dim*dim); temp2.reinit(dim*dim, dim*dim);
      right(temp1, temp3); left(temp2, temp);
      temp1.add(1.0, temp2);
      TM.mmult(temp2, temp1); temp2.Tvmult(temps2, P_as_vec3);

      if ((resolved_shear_tau_trial(iter) < 0.0) ^ (resolved_shear_tau_trial(iter3) < 0.0))
      temps2.equ(-1.0*x_beta_old(iter3), temps2);
      else
      temps2.equ(x_beta_old(iter3), temps2);

      s2 += temps2;

    }

    for (unsigned int index1 = 0;index1 < dim*dim;index1++) {
      p(alpha1, index1) = s1(index1) - s2(index1);
    }

  }

  A_PA.reinit(n_PA, n_PA);

  for (unsigned int i = 0;i < (n_PA);i++) {
    for (unsigned int j = 0;j < n_PA;j++) {
      A_PA[i][j] = A[PA(i)][PA(j)];

    }
  }

  temp.reinit(n_PA,n_PA);
  temp.invert(A_PA);


  temp.mmult(dgammadEtrial, p); dgammadEtrial.mmult(dgammadF, EtF);

}


std::vector<unsigned int> myPA(n_PA);
for (unsigned int i = 0;i < n_PA;i++) {
  myPA[i] = PA[i];
}
for (unsigned int i = 0;i < (n_twin_systems / 2);i++) {

  int num_items1 = std::count(myPA.begin(), myPA.end(), (i + n_slip_systemsWOtwin));
  int num_items2 = std::count(myPA.begin(), myPA.end(), (i + n_twin_systems / 2 + n_slip_systemsWOtwin));
  if (num_items1 == 1) {
    std::vector<unsigned int>::iterator result1;
    result1=std::find(myPA.begin(), myPA.end(), (i + n_slip_systemsWOtwin));
    unsigned int pos1 = std::distance(myPA.begin(), result1);
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        dgammadEmat1[Region][i][k][j] = dgammadF[pos1][dim*j + k];
      }
    }
  }
  if (num_items2 == 1) {
    std::vector<unsigned int>::iterator result2;
    result2=std::find(myPA.begin(), myPA.end(), (i + n_twin_systems / 2 + n_slip_systemsWOtwin));
    unsigned int pos2 = std::distance(myPA.begin(), result2);
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        dgammadEmat1[Region][i][k][j] = -dgammadF[pos2][dim*j + k];
      }
    }
  }
}




s.reinit(dim*dim, dim*dim); s = 0.0;

int alpha1 = 0;
if (n_PA > 0) {
  for (unsigned int i = 0;i < n_PA;i++) {
    alpha1 = active(i);
    iter = 0;
    temp3.reinit(dim, dim);temp.reinit(dim, dim);
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        dgammadEmat[k][j] = dgammadEtrial[i][dim*j + k];
        temp[j][k] = B[dim*alpha1 + j][k];

      }
    }

    temp1.reinit(dim*dim, dim*dim); right(temp1, dgammadEmat);
    ElasticProd(temp3, temp, ElasticityTensor);
    temp2.reinit(dim*dim, dim*dim); tracev(temp2, temp1, temp3);

    if (resolved_shear_tau_trial(alpha1) > 0) {
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

if (n_PA > 0) {

  for (unsigned int i = 0;i < n_PA;i++) {
    alpha1 = active(i);
    F_temp.reinit(dim, dim);
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        F_temp[j][k] = SCHMID_TENSOR1[dim*alpha1 + j][k];
        temp3[k][j] = F_temp[j][k];
      }
    }

    right(temp1, F_temp);
    left(temp2, temp3);
    temp1.add(1.0, temp2);
    TM.mmult(temp2, temp1);
    if (resolved_shear_tau_trial(alpha1) > 0) {
      temp2.equ(-1.0*x_beta_old(alpha1), temp2);
    }
    else {
      temp2.equ(1.0*x_beta_old(alpha1), temp2);
    }
    smat3.add(1.0, temp2);
  }
}


FullMatrix<double> tangent_moduli(dim*dim, dim*dim), dgammadFmat(dim, dim);

tangent_moduli = 0.0; tangent_moduli.add(1.0, smat1); tangent_moduli.add(1.0, smat3);
smat2 = 0.0;

if (n_PA > 0) {

  for (unsigned int i = 0;i < n_PA;i++) {
    alpha1 = active(i);
    temp2.reinit(dim, dim);
    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        dgammadFmat[k][j] = dgammadF[i][dim*j + k];
        temp2[j][k] = SCHMID_TENSOR1[dim*alpha1 + j][k];
      }
    }

    right(temp1, dgammadFmat);
    F_trial.mmult(temp3, temp2);
    tracev(temp2, temp1, temp3);
    if (resolved_shear_tau_trial(alpha1) > 0) {
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


for (unsigned int k = 0;k < dim*dim;k++) {
  for (unsigned int l = 0;l < dim*dim;l++) {
    PK1_Stiff_Region[Region][k][l] = PK1_Stiff[k][l];
  }
}



for (unsigned int i = 0;i < dim;i++) {
  for (unsigned int j = 0;j < dim;j++) {
    Fe_iter[cellID][quadPtID][i][j + alpha*dim] = FE_tau[i][j];
    Fp_iter[cellID][quadPtID][i][j + alpha*dim] = FP_tau[i][j];
  }
}
// Update the history variables

sres_tau.reinit(n_slip_systems);
sres_tau = s_alpha_tau;
for (unsigned int i = 0;i < n_slip_systems;i++) {
  s_alpha_iter[cellID][quadPtID][alpha*n_Tslip_systems + i] = sres_tau[i];
}
for (unsigned int i = 0;i < n_slip_systems;i++) {
  tslipvfsys[alpha][i] = x_beta_old[i];
}
Vector<double> rold(dim), rnew(dim);

for (unsigned int i = 0;i < dim;i++) {
  rold[i] = rotnew_conv[cellID][quadPtID][alpha*dim + i];
}

reorient2(rnew,rold,FE_tau,FE_t);
for (unsigned int i = 0;i < dim;i++) {
  rotnew_iter[cellID][quadPtID][alpha*dim + i] = rnew[i];
}


}

std::vector<unsigned int> tActiveTwinSystems(NumberOfTwinnedRegionK), DeactiveTwinSystems(n_twin_systems_Size - NumberOfTwinnedRegionK);
unsigned int jj = 0;
unsigned int kk = 0;
for (unsigned int i = 0;i < n_twin_systems_Size;i++) {
  if (activeTwinSystems[i] > 0) {
    tActiveTwinSystems[jj] = activeTwinSystems[i];
    jj = jj + 1;
  }
}
for (unsigned int j = 1;j < n_twin_systems_Size + 1;j++) {
  unsigned int t = 1;
  for (unsigned int i = 0;i < NumberOfTwinnedRegionK;i++) {
    if (j == tActiveTwinSystems[i]) {
      t = 0;
      break;
    }
  }
  if (t == 1){
    DeactiveTwinSystems[kk] = j;
    kk = kk + 1;
  }
}
unsigned int NumberOfNonTwinnedRegionK = n_twin_systems_Size - NumberOfTwinnedRegionK;
Vector<double> twinvfNucleated(NumberOfTwinnedRegionK), RegionsTwinvf(1 + NumberOfTwinnedRegionK);
double Totaltwinvf = 0;
double twintemp = 0;
RegionsTwinvf = 0;
twinvfNucleated = 0;
for (unsigned int i = 0;i < NumberOfTwinnedRegionK;i++) {
  twintemp = ttwinvf[tActiveTwinSystems[i]-1];
  twinvfNucleated[i] = twintemp;
  Totaltwinvf = Totaltwinvf + twintemp;
  RegionsTwinvf[i + 1] = twintemp;
}
RegionsTwinvf[0] = 1 - Totaltwinvf;



P.reinit(dim, dim);
P = 0;
for (unsigned int Region = 0;Region < NumberOfTwinnedRegionK + 1;Region++) {
  for (unsigned int i = 0;i < dim;i++) {
    for (unsigned int j = 0;j < dim;j++) {
      P[i][j] = P[i][j] + RegionsTwinvf[Region] * P_tau[Region][i][j];
    }
  }
}

T.reinit(dim, dim);
T = 0;
temp2.reinit(dim, dim);
temp.reinit(dim, dim);
temp=F_tau; P.mTmult(T, temp);
T.equ(1/det_F_tau, T);

PK1_Stiff_RVE = 0;


for (unsigned int Region = 0;Region < NumberOfTwinnedRegionK + 1;Region++) {
  for (unsigned int i = 0;i < dim*dim;i++) {
    for (unsigned int j = 0;j < dim*dim;j++) {
      PK1_Stiff_RVE[i][j] = PK1_Stiff_RVE[i][j] + RegionsTwinvf[Region] * PK1_Stiff_Region[Region][i][j];
    }
  }
}

FullMatrix<double> UntwinnedRegionsTwinvfDiff(dim, dim), dgammadEmatRegion(dim,dim);
UntwinnedRegionsTwinvfDiff = 0;
dgammadEmatRegion = 0;
for (unsigned int ii = 0;ii < NumberOfTwinnedRegionK;ii++) {
  int alpha = activeTwinSystems[ii];
  for (unsigned int i = 0;i < dim;i++) {
    for (unsigned int j = 0;j < dim;j++) {
      dgammadEmatRegion[i][j] = ((1 - RegionsTwinvf[0])*dgammadEmat1[0][alpha-1][i][j] - RegionsTwinvf[1 + ii] * dgammadEmat1[1 + ii][0][i][j]) / this->userInputs.twinShear1;
      UntwinnedRegionsTwinvfDiff[i][j] = UntwinnedRegionsTwinvfDiff[i][j] + dgammadEmatRegion[i][j];
    }
  }
  temp1.reinit(dim*dim, dim*dim); right(temp1, dgammadEmatRegion);
  temp3.reinit(dim, dim);
  for (unsigned int i = 0;i < dim;i++) {
    for (unsigned int j = 0;j < dim;j++) {
      temp3[i][j] = P_tau[ii + 1][i][j];
    }
  }
  temp2.reinit(dim*dim, dim*dim); tracev(temp2, temp1, temp3);
  for (unsigned int i = 0;i < dim*dim;i++) {
    for (unsigned int j = 0;j < dim*dim;j++) {
      PK1_Stiff_RVE[i][j] = PK1_Stiff_RVE[i][j] + temp2[i][j];
    }
  }
}
temp1.reinit(dim*dim, dim*dim); right(temp1, UntwinnedRegionsTwinvfDiff);
temp3.reinit(dim, dim);
for (unsigned int i = 0;i < dim;i++) {
  for (unsigned int j = 0;j < dim;j++) {
    temp3[i][j] = P_tau[0][i][j];
  }
}
temp2.reinit(dim*dim, dim*dim); tracev(temp2, temp1, temp3);
for (unsigned int i = 0;i < dim*dim;i++) {
  for (unsigned int j = 0;j < dim*dim;j++) {
    PK1_Stiff_RVE[i][j] = PK1_Stiff_RVE[i][j] - temp2[i][j];
  }
}


dP_dF = 0.0;
FullMatrix<double> L(dim, dim);
L = IdentityMatrix(dim);
for (unsigned int m = 0;m < dim;m++) {
  for (unsigned int n = 0;n < dim;n++) {
    for (unsigned int o = 0;o < dim;o++) {
      for (unsigned int p = 0;p < dim;p++) {
        for (unsigned int i = 0;i < dim;i++) {
          for (unsigned int j = 0;j < dim;j++) {
            for (unsigned int k = 0;k < dim;k++) {
              for (unsigned int l = 0;l < dim;l++) {
                dP_dF[m][n][o][p] = dP_dF[m][n][o][p] + PK1_Stiff_RVE(dim*i + j, dim*k + l)*L(i, m)*L(j, n)*L(k, o)*L(l, p);;
              }
            }
          }
        }
      }
    }
  }
}







unsigned int j=0;
for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {
    slipfraction_iter[cellID][quadPtID][j*n_slip_systemsWOtwin+i] = slipfraction_conv[cellID][quadPtID][j*n_slip_systemsWOtwin + i] + tslipvfsys[j][i]*RegionsTwinvf[0];
}
for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {
  for (unsigned int j = 1;j < n_twin_systems_Size + 1;j++) {
    slipfraction_iter[cellID][quadPtID][j*n_slip_systemsWOtwin+i] = slipfraction_conv[cellID][quadPtID][j*n_slip_systemsWOtwin + i] + tslipvfsys[j][i]*ttwinvf[j-1];
  }
}

for (unsigned int i = 0;i < n_twin_systems_Size*2;i++) {
  TwinOutputfraction_iter[cellID][quadPtID][i] = TwinOutputfraction_conv[cellID][quadPtID][i] + tslipvfsys[0][n_slip_systemsWOtwin + i]*RegionsTwinvf[0];
}

for (unsigned int i = 0;i < 2;i++) {
  for (unsigned int j = 0;j < n_twin_systems_Size ;j++) {
    TwinOutputfraction_iter[cellID][quadPtID][n_twin_systems_Size*2+j*2+i] = TwinOutputfraction_conv[cellID][quadPtID][n_twin_systems_Size*2+j*2+i] + tslipvfsys[j+1][n_slip_systemsWOtwin+i]*ttwinvf[j];
  }
}


for (unsigned int i = 0;i < n_twin_systems_Size;i++) {
  ttwinvf[i] = ttwinvf[i] + (1 - Totaltwinvf)*(tslipvfsys[0][n_slip_systemsWOtwin + i] - tslipvfsys[0][n_slip_systemsWOtwin + n_twin_systems_Size + i]) / this->userInputs.twinShear1;
}

unsigned int numberOfTwinnedRegion = NumberOfTwinnedRegionK;
std::vector<unsigned int> ActiveTwinSystemsR = tActiveTwinSystems;

for (unsigned int i = 0;i < NumberOfTwinnedRegionK;i++) {
  unsigned int alpha = tActiveTwinSystems[i];
  ttwinvf[alpha - 1] = ttwinvf[alpha - 1] * (1 - (tslipvfsys[alpha][n_slip_systemsWOtwin] - tslipvfsys[alpha][n_slip_systemsWOtwin + 1]) / this->userInputs.twinShear1);
  if (ttwinvf[alpha - 1] < this->userInputs.twinThresholdFraction1) {
    if (numberOfTwinnedRegion > 1) {
      numberOfTwinnedRegion = numberOfTwinnedRegion - 1;
      ActiveTwinSystemsR.resize(numberOfTwinnedRegion);
      unsigned  int ii = 0;
      for (unsigned int j = 0;j < numberOfTwinnedRegion + 1;j++) {
        if (alpha == tActiveTwinSystems[j]) {
          continue;
        }
        ActiveTwinSystemsR[ii] = tActiveTwinSystems[j];
        ii = ii + 1;
      }
      tActiveTwinSystems.resize(numberOfTwinnedRegion);
      tActiveTwinSystems = ActiveTwinSystemsR;
    }

    else {
      numberOfTwinnedRegion = numberOfTwinnedRegion - 1;
      tActiveTwinSystems.resize(1);
      tActiveTwinSystems[0] = 0;
    }
    temp = IdentityMatrix(dim);
    //unsigned int alpha2 = alpha*dim;
    for (unsigned int k = 0;k < dim;k++) {
      for (unsigned int l = 0;l < dim;l++) {
        Fe_iter[cellID][quadPtID][k][l + alpha*dim] = temp[k][l];
        Fp_iter[cellID][quadPtID][k][l + alpha*dim] = temp[k][l];
      }
    }

    for (unsigned int k = 0;k < n_Tslip_systems;k++) {
      s_alpha_iter[cellID][quadPtID][alpha*n_Tslip_systems +i] = 0;
    }
    for (unsigned int i = 0;i < dim;i++) {
      rotnew_iter[cellID][quadPtID][alpha*dim+i] = rot[cellID][quadPtID][alpha*dim + i];
    }
    for (unsigned int i = 0;i < n_slip_systemsWOtwin;i++) {
      slipfraction_iter[cellID][quadPtID][alpha*n_slip_systemsWOtwin + i] = 0;
    }
    /*tslipvfsys[alpha] = 0;*/
    TwinFlag_iter[cellID][quadPtID][alpha - 1] = 0;

  }
}
NumberOfTwinnedRegionK = numberOfTwinnedRegion;

for (unsigned int i = 0;i < NumberOfNonTwinnedRegionK;i++) {
  if ((!this->userInputs.enableOneTwinSys_Reorien)||(NumberOfTwinnedRegionK==0)){
    unsigned int alpha = DeactiveTwinSystems[i];
    if (ttwinvf[alpha - 1] >= this->userInputs.twinThresholdFraction1) {
      if (NumberOfTwinnedRegionK > 0) {
        NumberOfTwinnedRegionK = NumberOfTwinnedRegionK + 1;
        if (NumberOfTwinnedRegionK>1){
          NumberOfTwinnedRegionK=NumberOfTwinnedRegionK;
        }
        ActiveTwinSystemsR.resize(NumberOfTwinnedRegionK);
        for (unsigned int j = 0;j < NumberOfTwinnedRegionK - 1;j++) {
          ActiveTwinSystemsR[j] = tActiveTwinSystems[j];
        }
        ActiveTwinSystemsR[NumberOfTwinnedRegionK - 1] = alpha;
      }
      else {
        NumberOfTwinnedRegionK = NumberOfTwinnedRegionK + 1;
        ActiveTwinSystemsR.resize(NumberOfTwinnedRegionK);
        ActiveTwinSystemsR[NumberOfTwinnedRegionK - 1] = alpha;
      }

      tActiveTwinSystems.resize(NumberOfTwinnedRegionK);
      tActiveTwinSystems = ActiveTwinSystemsR;

      //unsigned int alpha2 = alpha*dim;
      for (unsigned int k = 0;k < dim;k++) {
        for (unsigned int l = 0;l < dim;l++) {
          Fe_iter[cellID][quadPtID][k][l + alpha*dim] = Fe_iter[cellID][quadPtID][k][l];
          Fp_iter[cellID][quadPtID][k][l + alpha*dim] = Fp_iter[cellID][quadPtID][k][l];
        }
      }


      for (unsigned int k = 0;k < n_slip_systemsWOtwin;k++) {
        s_alpha_iter[cellID][quadPtID][n_Tslip_systems*alpha+k] = s_alpha_iter[cellID][quadPtID][k];
      }
      s_alpha_iter[cellID][quadPtID][alpha*n_Tslip_systems +n_slip_systemsWOtwin]= this->userInputs.initialSlipResistanceTwin1[n_twin_systems_Size * 2 - 1];
      s_alpha_iter[cellID][quadPtID][alpha*n_Tslip_systems +n_slip_systemsWOtwin+1]= this->userInputs.initialSlipResistanceTwin1[n_twin_systems_Size * 2 - 1];
      for (unsigned int k = 0;k < dim;k++) {
        rotnew_iter[cellID][quadPtID][alpha*dim + k] = rot[cellID][quadPtID][alpha*dim + k];
      }
      for (unsigned int k = 0;k < n_slip_systemsWOtwin;k++) {
        slipfraction_iter[cellID][quadPtID][alpha*n_slip_systemsWOtwin + k] = 0;
      }

      TwinFlag_iter[cellID][quadPtID][alpha - 1] = 1;


      if (this->userInputs.enableOneTwinSys_Reorien){
        for (unsigned int kk = 0;kk < n_twin_systems_Size;kk++) {
          if (kk!=(alpha - 1)){
            s_alpha_iter[cellID][quadPtID][kk+n_slip_systemsWOtwin] = s_alpha_iter[cellID][quadPtID][kk+n_slip_systemsWOtwin]*1000000000;
          }
        }
      }

    }
  }
}

NumberOfTwinnedRegion_iter[cellID][quadPtID]= NumberOfTwinnedRegionK;

twinvfNucleated.reinit(NumberOfTwinnedRegionK); RegionsTwinvf.reinit(1 + NumberOfTwinnedRegionK);
Totaltwinvf = 0;
twintemp = 0;
RegionsTwinvf = 0;
twinvfNucleated = 0;
for (unsigned int i = 0;i < NumberOfTwinnedRegionK;i++) {
  twintemp = ttwinvf[tActiveTwinSystems[i]-1];
  twinvfNucleated[i] = twintemp;
  Totaltwinvf = Totaltwinvf + twintemp;
  RegionsTwinvf[i + 1] = twintemp;
}
RegionsTwinvf[0] = 1 - Totaltwinvf;



///////The twin volume fraction should be updated first and then update TwinMaxFlag!!!!!!

TotaltwinvfK[cellID][quadPtID] = Totaltwinvf;

if (Totaltwinvf > this->userInputs.twinSaturationFactor1) {
  TwinMaxFlag_iter[cellID][quadPtID] = 0;
}
else{
  TwinMaxFlag_iter[cellID][quadPtID] = 1;
}
for (unsigned int i = 0;i < n_twin_systems_Size;i++) {
  ActiveTwinSystems_iter[cellID][quadPtID][i] = 0;
}
for (unsigned int i = 0;i < NumberOfTwinnedRegionK;i++) {
  ActiveTwinSystems_iter[cellID][quadPtID][i] = tActiveTwinSystems[i];
}

twinfraction_iter[cellID][quadPtID]= ttwinvf;




}

#include "../../../include/crystalPlasticity_template_instantiations.h"
