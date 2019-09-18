#include "../../../include/crystalPlasticity.h"
#include "userFunctions.cc"

template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
  unsigned int quadPtID,
  unsigned int StiffnessCalFlag)
  {

    multiphaseInit(cellID,quadPtID);
    F_tau=F; // Deformation Gradient
    FullMatrix<double> FE_t(dim,dim),FP_t(dim,dim);  //Elastic and Plastic deformation gradient
    Vector<double> s_alpha_t(n_Tslip_systems); // Slip resistance
    Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)

    // Tolerance

    double tol1=this->userInputs.modelStressTolerance;
    std::cout.precision(16);

    ////////////////////// The following parameters must be read in from the input file ///////////////////////////
    //double delgam_ref = 0.0001; // Reference slip increment
    //double strexp=1.0/50.0; // Strain rate sensitivity exponent ; the higher the less sensitive

    double delgam_ref = UserMatConstants(0); // Reference slip increment
    double strexp=UserMatConstants(1); // Strain rate sensitivity exponent ; the higher the less sensitive
    double sliptol = UserMatConstants(2);
    double tol2=UserMatConstants(3); // Slip system resistance tolerance for constitutive model loop
    double tol3=UserMatConstants(4); // Stress tensor tolerance
    // double corrfac=UserMatConstants(5) ; // Tolerance factor
    double tolstr=UserMatConstants(5); // Initial CRSS used in constitutive model to accept correction
    unsigned int nitr1=UserMatConstants(6),nitr2=UserMatConstants(7); // Maximum number of iterations for the Newton-Raphson scheme for the outer and inner loop
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    FE_t=Fe_conv[cellID][quadPtID];
    FP_t=Fp_conv[cellID][quadPtID];
    s_alpha_t=s_alpha_conv[cellID][quadPtID];
    rot1=rot_conv[cellID][quadPtID];

    // Rotation matrix of the crystal orientation
    FullMatrix<double> rotmat(dim,dim);
    rotmat=0.0;
    odfpoint(rotmat,rot1);

    FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim),temp4(dim,dim),temp5(dim,dim),temp6(dim,dim); // Temporary matrices
    FullMatrix<double> T_tau(dim,dim),P_tau(dim,dim);
    FullMatrix<double> FP_inv_t(dim,dim),FE_tau_trial(dim,dim),CE_tau_trial(dim,dim),Ee_tau_trial(dim,dim),FP_inv_tau(dim,dim),F_inv_tau(dim,dim);
    FullMatrix<double> SCHMID_TENSOR1(n_Tslip_systems*dim,dim),B(n_Tslip_systems*dim,dim),C(n_Tslip_systems*dim,dim);
    Vector<double> m1(dim),n1(dim);


    // Elastic Modulus

    FullMatrix<double> Dmat(2*dim,2*dim),Dmat2(2*dim,2*dim),TM(dim*dim,dim*dim);
    Vector<double> vec1(2*dim),vec2(dim*dim);

    elasticmoduli(Dmat2, rotmat, elasticStiffnessMatrix);

    //Elastic Stiffness Matrix Dmat
    Dmat.reinit(6,6) ; Dmat = 0.0;

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

    vec2(0)=0;vec2(1)=5;vec2(2)=4;vec2(3)=5;vec2(4)=1;vec2(5)=3;vec2(6)=4;vec2(7)=3;vec2(8)=2;



    for(unsigned int i=0;i<9;i++){
      for(unsigned int j=0;j<9;j++){
        TM[i][j]=Dmat2(vec2(i),vec2(j));
      }
    }

    Vector<double> s_alpha_tau;
    Vector<double> h_beta(n_Tslip_systems),h0(n_Tslip_systems),a_pow(n_Tslip_systems),s_s(n_Tslip_systems);
    FullMatrix<double> h_alpha_beta_t(n_Tslip_systems,n_Tslip_systems);
    Vector<double> resolved_shear_tau(n_Tslip_systems);

    FullMatrix<double> PK1_Stiff(dim*dim,dim*dim);
    FullMatrix<double> T_star_tau(dim,dim),T_star_tau_trial(dim,dim),mtemp(dim,dim);
    Vector<double> vtemp(6),tempv1(6);
    double det_FE_tau, det_F_tau, det_FP_tau;


    FP_inv_t = 0.0; FP_inv_t.invert(FP_t);
    FE_tau_trial=0.0;
    F_tau.mmult(FE_tau_trial,FP_inv_t) ;

    // CE_tau_trial is the same as A matrix - Kalidindi's thesis
    CE_tau_trial=0.0;
    FE_tau_trial.Tmmult(CE_tau_trial,FE_tau_trial);

    temp.reinit(dim,dim);
    Ee_tau_trial=CE_tau_trial;
    temp=IdentityMatrix(dim);
    for(unsigned int i=0;i<dim;i++){
      for(unsigned int j=0;j<dim;j++){
        Ee_tau_trial[i][j] = 0.5*(Ee_tau_trial[i][j]-temp[i][j]); // Compute the trial elastic Green-Lagrange strain tensor
      }
    }

    // Calculate the trial stress T_star_tau_trial

    tempv1=0.0;
    Dmat.vmult(tempv1, vecform(Ee_tau_trial));
    matform(T_star_tau_trial,tempv1);

    // Loop over slip systems to construct relevant matrices - Includes both slip and twin(considered as pseudo-slip) systems
    for (unsigned int i = 0;i<n_Tslip_systems;i++) {

      for (unsigned int j = 0;j<dim;j++) {
        m1(j) = m_alpha[i][j];
        n1(j) = n_alpha[i][j];
      }

      temp = 0.0;
      temp2 = 0.0;
      for (unsigned int j = 0;j<dim;j++) {
        for (unsigned int k = 0;k<dim;k++) {
          temp[j][k] = m1(j)*n1(k);         // Construct Schmid tensor matrix relative to the crystal reference frame
        }
      }
      // Transform the Schmid tensor matrix to sample coordinates
      rotmat.mmult(temp2, temp);
      temp2.mTmult(temp, rotmat);

      for (unsigned int j = 0;j<dim;j++) {
        for (unsigned int k = 0;k<dim;k++) {
          SCHMID_TENSOR1[dim*i + j][k] = temp[j][k]; // Store rotated Schmid tensor matrix
        }
      }

      // Construct B matrix - Kalidindi's thesis
      CE_tau_trial.mmult(temp2, temp);
      temp2.symmetrize();
      vtemp=0.0;
      mtemp=0.0;
      temp2.equ(2.0,temp2);
      // Construct C matrix - Kalidindi's thesis
      Dmat.vmult(vtemp, vecform(temp2));
      matform(mtemp,vtemp);
      mtemp.equ(0.5,mtemp);

      for (unsigned int j = 0;j<dim;j++) {
        for (unsigned int k = 0;k<dim;k++) {
          B[dim*i + j][k] = temp2[j][k];
          C[dim*i + j][k] = mtemp[j][k];
        }
      }
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////Start Nonlinear iteration for Slip increments////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    unsigned int itr1=0,itr2=0;
    double dffhrdn=1.0,dffstr=1.0,sgnm,dffslip=1.0;
    double sctmp1,sctmp2;
    FullMatrix<double> G_iter(dim,dim);
    FullMatrix<double> btemp1(dim*dim,dim*dim),J_iter(dim*dim,dim*dim),J_iter_inv(dim*dim,dim*dim);
    FullMatrix<double> T_star_iter(dim,dim),T_star_iterp(dim,dim);
    FullMatrix<double> CE_t(dim,dim);
    Vector<double> s_alpha_it(n_Tslip_systems),s_alpha_iterp(n_Tslip_systems),delgam_tau(n_Tslip_systems),delgam_tau_iter(n_Tslip_systems),diffgam(n_Tslip_systems);
    Vector<double> vtmp1(2*dim),vtmp2(2*dim),vtmp3(2*dim),vtmp4(n_Tslip_systems);
    Vector<double> nv1(dim*dim),nv2(dim*dim);


    FE_t.Tmmult(CE_t,FE_t);
    T_star_iter.equ(1.0,T_star_tau_trial);
    s_alpha_it=s_alpha_t;

    for(unsigned int i=0; i<n_Tslip_systems;i++){
      delgam_tau(i) = 1 ;
    }


    // Loop to check for the difference in CRSS in subsequent Newton-Raphson iterations
    while(dffhrdn>tol2 && dffslip>sliptol && itr1<nitr1){
      // Iterant 1
      itr1 = itr1+1 ;
      delgam_tau_iter.equ(1.0,delgam_tau);

      // Loop to check the difference in stress components in subsequent Newton-Raphson iterations
      while(dffstr>tol3 && itr2<nitr2){
        // Iterant 2
        itr2 = itr2+1;

        // Residual for the non-linear algebraic equation
        G_iter=0.0;
        // Jacobian for the Newton-Raphson iteration
        J_iter=IdentityMatrix(dim*dim);

        // Loop over slip systems
        for (unsigned int i = 0;i<n_Tslip_systems;i++){

          for (unsigned int j = 0;j < dim;j++) {
            for (unsigned int k = 0;k < dim;k++) {
              temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
              temp4[j][k]=C[dim*i + j][k];
            }
          }

          T_star_iter.mTmult(temp2,temp);
          sctmp1=temp2.trace();


          vecform9(nv1,temp4);
          vecform9(nv2,temp);

          if(i<n_slip_systems){ // For slip systems due to symmetry of slip
            if(sctmp1<0)
            sgnm=-1;
            else
            sgnm=1 ;
          }
          else               // For twin systems due to asymmetry of slip
          {
            if(sctmp1<=0)
            sgnm=0;
            else
            sgnm=1 ;
          }
          delgam_tau(i) = delgam_ref*pow(fabs(sctmp1/s_alpha_it(i)),1.0/strexp)*sgnm ;

          // Modification to the residual inside loop
          G_iter.add(delgam_tau(i),temp4);

          for (unsigned int j = 0;j < dim*dim;j++) {
            for (unsigned int k = 0;k < dim*dim;k++) {
              btemp1[j][k]=nv1(j)*nv2(k);
            }
          }

          if(i<n_slip_systems){
            sgnm = 1 ;
          }
          else{
            if(sctmp1<=0)
            sgnm=0;
            else
            sgnm=1 ;
          }
          sctmp2=delgam_ref/(strexp*s_alpha_it(i))*pow(fabs(sctmp1/s_alpha_it(i)),(1.0/strexp - 1.0))*sgnm;
          btemp1.equ(sctmp2,btemp1);
          // Modification to the Jacobian of the Newton-Raphson iteration
          J_iter.add(1.0,btemp1);
        }
        // Final modification to the residual outside the loop
        G_iter.add(1.0,T_star_iter,-1.0,T_star_tau_trial);




        // Invert Jacobian
        J_iter_inv.invert(J_iter);
        vecform9(nv1,G_iter);
        J_iter_inv.vmult(nv2,nv1);
        matform9(temp6,nv2);
        temp6.equ(-1.0,temp6);
        T_star_iterp=0.0;
        T_star_iterp.add(1.0,T_star_iter,1.0,temp6);

        // Criteria to accept or modify the Newton correction
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            if(temp6[j][k]<0)
            sgnm=-1;
            else
            sgnm=1 ;

            if(fabs(temp6[j][k])>tolstr){
              T_star_iterp[j][k]=T_star_iter[j][k]+tolstr*sgnm;
            }
          }


        }


        vtmp3=0.0;
        temp6=0.0;
        temp6.add(1.0,T_star_iterp,-1.0,T_star_iter);
        vtmp3.equ(1.0,vecform(temp6));
        dffstr=vtmp3.l2_norm();
        T_star_iter.equ(1.0,T_star_iterp);

      } // inner while

      // Single slip hardening rate
      for(unsigned int i=0;i<n_slip_systems;i++){
        h_beta(i)=initialHardeningModulus[i]*pow((1-s_alpha_it(i)/saturationStress[i]),powerLawExponent[i]);
      }


      for(unsigned int i=0;i<n_twin_systems;i++){
        h_beta(n_slip_systems+i)=initialHardeningModulusTwin[i]*pow((1-s_alpha_it(n_slip_systems+i)/saturationStressTwin[i]),powerLawExponentTwin[i]);
      }


      for(unsigned int i=0;i<n_Tslip_systems;i++){
        for(unsigned int j=0;j<n_Tslip_systems;j++){
          h_alpha_beta_t[i][j] = q[i][j]*h_beta(j);
        }
      }

      s_alpha_iterp=s_alpha_t;

      temp=0.0;
      temp1=0.0;
      temp2=0.0;
      temp3=0.0;
      temp4=0.0;
      temp5=0.0;
      temp6=0.0;
      sctmp1=0.0;
      vtmp1=0.0;
      vtmp2=0.0;

      for (unsigned int i = 0;i<n_Tslip_systems;i++){
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
          }
        }
        T_star_iter.mTmult(temp1,temp);
        sctmp1=temp1.trace();
        resolved_shear_tau(i)=sctmp1;

        if(i<n_slip_systems){ // For slip systems due to symmetry of slip
          if(sctmp1<0)
          sgnm=-1;
          else
          sgnm=1 ;
        }
        else               // For twin systems due to asymmetry of slip
        {
          if(sctmp1<=0)
          sgnm=0;
          else
          sgnm=1 ;
        }


        delgam_tau(i)=delgam_ref*pow(fabs(resolved_shear_tau(i)/s_alpha_it(i)),(1.0/strexp))*sgnm;

        for (unsigned int j = 0;j<n_Tslip_systems;j++){
          s_alpha_iterp(j)=s_alpha_iterp(j)+h_alpha_beta_t[j][i]*fabs(delgam_tau(i));
        }

        // Check if the slip system resistances exceed their corresponding saturation stress. If yes, set them equal to the saturation stress

        for(unsigned int i=0;i<n_slip_systems;i++){
          if(s_alpha_iterp[i] >= saturationStress[i])
          s_alpha_iterp[i] =  saturationStress[i] ;
        }


        for(unsigned int i=0;i<n_twin_systems;i++){
          if(s_alpha_iterp[n_slip_systems+i] >= saturationStressTwin[i])
          s_alpha_iterp[n_slip_systems+i] =  saturationStressTwin[i] ;
        }


      }

      diffgam.equ(1.0,delgam_tau);
      diffgam.add(-1.0,delgam_tau_iter);
      dffslip = diffgam.l2_norm();

      vtmp4=0.0;
      vtmp4.add(1.0,s_alpha_iterp,-1.0,s_alpha_it);
      dffhrdn=vtmp4.l2_norm();
      s_alpha_it=s_alpha_iterp;
    } // outer while



    ////////////////////////////////////End Nonlinear iteration for Slip increments////////////////////////////////////

    s_alpha_tau=s_alpha_it;
    FP_tau=IdentityMatrix(dim);
    temp=0.0;
    for (unsigned int i = 0;i<n_Tslip_systems;i++){
      for (unsigned int j = 0;j < dim;j++) {
        for (unsigned int k = 0;k < dim;k++) {
          temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
        }
      }

      if(i<n_slip_systems){ // For slip systems due to symmetry of slip
        if(resolved_shear_tau(i)<0)
        sgnm=-1;
        else
        sgnm=1 ;
      }
      else               // For twin systems due to asymmetry of slip
      {
        if(resolved_shear_tau(i)<=0)
        sgnm=0;
        else
        sgnm=1 ;
      }

      delgam_tau(i)=delgam_ref*pow(fabs(resolved_shear_tau(i)/s_alpha_tau(i)),(1.0/strexp))*sgnm;
      FP_tau.add(-1*delgam_tau(i),temp);
    }
    temp=0.0;
    temp.invert(FP_tau);
    temp.mmult(FP_tau,FP_t);

    det_FP_tau=FP_tau.determinant();
    FP_tau.equ(pow(det_FP_tau,-1.0/3),FP_tau);

    FP_inv_tau = 0.0; FP_inv_tau.invert(FP_tau);
    FE_tau = 0.0;
    F_tau.mmult(FE_tau, FP_inv_tau);
    temp.reinit(dim, dim);
    det_FE_tau = FE_tau.determinant();
    T_star_tau.equ(1.0,T_star_iter);
    FE_tau.mmult(temp, T_star_tau); temp.equ(1.0 / det_FE_tau, temp);
    temp.mTmult(T_tau, FE_tau);

    det_F_tau = F_tau.determinant();
    temp.invert(F_tau);
    F_inv_tau.equ(1.0,temp);
    T_tau.mTmult(P_tau, temp);
    P_tau.equ(det_F_tau, P_tau);

    for (unsigned int i=0;i<n_twin_systems;i++){
      twinfraction_iter[cellID][quadPtID][i]=twinfraction_conv[cellID][quadPtID][i]+delgam_tau[i+n_slip_systems]/twinShear;
    }

    for (unsigned int i=0;i<n_slip_systems;i++){
      slipfraction_iter[cellID][quadPtID][i]=slipfraction_conv[cellID][quadPtID][i]+fabs(delgam_tau[i]);
    }


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// Computing Algorithmic Tangent Modulus ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    FullMatrix<double>  cnt1(dim*dim,dim*dim),cnt2(dim*dim,dim*dim),cnt3(dim*dim,dim*dim),cnt4(dim*dim,dim*dim); // Variables to track individual contributions
    FullMatrix<double> dFedF(dim*dim,dim*dim); // Meaningful variables
    FullMatrix<double>  ntemp1(dim*dim,dim*dim),ntemp2(dim*dim,dim*dim),ntemp3(dim*dim,dim*dim),ntemp4(dim*dim,dim*dim),ntemp5(dim*dim,dim*dim),ntemp6(dim*dim,dim*dim); // Temporary variables
    double mulfac;

    if (StiffnessCalFlag==1){

      // Contribution 1 - Most straightforward because no need to invoke constitutive model
      FE_tau.mmult(temp,T_star_tau);
      temp.mTmult(temp1,FE_tau);
      temp1.mTmult(temp2,F_inv_tau);
      left(ntemp1,temp2);
      left(ntemp2,F_inv_tau);
      trpose(ntemp3,ntemp2);
      ntemp1.mmult(cnt4,ntemp3);
      cnt4.equ(-1.0,cnt4);

      // Compute dFedF
      ntemp4=0.0;
      for (unsigned int i = 0;i<n_Tslip_systems;i++){
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
          }
        }


        traceval(ntemp1,temp);

        temp1=0.0;
        temp1.add(0.5,temp);
        temp1.Tadd(0.5,temp);

        Dmat.vmult(vtmp1,vecform(temp1));
        matform(temp1,vtmp1);
        temp1.mTmult(temp2,FE_tau);
        left(ntemp2,temp2);
        ntemp1.mmult(ntemp3,ntemp2);


        if(i<n_slip_systems){
          sgnm = 1 ;
        }
        else{
          if(sctmp1<=0)
          sgnm=0;
          else
          sgnm=1 ;
        }

        mulfac = delgam_ref*1.0/strexp*1.0/s_alpha_tau(i)*pow(fabs(resolved_shear_tau(i)/s_alpha_tau(i)),1.0/strexp - 1.0)*sgnm;
        ntemp4.add(mulfac,ntemp3);
      }

      left(ntemp1,FE_tau_trial);
      ntemp1.mmult(ntemp2,ntemp4);
      temp1=IdentityMatrix(dim);
      left(ntemp3,temp1);
      ntemp2.add(1.0,ntemp3);
      right(ntemp3,FP_inv_tau);
      ntemp1.invert(ntemp2);
      ntemp1.mmult(dFedF,ntemp3);


      // Compute remaining contributions which depend solely on dFedF

      // Contribution 1
      T_star_tau.mTmult(temp1,FE_tau);
      temp1.mmult(temp2,F_inv_tau);
      right(ntemp1,temp2);
      ntemp1.mmult(cnt1,dFedF);

      // Contribution 2
      temp=0.0;
      temp.Tadd(1.0,FE_tau);
      left(ntemp1,temp);
      TM.mmult(ntemp2,ntemp1);
      ntemp1.equ(0.5,ntemp2);

      left(ntemp2,temp);
      trpose(ntemp3,ntemp2);
      TM.mmult(ntemp2,ntemp3);
      ntemp2.equ(0.5,ntemp2);

      ntemp3=0.0;
      ntemp3.add(1.0,ntemp1,1.0,ntemp2);
      ntemp3.mmult(ntemp4,dFedF);

      left(ntemp1,FE_tau);

      temp1=0.0;
      temp1.Tadd(1.0,FE_tau);
      temp1.mTmult(temp2,F_inv_tau);
      right(ntemp2,temp2);

      ntemp1.mmult(ntemp3,ntemp4);
      ntemp3.mmult(cnt2,ntemp2);

      // Contribution 3
      FE_tau.mmult(temp1,T_star_tau);
      left(ntemp1,temp1);

      left(ntemp2,F_inv_tau);
      ntemp2.mmult(ntemp3,dFedF);
      trpose(ntemp4,ntemp3);
      ntemp1.mmult(cnt3,ntemp4);


      // Assemble contributions to PK1_Stiff

      PK1_Stiff=0.0;
      PK1_Stiff.add(1.0,cnt1);
      PK1_Stiff.add(1.0,cnt2);
      PK1_Stiff.add(1.0,cnt3);
      PK1_Stiff.add(1.0,cnt4);

      ////////////////// End Computation ////////////////////////////////////////


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
                      dP_dF[m][n][o][p] = dP_dF[m][n][o][p] + PK1_Stiff(dim*i + j, dim*k + l)*L(i, m)*L(j, n)*L(k, o)*L(l, p);
                    }
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


    sres_tau.reinit(n_Tslip_systems);
    sres_tau = s_alpha_tau;

    // Update the history variables
    Fe_iter[cellID][quadPtID]=FE_tau;
    Fp_iter[cellID][quadPtID]=FP_tau;
    s_alpha_iter[cellID][quadPtID]=sres_tau;


    /////// EXTRA STUFF FOR REORIENTATION POST TWINNING ////////////////////

    std::vector<double> local_twin;
    std::vector<double>::iterator result;
    double twin_pos, twin_max;
    Vector<double> quat1(4), rod(3), quat2(4), quatprod(4);

    if (enableTwinning){
      if (!this->userInputs.enableMultiphase){
        if (F_r > 0) {
          F_T = twinThresholdFraction + (twinSaturationFactor*F_e / F_r);
        }
        else {
          F_T = twinThresholdFraction;
        }
      }
      else{
        F_T = twinThresholdFraction;
      }

      //////Eq. (13) in International Journal of Plasticity 65 (2015) 61–84
      if (F_T > 1.0) {
        F_T = 1.0;
      }
      local_twin.resize(n_twin_systems,0.0);
      local_twin=twinfraction_iter[cellID][quadPtID];
      result = std::max_element(local_twin.begin(), local_twin.end());
      twin_pos= std::distance(local_twin.begin(), result);
      twin_max=local_twin[twin_pos];
      if (twin_conv[cellID][quadPtID] != 1.0) {
        if(F_r>0){
          if(twin_max > F_T){

            rod(0) = rot_conv[cellID][quadPtID][0];rod(1) = rot_conv[cellID][quadPtID][1];rod(2) = rot_conv[cellID][quadPtID][2];
            odfpoint(rotmat, rod);
            rod2quat(quat2, rod);
            quat1(0) = 0;
            quat1(1) = n_alpha[n_slip_systems + twin_pos][0];
            quat1(2) = n_alpha[n_slip_systems + twin_pos][1];
            quat1(3) = n_alpha[n_slip_systems + twin_pos][2];

            quatproduct(quatprod, quat2, quat1);


            quat2rod(quatprod, rod);

            odfpoint(rotmat, rod);

            rot_iter[cellID][quadPtID][0] = rod(0);rot_iter[cellID][quadPtID][1] = rod(1);rot_iter[cellID][quadPtID][2] = rod(2);
            rotnew_iter[cellID][quadPtID][0] = rod(0);rotnew_iter[cellID][quadPtID][1] = rod(1);rotnew_iter[cellID][quadPtID][2] = rod(2);
            twin_iter[cellID][quadPtID] = 1.0;
            for (unsigned int i = 0;i < n_twin_systems;i++) {
              s_alpha_iter[cellID][quadPtID][n_slip_systems + i] =100000;
            }
          }
        }
      }
    }


  }
  #include "../../../include/crystalPlasticity_template_instantiations.h"
