#include "../../../include/crystalPlasticity.h"
#include "userFunctions.cc"
//////////////////////////////////////////////////////////////////////////
//////calculatePlasticity.cc numerically integrates the constitive model.
//////This calculatePlasticity.cc is based on the following rate-dependent crystal plasticity model:
//////SR Kalidindi, Polycrystal plasticity: constitutive modeling and deformation processing,
////// PhD thesis, MIT, 1992. In addition to isotropic hardening tehre is also kinematic hardening
////// where the backstress evolves based on the OW hardening model. The guess stress to start the
////// Newton-Raphson iteration is the previously converged stress. Exponential update of Fp is
//// implemented along with a line-search algorithm to solve the nonlinear system. The
//// tangent modulus computation is also more involved.
//////To use this file, one should copy this calculatePlasticity.cc into the following folder:
//////    plasticity/src/materialModels/crystalPlasticity/
////// (It should replace the original calculatePlasticity.cc inside that folder)
//////Next, they should copy the file userFunctions.cc into the following folder:
//////    plasticity/src/materialModels/crystalPlasticity/
//////
//////Finally, the PRISMS-Plasticity should be recompiled.
//////////////////////////////////////////////////////////////////////////////
////
//
template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
  unsigned int quadPtID,
  unsigned int StiffnessCalFlag)
  {

    multiphaseInit(cellID,quadPtID);

    /////////////////////////////////////////////////////////
    FullMatrix<double> FE_t(dim,dim),FP_t(dim,dim),F_t(dim,dim);  //Elastic, Plastic, and total deformation gradient
    Vector<double> s_alpha_t(n_Tslip_systems),slipfraction_t(n_slip_systems),twinfraction_t(n_twin_systems); // Slip resistance
    Vector<double> W_kh_t(n_Tslip_systems),W_kh_t1(n_Tslip_systems),W_kh_t2(n_Tslip_systems),signed_slip_t(n_Tslip_systems) ;
    Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)
    FullMatrix<double> Tinter_diff_guess(dim,dim), Ep_t (dim,dim), iMinusL (dim,dim);
    double Ep_eff_cum_t;
    double tol1, delgam_ref,strexp,locres_tol,locres_tol2,h_back1,h_back2,r_back1,r_back2,m_back1,m_back2,back_lim_1,back_lim_2,b_back1,b_back2;
    unsigned int nitr1,nitr2;
    unsigned int ii;
    FullMatrix<double> rotmat(dim,dim); // Rotation matrix of the crystal orientation
    FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim),temp4(dim,dim) ; // Temporary matrices
    FullMatrix<double> temp5(dim,dim),temp6(dim,dim), temp7(dim,dim), temp8(dim,dim); // Temporary matrices
    FullMatrix<double> T_tau(dim,dim),P_tau(dim,dim); // Stress measures for current timestep
    FullMatrix<double> FP_inv_t(dim,dim),FE_tau_trial(dim,dim),CE_tau_trial(dim,dim) ; // Kinematic descriptors
    FullMatrix<double>FP_inv_tau(dim,dim),F_inv_tau(dim,dim); // Kinematic descriptors
    FullMatrix<double> SCHMID_TENSOR1(n_Tslip_systems*dim,dim),Normal_SCHMID_TENSOR1(n_Tslip_systems*dim,dim); // Projection matrices
    Vector<double> m1(dim),n1(dim);

    // Elastic Modulus
    FullMatrix<double> Dmat(2*dim,2*dim),Dmat2(2*dim,2*dim),TM(dim*dim,dim*dim);
    Vector<double> vec2(dim*dim);
    Vector<double> s_alpha_tau(n_Tslip_systems),slipfraction_tau(n_slip_systems),twinfraction_tau(n_twin_systems) ;
    Vector<double> W_kh_tau(n_Tslip_systems),W_kh_tau1(n_Tslip_systems),W_kh_tau2(n_Tslip_systems); // Converged backstresses
    Vector<double> W_kh_tau_it(n_Tslip_systems),W_kh_tau1_it(n_Tslip_systems),W_kh_tau2_it(n_Tslip_systems); // Iterative backstresses
    Vector<double> h_beta(n_Tslip_systems),h0(n_Tslip_systems),a_pow(n_Tslip_systems),s_s(n_Tslip_systems); // Isotropic hardening parameters
    FullMatrix<double> h_alpha_beta_t(n_Tslip_systems,n_Tslip_systems),Ep_tau(dim,dim),del_Ep_tau(dim,dim);
    Vector<double> resolved_shear_tau(n_Tslip_systems), normal_stress_tau(n_Tslip_systems),signed_slip_tau(n_Tslip_systems);
    double Ep_eff_cum_tau,del_Ep_eff_cum_tau;
    FullMatrix<double> PK1_Stiff(dim*dim,dim*dim); // Tangent modulus
    FullMatrix<double> T_star_tau(dim,dim),mtemp(dim,dim);
    Vector<double> vtemp(2*dim);
    double det_FE_tau, det_F_tau, trny_op, trny_in;
    double m1_norm, n1_norm ;

    unsigned int itr1, itr2;
    double sctmp1, sctmp2, sctmp3, sctmp4; // Scalars used as temporary variables
    Vector<double> G_iter(2*dim),locres_vec(2*dim+2*n_Tslip_systems),stateVar_it(2*dim+2*n_Tslip_systems),stateVar_temp(2*dim+2*n_Tslip_systems),stateVar_diff(2*dim+2*n_Tslip_systems),T_star_iter_vec(2*dim);
    Vector<double> gradFold(2*dim+2*n_Tslip_systems) ;
    FullMatrix<double> btemp1(2*dim,2*dim),J_iter(2*dim+2*n_Tslip_systems,2*dim+2*n_Tslip_systems),J_iter_inv(2*dim+2*n_Tslip_systems,2*dim+2*n_Tslip_systems),LP_acc(dim,dim),LP_acc2(dim,dim) ;
    FullMatrix<double> J_iter_cp(2*dim+2*n_Tslip_systems,2*dim+2*n_Tslip_systems);
    FullMatrix<double> T_star_iter(dim,dim),T_star_iterp(dim,dim);
    FullMatrix<double> CE_t(dim,dim);
    Vector<double> s_alpha_it(n_Tslip_systems),s_alpha_iterp(n_Tslip_systems),delgam_tau(n_Tslip_systems),delgam_tau_iter(n_Tslip_systems);

    Vector<double> vtmp1(2*dim),vtmp2(2*dim),vtmp3(2*dim);
    Vector<double> vtmp4(n_Tslip_systems);
    Vector<double> nv1(2*dim),nv2(2*dim);
    double locres, locres2;
    double sgnm, sgnm2 , sgnm3 , sgnm4, Fold ;

    FullMatrix<double>  cnt1(dim*dim,dim*dim),cnt2(dim*dim,dim*dim),cnt3(dim*dim,dim*dim),cnt4(dim*dim,dim*dim); // Variables to track individual contributions
    FullMatrix<double> dFedF(dim*dim,dim*dim),dFpdFe(dim*dim,dim*dim); // Meaningful variables
    FullMatrix<double>  ntemp1(dim*dim,dim*dim),ntemp2(dim*dim,dim*dim),ntemp3(dim*dim,dim*dim),ntemp4(dim*dim,dim*dim),ntemp5(dim*dim,dim*dim),ntemp6(dim*dim,dim*dim); // Temporary variables
    FullMatrix<double> Pmat(n_Tslip_systems,n_Tslip_systems),Qmat(n_Tslip_systems,n_Tslip_systems),Pmat_inv(n_Tslip_systems,n_Tslip_systems),Rmat(n_Tslip_systems,n_Tslip_systems) ;
    FullMatrix<double> Smat(n_Tslip_systems,n_Tslip_systems) ;
    FullMatrix<double> Qmat2(n_Tslip_systems,dim*dim),dgamdFe(n_Tslip_systems,dim*dim) ;
    Vector<double> nvec1(dim*dim),nvec2(dim*dim) ;

    double modifier1_num, modifier1_den, modifier2_num, modifier2_den, modifier;

    double mulfac;
    FullMatrix<double> L(dim, dim);

    std::vector<double> local_twin;
    std::vector<double>::iterator result;
    double twin_pos, twin_max;
    Vector<double> quat1(4), rod(3), quat2(4), quatprod(4);

    double Criteria_Delta_F,inverseNumberOfCuts,Criteria_Delta_F_2;
    unsigned int  numberOfCuts;
    FullMatrix<double> Delta_F,div_Delta_F;


    //////////////////////////////////////////////////////////


    std::cout.precision(16);
    F_tau=F; // Deformation Gradient

    tol1=this->userInputs.modelStressTolerance;
    delgam_ref = UserMatConstants(0)*this->userInputs.delT ; // Reference slip increment
    strexp=UserMatConstants(1); // Strain rate sensitivity exponent ; the higher the less sensitive
    locres_tol=UserMatConstants(2); // Tolerance for innter residual of nonlinear constitutive model
    locres_tol2=UserMatConstants(3); // Tolerance for outer slip resistance loop
    nitr1=UserMatConstants(4); // Maximum number of iterations for inner loop convergence
    nitr2=UserMatConstants(5); // Maximum number of iterations for outer loop convergence
    ///////////////////////////////Backstress evolution parameters/////////
    h_back1=UserMatConstants(6);h_back2=UserMatConstants(7);
    r_back1=UserMatConstants(8);r_back2=UserMatConstants(9);
    m_back1=UserMatConstants(10);m_back2=UserMatConstants(11);



    // To deal with the scenario when there might be no backstress evolution
    if(fabs(r_back1)>1.0e-7)
    b_back1 = h_back1/r_back1 ;
    else
    b_back1 = 1 ;

    if(fabs(r_back2)>1.0e-7)
    b_back2 = h_back2/r_back2 ;
    else
    b_back2 = 1 ;

    // Limiting value of backstress
    back_lim_1 = b_back1 ;
    back_lim_2 = b_back2 ;
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////Reading the previous converged state Variables////////////////
    FE_t=Fe_conv[cellID][quadPtID];
    FP_t=Fp_conv[cellID][quadPtID];
    rot1=rot_conv[cellID][quadPtID];
    Tinter_diff_guess = TinterStress_diff[cellID][quadPtID] ; // Stress increment from previous increment

    for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
      W_kh_t1(i) = stateVar_conv[cellID][quadPtID][i];
      W_kh_t2(i) = stateVar_conv[cellID][quadPtID][i+n_Tslip_systems];
      signed_slip_t(i)=stateVar_conv[cellID][quadPtID][i+2*n_Tslip_systems];
      s_alpha_t(i)=s_alpha_conv[cellID][quadPtID][i];
    }

    ii=0;
    for(unsigned int i=0 ; i<dim ; i++){
      for(unsigned int j=0 ; j<dim ; j++){
        Ep_t[i][j]=stateVar_conv[cellID][quadPtID][ii+4*n_Tslip_systems];
        ii=ii+1;
      }
    }

    Ep_eff_cum_t=stateVar_conv[cellID][quadPtID][4*n_Tslip_systems+dim*dim];

    for(unsigned int i=0 ; i<n_slip_systems ; i++)
    slipfraction_t(i) = slipfraction_conv[cellID][quadPtID][i] ;

    for(unsigned int i=0 ; i<n_twin_systems ; i++)
    twinfraction_t(i) = twinfraction_conv[cellID][quadPtID][i] ;
    /////////////////////////////////////////////////////////////////////////////

    ////////////////////////Rotating the Elsatic modulus from crystal to sample//////////////////
    rotmat=0.0;
    odfpoint(rotmat,rot1);
    elasticmoduli(Dmat2, rotmat, elasticStiffnessMatrix);
    /////////////////////////////////////////////////////////////////////////////////////////////

    ///////////////////////////Calculation of elastic stiffness matrix in different shapes required in the implementation//////////////
    //Elastic Stiffness Matrix Dmat
    Dmat.reinit(6,6) ;
    Dmat = 0.0;
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
    //////////////////////////////////////////////////////////////////

    /////////////Building the Schmid Tensor//////////////////////////
    // Loop over slip systems to construct relevant matrices - Includes both slip and twin(considered as pseudo-slip) systems
    for (unsigned int i = 0;i<n_Tslip_systems;i++) {

      for (unsigned int j = 0;j<dim;j++) {
        m1(j) = m_alpha[i][j];
        n1(j) = n_alpha[i][j];
      }
      temp = 0.0;
      temp2 = 0.0;
      temp3 = 0.0;
      temp4 = 0.0 ;
      for (unsigned int j = 0;j<dim;j++) {
        for (unsigned int k = 0;k<dim;k++) {
          temp[j][k] = m1(j)*n1(k);         // Construct Schmid tensor matrix relative to the crystal reference frame
          temp3[j][k] = n1(j)*n1(k);
        }
      }
      // Transform the Schmid tensor matrix to sample coordinates
      rotmat.mmult(temp2, temp);
      temp2.mTmult(temp, rotmat);

      rotmat.mmult(temp4, temp3);
      temp4.mTmult(temp3, rotmat);

      for (unsigned int j = 0;j<dim;j++) {
        for (unsigned int k = 0;k<dim;k++) {
          SCHMID_TENSOR1[dim*i + j][k] = temp[j][k]; // Store rotated Schmid tensor matrix
          Normal_SCHMID_TENSOR1[dim*i + j][k] = temp3[j][k];
        }
      }

    }
    //////////////////////////////////////////////////////

    ////////////////////////Initialization of variables/////////////
    W_kh_tau1 = W_kh_t1 ;
    W_kh_tau2 = W_kh_t2 ;
    W_kh_tau = 0.0 ;
    W_kh_tau.add(1.0,W_kh_t1,1.0,W_kh_t2) ;
    W_kh_tau_it = W_kh_tau ;

    slipfraction_tau = slipfraction_t ;
    signed_slip_tau  = signed_slip_t;
    Ep_tau=Ep_t;
    Ep_eff_cum_tau=Ep_eff_cum_t;
    twinfraction_tau = twinfraction_t ;
    //////////////////////////////////////////////////////////////////

    /////////////Dividing the applied DeltaF to small increments///////
    FE_t.mmult(F_t,FP_t) ;

    Delta_F.reinit(dim,dim);div_Delta_F.reinit(dim,dim);
    Delta_F=0;
    div_Delta_F=0;

    if (!this->userInputs.flagTaylorModel){
      Delta_F.add(-1,F_t,1.0,F_tau) ;
      if (this->userInputs.numberTaylorSubsteps==1){
        Criteria_Delta_F=fabs(Delta_F[0][0])+fabs(Delta_F[1][1])+fabs(Delta_F[2][2])+2*fabs(Delta_F[0][1])+2*fabs(Delta_F[0][2])+2*fabs(Delta_F[1][2])+2*fabs(Delta_F[1][0])+2*fabs(Delta_F[2][0])+2*fabs(Delta_F[2][1]);
        numberOfCuts=std::floor(Criteria_Delta_F/this->userInputs.criticalDeltaFCriteria);
      }
      else{
        numberOfCuts=this->userInputs.numberTaylorSubsteps;
      }

      if (numberOfCuts==0) numberOfCuts=1;
      inverseNumberOfCuts=1.0/numberOfCuts;
      div_Delta_F.add(inverseNumberOfCuts,Delta_F);
    }
    else{
      numberOfCuts=this->userInputs.numberTaylorSubsteps;
      if (numberOfCuts==0) numberOfCuts=1;
      inverseNumberOfCuts=1.0/numberOfCuts;
    }

    ///Criteria_Delta_F_2 defines if the current F_t is similar to identity matrix (right after reorientation)
    temp=0;
    temp1=IdentityMatrix(dim);
    temp.add(-1,F_t,1.0,temp1) ;
    Criteria_Delta_F_2=0;
    for (unsigned int jj = 0;jj<dim;jj++) {
      for (unsigned int kk = 0;kk<dim;kk++) {
        Criteria_Delta_F_2=Criteria_Delta_F_2+temp[jj][kk]*temp[jj][kk];
      }
    }
    ///The following condition is set for the case of right after twin reorientation. In that case, twin_conv becomes 1 and
    ///F_t becomes the identity matrix (or the difference of F_t from identity matrix becomes very small). In that case, the dt
    ///should not be divided by numberOfCuts, because numberOfCuts is artificially small!!! We keep Dt as it is for this step.
    if ((twin_conv[cellID][quadPtID] != 1.0)||(Criteria_Delta_F_2>1e-8)) {
      ///The Dt for constitutive model integration (not the Algorithmic tangent modulus) should be divided by numberOfCuts because
      ///in the calculation within each mini-step, we update the variable_t to the previously updated variable_tau.
      ///Accordingly, we should use Dt/numberOfCuts as the new timestep in this part (not the Algorithmic tangent modulus).
      ///Updating Dt is refelected here as updating delgam_ref (check line 111 for original calculation of delgam_ref).
      delgam_ref=delgam_ref/numberOfCuts;
    }
    else{
      ///This is the condition when we have right after reorientation
      ///In the case of Taylor model, we should march all the way to the current F from the identity matrix.
      ///The number of Cuts should be defined in this case as below:
      if ((this->userInputs.flagTaylorModel)||(this->userInputs.numberTaylorSubsteps>1)){
        delgam_ref=delgam_ref/numberOfCuts;
        numberOfCuts=numberOfCuts*(this->currentIncrement+1);
        inverseNumberOfCuts=1.0/(this->currentIncrement+1.0);
        div_Delta_F*=inverseNumberOfCuts;
      }
    }

    F_tau=F_t;
    //////////////////////////////////////////////////////////////////

    //////////////Initialization of the varaibles for nonlinear loops which are stress tensor, backstress terms, and slip resistances///////////
    FE_t.Tmmult(CE_t,FE_t);
    temp =0 ;
    temp1= IdentityMatrix(dim) ;
    temp.add(0.5,CE_t,-0.5,temp1) ;
    Dmat.vmult(vtmp1,vecform(temp)) ;
    T_star_iter = 0 ;
    matform(T_star_iter,vtmp1) ;
    if(this->userInputs.enableAdvRateDepModel)
    T_star_iter.add(1.0,Tinter_diff_guess) ;
    nv1 = vecform(T_star_iter) ;
    for(unsigned int j=0;j<2*dim;j++){
      stateVar_it(j) = nv1(j) ;
    }

    for(unsigned int j=0;j<n_Tslip_systems;j++){
      stateVar_it(j+2*dim) = W_kh_tau1(j) ;
    }

    for(unsigned int j=0;j<n_Tslip_systems;j++){
      stateVar_it(j+2*dim+n_Tslip_systems) = W_kh_tau2(j) ;
    }
    s_alpha_it=s_alpha_t;
    ///////////////////////////////////////////////////////////////////////////////////////////////

    for(unsigned int Inc=0;Inc<numberOfCuts;Inc++){

      if (!this->userInputs.flagTaylorModel){
        F_tau.add(1,div_Delta_F) ;
        if (Inc==(numberOfCuts-1)) F_tau=F;
      }
      else{
        iMinusL=IdentityMatrix(dim);
        temp.reinit(dim,dim);
        temp1.reinit(dim,dim);
        temp1=this->targetVelGrad;
        ///We should use Dt/numberOfCuts as new Dt in this part (not the Algorithmic tangent modulus).
        temp1*=this->delT/this->userInputs.numberTaylorSubsteps;
        iMinusL.add(-1,temp1); //I-L
        temp.invert(iMinusL); //inverse(I-L)
        temp.mmult(F_tau,F_tau); // F=inverse(I-L)*Fprev
        if (numberOfCuts==1) F_tau=F;
      }

      del_Ep_tau=0;
      FP_inv_t = 0.0;
      FP_inv_t.invert(FP_t);
      FE_tau_trial=0.0;
      F_tau.mmult(FE_tau_trial,FP_inv_t) ;

      // CE_tau_trial is the same as A matrix - Kalidindi's thesis
      CE_tau_trial=0.0;
      FE_tau_trial.Tmmult(CE_tau_trial,FE_tau_trial);

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      ////////////////////////////////////Start Nonlinear iteration for Slip increments////////////////////////////////////
      ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
      itr1=0; itr2 = 0;
      locres = 1.0; locres2 = 1.0;

      while(locres2>locres_tol2 && itr2 < nitr2){
        itr2 = itr2 + 1  ;
        itr1 = 0 ;
        locres = 1.0 ;
        // Loop until the residual is greater than tolerance or the  number of iterations crosses a certain specified maximum
        while(locres>locres_tol && itr1<nitr1){
          itr1 = itr1+1;
          // Residual for the stress part of the non-linear algebraic equation
          G_iter=0.0;
          // Jacobian for the Newton-Raphson iteration
          J_iter=IdentityMatrix(2*dim+2*n_Tslip_systems);
          locres_vec = 0.0 ;

          for (unsigned int j = 0;j < 2*dim;j++) {
            nv1(j) = stateVar_it(j) ;
          }
          matform(T_star_iter,nv1) ;


          for (unsigned int j = 0;j < n_Tslip_systems;j++) {
            W_kh_tau1_it(j) = stateVar_it(j+2*dim) ;
          }

          for (unsigned int j = 0;j < n_Tslip_systems;j++) {
            W_kh_tau2_it(j) = stateVar_it(j+2*dim+n_Tslip_systems) ;
          }

          W_kh_tau_it = 0.0 ;
          W_kh_tau_it.add(1.0,W_kh_tau1_it,1.0,W_kh_tau2_it) ;

          LP_acc = 0 ;
          for(unsigned int i=0 ; i<n_Tslip_systems;i++){
            for (unsigned int j = 0;j < dim;j++) {
              for (unsigned int k = 0;k < dim;k++) {
                temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
              }
            }

            T_star_iter.mTmult(temp2,temp);
            sctmp1=temp2.trace();
            resolved_shear_tau(i) = sctmp1 ;

            if(i<n_slip_systems){ // For slip systems due to symmetry of slip
              if((sctmp1-W_kh_tau_it(i))<0)
              sgnm=-1;
              else
              sgnm=1 ;
            }
            else               // For twin systems due to asymmetry of slip
            {
              if((sctmp1-W_kh_tau_it(i))<=0)
              sgnm=0;
              else
              sgnm=1 ;
            }

            delgam_tau(i) = delgam_ref*pow(fabs((sctmp1 - W_kh_tau_it(i))/s_alpha_it(i)),1.0/strexp)*sgnm ;
            LP_acc.add(delgam_tau(i),temp) ;

          }

          LP_acc2.equ(1.0,LP_acc) ;
          LP_acc.equ(-1.0,LP_acc) ;
          temp1.equ(1.0,matrixExponential(LP_acc)) ;
          LP_acc.equ(1.0,temp1) ;


          // Loop over slip systems
          for (unsigned int i = 0;i<n_Tslip_systems;i++){

            for (unsigned int j = 0;j < dim;j++) {
              for (unsigned int k = 0;k < dim;k++) {
                temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];

              }
            }

            temp2 = 0.0 ;
            temp2.add(1.0,temp) ;
            temp2.Tadd(1.0,temp) ;


            CE_tau_trial.mmult(temp5,LP_acc) ;
            temp1 = 0;
            temp1.Tadd(1.0,LP_acc2) ;
            temp1.equ(-1.0,temp1) ;
            temp3 = 0 ;
            temp3.Tadd(1.0,temp) ;
            temp4.equ(1.0,matrixExponentialGateauxDerivative2(temp1,temp3)) ;
            temp4.mmult(temp3,temp5) ;
            temp4 = 0 ;
            temp4.add(0.5,temp3) ;
            temp4.Tadd(0.5,temp3) ;
            Dmat.vmult(vtmp1,vecform(temp4)) ;


            nv1.equ(1.0,vtmp1) ;




            nv2 = vecform(temp2);
            nv2(0) = nv2(0)/2.0 ;
            nv2(1) = nv2(1)/2.0 ;
            nv2(2) = nv2(2)/2.0 ;


            btemp1.outer_product(nv1,nv2) ;

            if(i<n_slip_systems){
              sgnm = 1 ;
            }
            else{
              if((resolved_shear_tau(i)-W_kh_tau_it(i))<=0)
              sgnm=0;
              else
              sgnm=1 ;
            }
            sctmp2=delgam_ref/(strexp*s_alpha_it(i))*pow(fabs((resolved_shear_tau(i)-W_kh_tau_it(i))/s_alpha_it(i)),(1.0/strexp - 1.0))*sgnm;

            btemp1.equ(sctmp2,btemp1);

            // Modification to the Jacobian of the Newton-Raphson iteration

            // Components 1:6
            for (unsigned int j = 0;j < 2*dim;j++) {
              for (unsigned int k = 0;k < 2*dim;k++) {
                J_iter[j][k] = J_iter[j][k] + btemp1(j,k);
              }
            }

            vtmp1.equ(-1.0*sctmp2,nv1) ;

            for (unsigned int j = 0;j < 2*dim;j++) {
              J_iter[j][i+2*dim] = J_iter[j][i+2*dim] + vtmp1(j);
              J_iter[j][i+2*dim+n_Tslip_systems] = J_iter[j][i+2*dim+n_Tslip_systems] + vtmp1(j);
            }

            // Components rest
            if(i<n_slip_systems){ // For slip systems due to symmetry of slip
              if((resolved_shear_tau(i)-W_kh_tau_it(i))<0)
              sgnm=-1;
              else
              sgnm=1 ;
            }
            else               // For twin systems due to asymmetry of slip
            {
              if((resolved_shear_tau(i)-W_kh_tau_it(i))<=0)
              sgnm=0;
              else
              sgnm=1 ;
            }

            // Backstress component 1

            if(W_kh_tau1_it(i)<0)
            sgnm2=-1;
            else
            sgnm2=1 ;

            sctmp3 = r_back1*pow(fabs(W_kh_tau1_it(i)/b_back1),m_back1)*W_kh_tau1_it(i)*sgnm*sctmp2 ;
            sctmp3 = sctmp3 - h_back1*sctmp2 ;


            for (unsigned int j = 0;j < 2*dim;j++) {
              J_iter(i+2*dim,j) = sctmp3*nv2(j) ;
            }

            J_iter(i+2*dim,i+2*dim) = J_iter(i+2*dim,i+2*dim) - sctmp3 ;
            J_iter(i+2*dim,i+2*dim+n_Tslip_systems) = J_iter(i+2*dim,i+2*dim+n_Tslip_systems) - sctmp3 ;


            J_iter(i+2*dim,i+2*dim) = J_iter(i+2*dim,i+2*dim) + abs(delgam_tau(i))*r_back1*pow(fabs(W_kh_tau1_it(i)/b_back1),m_back1) ;
            if (fabs(m_back1)>1e-10) {
              J_iter(i+2*dim,i+2*dim) = J_iter(i+2*dim,i+2*dim) + abs(delgam_tau(i))*W_kh_tau1_it(i)*r_back1*m_back1/b_back1*pow(fabs(W_kh_tau1_it(i)/b_back1),m_back1-1.0)*sgnm2 ;
            }
            locres_vec(i+2*dim) = W_kh_tau1_it(i) - W_kh_t1(i) - h_back1*delgam_tau(i) + r_back1*pow(fabs(W_kh_tau1_it(i)/b_back1),m_back1)*W_kh_tau1_it(i)*fabs(delgam_tau(i)) ;

            // Backstress component 2

            sctmp4 = r_back2*pow(fabs(W_kh_tau2_it(i)/b_back2),m_back2)*W_kh_tau2_it(i)*sgnm*sctmp2 ;
            sctmp4 = sctmp4 - h_back2*sctmp2 ;

            if(W_kh_tau2_it(i)<0)
            sgnm2=-1;
            else
            sgnm2=1 ;

            for (unsigned int j = 0;j < 2*dim;j++) {
              J_iter(i+2*dim+n_Tslip_systems,j) = sctmp4*nv2(j) ;
            }

            J_iter(i+2*dim+n_Tslip_systems,i+2*dim+n_Tslip_systems) = J_iter(i+2*dim+n_Tslip_systems,i+2*dim+n_Tslip_systems) - sctmp4 ;
            J_iter(i+2*dim+n_Tslip_systems,i+2*dim) = J_iter(i+2*dim+n_Tslip_systems,i+2*dim) - sctmp4 ;

            J_iter(i+2*dim+n_Tslip_systems,i+2*dim+n_Tslip_systems) = J_iter(i+2*dim+n_Tslip_systems,i+2*dim+n_Tslip_systems) + abs(delgam_tau(i))*r_back2*pow(fabs(W_kh_tau2_it(i)/b_back2),m_back2) ;
            if (fabs(m_back2)>1e-10) {
              J_iter(i+2*dim+n_Tslip_systems,i+2*dim+n_Tslip_systems) = J_iter(i+2*dim+n_Tslip_systems,i+2*dim+n_Tslip_systems) + abs(delgam_tau(i))*W_kh_tau2_it(i)*r_back2*m_back2/b_back2*pow(fabs(W_kh_tau2_it(i)/b_back2),m_back2-1.0)*sgnm2 ;
            }
            locres_vec(i+2*dim+n_Tslip_systems) = W_kh_tau2_it(i) - W_kh_t2(i) - h_back2*delgam_tau(i) + r_back2*pow(fabs(W_kh_tau2_it(i)/b_back2),m_back2)*W_kh_tau2_it(i)*fabs(delgam_tau(i)) ;

          }
          // Construct residual
          vtmp1 = vecform(T_star_iter) ;

          CE_tau_trial.mmult(temp1,LP_acc) ;
          LP_acc.Tmmult(temp2,temp1) ;
          temp1 = IdentityMatrix(dim) ;
          temp3 = 0 ;
          temp3.add(0.5,temp2,-0.5,temp1) ;
          vtmp2 = 0 ;
          Dmat.vmult(vtmp2,vecform(temp3)) ;
          G_iter = 0 ;
          G_iter.add(1.0,vtmp1,-1.0,vtmp2) ;

          for (unsigned int j = 0;j < 2*dim;j++) {
            locres_vec(j) = G_iter(j) ;
          }

          // Invert Jacobian
          J_iter_inv.invert(J_iter);

          J_iter_inv.vmult(stateVar_diff,locres_vec);
          stateVar_diff.equ(-1.0,stateVar_diff) ;
          Fold = locres_vec.l2_norm() ;
          Fold = 0.5*Fold*Fold ;
          J_iter_cp = 0 ;
          J_iter_cp.add(0.5,J_iter) ;
          J_iter_cp.Tadd(0.5,J_iter) ;
          J_iter_cp.vmult(gradFold,locres_vec) ;

          // Implement cubic line search
          lnsrch(stateVar_temp,2*dim+2*n_Tslip_systems,stateVar_it, Fold,gradFold,stateVar_diff,delgam_ref,strexp,SCHMID_TENSOR1,n_slip_systems,n_Tslip_systems,s_alpha_it,Dmat, CE_tau_trial,W_kh_t1, W_kh_t2, h_back1, h_back2,m_back1,m_back2,r_back1, r_back2,b_back1, b_back2) ;

          for (unsigned int j = 0;j < n_Tslip_systems;j++) {
            if(stateVar_temp(2*dim+j)<0)
            sgnm2=-1;
            else
            sgnm2=1 ;

            if(fabs(stateVar_temp(2*dim+j)) >= back_lim_1)
            stateVar_temp(2*dim+j) = back_lim_1*sgnm2 ;
          }


          for (unsigned int j = 0;j < n_Tslip_systems;j++) {
            if(stateVar_temp(2*dim+j+n_Tslip_systems)<0)
            sgnm2=-1;
            else
            sgnm2=1 ;

            if(fabs(stateVar_temp(2*dim+j+n_Tslip_systems)) >= back_lim_2)
            stateVar_temp(2*dim+j+n_Tslip_systems) = back_lim_2*sgnm2 ;
          }

          stateVar_diff = 0.0 ;
          stateVar_diff.add(1.0,stateVar_temp,-1.0,stateVar_it);

          locres = stateVar_diff.l2_norm() ;
          stateVar_it.equ(1.0,stateVar_temp) ;

        } // inner while

        for(unsigned int j=0;j<2*dim;j++){
          T_star_iter_vec(j) = stateVar_it(j) ;
        }

        matform(T_star_iter,T_star_iter_vec) ;

        for(unsigned int j=0;j<n_Tslip_systems;j++){
          W_kh_tau1(j) = stateVar_it(j+2*dim)  ;
          W_kh_tau2(j) = stateVar_it(j+2*dim+n_Tslip_systems) ;
          W_kh_tau(j) = W_kh_tau1(j) + W_kh_tau2(j) ;
        }

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

        s_alpha_iterp = s_alpha_t;

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
          sctmp1 = sctmp1 - W_kh_tau(i) ;
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
          delgam_tau(i)=delgam_ref*pow(fabs((resolved_shear_tau(i)-W_kh_tau(i))/s_alpha_it(i)),(1.0/strexp))*sgnm;
          for (unsigned int j = 0;j<n_Tslip_systems;j++){
            s_alpha_iterp(j)=s_alpha_iterp(j)+h_alpha_beta_t[j][i]*fabs(delgam_tau(i));
          }

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

        vtmp4=0.0;
        vtmp4.add(1.0,s_alpha_iterp,-1.0,s_alpha_it);
        locres2=vtmp4.l2_norm();
        s_alpha_it=s_alpha_iterp;

      } // Outer while
      ////////////////////////////////////End Nonlinear iteration for Slip increments////////////////////////////////////
    ////previous end}

    s_alpha_tau=s_alpha_it;
    s_alpha_t=s_alpha_tau;
    for(unsigned int j=0 ; j<2*dim ; j++)
    nv1(j) = stateVar_it(j) ;
    matform(T_star_iter,nv1) ;

    for(unsigned int j=0 ; j<n_Tslip_systems ; j++){
      W_kh_tau1(j) = stateVar_it(2*dim+j) ;
      W_kh_tau2(j) = stateVar_it(2*dim+j+n_Tslip_systems) ;
      W_kh_tau(j) = W_kh_tau1(j) + W_kh_tau2(j) ;
    }

    FP_tau = 0;
    temp=0.0;
    LP_acc = 0;

    for (unsigned int i = 0;i<n_Tslip_systems;i++){
      for (unsigned int j = 0;j < dim;j++) {
        for (unsigned int k = 0;k < dim;k++) {
          temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
        }
      }

      T_star_iter.mTmult(temp2,temp);
      resolved_shear_tau(i)=temp2.trace();

      if(i<n_slip_systems){ // For slip systems due to symmetry of slip
        if((resolved_shear_tau(i) - W_kh_tau(i))<0)
        sgnm=-1;
        else
        sgnm=1 ;
      }
      else               // For twin systems due to asymmetry of slip
      {
        if((resolved_shear_tau(i) - W_kh_tau(i))<=0)
        sgnm=0;
        else
        sgnm=1 ;
      }

      delgam_tau(i)=delgam_ref*pow(fabs((resolved_shear_tau(i) - W_kh_tau(i))/s_alpha_tau(i)),(1.0/strexp))*sgnm;
      LP_acc.add(delgam_tau(i),temp);
    }

    for(unsigned int i=0 ; i<n_slip_systems ; i++)
    slipfraction_tau(i) = slipfraction_t(i) + fabs(delgam_tau(i)) ;
    slipfraction_t=slipfraction_tau;

    for(unsigned int i=0 ; i<n_Tslip_systems ; i++)
    signed_slip_tau(i) = signed_slip_t(i) + delgam_tau(i) ;
    signed_slip_t=signed_slip_tau;

    for (unsigned int i = 0;i<n_Tslip_systems;i++){
      temp7=0.0;
      for (unsigned int j = 0;j < dim;j++) {
        for (unsigned int k = 0;k < dim;k++) {
          temp7[j][k]=SCHMID_TENSOR1[dim*i + j][k];
        }
      }
      for (unsigned int j = 0;j < dim;j++) {
        for (unsigned int k = 0;k < dim;k++) {
          del_Ep_tau[j][k]=del_Ep_tau[j][k]+ delgam_tau(i)*0.5*(temp7[j][k]+temp7[k][j]);
        }
      }
    }

    for (unsigned int j = 0;j < dim;j++) {
      for (unsigned int k = 0;k < dim;k++) {
        Ep_tau[j][k]=Ep_t[j][k]+del_Ep_tau[j][k];

      }
    }
    Ep_t=Ep_tau;

    del_Ep_eff_cum_tau=sqrt(2.0/3.0)*del_Ep_tau.frobenius_norm();
    Ep_eff_cum_tau=Ep_eff_cum_t+del_Ep_eff_cum_tau;
    Ep_eff_cum_t=Ep_eff_cum_tau;

    for(unsigned int i=0 ; i<n_twin_systems ; i++)
    twinfraction_tau(i) = twinfraction_t(i) + fabs(delgam_tau(i+n_slip_systems))/twinShear;
    twinfraction_t=twinfraction_tau;
    temp.equ(1.0,matrixExponential(LP_acc)) ;
    temp.mmult(FP_tau,FP_t);
    FP_t=FP_tau;
  }

    FP_inv_tau = 0.0; FP_inv_tau.invert(FP_tau);
    FE_tau = 0.0;
    F_tau.mmult(FE_tau, FP_inv_tau);

    ///////////////Calculation of Cuachy Stress Tensor////////////////////
    temp.reinit(dim, dim);
    det_FE_tau = FE_tau.determinant();
    T_star_tau.equ(1.0,T_star_iter);
    FE_tau.mmult(temp, T_star_tau);
    temp.equ(1.0 / det_FE_tau, temp);
    temp.mTmult(T_tau, FE_tau);

    det_F_tau = F_tau.determinant();
    temp.invert(F_tau);
    F_inv_tau.equ(1.0,temp);
    T_tau.mTmult(P_tau, temp);
    P_tau.equ(det_F_tau, P_tau);
    //////////////////////////////////////////////////////////////////////

    ///////////////Calculation of normal stress for each slip system///////
    for (unsigned int i = 0;i<n_Tslip_systems;i++){
      temp7=0.0;
      temp8=0.0;
      for (unsigned int j = 0;j < dim;j++) {
        for (unsigned int k = 0;k < dim;k++) {
          temp7[j][k]=Normal_SCHMID_TENSOR1[dim*i + j][k];
        }
      }

      T_star_tau.mTmult(temp8,temp7);
      normal_stress_tau(i)=temp8.trace();

    }
    //////////////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //////////////////////////////////// Computing Algorithmic Tangent Modulus ////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    if(StiffnessCalFlag==1){
      ///////////////Reading the previous converged state Variables for Algorithmic Tangent Modulus calculation////////////////
        FE_t=Fe_conv[cellID][quadPtID];
        FP_t=Fp_conv[cellID][quadPtID];

        for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
          W_kh_t1(i) = stateVar_conv[cellID][quadPtID][i];
          W_kh_t2(i) = stateVar_conv[cellID][quadPtID][i+n_Tslip_systems];
          signed_slip_t(i)=stateVar_conv[cellID][quadPtID][i+2*n_Tslip_systems];
          s_alpha_t(i)=s_alpha_conv[cellID][quadPtID][i];
        }

        ii=0;
        for(unsigned int i=0 ; i<dim ; i++){
          for(unsigned int j=0 ; j<dim ; j++){
            Ep_t[i][j]=stateVar_conv[cellID][quadPtID][ii+4*n_Tslip_systems];
            ii=ii+1;
          }
        }

        Ep_eff_cum_t=stateVar_conv[cellID][quadPtID][4*n_Tslip_systems+dim*dim];

        for(unsigned int i=0 ; i<n_slip_systems ; i++)
        slipfraction_t(i) = slipfraction_conv[cellID][quadPtID][i] ;

        for(unsigned int i=0 ; i<n_twin_systems ; i++)
        twinfraction_t(i) = twinfraction_conv[cellID][quadPtID][i] ;
      /////////////////////////////////////////////////////////////////////////////

      ///For the Algorithmic Tangent Modulus calculation, the Dt should be returned back to the increment Dt. The reason
      ///for this is the initial value are brought back to the previously converged values.
      ///Accordingly, the time between the variable_t and variable_tau is the original Dt and not the one divided by numberOfCuts
      ///for constitutive model integration. Hence, Dt should be brought back to the original value in line 111 as below:
      delgam_ref = UserMatConstants(0)*this->userInputs.delT ;
      // Contribution 1 - Most straightforward because no need to invoke constitutive model
      FE_tau.mmult(temp,T_star_tau);
      temp.mTmult(temp1,FE_tau);
      temp1.mTmult(temp2,F_inv_tau);
      left(ntemp1,temp2);
      left(ntemp2,F_inv_tau);
      trpose(ntemp3,ntemp2);
      ntemp1.mmult(cnt4,ntemp3);
      cnt4.equ(-1.0,cnt4);

      // Part from isotropic hardening
      Pmat = IdentityMatrix(n_Tslip_systems)  ;

      for(unsigned int j= 0 ; j< n_slip_systems ; j++) {
        sctmp1 = powerLawExponent[j]*initialHardeningModulus[j]/saturationStress[j]*pow((1-s_alpha_tau(j)/saturationStress[j]),powerLawExponent[j]-1)*fabs(delgam_tau(j)) ;
        for(unsigned int i= 0 ; i<n_Tslip_systems ; i++){
          Pmat[i][j] =  Pmat[i][j] + sctmp1*q[i][j] 	;
        }
      }

      for(unsigned int j= 0 ; j< n_twin_systems ; j++) {
        sctmp1 = powerLawExponentTwin[j]*initialHardeningModulusTwin[j]/saturationStressTwin[j]*pow((1-s_alpha_tau(j+n_slip_systems)/saturationStressTwin[j]),powerLawExponentTwin[j]-1)*fabs(delgam_tau(j+n_slip_systems)) ;
        for(unsigned int i= 0 ; i<n_Tslip_systems ; i++){
          Pmat[i][j+n_slip_systems] = Pmat[i][j+n_slip_systems] + sctmp1*q[i][j+n_slip_systems] 	;
        }
      }

      Qmat = 0  ;
      for(unsigned int j= 0 ; j< n_slip_systems ; j++) {
        if(delgam_tau(j) <=0)
        sgnm = -1 ;
        else
        sgnm = 1 ;

        sctmp1 = initialHardeningModulus[j]*pow((1-s_alpha_tau(j)/saturationStress[j]),powerLawExponent[j])*sgnm ;
        for(unsigned int i= 0 ; i<n_Tslip_systems ; i++){
          Qmat[i][j] =  Qmat[i][j] + sctmp1*q[i][j]	;
        }
      }

      for(unsigned int j= 0 ; j< n_twin_systems ; j++) {
        if(delgam_tau(j+n_slip_systems) <=0)
        sgnm = 0 ;
        else
        sgnm = 1 ;
        sctmp1 = initialHardeningModulusTwin[j]*pow((1-s_alpha_tau(j+n_slip_systems)/saturationStressTwin[j]),powerLawExponentTwin[j])*sgnm ;
        for(unsigned int i= 0 ; i<n_Tslip_systems ; i++){
          Qmat[i][j+n_slip_systems] =  Qmat[i][j+n_slip_systems] + sctmp1*q[i][j+n_slip_systems]	;
        }
      }

      Pmat_inv.invert(Pmat) ;
      Pmat_inv.mmult(Pmat,Qmat) ;

      for(unsigned int j= 0 ; j< n_Tslip_systems ; j++) {
        sctmp1 = delgam_tau(j)/strexp/s_alpha_tau(j) ;
        for(unsigned int i= 0 ; i<n_Tslip_systems ; i++){
          Pmat[i][j] = Pmat[i][j]*sctmp1 ;
        }
      }

      Qmat2 = 0 ;
      for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
          }
        }

        temp2.equ(1.0,temp) ;
        temp2.symmetrize() ;
        vtmp1 = vecform(temp2) ;


        Dmat.vmult(vtmp2,vtmp1) ;
        matform(temp2,vtmp2) ;
        FE_tau.mmult(temp3,temp2) ;
        temp2 = 0 ;
        temp2.Tadd(1.0,temp3) ;

        if(i<n_slip_systems){
          sgnm = 1 ;
        }
        else{
          if((resolved_shear_tau(i)-W_kh_tau(i))<=0)
          sgnm=0;
          else
          sgnm=1 ;
        }

        sctmp2=delgam_ref/(strexp*s_alpha_tau(i))*pow(fabs((resolved_shear_tau(i)-W_kh_tau(i))/s_alpha_tau(i)),(1.0/strexp - 1.0))*sgnm;

        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            Qmat2[i][3*j+k] = sctmp2*temp2[k][j] ;
          }
        }
      }

      // Part from backstress
      Rmat = 0 ;

      for(unsigned int i=0 ; i< n_Tslip_systems ; i++){

        if(i<n_slip_systems){
          sgnm = 1 ;
        }
        else{
          if((resolved_shear_tau(i)-W_kh_tau_it(i))<=0)
          sgnm=0;
          else
          sgnm=1 ;
        }
        sctmp2 = delgam_ref/(strexp*s_alpha_tau(i))*pow(fabs((resolved_shear_tau(i)-W_kh_tau(i))/s_alpha_tau(i)),(1.0/strexp - 1.0))*sgnm;

        if(i<n_slip_systems){ // For slip systems due to symmetry of slip
          if(delgam_tau(i)<0)
          sgnm=-1;
          else
          sgnm=1 ;
        }
        else               // For twin systems due to asymmetry of slip
        {
          if(delgam_tau(i)<=0)
          sgnm=0;
          else
          sgnm=1 ;
        }

        modifier1_num = h_back1  - r_back1*pow(fabs(W_kh_tau1(i))/b_back1,m_back1)*W_kh_tau1(i)*sgnm;
        modifier2_num = h_back2  - r_back2*pow(fabs(W_kh_tau2(i))/b_back2,m_back2)*W_kh_tau2(i)*sgnm;

        modifier1_den = 1 + (m_back1+1)*r_back1*pow(fabs(W_kh_tau1(i))/b_back1,m_back1)*fabs(delgam_tau(i)) ;
        modifier2_den = 1 + (m_back2+1)*r_back2*pow(fabs(W_kh_tau2(i))/b_back2,m_back2)*fabs(delgam_tau(i)) ;

        sctmp1 = modifier1_num/modifier1_den + modifier2_num/modifier2_den ;

        modifier = 1 + sctmp1*sctmp2 ;
        Rmat[i][i]	= modifier ;

      }

      Smat = 0 ;
      Smat.add(1.0,Pmat,1.0,Rmat) ;
      Pmat.invert(Smat) ;
      Pmat.mmult(dgamdFe,Qmat2) ;

      dFpdFe = 0 ;

      for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
        for (unsigned int j = 0;j < dim;j++) {
          for (unsigned int k = 0;k < dim;k++) {
            temp[j][k]=SCHMID_TENSOR1[dim*i + j][k];
          }
        }

        temp2.equ(1.0,matrixExponentialGateauxDerivative2(LP_acc,temp)) ;
        temp2.mmult(temp3,FP_t) ;
        vecform9(nvec1,temp3) ;

        for(unsigned int j=0 ; j<dim*dim ; j++){
          nvec2(j) = dgamdFe[i][j] ;
        }

        ntemp1.outer_product(nvec1,nvec2) ;
        dFpdFe.add(1.0,ntemp1) ;

      }
      left(ntemp1,FE_tau) ;
      ntemp1.mmult(ntemp2,dFpdFe) ;
      right(ntemp1,FP_tau) ;
      ntemp3 = 0 ;
      ntemp3.add(1.0,ntemp1,1.0,ntemp2) ;
      dFedF.invert(ntemp3) ;

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

    }// Stiffness cal flag end

    P.reinit(dim, dim);
    P = P_tau;
    T = T_tau;

    T_inter.reinit(dim,dim) ;
    T_inter = T_star_tau ;

    sres_tau.reinit(n_Tslip_systems);
    sres_tau = s_alpha_tau;

    // Update the history variables
    Fe_iter[cellID][quadPtID]=FE_tau;
    Fp_iter[cellID][quadPtID]=FP_tau;
    s_alpha_iter[cellID][quadPtID]=sres_tau;
    W_kh_iter[cellID][quadPtID]=W_kh_tau;

    for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
      stateVar_iter[cellID][quadPtID][i]=W_kh_tau1(i);
      stateVar_iter[cellID][quadPtID][i+n_Tslip_systems]=W_kh_tau2(i);
      stateVar_iter[cellID][quadPtID][i+2*n_Tslip_systems]=signed_slip_tau(i);
      stateVar_iter[cellID][quadPtID][i+3*n_Tslip_systems]=normal_stress_tau(i);
    }
    ii=0;
    for(unsigned int i=0 ; i<dim ; i++){
      for(unsigned int j=0 ; j<dim ; j++){
        stateVar_iter[cellID][quadPtID][ii+4*n_Tslip_systems]=Ep_tau[i][j];
        ii=ii+1;
      }
    }
    stateVar_iter[cellID][quadPtID][dim*dim+4*n_Tslip_systems]=Ep_eff_cum_tau;



    for(unsigned int i=0 ; i<n_twin_systems ; i++)
    twinfraction_iter[cellID][quadPtID][i]=twinfraction_tau(i);

    for(unsigned int i=0 ; i<n_slip_systems ; i++)
    slipfraction_iter[cellID][quadPtID][i]=slipfraction_tau(i);

    if (this->userInputs.flagTaylorModel){
      F=F_tau; // Updating Deformation Gradient if it is Taylor model
    }
    /////// REORIENTATION Due to TWINNING ////////////////////

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
            Fe_iter[cellID][quadPtID]=IdentityMatrix(dim);
            Fp_iter[cellID][quadPtID]=IdentityMatrix(dim);
          }
        }
      }
    }

  }
  #include "../../../include/crystalPlasticity_template_instantiations.h"
