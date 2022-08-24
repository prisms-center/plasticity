#include "../../../include/continuumPlasticity.h"

template <int dim>
void continuumPlasticity<dim>::calculatePlasticity(unsigned int cellID,
						   unsigned int quadPtID)
{
  //Actual and trial values of the equivalent plastic strain, respectively
  double alpha, alpha_TR;
  //Actual and trial values of the back stress vector, respectively
  Vector<double> xi(dim), xi_TR(dim);
  //Inverse deformation gradient and elastic left C-G tensor, repsectively
  FullMatrix<double> F_inv(dim,dim), b_e(dim,dim);
  //Trial and actual values of the inverse of plastic right C-G tensor, respectively
  FullMatrix<double> invCp_TR(dim,dim), invCp(dim,dim);
  //Isotropic tensor (identity matrix)
  FullMatrix<double> Identity(dim,dim);
  //Initialize tensors to zero
  F_inv=0; c=0; invCp=0; tau=0; b_e=0;
  //Initialize isotropic tensor to correct values.
  Identity = IdentityMatrix(dim);

  //Set the trial values of invCP, alpha, and xi for the current increment equal to the
  //converged values of invCP and alpha from the previous increment
  //see equations (23) and (24) of the formulation
  invCp_TR = histInvCP_conv[cellID][quadPtID];
  alpha_TR = histAlpha_conv[cellID][quadPtID];
  xi_TR = histXi_conv[cellID][quadPtID];

  //Initialize the "identity vector", i.e. the vector [1,1,1] in 3D
  Vector<double> Ident_vec(3); Ident_vec = 0.; Ident_vec.add(1.);
  //Find the inverse of the deformation gradient
  F_inv.invert(F);

  //Elastic trial left C-G tensor (b_eTR) at n+1
  FullMatrix<double> b_eTR(dim,dim), invCp_TR_FT(dim,dim), be_FinvT(dim,dim);
  //invCp_FT = invCp_TR*F^T (temporary variable)
  invCp_TR.mTmult(invCp_TR_FT,F);
  //b_eTR = F*invCp_TR*F^T (see equation (24) of the formulation)
  F.mmult(b_eTR,invCp_TR_FT);

  //Compute the eigenvalues/eigenvectors of b_eTR. These are used in the spectral
  //decomposition of b_e and tau. See equations (27)-(29).

  //Eigenvalues of b_TR^e are square of trial elastic stretches.
  //Compare to equation (28) of the formulation
  Vector<double> eigTR(dim);
  //Some temporary objects are used in the computation.
  LAPACKFullMatrix<double> b_eTReig(dim,dim), Ident(dim,dim);
  std::vector< Vector<double> > eigVec(dim);
  b_eTReig = b_eTR;
  Ident = Identity;
  /* if(this->currentIncrement == 1){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
	this->pcout << invCp_TR[i][j] << " ";
      }
      this->pcout << std::endl;
    }
    this->pcout << std::endl;
    }*/

  //Use LAPACK to solve for the eigenvalues/eigenvectors.
  b_eTReig.compute_generalized_eigenvalues_symmetric(Ident,eigVec);
  //The dyadic product of the eigenvectors is used in the spectral
  //decomposition. Again, see equations (28) and (29).
  std::vector< FullMatrix<double> > eigDyad(dim);
  for(unsigned int i=0; i<dim; i++){
    eigTR(i) = abs(b_eTReig.eigenvalue(i));
    eigDyad[i].outer_product(eigVec[i],eigVec[i]);
  }

  //These variables are used in the elastoplastic tangent. See equation (42).
  Vector<double> d_A(dim); d_A = 0.;
  for(unsigned int i=0; i<dim; i++){
    d_A(i) = 1.;
    for(unsigned int j=1; j<dim; j++){
      //Since we will divide by d_A, we need to be careful if it is almost zero.
      if(std::abs(eigTR(i) - eigTR((i+j)%3)) < 1.e-13){
	d_A(i) *= copysign(5.e-8,eigTR(i) - eigTR((i+j)%3));
      }
      else{
	d_A(i) *= eigTR(i) - eigTR((i+j)%3);
      }
    }
  }

  //The principal elastic stretches
  Vector<double> lambda_e(dim);
  //The principal stresses and the derivative of the yield function w.r.t. the principal stresses
  Vector<double> beta(dim), nuBar(dim);
  for(unsigned int i=0; i<dim; i++){
    //The principal elastic stretches, lambda_e(i)
    varStrainEnergy[2+i] = sqrt(eigTR(i));
    lambda_e[i] = sqrt(eigTR(i));
  }
  double gamma_Dt = 0.;
  for(unsigned int i=0; i<dim; i++){
    //Calculate the principal stresses, beta, using equation (32)
    varYield[i] = strain_energy.grad(varStrainEnergy,2+i)*varStrainEnergy[2+i]; //Beta
    beta(i) = varYield[i];
  }
  //Initially take alpha as alpha_TR
  alpha = alpha_TR;
  //Initially take xi as xi_TR
  xi = xi_TR;
  //alpha is the first variable needed in the hardening function
  varIsoHardening[0] = alpha;
  //The variable, q, needed in the yield function gets its value from the hardening function
  varYield[4] = harden(varIsoHardening);
	//The back stress is needed for the yield function
  for(unsigned int i=0; i<dim; i++){
    varYield[dim+2+i] = xi(i);
  }
  //Get the value of the yield function, based on the trial values
  double yield_TR = yield(varYield);
  if (pausePlastic)
      yield_TR = 0;

	//if(cellID == 0)
	//	std::cout << "Trial yield: " << yield(varYield) << std::endl;

	//Used as a fifth unknown in the nonlinear solve, zeta is to be equal to (df/dq)_{n+1}
	double zeta = yield.grad(varYield,4);

  //If the yield function is greater than zero, plastic flow has occured.
  //We will need to update the plastic variables (see section 3 of the formulation) using Newton-Raphson.
  //The independent variables that we solve for are the 3 principal elastic stretches (lambda_e(i))
  //and the equivalent plastic strain (alpha). The four equations used are equation (30) and a variation
  //of equation (26).
  if(yield_TR > 0){
    //Report onset of plasticity
    if(plasticOnset == false){
      this->pcout << "\ncontinuumPlasticity: Onset of plasticity\n\n";
      plasticOnset = true;
    }

    //Being careful dividing by zero...
    if(alpha_TR == 0){
      varIsoHardening[0] = 1.e-4; //So that the derivativate of the hardening function isn't undefined.
    }
    //Set the tolerance for convergence
    const double Tolerance = 1.e-12;
    //Declare the jacobian and residual to solve for the plastic variables.
    //(product = inv(jacobian)*residual)
    Vector<double> residual(2+dim), product(2+dim);
    FullMatrix<double> jacobian(2+dim,2+dim);
    unsigned int iter=0;
    //Iterate until convergence is met
    while(true){

      //Fill in Residual vector and Jacobian matrix for current time step
      //(Note: the derivation of the jacobian isn't included in the formulation)
      jacobian = 0.;
      for (unsigned int B=0; B<dim; ++B){
	for (unsigned int A=0; A<dim; ++A){
	  for (unsigned int C=0; C<dim; ++C){
	    jacobian[A][B] -= gamma_Dt*yield.hess(varYield,A,C)*(strain_energy.grad(varStrainEnergy,2+C)*(B==C) +
								 strain_energy.hess(varStrainEnergy,2+B,2+C)*varStrainEnergy[2+C]);
	  }
	  jacobian[A][B] -= 1./varStrainEnergy[2+A]*(A==B);
	  jacobian[3][A] += yield.grad(varYield,B)*(strain_energy.grad(varStrainEnergy,2+B)*(A==B) +
						    strain_energy.hess(varStrainEnergy,2+A,2+B)*varStrainEnergy[2+B]);
	}
	jacobian[B][3] = -(1./yield.grad(varYield,4) - (varIsoHardening[0] - alpha_TR)*pow(yield.grad(varYield,4),-2.)*
			   yield.hess(varYield,4,4)*harden.grad(varIsoHardening,0))*yield.grad(varYield,B);
	//This equation is half of the log of equation(26)
	residual(B) = 0.5*log(eigTR(B)) - gamma_Dt*yield.grad(varYield,B) - log(varStrainEnergy[2+B]);
      }
      jacobian[3][3] = yield.grad(varYield,4)*harden.grad(varIsoHardening,0);
			jacobian[4][3] = -yield.hess(varYield,4,4)*harden.grad(varIsoHardening,0);
			jacobian[4][4] = 1;

      //Equation (30)
      residual(3) = yield(varYield);
			residual(4) = zeta - yield.grad(varYield,4);

      //Invert jacobian matrix
      jacobian.gauss_jordan();
      //Multiply the jacobian inverse and the function, the result is "product"
      jacobian.vmult(product,residual);

      //Update Newton-Raphson variables
      for(unsigned int i = 0; i < dim; i++){
	//Update lambda_e(i)
	varStrainEnergy[2+i] -= product(i);
      }
      //Update alpha
      varIsoHardening[0] -= product(3);
			//Update zeta
			zeta -= product(4);

      //Update other variables
      for(unsigned int i = 0; i < dim; i++){
	//Update beta, equation (9)
	varYield[i] = strain_energy.grad(varStrainEnergy,2+i)*varStrainEnergy[2+i]; //beta(i)
      }
      //Update q, directly from the hardening function
      varYield[4] = harden(varIsoHardening);
      //Update gamma_Dt, based on equation (25). Pull in a little from the yield surface for stability.
      gamma_Dt = (varIsoHardening[0] - alpha_TR)/zeta;
			//Update xi
	  	for(unsigned int i=0; i<dim; i++){
				nuBar(i) = (log(varStrainEnergy[2+i]) - 0.5*log(eigTR(i)))/(-gamma_Dt);
    		varYield[dim+2+i] = xi_TR(i) + 2./3.*gamma_Dt*properties.H*nuBar(i);
  		}

      //Check for convergence
      if(product.l2_norm() < Tolerance || residual.l2_norm() < Tolerance){
	break;
      }
      if(iter>30){
	this->pcout << "  During update of plastic variables: Maximum number of iterations reached without convergence. \n";
	this->pcout <<  "  Consider using a smaller load or a higher number of increments. \n";
	if (properties.stopOnConvergenceFailure) {exit (1);}
	else {this->pcout << "   stopOnConvergenceFailure==false, so marching ahead\n";}
	break;
      }
      iter++;
    }
		//if(cellID == 0)
		//	std::cout << "Before: " << yield(varYield) << " " << gamma_Dt << std::endl;
//*
    //For stability, pull in a little from the yield surface. Recalculate the updated
    //variables with this adjusted value for gamma_Dt.
    gamma_Dt *= 1 - 1.e-8;

    //Update alpha, equation (25)
    alpha = alpha_TR + gamma_Dt*yield.grad(varYield,4);
    varIsoHardening[0] = alpha;
    //Update the isotropic hardening, q
    varYield[4] = harden(varIsoHardening);
    for(unsigned int i = 0; i < dim; i++){
      //Update the principal stresses
      beta(i) = varYield[i];
      //Update the derivative of the yield function, f w.r.t. beta
      nuBar(i) = yield.grad(varYield,i);
      //Update the principal elastic stretches, equation (26)
      lambda_e(i) = exp(-gamma_Dt*nuBar(i))*sqrt(eigTR(i));
			varStrainEnergy[2+i] = lambda_e(i);
			//Update the back stress, xi
    	xi(i) = xi_TR(i) + 2./3.*gamma_Dt*properties.H*nuBar(i);
			varYield[dim+2+i] = xi(i);
    }//*/

		//if(cellID == 0)
		//	std::cout << "After: " << yield(varYield) << " " << gamma_Dt << std::endl;

  }
  //Find tau, equation (29)
  tau.equ(beta(0),eigDyad[0],beta(1),eigDyad[1],beta(2),eigDyad[2]);
  //Find b_e, equation (28)
  b_e.equ(std::pow(lambda_e(0),2),eigDyad[0],std::pow(lambda_e(1),2),eigDyad[1],std::pow(lambda_e(2),2),eigDyad[2]);
  //Find invCP, based on equation (4)
  b_e.mTmult(be_FinvT,F_inv); //temporary variable
  F_inv.mmult(invCp,be_FinvT); //invCp = F_inv*b_e*F_inv^T

  //Store history variables for this iteration
  histInvCP_iter[cellID][quadPtID] = invCp;
  histAlpha_iter[cellID][quadPtID] = alpha;
  histXi_iter[cellID][quadPtID] = xi;

  //Store the von Mises stress to project to the nodes and include in output file
  projectVonMisesStress[cellID][quadPtID] = std::sqrt(0.5*(std::pow(beta(0) - beta(1),2.) +
							   std::pow(beta(1) - beta(2),2.) +
							   std::pow(beta(2) - beta(0),2.)));

  //Determine algorithmic elastoplastic tangent (section 4 of the formulation)
  FullMatrix<double> a_ep(dim,dim), a_e(dim,dim);
  for (unsigned int A=0; A<dim; ++A){
    for (unsigned int B=0; B<dim; ++B){
      //Equation (36)
      a_e[A][B] = strain_energy.hess(varStrainEnergy,2+A,2+B)*varStrainEnergy[2+A]*varStrainEnergy[2+B] +
	strain_energy.grad(varStrainEnergy,2+A)*varStrainEnergy[2+A]*(A==B);
    }
  }
  if(yield_TR <= 0){
    //If elastic...
    a_ep = a_e;
  }
  else{
    //If plastic...

    //The next several lines go into equation (35)
    FullMatrix<double> a_eInv(dim,dim), f2_B(dim,dim);
    a_eInv.invert(a_e);
    //Find the 2nd derivatives of f w.r.t. beta
    for (unsigned int A=0; A<dim; ++A){
      for (unsigned int B=0; B<dim; ++B){
	f2_B[A][B] = yield.hess(varYield,A,B);
      }
    }

    FullMatrix<double> h1(dim,dim), h2(dim,dim), h3(dim,dim), h3h2(dim,dim), h2h3h2(dim,dim), h2_etc(dim,dim);
		double h4;
		//Find h1
		h1.equ(1,Identity,2./3.*properties.H*gamma_Dt,f2_B);
		h1.gauss_jordan();
		//Find h2
		h1.mmult(h2,f2_B); //h1*f2_B
		h2.equ(1,Identity,-2./3.*properties.H*gamma_Dt,h2);
		//Find h3
		f2_B.mmult(h3,h2); //f2_B*h2
    h3.equ(1.,a_eInv,gamma_Dt,h3);
    h3.gauss_jordan();
		//Find h4
		h4 = 1. - gamma_Dt*harden.grad(varIsoHardening,0)*yield.hess(varYield,4,4);
		//Find h3*h2
		h3.mmult(h3h2,h2);
		//Find h2*h3*h2
		h2.mmult(h2h3h2,h3h2);
		//Find h2*h3*h2 + 2/3*H*h1
		h2_etc.equ(1,h2h3h2,2./3.*properties.H,h1);

    Vector<double> h3h2_nuBar(dim), h2etc_nuBar(dim);
    //Temporary variables, recall that nuBar is the partial of f w.r.t. beta
    h3h2.vmult(h3h2_nuBar,nuBar);
		h2_etc.vmult(h2etc_nuBar,nuBar);

    //The next two lines are equation (34)
    a_ep.outer_product(h3h2_nuBar,h3h2_nuBar);
    a_ep.equ(1.,h3,-h4/(h4*(nuBar*h2etc_nuBar) - pow(yield.grad(varYield,4),2.)*harden.grad(varIsoHardening,0)),a_ep);

  }
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	for (unsigned int l=0; l<dim; ++l){
	  if(this->currentIncrement == 0 && this->currentIteration==0){
	    //For readability, extract the Lame parameters
	    double lambda = varStrainEnergy[0], mu = varStrainEnergy[1];
	    //c_{ijkl}=C_{IJKL}, the standard modulus, in the limit as stretches go to one
	    //This is done to avoid dividing by zero.
	    c[i][j][k][l] = lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));
	  }
	  else{
	    //These two for loops are computing equations (33) and (38) of the formulation
	    for (unsigned int A=0; A<dim; ++A){
	      for (unsigned int B=0; B<dim; ++B){
		c[i][j][k][l] += a_ep[A][B]*eigDyad[A][i][j]*eigDyad[B][k][l];
	      }
	      c[i][j][k][l] += 2.*beta(A)/d_A(A)*(0.5*(b_eTR[i][k]*b_eTR[j][l] + b_eTR[i][l]*b_eTR[j][k])
						  - b_eTR[i][j]*b_eTR[k][l]
						  - b_eTR.determinant()/eigTR(A)*(0.5*((i==k)*(j==l) + (i==l)*(j==k))
										  - ((i==j) - eigDyad[A][i][j])*((k==l) - eigDyad[A][k][l]))
						  + eigTR(A)*(b_eTR[i][j]*eigDyad[A][k][l]
							      + eigDyad[A][i][j]*b_eTR[k][l]
							      + (b_eTR.trace() - 4.*eigTR(A))*eigDyad[A][i][j]*eigDyad[A][k][l]));
	    }
	  }
	}
      }
    }
  }
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
