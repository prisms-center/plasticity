//implementation of the continuum plasticity material model
#ifndef CONTINUUMPLASTICITY_H
#define CONTINUUMPLASTICITY_H

//dealii headers
#include "../../../include/ellipticBVP.h"
#include "../../../src/enrichmentModels/enhancedStrain.h"
#include "IntegrationTools/PFunction.hh"
#include "models/PLibrary.hh"


typedef struct {
  double lambda, mu, tau_y, K, H;
  std::string yieldModel, strainEnergyModel, isoHardeningModel;
} materialProperties;

//material model class for continuum plasticity
//derives from ellipticBVP base abstract class
template <int dim>
class continuumPlasticity : public ellipticBVP<dim>
{
 public:
  /**
   *continuumPlasticity class constructor.
   */
  continuumPlasticity();
  /**
   *Structure to hold the material parameters and model names.
   */
  materialProperties properties;
 private:
  /**
   *Initialize and resize class data structures.
   */
  void init(unsigned int num_quad_points);
  void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  //void mesh();
  //void markBoundaries();
  //void applyDirichletBCs();
  /**
   *Update the plastic variables. For a given element and quadrature point, 
   *this takes the proposed deformation gradient for the current iteration
   *and the plastic variables for the previous increment to calculate the 
   *current values for the plastic variables, as well as the stress and tangent modulus.
   */
  void calculatePlasticity(unsigned int cellID,
			   unsigned int quadPtID);
  void getElementalValues(FEValues<dim>& fe_values,
			  unsigned int dofs_per_cell,
			  unsigned int num_quad_points,
			  FullMatrix<double>& elementalJacobian,
			  Vector<double>&     elementalResidual);
  void updateAfterIteration();
  void updateAfterIncrement();

  /**
   *Deformation gradient tensor
   */
  FullMatrix<double> F;
  /**
   *Kirchhoff stress tensor
   */
  FullMatrix<double> tau;
  /**
   *Tangent modulus
   */
  Tensor<4,dim,double> c;
  /**
   *Instantian of the enhanced strain class, used to prevent element locking.
   */
  enhancedStrain<dim> enhStrain;

  /**
   *Converged results of invCP (the inverse of the plastic right Cauchy-Green tensor)
   *for the previous increment, indexed by element number and quadrature point number.
   */
  std::vector< std::vector< FullMatrix<double> > > histInvCP_conv;
  /**
   *Stores values of invCP (the inverse of the plastic right Cauchy-Green tensor)
   *for the most recent iteration, indexed by element number and quadrature point number.
   *Once the solution for the current increment converges, this information will be 
   *transferred to "histInvCP_conv".
   */
  std::vector< std::vector< FullMatrix<double> > > histInvCP_iter;
  /**
   *Converged results of xi (deviatoric back stress)
   *for the previous increment, indexed by element number and quadrature point number.
   */
  std::vector< std::vector< Vector<double> > > histXi_conv;
  /**
   *Stores values of xi (deviatoric back stress)
   *for the most recent iteration, indexed by element number and quadrature point number.
   *Once the solution to the current increment converges, this information will be 
   *transferred to "histXi_conv".
   */
  std::vector< std::vector< Vector<double> > > histXi_iter;
  /**
   *Converged results of alpha (the equivalent plastic strain)
   *for the previous increment, indexed by element number and quadrature point number.
   */
  std::vector< std::vector< double > > histAlpha_conv;
  /**
   *Stores values of alpha (the equivalent plastic strain)
   *for the most recent iteration, indexed by element number and quadrature point number.
   *Once the solution to the current increment converges, this information will be 
   *transferred to "histAlpha_conv".
   */
  std::vector< std::vector< double > > histAlpha_iter;
  /**
   *Store the von Mises stress for each element/quadrature point to project when the
   *solution for the current increment converges.
   */
  std::vector< std::vector< double > > projectVonMisesStress;
  /**
   *Marker to show when plasticity first occurs.
   */
  bool plasticOnset;
  /**
   *Marker to show if the function "init" has been called.
   */
  bool initCalled;
  /**
   *Vector of variables/parameters for the elastic strain energy density function (pfunction).
   */
  std::vector<double> varStrainEnergy;
  /**
   *Vector of variables/parameters for the yield function (pfunction).
   */
  std::vector<double> varYield;
  /**
   *Vector of variables/parameters for the isotropic hardening function (pfunction).
   */
  std::vector<double> varIsoHardening;
  /**
   *Pfunction structure for the elastic strain energy density function.
   */
  PRISMS::PFunction< std::vector<double>, double> strain_energy;
  /**
   *Pfunction structure for the yield function.
   */
  PRISMS::PFunction< std::vector<double>, double> yield;
  /**
   *Pfunction structure for the isotropic hardening function.
   */
  PRISMS::PFunction< std::vector<double>, double> harden;
};

//constructor
template <int dim>
continuumPlasticity<dim>::continuumPlasticity()
: 
ellipticBVP<dim>(),
  enhStrain(this->FE,
	    this->pcout)
{
  //initialize "initCalled"
  initCalled = false;

  //post processing (set up projection of von Mises stress and equivalent plastic strain
  ellipticBVP<dim>::numPostProcessedFields=2;
  ellipticBVP<dim>::postprocessed_solution_names.push_back("alpha");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("tau_vm");
}

template <int dim>
void continuumPlasticity<dim>::init(unsigned int num_quad_points)
{
  //Get the total numbers of elements for this processor (num_local_cells)
  unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();
  //Initiate the enhanced strain object with the number elements
  enhStrain.init_enh_dofs(num_local_cells);

  //Resize the deformation gradient and Kirchhoff stress tensors
  F.reinit(dim, dim);
  tau.reinit(dim, dim);

  //Checkout the elastic strain energy density function from the pfunction library.
  //These are located in the "models" folder. First check that an available
  //model was requested.
  if(properties.strainEnergyModel == "quadlog" || 
     properties.strainEnergyModel == "neohook" || 
     properties.strainEnergyModel == "stvenkir"){
    PRISMS::PLibrary::checkout(properties.strainEnergyModel, strain_energy);	
  }
  else{
    this->pcout << "Using user defined strain energy density function.\n";
    PRISMS::PLibrary::checkout(properties.strainEnergyModel, strain_energy);
  }

  /*strain energy function inputs:
    0: lambda, "First Lame parameter"
    1: mu, "Second Lame parameter"
    2: lambda1, "First principle stretch"
    3: lambda2, "Second principle stretch"
    4: lambda3, "Third principle stretch"*/

  //Resize the vector of parameters/variables used in the strain energy function
  varStrainEnergy.resize(5,1.);
  //Specify the material parameters used in the strain energy function
  varStrainEnergy[0] = properties.lambda;
  varStrainEnergy[1] = properties.mu;

  //For now, specify the hardening model here and check it out from the
  //pfunction library (NOTE: this calculates the value for "q", the conjugate
  // stress-like quantity, used in the yield function).
  //See, for example, the end of section 2.1 of the formulation.
  if(properties.isoHardeningModel == "linear_hardening"){
    PRISMS::PLibrary::checkout(properties.isoHardeningModel, harden);
  }
  else{
    this->pcout << "Using user defined isotropic function.\n";
    PRISMS::PLibrary::checkout(properties.isoHardeningModel, harden);
  }

  /*hardening function inputs:
    0: alpha, "Equivalent plastic strain"
    1: K, "Hardening parameter"*/

  //Resize the vector of parameters/variables used in the hardening function
  varIsoHardening.resize(2,0.);
  //Specify the material parameter used in the hardening function
  varIsoHardening[1] = properties.K;

  //Checkout the yield function from the pfunction library.
  //This is located in the "models" folder. First check that an available
  //model was requested.
  //See equation (16) of the formulation.
  if(properties.yieldModel == "von_mises"){
    PRISMS::PLibrary::checkout(properties.yieldModel, yield);
  }
  else{
    this->pcout << "Using user defined yield function.\n";
    PRISMS::PLibrary::checkout(properties.yieldModel, yield);
  }
		
  /*yield function inputs:
    0: beta1, "First principle stress"
    1: beta2, "Second principle stress"
    2: beta3, "Third principle stress"
    3: tau_y, "Yield stress"
    4: q, "Conjugate stress-like quantity"*/

  //Resize the vector of parameters/variables used in the yield function
  varYield.resize(8,0.);
  //Specify the material parameter used in the yield function
  varYield[3] = properties.tau_y;

  //Resize the vectors of history variables according to the number of elements and quadrature points
	Vector<double> zero_vec(dim); zero_vec = 0.;
  histInvCP_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
  histAlpha_conv.resize(num_local_cells,std::vector<double>(num_quad_points,0));
  histXi_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,zero_vec));
  histInvCP_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
  histAlpha_iter.resize(num_local_cells,std::vector<double>(num_quad_points,0));
  histXi_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,zero_vec));

  //Resize the vector of vector used to store the von Mises stress
  projectVonMisesStress.resize(num_local_cells,std::vector<double>(num_quad_points,0));

  //Initialize the "plasticOnset" marker
  plasticOnset = false;
  //Change "initCalled" to "true" to show that the "init" function has now been called
  initCalled = true;
}

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
	if (stopOnConvergenceFailure) {exit (1);}
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

//implementation of the getElementalValues method
template <int dim>
void continuumPlasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
						  unsigned int dofs_per_cell,
						  unsigned int num_quad_points,
						  FullMatrix<double>& elementalJacobian,
						  Vector<double>&     elementalResidual)
{

  //Initialized history variables and pfunction variables if unititialized
  if(initCalled == false){
    init(num_quad_points);
  }

  //Get the element number.
  unsigned int cellID = fe_values.get_cell()->user_index();
  //Reset the vectors and matrices in static condensation to zero.
  //This is necessary because we are using an enhanced strain (see the paper
  //cited in the formulation).
  enhStrain.staticCondensationData[cellID].reset();

  //Vector relating local to global degree of freedom numbers
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  //Local displacement vector
  Vector<double> Ulocal(dofs_per_cell);

  //Initiate element iterator
  typename DoFHandler<dim>::active_cell_iterator cell(& this->triangulation,
						      fe_values.get_cell()->level(),
						      fe_values.get_cell()->index(),
						      & this->dofHandler);
  //Pass in the element number
  cell->set_user_index(fe_values.get_cell()->user_index());
  //Populate the local to global dof vector
  cell->get_dof_indices (local_dof_indices);

  //Populate the locacal displacement vector
  for(unsigned int i=0; i<dofs_per_cell; i++){
    Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
  }
  //Initialize the enhanced strain object with this information.
  enhStrain.reinit(Ulocal, cell);

  //loop over quadrature points
  for (unsigned int q=0; q<num_quad_points; ++q){
    //Get enhanced deformation gradient
    enhStrain.get_F_enh(q, F);

    //Update strain, stress, and tangent for current time step/quadrature point
    calculatePlasticity(cellID, q);

    //Update block matrices and vectors in enhanced strain
    enhStrain.create_block_mat_vec(F, tau, c, q);

    //Pass local matrices and vectors to static condensation
    //Original dofs matrix
    enhStrain.staticCondensationData[cellID].K.add(enhStrain.fe_values.JxW(q),enhStrain.Klocal);
    //Cross terms matrix
    enhStrain.staticCondensationData[cellID].G.add(enhStrain.fe_values.JxW(q),enhStrain.Glocal);
    //Enhanced dofs matrix
    enhStrain.staticCondensationData[cellID].M.add(enhStrain.fe_values.JxW(q),enhStrain.Mlocal);
    //Original dofs vector
    enhStrain.staticCondensationData[cellID].F.add(enhStrain.fe_values.JxW(q),enhStrain.Flocal);
    //Enhanced dofs vector
    enhStrain.staticCondensationData[cellID].H.add(enhStrain.fe_values.JxW(q),enhStrain.Hlocal);
  }
  //Perform static condensation
  enhStrain.staticCondensationData[cellID].staticCondense();
  //Retrieve the element level jacobian
  elementalJacobian = enhStrain.staticCondensationData[cellID].K2;
  //Retrieve the element level residual
  elementalResidual.equ(-1.,enhStrain.staticCondensationData[cellID].F2);
}

//implementation of the getElementalValues method
template <int dim>
void continuumPlasticity<dim>::updateAfterIteration()
{
  //After solving for the nodal values, calculate the enhanced degrees of freedom.
  std::vector<types::global_dof_index> local_dof_indices(this->FE.dofs_per_cell);
  Vector<double> dUlocal(this->FE.dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(),
    endc = this->dofHandler.end();
  //Loop over elements
  for (; cell!=endc; ++cell) {
    //For the current processor
    if (cell->is_locally_owned()){
      //Update local to global dof map
      cell->get_dof_indices(local_dof_indices);
      //Update local displacement vector
      for(unsigned int i=0; i<dUlocal.size(); i++){
	dUlocal[i] = this->solutionIncWithGhosts[local_dof_indices[i]];
      }
      //Update the enhanced degrees of freedom at each iteration
      enhStrain.updateAlpha(dUlocal, cell->user_index());
    }
  }
}

//implementation of the getElementalValues method
template <int dim>
void continuumPlasticity<dim>::updateAfterIncrement()
{
  //Update the history variables when convergence is reached for the current increment
  histInvCP_conv = histInvCP_iter;
  histAlpha_conv = histAlpha_iter;
  histXi_conv = histXi_iter;

  //fill in post processing field values
  unsigned int cellID=0;
  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), 
    endc = this->dofHandler.end();
  for (; cell!=endc; ++cell) {
    if (cell->is_locally_owned()){
      //loop over quadrature points
      for (unsigned int q=0; q<histAlpha_conv[cellID].size(); ++q){
	//Add the equivalent plastic strain
	this->postprocessValues(cellID, q, 0, 0)=histAlpha_conv[cellID][q];
	//Add the von Mises stress
	this->postprocessValues(cellID, q, 1, 0)=projectVonMisesStress[cellID][q];
      }
      cellID++;
    }
  }

  //call base class project() function to project post processed fields
  ellipticBVP<dim>::project();
}

#endif
