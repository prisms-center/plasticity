//implementation of the continuum plasticity material model
#ifndef CONTINUUMPLASTICITY_H
#define CONTINUUMPLASTICITY_H

//dealii headers
#include "../../../include/ellipticBVP.h"
#include "../../../src/enrichmentModels/enhancedStrain.cc"

typedef struct {
double lambda, mu, tau_y;
std::string yieldModel, strainEnergyModel;
} materialProperties;

//
//material model class for continuum plasticity
//derives from ellipticBVP base abstract class
template <int dim>
class ContinuumPlasticity : public ellipticBVP<dim>
{
public:
  ContinuumPlasticity();
private:
  void markBoundaries();
  void applyDirichletBCs();
	void calculatePlasticity(unsigned int elemNum,
				unsigned int quadNum);
  void getElementalValues(FEValues<dim>& fe_values,
			  unsigned int dofs_per_cell,
			  unsigned int num_quad_points,
			  FullMatrix<double>& elementalJacobian,
			  Vector<double>&     elementalResidual);
  void updateAfterIteration();
  void updateAfterIncrement();

  //quadrature level data structures
  FullMatrix<double> F, tau;
  Tensor<4,dim,double> c;
  enhancedStrain<dim> enhStrain;

  //Store history variables
  std::vector< std::vector< FullMatrix<double> > > histInvCP_conv;
  std::vector< std::vector< FullMatrix<double> > > histInvCP_iter;
  std::vector< std::vector< double > > histAlpha_conv;
  std::vector< std::vector< double > > histAlpha_iter;
  unsigned int plasticOnset; //Marker to show when plasticity first occurs
  std::vector<double> varSE, varYF, varIH; //Vector of variables/parameter for the strain energy, yield function, and isotropic hardening
  PRISMS::PFunction< std::vector<double>, double> strain_energy;
  PRISMS::PFunction< std::vector<double>, double> yield;
  PRISMS::PFunction< std::vector<double>, double> harden;
};

//constructor
template <int dim>
ContinuumPlasticity<dim>::ContinuumPlasticity() : 
ellipticBVP<dim>()
{
	F.reinit(dim, dim);
	tau.reinit(dim, dim);

	if(constLaw == "quadlog" || constLaw == "neohook" || constLaw == "stvenkir"){
  	PRISMS::PLibrary::checkout(constLaw, strain_energy);	
	}
	else{
		pcout << "Error: Constitutive law not recognized.\n";
		exit(1);
	}

  /*strain energy function inputs:
	lambda mu lambda1 lambda2 lambda3
  0: "First Lame parameter"
  1: "Second Lame parameter"
  2: "First principle stretch"
  3: "Second principle stretch"
  4: "Third principle stretch"
  */
	varSE.resize(5,1.);
	varSE[0] = lambda;
	varSE[1] = mu;

  PRISMS::PLibrary::checkout("hardening", harden);

	varIH.resize(1,0.);
//	varIH[1] = kay;		
		
	/*yield function inputs:
	beta1 beta2 beta3 tau_y q
	0: "First principle stress"
	1: "Second principle stress"
	2: "Third principle stress"
	3: "Yield stress"
	4: "Conjugate stress-like quantity"*/

	if(yieldName == "von_mises"){
	  PRISMS::PLibrary::checkout(yieldName, yield);
	}
	else{
		pcout << "Error: Yield function not recognized.\n";
		exit(1);
	}

	varYF.resize(5,0.);
	varYF[3] = tau_y;

	//Resize the vectors of history variables
	histInvCP_conv.resize(numElem,std::vector<FullMatrix<double> >(n_q_pts,IdentityMatrix(dim)));
	histAlpha_conv.resize(numElem,std::vector<double>(n_q_pts,0));
	histInvCP_iter.resize(numElem,std::vector<FullMatrix<double> >(n_q_pts,IdentityMatrix(dim)));
	histAlpha_iter.resize(numElem,std::vector<double>(n_q_pts,0));

	plasticOnset = 0;
}

template <int dim>
void ContinuumPlasticity<dim>::calculatePlasticity(unsigned int elemNum,
							unsigned int quadNum)
{
	double alpha, alpha_TR;
  FullMatrix<double> F_inv, b_e, invCp_TR, invCp, Identity;
	F_inv=0; c=0; invCp=0; tau=0; b_e=0;
	Identity = IdentityMatrix(dim);

	invCp_TR = histInvCP_conv[elemNum][quadNum];
	alpha_TR = histAlpha_conv[elemNum][quadNum];

	Vector<double> Ident_vec(3); Ident_vec = 0.; Ident_vec.add(1.); //"Identity" vector

	F_inv.invert(F);

	//Elastic trial left C-G tensor (b_eTR) at n+1
	FullMatrix<double> b_eTR(dim,dim), invCp_TR_FT(dim,dim), be_FinvT(dim,dim);
	invCp_TR.mTmult(invCp_TR_FT,F); //invCp_FT = invCp_TR*F^T
	F.mmult(b_eTR,invCp_TR_FT); //b_eTR = F*invCp_TR*F^T

	//Eigenvalues of b_eTR with LAPACK
	Vector<double> eigTR(dim);
	LAPACKFullMatrix<double> b_eTReig(dim,dim), Ident(dim,dim);
	std::vector< Vector<double> > eigVec(dim);
	b_eTReig = b_eTR;
	Ident = Identity;

	b_eTReig.compute_generalized_eigenvalues_symmetric(Ident,eigVec);
	std::vector< SymmetricTensor<2,dim,double> > eigDyad_tens(dim);
	std::vector< FullMatrix<double> > eigDyad(dim);
	Tensor<2,dim,double> eigTens;
	for(unsigned int i=0; i<dim; i++){
		//Eigenvalues of b_TR^e are square of trial elastic stretches
		eigTR(i) = abs(b_eTReig.eigenvalue(i));
		eigDyad[i].outer_product(eigVec[i],eigVec[i]);
	}

	Vector<double> d_A(dim); d_A = 0.;
	for(unsigned int i=0; i<dim; i++){
		d_A(i) = 1.;
		for(unsigned int j=1; j<dim; j++){
			if(std::abs(eigTR(i) - eigTR((i+j)%3)) < 1.e-13){
				d_A(i) *= copysign(5.e-8,eigTR(i) - eigTR((i+j)%3));
			}
			else{
				d_A(i) *= eigTR(i) - eigTR((i+j)%3);
			}
		}
	}

	//Logarithmic stretches (n+1)
	//If elastic, trial values hold
	Vector<double> epsTRe(dim), eps_e(dim), beta(dim), nu_bar(dim);
	for(unsigned int i=0; i<dim; i++){
		varSE[2+i] = sqrt(eigTR(i));
		epsTRe(i) = 0.5*log(eigTR(i));
	}
	double gamma_Dt = 0.;
	eps_e = epsTRe;
	for(unsigned int i=0; i<dim; i++){
		varYF[i] = strain_energy.grad(varSE,2+i)*varSE[2+i]; //Beta
		beta(i) = varYF[i];
	}
	alpha = alpha_TR;
	varIH[0] = alpha;
	varYF[4] = harden(varIH); //q
	double yield_TR = yield(varYF);
	if(yield_TR > 0){ //Plastic
		//Report onset of plasticity
		if(plasticOnset == 0){
			pcout << "Onset of plasticity\n";
			plasticOnset = 1;
		}

		//Newton-Raphson routine to solve for alpha, Beta => eps_e
		if(alpha_TR == 0){
			varIH[0] = 1.e-4; //So that the derivativate of the hardening function isn't undefined.
		}
		const double Tolerance = 1.e-12;
		Vector<double> Function(1+dim), product(1+dim);
		FullMatrix<double> jacobian(1+dim,1+dim);
		unsigned int iter=0;
		while(true){
			
			//Fill in Function vector and Jacobian matrix for current time step
 			jacobian = 0.;
			for (unsigned int B=0; B<dim; ++B){
				for (unsigned int A=0; A<dim; ++A){
					for (unsigned int C=0; C<dim; ++C){
						jacobian[A][B] -= gamma_Dt*yield.hess(varYF,A,C)*(strain_energy.grad(varSE,2+C)*(B==C) + 
															strain_energy.hess(varSE,2+B,2+C)*varSE[2+C]);
					}
					jacobian[A][B] -= 1./varSE[2+A]*(A==B);
					jacobian[3][A] += yield.grad(varYF,B)*(strain_energy.grad(varSE,2+B)*(A==B) + 
														strain_energy.hess(varSE,2+A,2+B)*varSE[2+B]);	
				}
				jacobian[B][3] = -(1./yield.grad(varYF,4) - (varIH[0] - alpha_TR)*pow(yield.grad(varYF,4),-2.)*
													yield.hess(varYF,4,4)*harden.grad(varIH,0))*yield.grad(varYF,B);
				Function(B) = epsTRe(B) - gamma_Dt*yield.grad(varYF,B) - log(varSE[2+B]);
			}
			jacobian[3][3] = yield.grad(varYF,4)*harden.grad(varIH,0);
			Function(3) = yield(varYF);

			//Invert jacobian matrix
			jacobian.gauss_jordan();
			jacobian.vmult(product,Function);

			//Update Newton-Raphson variables
			for(unsigned int i = 0; i < dim; i++){
				varSE[2+i] -= product(i); //lambda_e(i)
			}

			//Update other variables
			for(unsigned int i = 0; i < dim; i++){
				varYF[i] = strain_energy.grad(varSE,2+i)*varSE[2+i]; //beta(i)
			}
			varIH[0] -= product(3); //alpha
			varYF[4] = harden(varIH); //q
			gamma_Dt = (varIH[0] - alpha_TR)/yield.grad(varYF,4);

			if(product.l2_norm() < Tolerance){break;}
			if(iter>30){
				pcout << "  Maximum number of iterations reached without convergence. \n";
				break;
			}
			iter++;
		}
		gamma_Dt *= 1 - 1.e-8;
		alpha = alpha_TR + gamma_Dt*yield.grad(varYF,4);
		varIH[0] = alpha;
		varYF[4] = harden(varIH);
		for(unsigned int i = 0; i < dim; i++){
			beta(i) = varYF[i];
			nu_bar(i) = yield.grad(varYF,i);
		}
		eps_e.equ(1.,epsTRe,-gamma_Dt,nu_bar);
	}
	tau.equ(beta(0),eigDyad[0],beta(1),eigDyad[1],beta(2),eigDyad[2]);
	b_e.equ(exp(2.*eps_e(0)),eigDyad[0],exp(2.*eps_e(1)),eigDyad[1],exp(2.*eps_e(2)),eigDyad[2]);
	b_e.mTmult(be_FinvT,F_inv);
	F_inv.mmult(invCp,be_FinvT);

	//Store history variables for this iteration
	histInvCP_iter[elemNum][quadNum] = invCp;
	histAlpha_iter[elemNum][quadNum] = alpha;

	//Determine algorithmic elastoplastic tangent
	FullMatrix<double> a_ep(dim,dim), a_e(dim,dim);
	for (unsigned int A=0; A<dim; ++A){
		for (unsigned int B=0; B<dim; ++B){
			a_e[A][B] = strain_energy.hess(varSE,2+A,2+B)*varSE[2+A]*varSE[2+B] + 
									strain_energy.grad(varSE,2+A)*varSE[2+A]*(A==B);
		}
	}
	if(yield_TR < 0){ //Elastic
		a_ep = a_e;
	}
	else{ //Plastic
		FullMatrix<double> a_eInv(dim,dim), f2_B(dim,dim);
		a_eInv.invert(a_e);
		for (unsigned int A=0; A<dim; ++A){
			for (unsigned int B=0; B<dim; ++B){
				f2_B[A][B] = yield.hess(varYF,A,B);
			}
		}

		FullMatrix<double> h(dim,dim);
		h.equ(1.,a_eInv,gamma_Dt,f2_B);
		h.gauss_jordan();

		Vector<double> h_nu_bar(dim);
		h.vmult(h_nu_bar,nu_bar);
	
		a_ep.outer_product(h_nu_bar,h_nu_bar);
		a_ep.equ(1.,h,-1./(nu_bar*h_nu_bar - pow(yield.grad(varYF,4),2.)/(1./harden.grad(varIH,0) - gamma_Dt*yield.hess(varYF,4,4))),a_ep);

	}
	for (unsigned int i=0; i<dim; ++i){
		for (unsigned int j=0; j<dim; ++j){
			for (unsigned int k=0; k<dim; ++k){
				for (unsigned int l=0; l<dim; ++l){
					if(currentIncrement == 0 && currentIteration==0){
						double lambda = varSE[0], mu = varSE[1];
						c[i][j][k][l] = lambda*(i==j)*(k==l) + mu*((i==k)*(j==l) + (i==l)*(j==k));
						//c_{ijkl}=C_{IJKL} in the limit as stretches go to one
					}
					else{
						for (unsigned int A=0; A<dim; ++A){
							for (unsigned int B=0; B<dim; ++B){
								c[i][j][k][l] += a_ep[A][B]*eigDyad[A][i][j]*eigDyad[B][k][l];
							}
							c[i][j][k][l] += 2.*beta(A)/d_A(A)*(0.5*(b_eTR[i][k]*b_eTR[j][l] + b_eTR[i][l]*b_eTR[j][k])
									- b_eTR[i][j]*b_eTR[k][l]	- b_eTR.determinant()/eigTR(A)*(0.5*((i==k)*(j==l) + (i==l)*(j==k))
									- ((i==j) - eigDyad[A][i][j])*((k==l) - eigDyad[A][k][l]))	+ eigTR(A)*(b_eTR[i][j]*eigDyad[A][k][l]
									+ eigDyad[A][i][j]*b_eTR[k][l] + (b_eTR.trace() - 4.*eigTR(A))*eigDyad[A][i][j]*eigDyad[A][k][l]));
						}
					}
				}
			}
		}
	}
}

//implementation of the getElementalValues method
template <int dim>
void ContinuumPlasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
							unsigned int dofs_per_cell,
							unsigned int num_quad_points,
							FullMatrix<double>& elementalJacobian,
							Vector<double>&     elementalResidual)
{
	unsigned int cellID = fe_values.get_cell()->user_index();

	//Reset the vectors and matrices in static condensation to zero
	enhStrain.staticCondensationData[cellID].reset();

	std::vector<double> local_dof_indices(dofs_per_elem);
	Vector<double> Ulocal(dofs_per_elem);

	fe_values.get_cell()->get_dof_indices (local_dof_indices);
	for(unsigned int i=0; i<dofs_per_elem; i++){
		Ulocal[i] = solutionWithGhosts[local_dof_indices[i]];
	}
	enhStrain.reinit(Ulocal, fe_values.get_cell());

	//loop over quadrature points
	for (unsigned int q=0; q<num_quad_points; ++q){
		//Get enhanced deformation gradient
		enhStrain.get_F_enh(q, F);

		//Update strain, stress, and tangent for current time step/quadrature point
		calculatePlasticity(cellID, q);

		//Update block matrices and vectors in enhanced strain
		enhStrain.create_block_mat_vec(F, tau, c, q);

		//Pass local matrices and vectors to static condensation
		enhStrain.staticCondensationData[cellID].K.add(enhStrain.fe_values.JxW(q),enhStrain.Klocal);
		enhStrain.staticCondensationData[cellID].G.add(enhStrain.fe_values.JxW(q),enhStrain.Glocal);
		enhStrain.staticCondensationData[cellID].M.add(enhStrain.fe_values.JxW(q),enhStrain.Mlocal);
		enhStrain.staticCondensationData[cellID].F.add(enhStrain.fe_values.JxW(q),enhStrain.Flocal);
		enhStrain.staticCondensationData[cellID].H.add(enhStrain.fe_values.JxW(q),enhStrain.Hlocal);
	}
	enhStrain.staticCondensationData[elem->user_index()].staticCondense();

	elementalJacobian = enhStrain.staticCondensationData[cellID].K2;
	elementalResidual.equ(-1.,enhStrain.staticCondensationData[cellID].F2);
}

//implementation of the getElementalValues method
template <int dim>
void ContinuumPlasticity<dim>::updateAfterIteration()
{
	//Update the enhanced degrees of freedom at each iteration
	enhStrain.updateAlpha(dUlocal, ElemNum);
}

//implementation of the getElementalValues method
template <int dim>template <int dim>
void ContinuumPlasticity<dim>::updateAfterIncrement()
{
	//Update the history variables when convergence is reached for the current increment
	histInvCP_conv = histInvCP_iter;
	histAlpha_conv = histAlpha_iter;
}

#endif
