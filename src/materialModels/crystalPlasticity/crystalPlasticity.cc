//implementation of the crystal plasticity material model
#ifndef CRYSTALPLASTICITY_H
#define CRYSTALPLASTICITY_H

//dealii headers
#include "../../../include/ellipticBVP.h"
#include "../../../src/utilityObjects/crystalOrientationsIO.cc"

typedef struct {
	unsigned int n_slip_systems; //No. of slip systems
	double q1,q2,a,h0,s_s,s0,C11,C12,C44;
	FullMatrix<double> m_alpha,n_alpha;	
} materialProperties;

//material model class for crystal plasticity
//derives from ellipticBVP base abstract class
template <int dim>
class crystalPlasticity : public ellipticBVP<dim>
{
public:
	crystalPlasticity();
	void mesh();
	void reorient();
	void tangent_modulus(FullMatrix<double> &F_trial, FullMatrix<double> &Fpn_inv, FullMatrix<double> &SCHMID_TENSOR1, FullMatrix<double> &A,FullMatrix<double> &A_PA,FullMatrix<double> &B,FullMatrix<double> &T_tau, FullMatrix<double> &PK1_Stiff, Vector<double> &active, Vector<double> &resolved_shear_tau_trial, Vector<double> &x_beta, Vector<double> &PA, int &n_PA, double &det_F_tau, double &det_FE_tau );
	void inactive_slip_removal(Vector<double> &inactive,Vector<double> &active, Vector<double> &x_beta, int &n_PA, Vector<double> &PA, Vector<double> b,FullMatrix<double> A,FullMatrix<double> A_PA);
	//material properties
	materialProperties properties;
	//orientation maps
	crystalOrientationsIO<dim> orientations;  
private:
	void init(unsigned int num_quad_points);
	void markBoundaries();
	void applyDirichletBCs();
	void calculatePlasticity(unsigned int cellID,
		unsigned int quadPtID);
	void getElementalValues(FEValues<dim>& fe_values,
		unsigned int dofs_per_cell,
		unsigned int num_quad_points,
		FullMatrix<double>& elementalJacobian,
		Vector<double>&     elementalResidual);
	void updateAfterIncrement();


	void odfpoint(FullMatrix <double> &OrientationMatrix,Vector<double> r);
	Vector<double> vecform(FullMatrix<double> A);
	void matform(FullMatrix<double> &A, Vector<double> Av); 
	void right(FullMatrix<double> &Aright,FullMatrix<double> elm);
	void symmf(FullMatrix<double> &A,FullMatrix<double> elm); 
	void left(FullMatrix<double> &Aleft,FullMatrix<double> elm);
	void ElasticProd(FullMatrix<double> &stress,FullMatrix<double> elm, FullMatrix<double> ElasticityTensor);
	void tracev(FullMatrix<double> &Atrace, FullMatrix<double> elm, FullMatrix<double> B);



	FullMatrix<double> F,F_tau,FP_tau,FE_tau,T,P;
	Tensor<4,dim,double> dP_dF;
	double No_Elem, N_qpts;

	//Store crystal orientations
	std::vector<std::vector<  Vector<double> > >  rot;
	std::vector<std::vector<  Vector<double> > >  rotnew;

	//Store history variables
	std::vector< std::vector< FullMatrix<double> > >   Fp_iter;
	std::vector< std::vector< FullMatrix<double> > > Fp_conv;
	std::vector< std::vector< FullMatrix<double> > >   Fe_iter;
	std::vector< std::vector< FullMatrix<double> > > Fe_conv;
	std::vector<std::vector<  Vector<double> > >  s_alpha_iter;
	std::vector<std::vector<  Vector<double> > >  s_alpha_conv;

	unsigned int n_slip_systems; //No. of slip systems
	FullMatrix<double> m_alpha,n_alpha,q,sres,Dmat;
	Vector<double> sres_tau;
	bool initCalled;

	//orientatations data for each quadrature point
	std::vector<std::vector<unsigned int> > quadratureOrientationsMap;  
	void loadOrientations();
};

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
	ellipticBVP<dim>::numPostProcessedFields=1;
	ellipticBVP<dim>::postprocessed_solution_names.push_back("post");
}

template <int dim>
void crystalPlasticity<dim>::loadOrientations(){
	QGauss<dim>  quadrature(quadOrder);
	const unsigned int num_quad_points = quadrature.size();
	FEValues<dim> fe_values (this->FE, quadrature, update_quadrature_points);
	//loop over elements
	typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
	for (; cell!=endc; ++cell) {
		if (cell->is_locally_owned()){
			quadratureOrientationsMap.push_back(std::vector<unsigned int>(num_quad_points,0));
			fe_values.reinit(cell);
			//loop over quadrature points
			for (unsigned int q=0; q<num_quad_points; ++q){
				double pnt[3];
				pnt[0]=fe_values.get_quadrature_points()[q][0];
				pnt[1]=fe_values.get_quadrature_points()[q][1];
				pnt[2]=fe_values.get_quadrature_points()[q][2];
				//get orientation ID and store it in quadratureOrientationsMap
				quadratureOrientationsMap.back()[q]=orientations.getMaterialID(pnt);
				//now one can ac orientations.getMaterialID(pnt) << std::endlcess the oreintation id for each quadrature point using quadratureOrientationsMap[cellID][q]
			}      
		}
	} 
}


template <int dim>
void crystalPlasticity<dim>::init(unsigned int num_quad_points)
{
	//call loadOrientations to load material orientations
	loadOrientations();

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

template <int dim>
void crystalPlasticity<dim>::calculatePlasticity(unsigned int cellID,
												 unsigned int quadPtID)
{

	F_tau=F; // Deformation Gradient
	FullMatrix<double> FE_t(dim,dim),FP_t(dim,dim);  //Elastic and Plastic deformation gradient 
	Vector<double> s_alpha_t(n_slip_systems); // Slip resistance
	Vector<double> rot1(dim);// Crystal orientation (Rodrigues representation)

	FE_t=Fe_conv[cellID][quadPtID];
	FP_t=Fp_conv[cellID][quadPtID];
	s_alpha_t=s_alpha_conv[cellID][quadPtID];
	rot1=rot[cellID][quadPtID];

	// Rotation matrix of the crystal orientation
	FullMatrix<double> rotmat(dim,dim);
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
	Vector<double> s_beta(n_slip_systems),h_beta(n_slip_systems);
	FullMatrix<double> h_alpha_beta_t(n_slip_systems,n_slip_systems),A(n_slip_systems,n_slip_systems);
	s_beta=s_alpha_t;

	// Single slip hardening rate
	for(unsigned int i=0;i<n_slip_systems;i++){
		h_beta=properties.h0*pow((1-s_beta(i)/properties.s_s),properties.a);
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

		T_star_tau.add(1.0,T_star_tau_trial);

		// % % % % % STEP 9 % % % % %

		temp.reinit(dim,dim);
		det_FE_tau=FE_tau.determinant();
		FE_tau.mmult(temp,T_star_tau); temp.equ(1.0/det_FE_tau,temp);
		temp.mTmult(T_tau,FE_tau);

		det_F_tau=F_tau.determinant();
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

	temp.reinit(dim,dim); T_tau.mTmult(temp,rotmat);
	rotmat.mmult(T_tau,temp);

	temp.reinit(dim,dim); P_tau.mTmult(temp,rotmat);
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
	FullMatrix<double> K_local(dofs_per_cell,dofs_per_cell);
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
		}

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
			}
		}
	}
	elementalJacobian = K_local;
	elementalResidual = Rlocal;
}


//implementation of the getElementalValues method
template <int dim>
void crystalPlasticity<dim>::updateAfterIncrement()
{
	reorient();

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
				orientations.addToOutputOrientations(temp);
			}      
			cellID++;
		}
	}
	orientations.writeOutputOrientations();

	//Update the history variables when convergence is reached for the current increment
	Fe_conv=Fe_iter;
	Fp_conv=Fp_iter;
	s_alpha_conv=s_alpha_iter;

	//fill in post processing field values
	cellID=0;
	cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
	for (; cell!=endc; ++cell) {
	  if (cell->is_locally_owned()){
	    //loop over quadrature points
	    for (unsigned int q=0; q<num_quad_points; ++q){
	      this->postprocessValues(cellID, q, 0, 0)=-1.9;
	    }
	    cellID++;
	  }
	}

	//call base class project() function to project post processed fields
	ellipticBVP<dim>::project();
}


template <int dim>
void crystalPlasticity<dim>::tracev(FullMatrix<double> &Atrace, FullMatrix<double> elm, FullMatrix<double> B) {

	Atrace.reinit(dim*dim,dim*dim);
	Vector<double> C(dim*dim);

	for(unsigned int i=0;i<9;i++){
		C(i)=elm(0,i)+elm(4,i)+elm(8,i);
	}

	for(unsigned int i=0;i<9;i++){
		for(unsigned int j=0;j<dim;j++){
			for(unsigned int k=0;k<dim;k++){
				Atrace(dim*j+k, i) =  B(j, k) * C(i);
			}
		}
	}

}  




template <int dim>
void crystalPlasticity<dim>::ElasticProd(FullMatrix<double> &stress,FullMatrix<double> elm, FullMatrix<double> ElasticityTensor) {  

	int index=0;
	stress.reinit(dim,dim);
	Vector<double> temp1(2*dim),temp2(2*dim);

	for(unsigned int i=0;i<dim;i++){

		for(unsigned int j=i;j<dim;j++){

			temp1(index) = elm(i, j);
			index = index+1;
		}
	}

	ElasticityTensor.vmult(temp2,temp1);

	index = 0;

	for(unsigned int i=0;i<dim;i++){

		for(unsigned int j=i;j<dim;j++){

			stress(i, j) = temp2(index);
			index = index+1;
		}
	}


	stress(1, 0) = stress(0, 1);
	stress(2, 0) = stress(0, 2);
	stress(2, 1) = stress(1, 2);

}






template <int dim>
void crystalPlasticity<dim>::right(FullMatrix<double> &Aright,FullMatrix<double> elm) {

	Aright.reinit(dim*dim,dim*dim);

	Aright=0.0;

	for(unsigned int i=0;i<dim;i++){

		for(unsigned int j=0;j<dim;j++){

			Aright[i][j]=elm(j,i) ;
			Aright[i+3][j+3]=elm(j,i) ;
			Aright[i+6][j+6]=elm(j,i) ;
		}
	}

}


template <int dim>
void crystalPlasticity<dim>::left(FullMatrix<double> &Aleft,FullMatrix<double> elm) {

	Aleft.reinit(dim*dim,dim*dim);

	Aleft=0.0;

	for(unsigned int i=0;i<dim;i++){

		for(unsigned int j=0;j<dim;j++){

			Aleft[dim*i][dim*j]=elm(i,j) ;
			Aleft[dim*i+1][dim*j+1]=elm(i,j) ;
			Aleft[dim*i+2][dim*j+2]=elm(i,j) ;
		}
	}

}


template <int dim>
void crystalPlasticity<dim>::symmf(FullMatrix<double> &A, FullMatrix<double> elm) {

	A.reinit(dim*dim,dim*dim);

	for(unsigned int i=0;i<9;i++){

		for(unsigned int j=0;j<9;j++){

			if (i == 1 || i == 3)
				A(i, j) = 0.5* (elm(1, j) + elm(3, j));
			else
				if (i == 2 || i == 6)
					A(i, j) = 0.5* (elm(2, j) + elm(6, j));
				else
					if( i == 5 || i == 7)
						A(i, j) = 0.5* (elm(5, j) + elm(7, j));
					else
						A(i, j) = elm(i, j);
		}
	}

}



template <int dim>
Vector<double> crystalPlasticity<dim>::vecform(FullMatrix<double> A) {

	Vector<double> Av(6);

	Av(0) =A[0][0];
	Av(1) =A[1][1];
	Av(2) =A[2][2];
	Av(3) =A[1][2];
	Av(4) =A[0][2];
	Av(5) =A[0][1];

	return Av;
}

template <int dim>
void crystalPlasticity<dim>::matform(FullMatrix<double> &A, Vector<double> Av) {

	A.reinit(dim,dim);

	A[0][0]=Av(0) ;
	A[1][1]=Av(1) ;
	A[2][2]=Av(2) ;
	A[1][2]=Av(3) ;
	A[0][2]=Av(4) ;
	A[0][1]=Av(5) ;
	A[2][1]=Av(3) ;
	A[2][0]=Av(4) ;
	A[1][0]=Av(5) ;

}

template <int dim>
void crystalPlasticity<dim>::odfpoint(FullMatrix <double> &OrientationMatrix,Vector<double> r) {


	//function OrientationMatrix = odfpoint(r)
	//%USAGE: [C] = odfpoint(R)
	//%TO OBTAIN ORIENTATION MATRICES FROM RODRIGUES FORM

	double rdotr = 0.0;

	for(unsigned int i=0;i<dim;i++){
		rdotr = rdotr + r(i)*r(i);
	}

	double term1 = 1.0 - (rdotr);
	double term2 = 1.0 + (rdotr);

	OrientationMatrix.reinit(dim,dim);OrientationMatrix = IdentityMatrix(dim);

	for(unsigned int i=0;i<dim;i++){

		OrientationMatrix[i][i]=OrientationMatrix[i][i]*term1;
	}

	for(unsigned int i=0;i<dim;i++){

		for(unsigned int j=0;j<dim;j++){
			OrientationMatrix[i][j] = OrientationMatrix[i][j] + 2*(r(i)*r(j));	
		}
	}

	OrientationMatrix[0][1] = OrientationMatrix[0][1]-2.0*r(2);
	OrientationMatrix[0][2] =  OrientationMatrix[0][2]+2.0*r(1);
	OrientationMatrix[1][2] = OrientationMatrix[1][2]-2.0*r(0);
	OrientationMatrix[1][0] =  OrientationMatrix[1][0]+2.0*r(2);
	OrientationMatrix[2][0] = OrientationMatrix[2][0]-2.0*r(1);
	OrientationMatrix[2][1] =  OrientationMatrix[2][1]+2.0*r(0);


	for(unsigned int i=0;i<dim;i++){

		for(unsigned int j=0;j<dim;j++){
			OrientationMatrix[i][j] = OrientationMatrix[i][j]/term2;	
		}
	}



}



template <int dim>
void crystalPlasticity<dim>::tangent_modulus(FullMatrix<double> &F_trial, FullMatrix<double> &Fpn_inv, FullMatrix<double> &SCHMID_TENSOR1, FullMatrix<double> &A,FullMatrix<double> &A_PA,FullMatrix<double> &B,FullMatrix<double> &T_tau, FullMatrix<double> &PK1_Stiff, Vector<double> &active, Vector<double> &resolved_shear_tau_trial, Vector<double> &x_beta, Vector<double> &PA, int &n_PA, double &det_F_tau, double &det_FE_tau ) {

	int iter=0,iter1=0,iter2=0,iter3=0;
	FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim);
	FullMatrix<double> Dmat2=Dmat;
	FullMatrix<double> ElasticityTensor(6,6), TM(9,9);
	Vector<double> vec1(6), vec2(9);

	Dmat2(3,3)=properties.C44; 	Dmat2(4,4)=properties.C44; 	Dmat2(5,5)=properties.C44;
	vec1(0)=0;vec1(1)=5;vec1(2)=4;vec1(3)=1;vec1(4)=3;vec1(5)=2;
	vec2(0)=0;vec2(1)=5;vec2(2)=4;vec2(3)=5;vec2(4)=1;vec2(5)=3;vec2(6)=4;vec2(7)=3;vec2(8)=2;

	for(unsigned int i=0;i<6;i++){
		for(unsigned int j=0;j<6;j++){
			ElasticityTensor[i][j]=Dmat2(vec1(i),vec1(j));
		}
	}


	for(unsigned int i=0;i<9;i++){
		for(unsigned int j=0;j<9;j++){
			TM[i][j]=Dmat2(vec2(i),vec2(j));
		}
	}

	FullMatrix<double> F_temp(dim,dim);
	F_temp=Fpn_inv;
	FullMatrix<double> scratch_1,scratch_2,EtF;
	right(scratch_1, F_temp);
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=0;j<dim;j++){
			F_temp[i][j]=F_trial[j][i];
		}
	}

	left(scratch_2,F_temp);
	F_temp.reinit(9,9);
	scratch_1.mmult(F_temp,scratch_2);
	symmf(EtF,F_temp);

	FullMatrix<double> s(dim,dim),p(n_PA,dim*dim),dgammadEmat(dim,dim), dgammadEtrial(n_PA,dim*dim),dgammadF(n_PA,dim*dim);
	Vector<double> P_as_vec3(dim*dim),s1(dim*dim),s2(dim*dim),temps2(dim*dim);

	if (n_PA > 0){

		s2=0.0;
		for(unsigned int alpha=0;alpha<n_PA;alpha++){
			iter=active(alpha);
			iter2=0;	
			for(unsigned int j=0;j<dim;j++){
				for(unsigned int k=0;k<dim;k++){
					s(j,k)=SCHMID_TENSOR1(dim*iter+j,k);
					P_as_vec3(iter2)=s(j,k);
					iter2++;
				}
			}

			TM.Tvmult(s1,P_as_vec3);
			if(resolved_shear_tau_trial(iter)<0){
				s1.equ(-1.0,s1);
			}
			s2=0.0;

			for(unsigned int beta=0;beta<n_PA;beta++){
				iter3=active(beta);

				temp3.reinit(dim,dim);
				temp.reinit(dim,dim);
				for(unsigned int k=0;k<dim;k++){
					for(unsigned int l=0;l<dim;l++){
						temp3(k,l)=SCHMID_TENSOR1(dim*iter3+k,l);
						temp(l,k)=temp3(k,l);
					}
				}

				temp1.reinit(dim*dim,dim*dim); temp2.reinit(dim*dim,dim*dim);
				right(temp1,temp3); left(temp2,temp);
				temp1.add(1.0,temp2);
				TM.mmult(temp2,temp1); temp2.Tvmult(temps2,P_as_vec3);

				if((resolved_shear_tau_trial(iter)<0.0)^(resolved_shear_tau_trial(iter3)<0.0))
					temps2.equ(-1.0*x_beta(iter3),temps2);
				else
					temps2.equ(x_beta(iter3),temps2);

				s2.add(temps2);

			}

			for(unsigned int index1=0;index1<dim*dim;index1++){
				p(alpha,index1)=s1(index1)-s2(index1);
			}

		}

		A_PA.reinit(n_PA,n_PA);
		for(unsigned int i=0;i<(n_PA);i++){

			for(unsigned int j=0;j<n_PA;j++){
				A_PA[i][j]=A[PA(i)][PA(j)];
			}
		}

		temp.reinit(n_PA,n_PA); temp.invert(A_PA);
		temp.mmult(dgammadEtrial,p); dgammadEtrial.mmult(dgammadF,EtF);

	}

	s.reinit(dim*dim,dim*dim); s=0.0;

	int alpha=0;
	if(n_PA>0){
		for(unsigned int i=0;i<n_PA;i++){
			alpha = active(i);
			iter=0;
			temp3.reinit(dim,dim);
			for(unsigned int j=0;j<dim;j++){
				for(unsigned int k=0;k<dim;k++){
					dgammadEmat[k][j]=dgammadEtrial[i][dim*j+k];
					temp3[j][k]=B[dim*alpha+j][k];

				}
			}

			temp1.reinit(dim*dim,dim*dim); right(temp1,dgammadEmat);
			temp.reinit(dim,dim);
			for(unsigned int j=0;j<dim;j++){
				for(unsigned int k=0;k<dim;k++){
					temp[j][k]=B[dim*alpha+j][k];
				}
			}

			temp3.reinit(dim,dim); temp3=0.0;	
			ElasticProd(temp3,temp,ElasticityTensor);	
			temp2.reinit(dim*dim,dim*dim); tracev(temp2,temp1,temp3);

			if(resolved_shear_tau_trial(alpha)>0){
				temp2.add(-1.5,temp2);	 
			}
			else{
				temp2.add(-0.5,temp2);	 
			}
			s.add(1.0,temp2);
		}
	} 

	FullMatrix<double> smat1(dim*dim,dim*dim),smat2(dim*dim,dim*dim), smat3(dim*dim,dim*dim);

	smat1=0.0; smat1.add(1.0,s); smat1.add(1.0,TM);
	smat3=0;

	if(n_PA>0){

		for(unsigned int i=0;i<n_PA;i++){
			alpha = active(i);
			F_temp.reinit(dim,dim);
			for(unsigned int j=0;j<dim;j++){
				for(unsigned int k=0;k<dim;k++){
					F_temp[j][k]=SCHMID_TENSOR1[dim*alpha+j][k];
					temp3[k][j]=F_temp[j][k];
				}
			}	

			right(temp1,F_temp);
			left(temp2,temp3);
			temp1.add(1.0,temp2);
			TM.mmult(temp2,temp1);
			if(resolved_shear_tau_trial(alpha)>0){
				temp2.equ(-1.0*x_beta(alpha),temp2); 
			}
			else{
				temp2.equ(1.0*x_beta(alpha),temp2);
			}
			smat3.add(1.0,temp2);
		}
	}


	FullMatrix<double> tangent_moduli(dim*dim,dim*dim),dgammadFmat(dim,dim);

	tangent_moduli=0.0; tangent_moduli.add(1.0,smat1); tangent_moduli.add(1.0,smat3);
	smat2=0.0;	

	if(n_PA>0){

		for(unsigned int i=0;i<n_PA;i++){
			alpha = active(i);
			temp2.reinit(dim,dim);
			for(unsigned int j=0;j<dim;j++){
				for(unsigned int k=0;k<dim;k++){
					dgammadFmat[k][j]=dgammadF[i][dim*j+k];
					temp2[j][k]=SCHMID_TENSOR1[dim*alpha+j][k];
				}
			}	

			right(temp1,dgammadFmat);
			F_trial.mmult(temp3,temp2);
			tracev(temp2,temp1,temp3);
			if(resolved_shear_tau_trial(alpha)>0){
				temp2.equ(-1.0,temp2); 
			}
			smat2.add(1.0,temp2);
		}		
	}

	temp2.reinit(dim,dim);
	smat1.reinit(dim*dim,dim*dim); 
	temp2.invert(FP_tau);
	right(smat1,temp2);


	FullMatrix<double> FeF(dim*dim,dim*dim);
	FeF=0.0; FeF.add(1.0,smat1); FeF.add(1.0,smat2);
	temp1.reinit(dim,dim); temp1.invert(FE_tau);
	F_temp.reinit(dim,dim);


	//term 1; -trace(dFe Fe_inv) T * (1/det(Fe)) = [Term1]{dF}
	F_temp=temp1;
	FullMatrix<double> C4(dim*dim,dim*dim);

	right(C4,F_temp);

	FullMatrix<double> Term1(dim*dim,dim*dim),Term2(dim*dim,dim*dim),Term3(dim*dim,dim*dim),Term24(dim*dim,dim*dim);
	tracev(Term3,C4,T_tau);
	Term3.mmult(Term1,FeF); Term1.equ(-1.0, Term1);

	//term 3; Fe * dT_bar * FeT * (1/det(Fe)) = [Term3]{dF}

	F_temp.reinit(dim,dim);
	F_temp=FE_tau;	

	FullMatrix<double> cauchy_Stiff(dim*dim,dim*dim);
	scratch_1.reinit(dim*dim,dim*dim); left(scratch_1,F_temp);
	temp1.reinit(dim,dim); temp1=IdentityMatrix(dim);

	temp2.reinit(dim,dim); temp2=F_temp;
	temp2.Tmmult(F_temp,temp1);

	right(C4,F_temp); C4.mmult(scratch_2,scratch_1);
	scratch_1.reinit(dim*dim,dim*dim); scratch_2.mmult(scratch_1,tangent_moduli);
	Term3.reinit(dim*dim,dim*dim); scratch_1.mmult(Term3,EtF);
	Term3.equ(1.0/det_FE_tau, Term3); 

	//term 2&4: (dFe * Fe_inv) * T + T * (dFe * Fe_inv)^T = 2Symm(dFe * Fe_inv * T) = [Term24]{dF}

	temp2.invert(FE_tau); temp2.mmult(F_temp,T_tau);
	right(C4,F_temp);
	temp3.reinit(dim*dim,dim*dim); temp3=C4;
	symmf(C4,temp3);
	Term24=0.0; C4.mmult(Term24,FeF); Term24.equ(2.0,Term24);
	cauchy_Stiff=0.0; cauchy_Stiff.add(1.0, Term1); cauchy_Stiff.add(1.0, Term24); cauchy_Stiff.add(1.0, Term3);
	temp1.reinit(dim,dim); temp1.invert(F_tau);

	C4.reinit(dim*dim,dim*dim); right(C4,temp1);
	Term1.reinit(dim*dim,dim*dim); tracev(Term1,C4,T_tau);

	scratch_1.reinit(dim*dim,dim*dim); right(scratch_1,F_temp);
	Term2.reinit(dim*dim,dim*dim); symmf(Term2,scratch_1);
	Term2.equ(2.0,Term2); Term2.add(-1.0,scratch_1); Term2.equ(-1.0,Term2);

	Term3.reinit(dim*dim,dim*dim); Term3=cauchy_Stiff;

	PK1_Stiff=0.0; PK1_Stiff.add(1.0,Term1); PK1_Stiff.add(1.0,Term2); PK1_Stiff.add(1.0,Term3);
	temp2.reinit(dim,dim); temp2=IdentityMatrix(dim);
	temp3.reinit(dim,dim); temp3=temp1; temp3.Tmmult(temp1,temp2);
	temp2.reinit(dim*dim,dim*dim); right(temp2,temp1);
	temp3.reinit(dim*dim,dim*dim); temp3=PK1_Stiff;
	temp3.mmult(PK1_Stiff,temp2);
	PK1_Stiff.equ(det_F_tau,PK1_Stiff);	

}





template <int dim>
void crystalPlasticity<dim>::inactive_slip_removal(Vector<double> &inactive,Vector<double> &active, Vector<double> &x_beta, int &n_PA, Vector<double> &PA, Vector<double> b,FullMatrix<double> A,FullMatrix<double> A_PA){

	FullMatrix<double> temp;
	int iter=0,iter1=0,iter2=0,iter3=0;
	Vector<double> x_beta1(n_PA), x_beta2(n_PA),b_PA(n_PA);
	double flag1=0;
	A_PA.reinit(n_PA,n_PA);
	for(unsigned int i=0;i<n_PA;i++){
		b_PA(i)=b(PA(i));      
		for(unsigned int j=0;j<n_PA;j++){
			A_PA[i][j]=A[PA(i)][PA(j)];
		}
	}    
	temp.reinit(n_PA,n_PA); temp.invert(A_PA);
	Vector<double> tempv3;	
	temp.vmult(x_beta1,b_PA);	

	// Determine the inactive slip ssytems
	for(unsigned int i=0;i<n_slip_systems;i++){
		for(unsigned int j=0;j<n_PA;j++){
			if((i) != PA(j)){
				iter++;
			}
		}
		if(iter==n_PA){
			inactive(iter2)=i;
			iter2++;
		}
		iter=0;

	}    

	x_beta2=x_beta1; x_beta1.reinit(n_slip_systems);
	x_beta1=0;
	for(unsigned int i=0;i<n_PA;i++){
		x_beta1(PA(i))=x_beta2(i);
	}

	Vector<double> row,PA_new,inactive2;
	double n_IA_new;

	// Continue the process till removal of all inactive slip systems
	while (flag1==0){

		iter=0;
		for(unsigned int i=0;i<n_slip_systems;i++){
			if(x_beta1(i)<0){
				iter++;
			}
		}

		if(iter>0){

			row.reinit(iter);
			iter=0;
			for(unsigned int i=0;i<n_slip_systems;i++){
				if(x_beta1(i)<0){
					row(iter)=i;
					iter++;
				}
			}

			n_IA_new=iter;
			PA_new.reinit(n_PA-n_IA_new);
			iter2=0; iter=0;
			for(unsigned int i=0;i<n_PA;i++){
				for(unsigned int j=0;j<n_IA_new;j++){
					if(PA(i) != row(j)){
						iter++;
					}
				}

				if(iter==n_IA_new){
					PA_new(iter2)=PA(i);
					iter2++;
				}
				iter=0;

			}

			PA.reinit(n_PA-n_IA_new); PA=PA_new;
			inactive2.reinit(n_slip_systems-n_PA); inactive2=inactive;
			inactive.reinit(n_slip_systems-n_PA+n_IA_new);

			for(unsigned int i=0;i<(n_slip_systems-n_PA+n_IA_new);i++){
				if(i<n_slip_systems-n_PA){
					inactive(i)=inactive2(i);

				}
				else{
					inactive(i)=row(i-(n_slip_systems-n_PA));

				}
			}

			A_PA.reinit(n_PA-n_IA_new,n_PA-n_IA_new); A_PA=0.0;
			b_PA.reinit(n_PA-n_IA_new); b_PA=0.0;


			for(unsigned int i=0;i<(n_PA-n_IA_new);i++){
				b_PA(i)=b(PA(i));      
				for(unsigned int j=0;j<(n_PA-n_IA_new);j++){
					A_PA[i][j]=A[PA(i)][PA(j)];
				}
			}
			temp.reinit(n_PA-n_IA_new,n_PA-n_IA_new); temp.invert(A_PA);


			x_beta1.reinit(n_PA-n_IA_new); temp.vmult(x_beta1,b_PA); 
			x_beta2.reinit(n_PA-n_IA_new); x_beta2=x_beta1;
			x_beta1.reinit(n_slip_systems);x_beta1=0.0;
			n_PA=n_PA-n_IA_new;

			for(unsigned int i=0;i<n_PA;i++){
				x_beta1(PA(i))=x_beta2(i);
			}

		}

		else
			flag1=1;

	}


	for(unsigned int i=0;i<n_slip_systems-n_PA;i++){
		for(unsigned int j=i+1;j<n_slip_systems-n_PA;j++){ 
			if(inactive(i)>inactive(j)){
				iter=inactive(i);  
				inactive(i)=inactive(j);	
				inactive(j)=iter;		
			}
		}
	}  


	for(unsigned int i=0;i<n_slip_systems;i++){
		if(x_beta1(i)<0)
			x_beta1(i)=0;

	}


	active.reinit(n_PA); active=PA;
	x_beta.reinit(n_slip_systems); x_beta=x_beta1;


}


template <int dim>
void crystalPlasticity<dim>::reorient() {
	//Update the history variables 

	LAPACKFullMatrix<double> C_old(dim,dim),C_new(dim,dim);

	Vector<double> eigenvalues(dim);
	FullMatrix<double> C_old_temp(dim,dim),C_new_temp(dim,dim),Fe_old(dim,dim), Fe_new(dim,dim), eigenvectors(dim,dim),Lambda(dim,dim),U_old(dim,dim),U_new(dim,dim),R_old(dim,dim),R_new(dim,dim),Omega(dim,dim),temp(dim,dim); 
	Lambda=IdentityMatrix(dim);
	Omega=0.0;
	Vector<double> rot1(dim),Omega_vec(dim),rold(dim),dr(dim),rnew(dim);
	FullMatrix<double> rotmat(dim,dim);
	//int itgno;
	unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();

	for (unsigned int i=0; i<num_local_cells; ++i) {
		for(unsigned int j=0;j<N_qpts;j++){

			C_old_temp=0.0;
			C_old=0.0;

			Fe_old=Fe_conv[i][j];
			Fe_new=Fe_iter[i][j];
			Fe_old.Tmmult(C_old_temp,Fe_old);  
			C_old=C_old_temp;
			C_old.compute_eigenvalues_symmetric(0.0,200000.0,1e-15,eigenvalues, eigenvectors);

			for(unsigned int k=0;k<dim;k++){
				Lambda(k,k)=sqrt(eigenvalues(k));
			}

			eigenvectors.mmult(U_old,Lambda);
			temp=U_old; temp.mTmult(U_old,eigenvectors);
			R_old.invert(U_old);
			temp=R_old; Fe_old.mmult(R_old,temp);
			Fe_new.Tmmult(C_new_temp,Fe_new);
			C_new=C_new_temp;
			C_new.compute_eigenvalues_symmetric(0.0,200000.0,1e-15,eigenvalues, eigenvectors);

			for(unsigned int k=0;k<dim;k++){
				Lambda(k,k)=sqrt(eigenvalues(k));
			}

			eigenvectors.mmult(U_new,Lambda);
			temp=U_new; temp.mTmult(U_new,eigenvectors);
			R_new.invert(U_new);
			temp=R_new; Fe_new.mmult(R_new,temp);

			Omega=0.0; Omega.add(1.0,R_new); Omega.add(-1.0,R_old);
			temp=Omega; temp.mTmult(Omega,R_new);


			rot1=rot[i][j];
			rold=rotnew[i][j];
			rotmat=0.0;
			odfpoint(rotmat,rot1);

			temp=Omega;
			temp.mTmult(Omega,rotmat);
			temp=Omega;
			rotmat.mmult(Omega,temp);


			Omega_vec(0)=-0.5*(Omega(1,2)-Omega(2,1));Omega_vec(1)=0.5*(Omega(0,2)-Omega(2,0));Omega_vec(2)=-0.5*(Omega(0,1)-Omega(1,0));

			double dot;
			dot=Omega_vec(0)*rold(0)+Omega_vec(1)*rold(1)+Omega_vec(2)*rold(2);
			Vector<double> cross(dim),dot_term(dim);

			dot_term(0)=dot*rold(0);dot_term(1)=dot*rold(1);dot_term(2)=dot*rold(2);

			cross(0)=Omega_vec(1)*rold(2)-Omega_vec(2)*rold(1);
			cross(1)=Omega_vec(2)*rold(0)-Omega_vec(0)*rold(2);
			cross(2)=Omega_vec(0)*rold(1)-Omega_vec(1)*rold(0);

			dr=0.0;	dr.add(1.0, Omega_vec); dr.add(1.0,dot_term); dr.add(1.0,cross); dr.equ(0.5,dr);

			rnew=0.0; rnew.add(1.0,rold); rnew.add(1.0,dr);


			rotnew[i][j]=rnew;


		}
	}


}


#endif
