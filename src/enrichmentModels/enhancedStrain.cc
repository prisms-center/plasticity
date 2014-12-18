//implementation of the enhanced strain enrichment model
/*Class to handle basis functions, etc. for enhanced strain
See "Improved versions of assumed enhanced strain," Simo & Armero*/
/*Although templated with dim, this is really only good for 3D at this point*/

#ifndef ENHANCEDSTRAIN_H
#define ENHANCEDSTRAIN_H


//dealii headers
#include "../../include/dealIIheaders.h"
//Static condensation
#include "../utilityObjects/staticCondensation.cc"

using namespace dealii;

template <int dim>
class  enhancedStrain{

template <int _dim>
friend class continuumPlasticity;

public:
  enhancedStrain(const FiniteElement<dim, dim> &fe, ConditionalOStream  pcout_temp);
  ~enhancedStrain();
 
private: 
  void reinit(Vector<double> Ulocal_temp, typename DoFHandler<dim>::active_cell_iterator elem);
  Tensor<1,dim,double> tilde_grad(unsigned int function_no, unsigned int point_no);
  Tensor<1,dim,double> bar_grad(unsigned int function_no, unsigned int point_no, Tensor<1,dim,double> alpha_4);
  Tensor<1,dim,double> center_grad(unsigned int function_no);
  void get_F_enh(unsigned int q, FullMatrix<double> &F);
  void create_block_mat_vec(FullMatrix<double> F, FullMatrix<double> tau, Tensor<4,dim,double> c_ep, unsigned int q);
		
  void updateAlpha(Vector<double> dUlocal, unsigned int cellID);
  unsigned int enhanced_system_to_component_index(unsigned int enhanced_elem_dof);
  unsigned int enhanced_dofs_per_cell();	
  unsigned int primary_enhanced_dofs_per_cell();

  void init_enh_dofs(unsigned int n_local_elems);
  Vector<double> Ulocal;
  typename DoFHandler<dim>::active_cell_iterator currentElem;
  QGauss<dim> center_quad;
  QGauss<dim>	quad_formula;
  FEValues<dim> fe_values;
  FEValues<dim> center_values;
  Vector<double> Alpha; //Global vectors
  std::vector<staticCondensation<8*dim,4*dim> > staticCondensationData; //Static condensation
  
  FullMatrix<double> Klocal, Glocal, Mlocal;
  Vector<double> Flocal, Hlocal;
  
  Tensor<2,dim> J_0, J_0inv;
  double j_0;
  bool enh_dofs_initialized;

  ConditionalOStream  *pcout;
};

template <int dim>
enhancedStrain<dim>::enhancedStrain(const FiniteElement<dim, dim> &fe, ConditionalOStream  pcout_temp)
	:
	center_quad(1),
	quad_formula(quadOrder),
	fe_values (fe, quad_formula, update_values   | update_gradients | update_quadrature_points | update_jacobians | update_JxW_values),
	center_values (fe, center_quad, update_values   | update_gradients | update_jacobians | update_quadrature_points | update_JxW_values)
	{
		pcout = &pcout_temp;

		unsigned int dofs_per_elem = fe.dofs_per_cell;
		unsigned int enh_dofs_per_elem = enhanced_dofs_per_cell();

		Klocal.reinit(dofs_per_elem,dofs_per_elem);
		Glocal.reinit(dofs_per_elem,enh_dofs_per_elem);
		Mlocal.reinit(enh_dofs_per_elem,enh_dofs_per_elem);
		Flocal.reinit(dofs_per_elem);
		Hlocal.reinit(enh_dofs_per_elem);

		enh_dofs_initialized = false;
}

template <int dim>
enhancedStrain<dim>::~enhancedStrain(){
	
}

template <int dim>
void enhancedStrain<dim>::init_enh_dofs(unsigned int n_local_elems){

		staticCondensationData.resize(n_local_elems);
		Alpha.reinit (4*dim*n_local_elems);
		Alpha=0;

		enh_dofs_initialized = true;
}

template <int dim>
void enhancedStrain<dim>::reinit(Vector<double> Ulocal_temp, typename DoFHandler<dim>::active_cell_iterator elem) {
	if(enh_dofs_initialized == false){
		*pcout << "Error: enhanced dofs not yet initialized.\n";
		exit(1);
	}

	currentElem = elem;
	Ulocal = Ulocal_temp; 

	fe_values.reinit (elem);
	center_values.reinit (elem);

	J_0 = center_values.jacobian(0);
	j_0 = determinant(J_0);
	J_0inv = invert(J_0);
}

template <int dim>
Tensor<1,dim,double> enhancedStrain<dim>::tilde_grad(unsigned int function_no, unsigned int point_no){
	//Enhanced strain gradient functions for the four incompatible shape functions

	if(enh_dofs_initialized == false){
		*pcout << "Error: enhanced dofs not yet initialized.\n";
		exit(1);
	}

	Tensor<1,dim,double> grad; grad = 0.;

	Tensor<2,dim,double> J_xi = fe_values.jacobian(point_no);
	double j_xi = determinant(J_xi);

	//Get the coordinates of the quadrature point in the unit domain
	Point<dim,double> xi = fe_values.get_quadrature().point(point_no);
	//Convert xi to the bi-unit domain
	for(unsigned int i=0; i<dim; i++){
		xi(i) = 2.*xi(i) - 1.;
	}

	if(function_no < 3*dim){	//	\tilde{N^i} = 0.5*(xi(i)^2 - 1), i = 0,1,2
		grad[function_no/dim] = xi(function_no/dim);
	}
	else if(function_no >= 3*dim && function_no < 4*dim){	//	\tilde{N^3} = xi(0)*xi(1)*xi(2)
		grad[0] = xi(1)*xi(2);
		grad[1] = xi(0)*xi(2);
		grad[2] = xi(0)*xi(1);
	}
	else{
		*pcout << "Error: invalid function number.\n";
		exit(1);
	}

	grad = j_0/j_xi*grad*J_0inv;
	
	return grad;
}

template <int dim>
Tensor<1,dim,double> enhancedStrain<dim>::bar_grad(unsigned int function_no, unsigned int point_no,Tensor<1,dim,double> alpha_4){
	//Enhanced strain gradient functions to replace the standard basis function gradients

	if(enh_dofs_initialized == false){
		*pcout << "Error: enhanced dofs not yet initialized.\n";
		exit(1);
	}

	Tensor<1,dim,double> grad;

	Tensor<2,dim> J_xi = fe_values.jacobian(point_no);
	double j_xi = determinant(J_xi);

	grad = (1. + alpha_4*tilde_grad(3*dim,point_no))*center_values.shape_grad(function_no,0);
	grad += j_0/j_xi*(fe_values.shape_grad(function_no,point_no) - center_values.shape_grad(function_no,0))*J_xi*J_0inv;

	return grad;
}

template <int dim>
Tensor<1,dim,double> enhancedStrain<dim>::center_grad(unsigned int function_no){
	return center_values.shape_grad(function_no,0);
}

template <int dim>
unsigned int enhancedStrain<dim>::enhanced_system_to_component_index(unsigned int enhanced_elem_dof){
	return enhanced_elem_dof % dim;
}

template <int dim>
unsigned int enhancedStrain<dim>::enhanced_dofs_per_cell(){
	return 4*dim;
}

template <int dim>
unsigned int enhancedStrain<dim>::primary_enhanced_dofs_per_cell(){
	return 3*dim;
}

template <int dim>
void enhancedStrain<dim>::get_F_enh(unsigned int q, FullMatrix<double> &F){

	if(enh_dofs_initialized == false){
		*pcout << "Error: enhanced dofs not yet initialized.\n";
		exit(1);
	}

	unsigned int dofs_per_elem = fe_values.get_fe().dofs_per_cell;
	unsigned int prim_enh_dofs_per_elem = primary_enhanced_dofs_per_cell();

	//Create vector of original location of each dof
	MappingQ1<dim,dim> mapping;
	std::vector<double>	nodalCoords (dofs_per_elem);
	for (unsigned int d=0; d<dofs_per_elem; ++d){
		unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
		nodalCoords[d] = mapping.transform_unit_to_real_cell(currentElem,fe_values.get_fe().unit_support_point(d))[i];
	}

	//define F_tilde
	Tensor<2,dim,double> F_tilde; F_tilde = 0.0;
	for (unsigned int A=0; A<prim_enh_dofs_per_elem; ++A){
		unsigned int i = enhanced_system_to_component_index(A);
		F_tilde[i] += Alpha(4*dim*currentElem->user_index()+A)*tilde_grad(A,q);
	}

	Tensor<1,dim,double> alpha_4;
	for (unsigned int i=0; i<dim; ++i){
		alpha_4[i] = Alpha(4*dim*currentElem->user_index()+9+i);
	}

	//define F_bar, F_0
	Tensor<2,dim,double> F_bar; F_bar = 0.0;
	for (unsigned int d=0; d<dofs_per_elem; ++d){
		unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
		F_bar[i] += (Ulocal[d] + nodalCoords[d])*bar_grad(d,q,alpha_4);
	}

	//define F_enhanced = GRADbar_{X}[phi] + GRADtilde_{X}[phi]
	Tensor<2,dim,double> F_enh; F_enh = 0.0;
	F_enh = F_bar + F_tilde;

	F.copy_from(F_enh);
}

template <int dim>
void enhancedStrain<dim>::create_block_mat_vec(FullMatrix<double> F, FullMatrix<double> tau, Tensor<4,dim,double> c_ep, unsigned int q){

	if(enh_dofs_initialized == false){
		*pcout << "Error: enhanced dofs not yet initialized.\n";
		exit(1);
	}

	Klocal = 0.; Glocal = 0.; Mlocal = 0.;
	Flocal = 0.; Hlocal = 0.;

	Tensor<2,dim,double> F_enh, Fenh_inv, tau_enh;
	unsigned int dofs_per_elem = fe_values.get_fe().dofs_per_cell;
	unsigned int prim_enh_dofs_per_elem = primary_enhanced_dofs_per_cell();
	unsigned int enh_dofs_per_elem = enhanced_dofs_per_cell();

	F.copy_to(F_enh);
	tau.copy_to(tau_enh);
	Fenh_inv = invert(F_enh);

	Tensor<1,dim,double> alpha_4;
	for (unsigned int i=0; i<dim; ++i){
		alpha_4[i] = Alpha(4*dim*currentElem->user_index()+9+i);
	}

	//Create vector of original location of each dof
	MappingQ1<dim,dim> mapping;
	std::vector<double>	nodalCoords (dofs_per_elem);
	Tensor<2,dim,double> F_0; F_0 = 0.0;
	for (unsigned int d=0; d<dofs_per_elem; ++d){
		unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
		nodalCoords[d] = mapping.transform_unit_to_real_cell(currentElem,fe_values.get_fe().unit_support_point(d))[i];
		F_0[i] += (Ulocal[d] + nodalCoords[d])*center_grad(d);
	}

	//Transform enhanced gradients to current configuration
	std::vector<Tensor<1,dim,double> > bar_del(dofs_per_elem),
																		 center_del(dofs_per_elem),
																		 tilde_del(prim_enh_dofs_per_elem);
	for (unsigned int d=0; d<dofs_per_elem; ++d){
		bar_del[d] = bar_grad(d,q,alpha_4)*Fenh_inv;
		center_del[d] = center_grad(d)*Fenh_inv;
	}
	for (unsigned int d=0; d<prim_enh_dofs_per_elem; d++){
		tilde_del[d] = tilde_grad(d,q)*Fenh_inv;
	}
	Tensor<2,dim,double> F0_Finv;
	F0_Finv = F_0*Fenh_inv;

	//Evaluate H = F_intTilde_elem
	for (unsigned int d=0; d<enh_dofs_per_elem; d++){
		unsigned int i = enhanced_system_to_component_index(d);
		if(d<prim_enh_dofs_per_elem){
			Hlocal(d) += tau_enh[i]*tilde_del[d];
		}
		else{
			Hlocal(d) += double_contract(F0_Finv,tau_enh)*tilde_grad(d,q)[i];
		}
	}

	//Evaluate F = F_intBar_elem
	for (unsigned int d=0; d<dofs_per_elem; ++d){
		unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
		Flocal(d) += tau_enh[i]*bar_del[d];
	}

	//Calculate M (material and geometric parts)
	for (unsigned int I=0; I<enh_dofs_per_elem; ++I){
		unsigned int j = enhanced_system_to_component_index(I);
		for (unsigned int J=0; J<enh_dofs_per_elem; ++J){
			unsigned int k = enhanced_system_to_component_index(J);
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int l=0; l<dim; ++l){
					if(I < prim_enh_dofs_per_elem && J < prim_enh_dofs_per_elem){
						Mlocal[I][J] += tilde_del[I][i]*(tau_enh[i][l]*(j==k) + c_ep[i][j][k][l])*tilde_del[J][l];
					}
					else if(I < prim_enh_dofs_per_elem){
						Mlocal[I][J] += (tilde_del[I][i]*tau_enh[i][l]*F0_Finv[j][l] + F0_Finv[i][l]*c_ep[i][l][j]*tilde_del[I])*
														tilde_grad(J,q)[k];
					}
				}
				if(I < prim_enh_dofs_per_elem && J >= prim_enh_dofs_per_elem){
					Mlocal[J][I] = Mlocal[I][J];
				}
			}
			if(I >= prim_enh_dofs_per_elem && J >= prim_enh_dofs_per_elem){
				Mlocal[I][J] += (contract3(F0_Finv,c_ep,F0_Finv) + double_contract(F0_Finv,tau_enh*F0_Finv))*
												tilde_grad(I,q)[j]*tilde_grad(J,q)[k];
			}
		}
	}

	//Fill in local K matrix (both material and geometric parts)
	for (unsigned int I=0; I<dofs_per_elem; ++I){
		unsigned int j = fe_values.get_fe().system_to_component_index(I).first;
		for (unsigned int J=0; J<dofs_per_elem; ++J){
			unsigned int k = fe_values.get_fe().system_to_component_index(J).first;
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int l=0; l<dim; ++l){
					Klocal[I][J] +=  bar_del[I][i]*(c_ep[i][j][k][l] + tau_enh[i][l]*(j==k))*bar_del[J][l];
				}
			}
		}
	}

	//Fill in local G matrix (both material and geometric parts)
	for (unsigned int I=0; I<dofs_per_elem; ++I){
		unsigned int j = fe_values.get_fe().system_to_component_index(I).first;
		for (unsigned int J=0; J<enh_dofs_per_elem; ++J){
			unsigned int k = enhanced_system_to_component_index(J);
			for (unsigned int i=0; i<dim; ++i){
				for (unsigned int l=0; l<dim; ++l){
					if(J<prim_enh_dofs_per_elem){
						Glocal[I][J] += bar_del[I][i]*(c_ep[i][j][k][l] + tau_enh[i][l]*(j==k))*tilde_del[J][l];
					}
					else{
						Glocal[I][J] += (F0_Finv[i][l]*c_ep[i][l][j]*bar_del[I] + (F0_Finv[j][i]*bar_del[I][l] + (j==i)*center_del[I][l])*tau_enh[i][l])*
										tilde_grad(J,q)[k];
					}
				}
			}
		}
	}

}

//Update the Alpha vector
template <int dim>
void enhancedStrain<dim>::updateAlpha(Vector<double> dUlocal, unsigned int cellID){

	if(enh_dofs_initialized == false){
		*pcout << "Error: enhanced dofs not yet initialized.\n";
		exit(1);
	}

	//loop over elements
	const unsigned int dofs_per_elem = fe_values.get_fe().dofs_per_cell;
	for (unsigned int n=0; n<dofs_per_elem; ++n){
		staticCondensationData[cellID].d(n) = dUlocal[n];
	}

	//Compute delta_Alpha
	staticCondensationData[cellID].recover();

	//Add delta_Alpha to Alpha
	for(unsigned int n=0; n<4*dim; n++){
		Alpha(4*dim*cellID+n) += staticCondensationData[cellID].s(n);
	}
}


#endif
