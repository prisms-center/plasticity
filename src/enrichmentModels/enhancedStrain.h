//implementation of the enhanced strain enrichment model
/*Class to handle basis functions, etc. for enhanced strain
  See "Improved versions of assumed enhanced strain," Simo & Armero*/
/*Although templated with dim, this is really only good for 3D at this point*/

#ifndef ENHANCEDSTRAIN_H
#define ENHANCEDSTRAIN_H


//dealii headers
#include "../../include/dealIIheaders.h"
#include "../../include/userInputParameters.h"

//Static condensation
#include "../utilityObjects/staticCondensation.cc"

using namespace dealii;

template <int dim>
class  enhancedStrain{

  template <int _dim>
  friend class continuumPlasticity;

public:
  /**
   *Class constructor. Takes as inputs the deal.II FiniteElement object used in the current problem,
   *and the parallel cout object.
   */
  enhancedStrain(const FiniteElement<dim, dim> &fe, ConditionalOStream  pcout_temp, userInputParameters _userInputs);
  /**
   *Class destructor.
   */
  ~enhancedStrain();
 
private: 
  /**
   *A one-time initialization based on the number of elements covered by the current processor.
   *Initializes the Alpha vector (global enhanced dof vector) and the static condensation object.
   */
  void init_enh_dofs(unsigned int n_local_elems);
  /**
   *Reinitilize the enhanced strain object for the current element. Takes the local displacment
   *vector and the deal.II element iterator for the current element.
   */
  void reinit(Vector<double> Ulocal_temp, typename DoFHandler<dim>::active_cell_iterator elem);
  /**
   *Enhanced strain gradient functions for the four incompatible shape functions
   */
  Tensor<1,dim,double> tilde_grad(unsigned int function_no, unsigned int point_no);
  /**
   *Enhanced strain gradient functions to replace the standard basis function gradients
   */
  Tensor<1,dim,double> bar_grad(unsigned int function_no, unsigned int point_no, Tensor<1,dim,double> alpha_4);
  /**
   *Standard strain gradient functions evaluated at the center of the element.
   */
  Tensor<1,dim,double> center_grad(unsigned int function_no);
  /**
   *Function to retrieve the enhanced deformation gradient, for a given quadrature point
   *and standard deformation gradient.
   */
  void get_F_enh(unsigned int q, FullMatrix<double> &F);
  /**
   *Function that creates the contribution of the specified quadrature point to the local K, M, and G matrices and
   *F and H vectors. K and F are associated with the standard (nodal) dofs, M and H are associated with the 
   *enhanced (interior) dofs, and G containes the cross terms. The deformation gradient, Kirchhoff stress, and 
   *elastoplastic tangent are needed, as well as the quadrature point number.
   */
  void create_block_mat_vec(FullMatrix<double> F, FullMatrix<double> tau, Tensor<4,dim,double> c_ep, unsigned int q);
  /**
   *Update the enhanced dofs (Alpha) after the standard dofs have been solved.
   */
  void updateAlpha(Vector<double> dUlocal, unsigned int cellID);
  /**
   *Function to convert from element level enhanced dof number (0-11) to the vector component (0-2).
   *Remember, the 12 enhanced dofs in an element are organized as 4 vectors with 3 components each (in 3D).
   */
  unsigned int enhanced_system_to_component_index(unsigned int enhanced_elem_dof);
  /**
   *Returns the number of enhanced dofs per element (4 interior dof (Alpha) vectors,
   *for a total 12 enhanced dofs in 3D).
   */
  unsigned int enhanced_dofs_per_cell();
  /**
   *Retrns the number of primary enhanced dofs per element (3 Alpha vectors, for total of
   *9 primary enhanced dofs in 3D). The fourth Alpha vector is treated a little differently.
   */
  unsigned int primary_enhanced_dofs_per_cell();

  /**
   *Local displacement vector.
   */
  Vector<double> Ulocal;
  /**
   *Deal.II element iterator with information about the current element.
   */
  typename DoFHandler<dim>::active_cell_iterator currentElem;
  /**
   *Deal.II quadrature object, to be used with a single quardrature point at the center of the element.
   */
  QGauss<dim> center_quad;
  /**
   *Deal.II quadrature object, to be used with standard Gaussian quadrature.
   */
  QGauss<dim>	quad_formula;
  /**
   *Deal.II object with information about the basis functions, quadrature points, etc. for the current element.
   */
  FEValues<dim> fe_values;
  /**
   *Similar to fe_values, but for the single quadrature point at the center of the element.
   */
  FEValues<dim> center_values;
  /**
   *Global vector of the enhanced degrees of freedom.
   */
  Vector<double> Alpha;
  /**
   *Static condensation object used to condense out the enhanced dofs (which are interior to the element).
   */
  std::vector<staticCondensation<8*dim,4*dim> > staticCondensationData;
  /**
   *Local matrices used to form the local jacobian matrix. K is associated with the standard (nodal) dofs,
   *M is associated with the enhanced (interior) dofs, and G containes the cross terms.
   */
  FullMatrix<double> Klocal, Glocal, Mlocal;
  /**
   *Local vectors used to form the local residual vector. F is associated with the standard (nodal) dofs
   *and H is associated with the enhanced (interior) dofs.
   */
  Vector<double> Flocal, Hlocal;
  /**
   *The jacobian matrix and its inverse at the center of the element (jacobian of the isoparametric mapping).
   */
  Tensor<2,dim> J_0, J_0inv;
  Tensor<1,dim,double> HlocalT;
  /**
   *The determinant of the jacobian matrix at the center of the element (jacobian of the isoparametric mapping).
   */
  double j_0;
  /**
   *Marker to track if the init_enh_dofs function was called.
   */
  bool enh_dofs_initialized;
  /**
   *Parallel cout object
   */
  ConditionalOStream  *pcout;
};

template <int dim>
enhancedStrain<dim>::enhancedStrain(const FiniteElement<dim, dim> &fe, ConditionalOStream  pcout_temp, userInputParameters _userInputs)
  :
  center_quad(1),
  quad_formula(_userInputs.quadOrder),
  fe_values (fe, quad_formula, update_values   | update_gradients | update_quadrature_points | update_jacobians | update_JxW_values),
  center_values (fe, center_quad, update_values   | update_gradients | update_jacobians | update_quadrature_points | update_JxW_values)
{
  pcout = &pcout_temp;

  //Retrieve number of standard and enhanced dofs per element
  unsigned int dofs_per_elem = fe.dofs_per_cell;
  unsigned int enh_dofs_per_elem = enhanced_dofs_per_cell();

  //Resize matrices and vectors used to form the local jacobian and residual 
  Klocal.reinit(dofs_per_elem,dofs_per_elem);
  Glocal.reinit(dofs_per_elem,enh_dofs_per_elem);
  Mlocal.reinit(enh_dofs_per_elem,enh_dofs_per_elem);
  Flocal.reinit(dofs_per_elem);
  Hlocal.reinit(enh_dofs_per_elem);

  //Initialize the marker to false
  enh_dofs_initialized = false;
}

template <int dim>
enhancedStrain<dim>::~enhancedStrain(){
	
}

template <int dim>
void enhancedStrain<dim>::init_enh_dofs(unsigned int n_local_elems){

  //After the mesh has been created in the ellipticBVP class, we take the number of elements
  //local to this processor to resize the static condensation object and the global
  //enhanced degree of freedom vector.
  staticCondensationData.resize(n_local_elems);
  Alpha.reinit (4*dim*n_local_elems);
  Alpha=0;

  //Mark that this function has been called.
  enh_dofs_initialized = true;
}

template <int dim>
void enhancedStrain<dim>::reinit(Vector<double> Ulocal_temp, typename DoFHandler<dim>::active_cell_iterator elem) {
  if(enh_dofs_initialized == false){
    *pcout << "Error: enhanced dofs not yet initialized.\n";
    exit(1);
  }

  //Update deal.II element iterator to the specified element
  currentElem = elem;
  //Update local displacement vector for current element
  Ulocal = Ulocal_temp; 

  //Update basis functions and quadrature points for current element
  fe_values.reinit (elem);
  center_values.reinit (elem);

  //Update center point jacobian of the isoparametric mapping
  J_0 = center_values.jacobian(0);
  j_0 = determinant(J_0);
  J_0inv = invert(J_0);
}

template <int dim>
Tensor<1,dim,double> enhancedStrain<dim>::tilde_grad(unsigned int function_no, unsigned int point_no){
  //Enhanced strain gradient functions for the four incompatible shape functions

  //Ensure the initiliaziation function has been called
  if(enh_dofs_initialized == false){
    *pcout << "Error: enhanced dofs not yet initialized.\n";
    exit(1);
  }

  //Tensor holding values of the specified enhanced strain gradient function at the given quad. point
  Tensor<1,dim,double> grad; grad = 0.;
  //Jacobian at quadrature point
  Tensor<2,dim,double> J_xi = fe_values.jacobian(point_no);
  double j_xi = determinant(J_xi);

  //Get the coordinates of the quadrature point in the unit domain
  Point<dim,double> xi = fe_values.get_quadrature().point(point_no);
  //Convert xi to the bi-unit domain (deal.II uses a unit domain)
  for(unsigned int i=0; i<dim; i++){
    xi(i) = 2.*xi(i) - 1.;
  }

  //Define enhanced shape function gradients (see paper cited at top)
  if(function_no < 3*dim){
    //\tilde{N^i} = 0.5*(xi(i)^2 - 1), i = 0,1,2, from functions listed 
    //in eqn (3.8) of the cited paper
    grad[function_no/dim] = xi(function_no/dim);
  }
  else if(function_no >= 3*dim && function_no < 4*dim){
    //\tilde{N^3} = xi(0)*xi(1)*xi(2), see equation (3.17) from cited paper
    grad[0] = xi(1)*xi(2);
    grad[1] = xi(0)*xi(2);
    grad[2] = xi(0)*xi(1);
  }
  else{
    *pcout << "Error: invalid function number.\n";
    exit(1);
  }

  //See the 2nd line of equation (4.1) of the cited paper
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

  //Hold enhanced strain gradient function value
  Tensor<1,dim,double> grad;

  //Jacobian at given quad. point
  Tensor<2,dim> J_xi = fe_values.jacobian(point_no);
  double j_xi = determinant(J_xi);

  //Modified basis function gradient to replace standard basis function gradient
  //See line 3 of equation (4.1) of paper cited at top, as well as (2.17).
  //Note that deal.II's basis function gradients are all w.r.t. the actual domain,
  //not the bi-unit or unit domain.
  grad = (1. + alpha_4*tilde_grad(3*dim,point_no))*center_grad(function_no);
  grad += j_0/j_xi*(fe_values.shape_grad(function_no,point_no) - center_grad(function_no))*J_xi*J_0inv;

  return grad;
}

template <int dim>
Tensor<1,dim,double> enhancedStrain<dim>::center_grad(unsigned int function_no){
  //Compute gradient at the center quadrature point
  return center_values.shape_grad(function_no,0);
}

template <int dim>
unsigned int enhancedStrain<dim>::enhanced_system_to_component_index(unsigned int enhanced_elem_dof){
  //Go from element dof number (0-11) to vector component (0-2)
  return enhanced_elem_dof % dim;
}

template <int dim>
unsigned int enhancedStrain<dim>::enhanced_dofs_per_cell(){
  //12 total enhanced dofs per element in 3D
  return 4*dim;
}

template <int dim>
unsigned int enhancedStrain<dim>::primary_enhanced_dofs_per_cell(){
  //The fourth Alpha vector is treated a little differently, so there are 9
  //"primary" enhancd dofs in 3D.
  return 3*dim;
}

template <int dim>
void enhancedStrain<dim>::get_F_enh(unsigned int q, FullMatrix<double> &F){

  if(enh_dofs_initialized == false){
    *pcout << "Error: enhanced dofs not yet initialized.\n";
    exit(1);
  }

  //Get the number of standard dofs and primary enhanced dofs per element
  unsigned int dofs_per_elem = fe_values.get_fe().dofs_per_cell;
  unsigned int prim_enh_dofs_per_elem = primary_enhanced_dofs_per_cell();

  //Create vector of original location of each standard dof in the element
  //This location is the x, y, or z location based on the orientation of the dof.
  MappingQ1<dim,dim> mapping;
  std::vector<double>	nodalCoords (dofs_per_elem);
  for (unsigned int d=0; d<dofs_per_elem; ++d){
    //Check the orientation of the standard dof (displacement in the x, y, or z direction)
    unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
    //Get the corresponding component of the position of the node in the initial configuration,
    //i.e. if the dof is displacement in the y direction at a particular node, store the y component
    //of that node's location. 
    nodalCoords[d] = mapping.transform_unit_to_real_cell(currentElem,fe_values.get_fe().unit_support_point(d))[i];
  }

  //define F_tilde (based on enhanced dofs)
  //see the 2nd half of equation (4.2) of the paper
  Tensor<2,dim,double> F_tilde; F_tilde = 0.0;
  for (unsigned int A=0; A<prim_enh_dofs_per_elem; ++A){
    unsigned int i = enhanced_system_to_component_index(A);
    F_tilde[i] += Alpha(4*dim*currentElem->user_index()+A)*tilde_grad(A,q);
  }

  //The fourth alpha vector
  Tensor<1,dim,double> alpha_4;
  for (unsigned int i=0; i<dim; ++i){
    alpha_4[i] = Alpha(4*dim*currentElem->user_index()+9+i);
  }

  //define F_bar, F_0 (based on standard dofs)
  //see the 1st half of equation (4.2) of the paper
  Tensor<2,dim,double> F_bar; F_bar = 0.0;
  for (unsigned int d=0; d<dofs_per_elem; ++d){
    unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
    //F is the derivative of (displacement + original position = deformed position)
    F_bar[i] += (Ulocal[d] + nodalCoords[d])*bar_grad(d,q,alpha_4);
  }

  //define F_enhanced = GRADbar_{X}[phi] + GRADtilde_{X}[phi]
  //equation (4.2) of the paper
  Tensor<2,dim,double> F_enh; F_enh = 0.0;
  F_enh = F_bar + F_tilde;

  //Copy F_enh into the user's F matrix 
  F.copy_from(F_enh);
}

template <int dim>
void enhancedStrain<dim>::create_block_mat_vec(FullMatrix<double> F, FullMatrix<double> tau, Tensor<4,dim,double> c_ep, unsigned int q){
  //Refer to the paper cited at the top

  if(enh_dofs_initialized == false){
    *pcout << "Error: enhanced dofs not yet initialized.\n";
    exit(1);
  }

  Klocal = 0.; Glocal = 0.; Mlocal = 0.;
  Flocal = 0.; Hlocal = 0.; HlocalT = 0.;

  Tensor<2,dim,double> F_enh, Fenh_inv, tau_enh;
  unsigned int dofs_per_elem = fe_values.get_fe().dofs_per_cell;
  unsigned int prim_enh_dofs_per_elem = primary_enhanced_dofs_per_cell();
  unsigned int enh_dofs_per_elem = enhanced_dofs_per_cell();

  //Copy inputs to the class variables
  F.copy_to(F_enh);
  tau.copy_to(tau_enh);
  Fenh_inv = invert(F_enh);

  //Retrieve the 4th Alpha vector for this element
  Tensor<1,dim,double> alpha_4;
  for (unsigned int i=0; i<dim; ++i){
    alpha_4[i] = Alpha(4*dim*currentElem->user_index()+9+i);
  }

  //Create vector of original location of each standard dof in the element
  //This location is the x, y, or z location based on the orientation of the dof.
  MappingQ1<dim,dim> mapping;
  std::vector<double>	nodalCoords (dofs_per_elem);
  Tensor<2,dim,double> F_0; F_0 = 0.0;
  for (unsigned int d=0; d<dofs_per_elem; ++d){
    //Check the orientation of the standard dof (displacement in the x, y, or z direction)
    unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
    //Get the corresponding component of the position of the node in the initial configuration,
    //i.e. if the dof is displacement in the y direction at a particular node, store the y component
    //of that node's location. 
    nodalCoords[d] = mapping.transform_unit_to_real_cell(currentElem,fe_values.get_fe().unit_support_point(d))[i];
    //F is the derivative of (displacement + original position = deformed position)
    F_0[i] += (Ulocal[d] + nodalCoords[d])*center_grad(d);
  }

  //Transform enhanced gradients to current configuration
  //see equation (4.3) of the cited paper
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
      //see line 2 of equation (4.9)
      Hlocal(d) += tau_enh[i]*tilde_del[d];
    }
    else{
      //see line 2 of equation (4.9)
      Hlocal(d) += scalar_product(F0_Finv,tau_enh)*tilde_grad(d,q)[i]; //double_contract wants symmetricTensor
    }
  }

  //Evaluate F = F_intBar_elem
  for (unsigned int d=0; d<dofs_per_elem; ++d){
    unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
    //see line 1 of equation (4.9)
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
	    //See line 3 of equations (4.14) and (4.18)
	    Mlocal[I][J] += tilde_del[I][i]*(tau_enh[i][l]*(j==k) + c_ep[i][j][k][l])*tilde_del[J][l];
	  }
	  else if(I < prim_enh_dofs_per_elem){
	    //See also line 2 of equation (4.19)
	    Mlocal[I][J] += (tilde_del[I][i]*tau_enh[i][l]*F0_Finv[j][l] + F0_Finv[i][l]*c_ep[i][l][j]*tilde_del[I])*
	      tilde_grad(J,q)[k];
	  }
	}
	if(I < prim_enh_dofs_per_elem && J >= prim_enh_dofs_per_elem){
	  Mlocal[J][I] = Mlocal[I][J];
	}
      }
      //See also line 3 of equation (4.19)
      if(I >= prim_enh_dofs_per_elem && J >= prim_enh_dofs_per_elem){
	Mlocal[I][J] += (contract3(F0_Finv,c_ep,F0_Finv) + scalar_product(F0_Finv,tau_enh*F0_Finv))*
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
	  //See line 1 of equations (4.14) and (4.18)
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
	    //See line 2 of equations (4.14) and (4.18)
	    Glocal[I][J] += bar_del[I][i]*(c_ep[i][j][k][l] + tau_enh[i][l]*(j==k))*tilde_del[J][l];
	  }
	  else{
	    //See also line 1 of equation (4.19)
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

  //loop over elements, pass the standard dofs (local displacment) to static condensation
  const unsigned int dofs_per_elem = fe_values.get_fe().dofs_per_cell;
  for (unsigned int n=0; n<dofs_per_elem; ++n){
    staticCondensationData[cellID].d(n) = dUlocal[n];
  }

  //Compute delta_Alpha useing static condensation object
  staticCondensationData[cellID].recover();

  //Add delta_Alpha to Alpha to update the enhanced dofs
  for(unsigned int n=0; n<4*dim; n++){
    Alpha(4*dim*cellID+n) += staticCondensationData[cellID].s(n);
  }
}


#endif
