//implementation of the continuum plasticity material model
#ifndef CONTINUUMPLASTICITY_H
#define CONTINUUMPLASTICITY_H

//dealii headers
#include "../../../include/ellipticBVP.h"
#include "../../../src/enrichmentModels/enhancedStrain.cc"

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
  void getElementalValues(FEValues<dim>& fe_values,
			  unsigned int dofs_per_cell,
			  unsigned int num_quad_points,
			  FullMatrix<double>& elementalJacobian,
			  Vector<double>&     elementalResidual);
  void updateAfterIteration();
  void updateAfterIncrement();

  //quadrature level data structures
  FullMatrix<double> F, F_inv, b_e, invCp_TR, invCp, S, tau, Identity;
  Tensor<4,dim,double> c;
  double alpha, alpha_TR;
  //bool init_called;
  enhancedStrain<dim> enhStrain;

  //Store history variables
  std::vector< std::vector< FullMatrix<double> > > histInvCP_conv;
  std::vector< std::vector< FullMatrix<double> > > histInvCP_iter;
  std::vector< std::vector< double > > histAlpha_conv;
  std::vector< std::vector< double > > histAlpha_iter;
  unsigned int plasticOnset; //Marker to show when plasticity first occurs
  bool firstRun; //Marker for the first iteration on the first step
  std::vector<double> varSE, varYF, varIH; //Vector of variables/parameter for the strain energy, yield function, and isotropic hardening
  PRISMS::PFunction< std::vector<double>, double> strain_energy;
  PRISMS::PFunction< std::vector<double>, double> yield;
  PRISMS::PFunction< std::vector<double>, double> harden;
};

//constructor
template <int dim>
ContinuumPlasticity<dim>::ContinuumPlasticity() : 
ellipticBVP<dim>()
{}

//implementation of the getElementalValues method
template <int dim>
void ContinuumPlasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
							unsigned int dofs_per_cell,
							unsigned int num_quad_points,
							FullMatrix<double>& elementalJacobian,
							Vector<double>&     elementalResidual)
{
  //loop over quadrature points
  for (unsigned int q_point=0; q_point<num_quad_points;++q_point){
    for (unsigned int i=0; i<dofs_per_cell; ++i){
      const unsigned int component_i = this->FE.system_to_component_index(i).first;
      for (unsigned int j=0; j<dofs_per_cell; ++j){
	const unsigned int component_j = this->FE.system_to_component_index(j).first;

	// klocal= grad(N)*C_ijkl*grad(N)                                                                                                                          
	elementalJacobian(i,j) +=
	  (
	   (fe_values.shape_grad(i,q_point)[component_i] *
	    fe_values.shape_grad(j,q_point)[component_j] *
	    lambda_value)
	   +
	   (fe_values.shape_grad(i,q_point)[component_j] *
	    fe_values.shape_grad(j,q_point)[component_i] *
	    mu_value)
	   +
	   ((component_i == component_j) ?
	    (fe_values.shape_grad(i,q_point) *
	     fe_values.shape_grad(j,q_point) *
	     mu_value)  :
	    0)
	   )
	  *
	  fe_values.JxW(q_point);
      }
      elementalResidual(i) +=0.0;
    }
  }
}

#endif
