//incomplete implementation of the elastic (Saint Venant-Kirchhoff) material model
//more elastic material models on their way
#ifndef ELASTIC_H
#define ELASTIC_H

//dealii headers
#include "../../../include/ellipticBVP.h"

typedef struct {
double lambda, mu;
std::string strainEnergyModel;
} materialProperties;


//
//material model class for elastic (Saint Venant-Kirchhoff) material model
//derives from ellipticBVP base abstract class
template <int dim>
class elasticity : public ellipticBVP<dim>
{
public:
  elasticity();
	materialProperties properties;
private:
  void markBoundaries();
  void applyDirichletBCs();
  void getElementalValues(FEValues<dim>& fe_values,
			  unsigned int dofs_per_cell,
			  unsigned int num_quad_points,
			  FullMatrix<double>& elementalJacobian,
			  Vector<double>&     elementalResidual);
};

//constructor
template <int dim>
elasticity<dim>::elasticity() : 
ellipticBVP<dim>()
{}

//implementation of the getElementalValues virtual method
template <int dim>
void elasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
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
	    properties.lambda)
	   +
	   (fe_values.shape_grad(i,q_point)[component_j] *
	    fe_values.shape_grad(j,q_point)[component_i] *
	    properties.mu)
	   +
	   ((component_i == component_j) ?
	    (fe_values.shape_grad(i,q_point) *
	     fe_values.shape_grad(j,q_point) *
	     properties.mu)  :
	    0)
	   )
	  *
	  fe_values.JxW(q_point);
      }
      elementalResidual(i) += - (0);
    }
  }
}

#endif
