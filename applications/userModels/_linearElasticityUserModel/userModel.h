//implementation of the Saint Venant–Kirchhoff elastic material model as a user model
#ifndef USERMODEL_H
#define USERMODEL_H

//headers
#include <deal.II/base/table.h>
#include "../../../include/ellipticBVP.h"

//Saint Venant–Kirchhoff elastic material model as a user model
//Derives from ellipticBVP base abstract class
template <int dim>
class userModel : public ellipticBVP<dim>
{
 public:
  //class constructor.
  userModel();
 private:
  //method to provide quadrature level residual and jacobian
  void getQuadratureValues(unsigned int elementID,
			   unsigned int numElemDofs,
			   unsigned int* componentIndices,
			   double* shapeValues,
			   double* shapeGrads,
			   double* F,
			   double* residual,
			   double* jacobian,
			   double* history); 
  //method to specify Dirichlet boundary conditions
  void setBoundaryValues(const Point<dim>& X, const unsigned int component, bool& bcFlag, double& bcValue);

  //define material properties
  double lambda, mu;
  Table<4,double> C;
};

//constructor
template <int dim>
userModel<dim>::userModel(): ellipticBVP<dim>(),  C(dim,dim,dim,dim){
  //set the number of quadrature point history variables required
  this->numQuadHistoryVariables=3; 

  //Elasticity tensor C_{ijkl}
  lambda=lame_lambda; mu=lame_mu;
  for (unsigned int i=0; i<dim; ++i){
    for (unsigned int j=0; j<dim; ++j){
      for (unsigned int k=0; k<dim; ++k){
	for (unsigned int l=0; l<dim; ++l){
	  C[i][j][k][l] = lambda*(i==j)*(k==l) + mu*((i==k)*(j==l)+(i==l)*(j==k));
	}
      }
    }
  }
}

//implementation of the user model
template <int dim>
void userModel<dim>::getQuadratureValues(unsigned int elementID,
					 unsigned int numElemDofs,
					 unsigned int* componentIndices,
					 double* shapeValues,
					 double* shapeGrads,
					 double* F,
					 double* residual,
					 double* jacobian,
					 double* history){
  //evaluate quadrature point residual
  for (unsigned int d1=0; d1<numElemDofs; ++d1) {
    residual[d1]=0.0; //no forcing terms
  }
  
  //evaluate quadrature point jacobian: K_{ij} = N_{i,k}*C_{ikjl}*N_{j,l}
  for (unsigned int d1=0; d1<numElemDofs; ++d1) {
    unsigned int i = componentIndices[d1];
    for (unsigned int d2=0; d2<numElemDofs; ++d2) {
      jacobian[d1*numElemDofs+d2]=0.0;
      unsigned int j = componentIndices[d2];
      for (unsigned int k = 0; k < dim; k++){
	for (unsigned int l= 0; l< dim; l++){
	  jacobian[d1*numElemDofs+d2] += shapeGrads[d1*dim+k]*C(i,k,j,l)*shapeGrads[d2*dim+l];
	}
      }
    }
  }

  //To store the quadrature point history variables, 
  //use the history[] array whose size is set by 
  //numQuadHistoryVariables in the constructor above.
}

#endif
