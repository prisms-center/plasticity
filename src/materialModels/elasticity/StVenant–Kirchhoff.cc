//implementation of the Saint Venant-Kirchhoff elastic material model
#ifndef STVENANTKIRSHHOFF_H
#define STVENANTKIRCHHOFF_H

//dealii headers
#include "../../../include/ellipticBVP.h"

//
//material model class for Saint Venant-Kirchhoff elastic model
//derives from ellipticBVP base abstract class
template <int dim>
class StVenantKirchhoff_Elastic : public ellipticBVP<dim>
{
public:
  StVenantKirchhoff_Elastic();
private:
  void getElementalValues(FullMatrix<double>& elementalJacobian,
			  Vector<double>&     elementalResidual);
};

//constructor
template <int dim>
StVenantKirchhoff_Elastic<dim>::StVenantKirchhoff_Elastic() : 
ellipticBVP<dim>()
{}

//implementation of the getElementalValues virtual method
template <int dim>
void StVenantKirchhoff_Elastic<dim>::getElementalValues(FullMatrix<double>& elementalJacobian,
							Vector<double>&     elementalResidual)
{
}

#endif
