#include "../../../../include/crystalPlasticity.h"
//constructor
template <int dim>
crystalPlasticity<dim>::crystalPlasticity(userInputParameters _userInputs) :
ellipticBVP<dim>(_userInputs),
F(dim,dim),
F_tau(dim,dim),
FP_tau(dim,dim),
FE_tau(dim,dim),
T(dim,dim),
P(dim,dim)
{
    initCalled = false;

    //post processing
    ellipticBVP<dim>::numPostProcessedFields=3;
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Grain_ID");
}
