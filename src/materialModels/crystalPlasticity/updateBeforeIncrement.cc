#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::updateBeforeIncrement()
{
    microvol=0.0;
    //call base class project() function to project post processed fields
    //ellipticBVP<dim>::project();
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
