#include "../../../include/continuumPlasticity.h"

template <int dim>
void continuumPlasticity<dim>::updateBeforeIncrement()
{
  microvol=0.0;
  local_strain = 0.0;
  local_stress = 0.0;
  local_microvol = 0.0;
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
