#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::addToQuadratureOutput(std::vector<double>& _QuadOutputs)
{
	outputQuadrature.push_back(_QuadOutputs);
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
