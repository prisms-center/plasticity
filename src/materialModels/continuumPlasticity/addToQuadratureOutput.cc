#include "../../../include/continuumPlasticity.h"

template <int dim>
void continuumPlasticity<dim>::addToQuadratureOutput(std::vector<double>& _QuadOutputs){
	outputQuadrature.push_back(_QuadOutputs);
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
