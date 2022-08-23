#include "../../../include/continuumPlasticity.h"

template <int dim>
void continuumPlasticity<dim>::updateBeforeIteration()
{
    //Currently, the initial iterations are hardcoded to elasticity for indentation problems
    if (this->userInputs.enableIndentationBCs) {
        if (this->currentIteration < 4) {
            pausePlastic = true;
            this->pcout
                    << " Plasticity disabled for iteration " << this->currentIteration << " \n";
        } else {
            pausePlastic = false;
            this->pcout
                    << " Plasticity enabled for iteration " << this->currentIteration << " \n";
        }
    }
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
