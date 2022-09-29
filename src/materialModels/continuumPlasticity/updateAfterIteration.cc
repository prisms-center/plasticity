#include "../../../include/continuumPlasticity.h"
//implementation of the getElementalValues method
template <int dim>
void continuumPlasticity<dim>::updateAfterIteration()
{
  //After solving for the nodal values, calculate the enhanced degrees of freedom.
  std::vector<unsigned int> local_dof_indices(this->FE.dofs_per_cell);
  Vector<double> dUlocal(this->FE.dofs_per_cell);
  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(),
    endc = this->dofHandler.end();
  //this->pcout<<"continuum::update after iteration\n";
  //Loop over elements
  for (; cell!=endc; ++cell) {
    //For the current processor
    if (cell->is_locally_owned()){
      //Update local to global dof map
      cell->get_dof_indices(local_dof_indices);
      //Update local displacement vector
      for(unsigned int i=0; i<dUlocal.size(); i++){
	dUlocal[i] = this->solutionIncWithGhosts[local_dof_indices[i]];
      }
      //Update the enhanced degrees of freedom at each iteration
      enhStrain.updateAlpha(dUlocal, cell->user_index());
    }
  }

}

#include "../../../include/continuumPlasticity_template_instantiations.h"
