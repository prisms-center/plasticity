#include "../../../include/continuumPlasticity.h"
//implementation of the getElementalValues method
template <int dim>
void continuumPlasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
						  unsigned int dofs_per_cell,
						  unsigned int num_quad_points,
						  FullMatrix<double>& elementalJacobian,
						  Vector<double>&     elementalResidual)
{

  //Initialized history variables and pfunction variables if unititialized
  if(initCalled == false){
    initcont(num_quad_points);
  }
  //std::cout<<"sanity check in continuumPlasticity<dim>::getElementalValues\n";
  //Get the element number.
  unsigned int cellID = fe_values.get_cell()->user_index();
  //Reset the vectors and matrices in static condensation to zero.
  //This is necessary because we are using an enhanced strain (see the paper
  //cited in the formulation).
  enhStrain.staticCondensationData[cellID].reset();

  //Vector relating local to global degree of freedom numbers
  std::vector<unsigned int> local_dof_indices(dofs_per_cell);
  //Local displacement vector
  Vector<double> Ulocal(dofs_per_cell);

  //Initiate element iterator
  typename DoFHandler<dim>::active_cell_iterator cell(& this->triangulation,
						      fe_values.get_cell()->level(),
						      fe_values.get_cell()->index(),
						      & this->dofHandler);
  //Pass in the element number
  cell->set_user_index(fe_values.get_cell()->user_index());
  //Populate the local to global dof vector
  cell->get_dof_indices (local_dof_indices);

  //Populate the locacal displacement vector
  for(unsigned int i=0; i<dofs_per_cell; i++){
    Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
  }
  //Initialize the enhanced strain object with this information.
  enhStrain.reinit(Ulocal, cell);

  //loop over quadrature points
  for (unsigned int q=0; q<num_quad_points; ++q){
    //Get enhanced deformation gradient
    enhStrain.get_F_enh(q, F);

    //Update strain, stress, and tangent for current time step/quadrature point
    calculatePlasticity(cellID, q);

    //Update block matrices and vectors in enhanced strain
    enhStrain.create_block_mat_vec(F, tau, c, q);

    //Pass local matrices and vectors to static condensation
    //Original dofs matrix
    enhStrain.staticCondensationData[cellID].K.add(enhStrain.fe_values.JxW(q),enhStrain.Klocal);
    //Cross terms matrix
    enhStrain.staticCondensationData[cellID].G.add(enhStrain.fe_values.JxW(q),enhStrain.Glocal);
    //Enhanced dofs matrix
    enhStrain.staticCondensationData[cellID].M.add(enhStrain.fe_values.JxW(q),enhStrain.Mlocal);
    //Original dofs vector
    enhStrain.staticCondensationData[cellID].F.add(enhStrain.fe_values.JxW(q),enhStrain.Flocal);
    //Enhanced dofs vector
    enhStrain.staticCondensationData[cellID].H.add(enhStrain.fe_values.JxW(q),enhStrain.Hlocal);
  }
  //Perform static condensation
  enhStrain.staticCondensationData[cellID].staticCondense();
  //Retrieve the element level jacobian
  elementalJacobian = enhStrain.staticCondensationData[cellID].K2;
  //Retrieve the element level residual
  elementalResidual.equ(-1.,enhStrain.staticCondensationData[cellID].F2);
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
