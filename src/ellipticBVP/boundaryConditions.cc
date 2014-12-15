//methods to apply Dirichlet boundary conditons 

#ifndef BOUNDARYCONDITIONS_H
#define BOUNDARYCONDITIONS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//methods to apply dirichlet BC's
template <int dim>
void ellipticBVP<dim>::applyDirichletBCs(){
  //default method to apply zero Dirichlet BC's on all components
  //of this field, given by implicitFieldIndex, on all boundary faces
  constraints.clear ();
  constraints.reinit (locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (dofHandler, constraints);
  std::vector<bool> zeroBoundaries (dim, true);
  VectorTools::interpolate_boundary_values (dofHandler,
					    0, 
					    ZeroFunction<dim>(dim),
					    constraints,
					    zeroBoundaries);
  constraints.close ();
}

#endif
