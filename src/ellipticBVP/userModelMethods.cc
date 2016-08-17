//Methods related to the user model functionality
#ifndef USERMODEL_SRC_H
#define USERMODEL_SRC_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

#ifdef enableUserModel
template <int dim>
void ellipticBVP<dim>::initQuadHistory(){
  QGauss<dim>  quadrature(quadOrder);
  const unsigned int   num_quad_points = quadrature.size();
  //initialize the quadHistory table
  quadHistory.reinit(TableIndices<3> (triangulation.n_locally_owned_active_cells(), num_quad_points, numQuadHistoryVariables));
}
#endif

#endif
