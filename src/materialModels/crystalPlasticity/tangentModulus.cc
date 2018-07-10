#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::tangent_modulus(FullMatrix<double> &F_trial, FullMatrix<double> &Fpn_inv, FullMatrix<double> &SCHMID_TENSOR1, FullMatrix<double> &A,FullMatrix<double> &A_PA,FullMatrix<double> &B,FullMatrix<double> &T_tau, FullMatrix<double> &PK1_Stiff, Vector<double> &active, Vector<double> &resolved_shear_tau_trial, Vector<double> &x_beta, Vector<double> &PA, unsigned int &n_PA, double &det_F_tau, double &det_FE_tau ) {

    
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
