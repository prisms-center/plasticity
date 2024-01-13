#include "../../../include/continuumPlasticity.h"

template <int dim>
void continuumPlasticity<dim>::initcont(unsigned int num_quad_points)
{
  // this->pcout << "continuumPlasticity<dim>::initcont "<<"\n";
  //Get the total numbers of elements for this processor (num_local_cells)
  unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();
  //Initiate the enhanced strain object with the number elements
  enhStrain.init_enh_dofs(num_local_cells);

  //Resize the deformation gradient and Kirchhoff stress tensors
  F.reinit(dim, dim);
  tau.reinit(dim, dim);
  T.reinit(dim, dim);

  local_strain.reinit(dim,dim);
  local_stress.reinit(dim,dim);
  global_strain.reinit(dim,dim);
  global_stress.reinit(dim,dim);
  local_strain=0.0;
  local_stress=0.0;
  global_strain=0.0;
  global_stress=0.0;

  FullMatrix<double> CauchyStress_init(dim,dim);
  CauchyStress_init=0;

  //Checkout the elastic strain energy density function from the pfunction library.
  //These are located in the "models" folder. First check that an available
  //model was requested.
  if(properties.strainEnergyModel == "quadlog" ||
     properties.strainEnergyModel == "neohook" ||
     properties.strainEnergyModel == "stvenkir"){
    PRISMS::PLibrary::checkout(properties.strainEnergyModel, strain_energy);
  }
  else{
    this->pcout << "Using user defined strain energy density function. "<<properties.strainEnergyModel<<"\n";
    PRISMS::PLibrary::checkout(properties.strainEnergyModel, strain_energy);
  }

  /*strain energy function inputs:
    0: lambda, "First Lame parameter"
    1: mu, "Second Lame parameter"
    2: lambda1, "First principle stretch"
    3: lambda2, "Second principle stretch"
    4: lambda3, "Third principle stretch"*/

  //Resize the vector of parameters/variables used in the strain energy function
  varStrainEnergy.resize(5,1.);
  //Specify the material parameters used in the strain energy function
  varStrainEnergy[0] = properties.lambda;
  varStrainEnergy[1] = properties.mu;

  //For now, specify the hardening model here and check it out from the
  //pfunction library (NOTE: this calculates the value for "q", the conjugate
  // stress-like quantity, used in the yield function).
  //See, for example, the end of section 2.1 of the formulation.
  if(properties.isoHardeningModel == "linear_hardening"){
    PRISMS::PLibrary::checkout(properties.isoHardeningModel, harden);
  }
  else{
    this->pcout << "Using user defined isotropic function.\n";
    PRISMS::PLibrary::checkout(properties.isoHardeningModel, harden);
  }

  /*hardening function inputs:
    0: alpha, "Equivalent plastic strain"
    1: K, "Hardening parameter"*/

  //Resize the vector of parameters/variables used in the hardening function
  varIsoHardening.resize(2,0.);
  //Specify the material parameter used in the hardening function
  varIsoHardening[1] = properties.K;

  //Checkout the yield function from the pfunction library.
  //This is located in the "models" folder. First check that an available
  //model was requested.
  //See equation (16) of the formulation.
  if(properties.yieldModel == "von_mises"){
    PRISMS::PLibrary::checkout(properties.yieldModel, yield);
  }
  else{
    this->pcout << "Using user defined yield function.\n";
    PRISMS::PLibrary::checkout(properties.yieldModel, yield);
  }

  /*yield function inputs:
    0: beta1, "First principle stress"
    1: beta2, "Second principle stress"
    2: beta3, "Third principle stress"
    3: tau_y, "Yield stress"
    4: q, "Conjugate stress-like quantity"*/

  //Resize the vector of parameters/variables used in the yield function
  varYield.resize(8,0.);
  //Specify the material parameter used in the yield function
  varYield[3] = properties.tau_y;

  //Resize the vectors of history variables according to the number of elements and quadrature points
	Vector<double> zero_vec(dim); zero_vec = 0.;
  histInvCP_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
  histAlpha_conv.resize(num_local_cells,std::vector<double>(num_quad_points,0));
  histXi_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,zero_vec));
  histInvCP_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
  histAlpha_iter.resize(num_local_cells,std::vector<double>(num_quad_points,0));
  histXi_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,zero_vec));
  CauchyStress.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,CauchyStress_init));

  //Resize the vector of vector used to store the von Mises stress
  projectVonMisesStress.resize(num_local_cells,std::vector<double>(num_quad_points,0));

  //Initialize projection of post process parameters???
  //ellipticBVP<dim>::initProjection();
  //Initialize the "plasticOnset" marker
  plasticOnset = false;
  //Change "initCalled" to "true" to show that the "init" function has now been called
  initCalled = true;
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
