//implementation of the continuum plasticity material model
#ifndef CONTINUUMPLASTICITY_H
#define CONTINUUMPLASTICITY_H


#include "../src/enrichmentModels/enhancedStrain.h"
#include "../utils/IntegrationTools/PFunction.hh"
#include "ellipticBVP.h"
#include "../src/materialModels/continuumPlasticity/models/PLibrary.hh"


typedef struct {
  double lambda, mu, tau_y, K, H;
  bool stopOnConvergenceFailure;
  std::string yieldModel, strainEnergyModel, isoHardeningModel;
} materialProperties;

//material model class for continuum plasticity
//derives from ellipticBVP base abstract class
template <int dim>
class continuumPlasticity : public ellipticBVP<dim>
{
 public:
  /**
   *continuumPlasticity class constructor.
   */
  continuumPlasticity(userInputParameters & _userInputs);
  /**
   *Structure to hold the material parameters and model names.
   */
  materialProperties properties;

 private:
  /**
  *Initialize and resize class data structures.
  */
  void initcont(unsigned int num_quad_points);
  //void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  //void mesh();
  //void markBoundaries();
  //void applyDirichletBCs();
  /**
   *Update the plastic variables. For a given element and quadrature point,
   *this takes the proposed deformation gradient for the current iteration
   *and the plastic variables for the previous increment to calculate the
   *current values for the plastic variables, as well as the stress and tangent modulus.
   */
  void calculatePlasticity(unsigned int cellID, unsigned int quadPtID);
  void getElementalValues(FEValues<dim>& fe_values, unsigned int dofs_per_cell, unsigned int num_quad_points, FullMatrix<double>& elementalJacobian, Vector<double>&elementalResidual);
  void updateAfterIteration();
  void updateAfterIncrement();
  void updateBeforeIteration();
  void writeQuadratureOutput(std::string _outputDirectory, unsigned int _currentIncrement);
  void addToQuadratureOutput(std::vector<double>& _QuadOutputs);
  void updateBeforeIncrement();

  /**
   *Deformation gradient tensor
   */
  FullMatrix<double> F;
  /**
   *Kirchhoff stress tensor
   */
  FullMatrix<double> tau;
  /**
   *Cauchy stress tensor
   */
  FullMatrix<double> T;

  /**
   *Tangent modulus
   */
  Tensor<4,dim,double> c;
  /**
   *Instantian of the enhanced strain class, used to prevent element locking.
   */
  enhancedStrain<dim> enhStrain;

  /**
  * volume weighted Cauchy stress per core
  */
  FullMatrix<double> local_stress;

  /**
  * volume weighted Lagrangian strain per core
  */
  FullMatrix<double> local_strain;


  /**
  * volume averaged global Cauchy stress
  */
  FullMatrix<double> global_stress;

  /**
  * volume averaged global Lagrangian strain
  */
  FullMatrix<double> global_strain;

  double local_microvol,microvol;

  std::vector< std::vector< FullMatrix<double> > > CauchyStress;

  std::vector<std::vector<double> > outputQuadrature;


  /**
   *Converged results of invCP (the inverse of the plastic right Cauchy-Green tensor)
   *for the previous increment, indexed by element number and quadrature point number.
   */
  std::vector< std::vector< FullMatrix<double> > > histInvCP_conv;
  /**
   *Stores values of invCP (the inverse of the plastic right Cauchy-Green tensor)
   *for the most recent iteration, indexed by element number and quadrature point number.
   *Once the solution for the current increment converges, this information will be
   *transferred to "histInvCP_conv".
   */
  std::vector< std::vector< FullMatrix<double> > > histInvCP_iter;
  /**
   *Converged results of xi (deviatoric back stress)
   *for the previous increment, indexed by element number and quadrature point number.
   */
  std::vector< std::vector< Vector<double> > > histXi_conv;
  /**
   *Stores values of xi (deviatoric back stress)
   *for the most recent iteration, indexed by element number and quadrature point number.
   *Once the solution to the current increment converges, this information will be
   *transferred to "histXi_conv".
   */
  std::vector< std::vector< Vector<double> > > histXi_iter;
  /**
   *Converged results of alpha (the equivalent plastic strain)
   *for the previous increment, indexed by element number and quadrature point number.
   */
  std::vector< std::vector< double > > histAlpha_conv;
  /**
   *Stores values of alpha (the equivalent plastic strain)
   *for the most recent iteration, indexed by element number and quadrature point number.
   *Once the solution to the current increment converges, this information will be
   *transferred to "histAlpha_conv".
   */
  std::vector< std::vector< double > > histAlpha_iter;
  /**
   *Store the von Mises stress for each element/quadrature point to project when the
   *solution for the current increment converges.
   */
  std::vector< std::vector< double > > projectVonMisesStress;
  /**
   *Marker to show when plasticity first occurs.
   */
  bool plasticOnset;
  /**
   *Marker to show if the function "init" has been called.
   */
  bool initCalled;
  /**
  *Marker to show if the plasticity should be neglected in the current iteration.
  */
  bool pausePlastic;
  /**
   *Vector of variables/parameters for the elastic strain energy density function (pfunction).
   */
  std::vector<double> varStrainEnergy;
  /**
   *Vector of variables/parameters for the yield function (pfunction).
   */
  std::vector<double> varYield;
  /**
   *Vector of variables/parameters for the isotropic hardening function (pfunction).
   */
  std::vector<double> varIsoHardening;
  /**
   *Pfunction structure for the elastic strain energy density function.
   */
  PRISMS::PFunction< std::vector<double>, double> strain_energy;
  /**
   *Pfunction structure for the yield function.
   */
  PRISMS::PFunction< std::vector<double>, double> yield;
  /**
   *Pfunction structure for the isotropic hardening function.
   */
  PRISMS::PFunction< std::vector<double>, double> harden;
};

#endif
