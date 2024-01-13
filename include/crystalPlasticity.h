//implementation of the crystal plasticity material model for FCC crystal structure
#ifndef CRYSTALPLASTICITY_H
#define CRYSTALPLASTICITY_H

//dealii headers
#include "ellipticBVP.h"

typedef struct {
  FullMatrix<double> m_alpha,n_alpha, eulerAngles2;
} materialProperties;
//material model class for crystal plasticity
//derives from ellipticBVP base abstract class
template <int dim>
class crystalPlasticity : public ellipticBVP<dim>
{
public:
  /**
  *crystalPlasticity class constructor.
  */
  crystalPlasticity(userInputParameters & _userInputs);
  /**
  *calculates the texture of the deformed polycrystal
  */
  void reorient();
  void reorient2(Vector<double> &rnew, Vector<double> rold, FullMatrix<double> FE_tau, FullMatrix<double> FE_t);
  /**
  *Initiation of Multiphase in calculatePlasticity.cc
  */
  void multiphaseInit(unsigned int cellID, unsigned int quadPtID);

  /**
  *calculates the material tangent modulus dPK1/dF at the quadrature point
  F_trial-trial Elastic strain (Fe_trial)
  Fpn_inv -inverse of plastic strain (Fp^-1)
  SCHMID_TENSOR1 (Schmid Tensor S)
  A->Stiffness Matrix A
  A-PA-> Stiffness matrix of active slip systems
  B->Refer to Eqn (5)
  T_tau-Cauchy Stress
  PK1Stiff-dP/dF -Tangent Modulus
  active-set of active slip systems
  resolved_shear_tau_trial- trial resolved shear stress
  x_beta-incremental shear strain delta_gamma
  PA-active slip systems
  det_F_tau-det(F_tau)
  det_FE_tau-det(FE_tau)
  */
  void tangent_modulus(FullMatrix<double> &F_trial, FullMatrix<double> &Fpn_inv, FullMatrix<double> &SCHMID_TENSOR1, FullMatrix<double> &A,FullMatrix<double> &A_PA,FullMatrix<double> &B,FullMatrix<double> &T_tau, FullMatrix<double> &PK1_Stiff, Vector<double> &active, Vector<double> &resolved_shear_tau_trial, Vector<double> &x_beta, Vector<double> &PA, unsigned int &n_PA, double &det_F_tau, double &det_FE_tau);
  /**
  *calculates the incremental shear strain in the slip systems by active set search and removal of inactive slip systems
  A->Stiffness Matrix A
  A-PA-> Stiffness matrix of active slip systems
  active-set of active slip systems
  inactive-set of inactive slip systems
  n_PA-No. of active slip systems
  x_beta-incremental shear strain delta_gamma
  PA-active slip systems
  */
  void inactive_slip_removal(Vector<double> &active,Vector<double> &x_beta_old, Vector<double> &x_beta, unsigned int &n_PA, unsigned int &n_Tslip_systems_Region, Vector<double> &PA, Vector<double> b,FullMatrix<double> A,FullMatrix<double> &A_PA);
  /**
  * Structure to hold material parameters
  */
  materialProperties properties;
  //orientation maps
  crystalOrientationsIO<dim> orientations;

private:
  void init(unsigned int num_quad_points);
  void init2(unsigned int num_quad_points);
  void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value);
  /**
  * Updates the stress and tangent modulus at a given quadrature point in a element for
  * the given constitutive model. Takes the deformation gradient at the current nonlinear
  * iteration and the elastic and plastic deformation gradient of the previous time step to
  *calculate the stress and tangent modulus
  */
  void calculatePlasticity(unsigned int cellID, unsigned int quadPtID, unsigned int StiffnessCalFlag);
  void getElementalValues(FEValues<dim>&fe_values, unsigned int dofs_per_cell, unsigned int num_quad_points, FullMatrix<double>&elementalJacobian, Vector<double>&elementalResidual);
  void updateAfterIncrement();
  void updateBeforeIteration();
  void updateBeforeIncrement();
  void writeQuadratureOutput(std::string _outputDirectory, unsigned int _currentIncrement);
  void addToQuadratureOutput(std::vector<double>& _QuadOutputs);
  /**
  *calculates the rotation matrix (OrientationMatrix) from the rodrigues vector (r)
  */
  void odfpoint(FullMatrix <double> &OrientationMatrix,Vector<double> r);
  /**
  *calculates the vector form (Voigt Notation) of the symmetric matrix A
  */
  Vector<double> vecform(FullMatrix<double> A);
  /**
  *calculates the symmetric matrix (A) from the vector form Av (Voigt Notation)
  */
  void matform(FullMatrix<double> &A, Vector<double> Av);
  /**
  *calculates the equivalent matrix Aright for the second order tensorial operation XA=B => A_r*{x}={b}
  */
  void right(FullMatrix<double> &Aright,FullMatrix<double> elm);
  /**
  *calculates the equivalent matrix A for the second order tensorial operation symm(AX)=B => A_r*{x}={b}
  */
  void symmf(FullMatrix<double> &A,FullMatrix<double> elm);
  /**
  *calculates the equivalent matrix Aleft for the second order tensorial operation AX=B => A_r*{x}={b}
  */
  void left(FullMatrix<double> &Aleft,FullMatrix<double> elm);
  /**
  *calculates the product of a fourth-order tensor and second-order Tensor to calculate stress
  */
  void ElasticProd(FullMatrix<double> &stress,FullMatrix<double> elm, FullMatrix<double> ElasticityTensor);
  /**
  *calculates the equivalent matrix Aleft for the second order tensorial operation trace(AX)=B => A_r*{x}={b}
  */
  void tracev(FullMatrix<double> &Atrace, FullMatrix<double> elm, FullMatrix<double> B);
  void rod2quat(Vector<double> &quat,Vector<double> &rod);
  void quatproduct(Vector<double> &quatp,Vector<double> &quat2,Vector<double> &quat1);
  void quat2rod(Vector<double> &quat,Vector<double> &rod);
  void elasticmoduli(FullMatrix<double> &Ar, FullMatrix<double> R, FullMatrix<double> Av);
  /**
  *calculates the matrix exponential of 3x3 matrix A
  */
  FullMatrix<double> matrixExponential(FullMatrix<double> A);
  /**
  *calculates the matrix exponential of 3x3 matrix A
  */
  FullMatrix<double> matrixExponential6(FullMatrix<double> A);
        /**
  *calculates the matrix exponential of 3x3 matrix A
  */
  FullMatrix<double> matrixExponentialGateauxDerivative(FullMatrix<double> A, FullMatrix<double> B);
        /**
   *calculates the matrix exponential of 3x3 matrix A
   */
  FullMatrix<double> matrixExponentialGateauxDerivative2(FullMatrix<double> A, FullMatrix<double> B);
  /**
   * Implements line search to solve the optimization problem
   */
  void lnsrch(Vector<double> &statenew, unsigned int n, Vector<double> stateold, double Fold, Vector<double> gradFold, Vector<double> srchdir,double delgam_ref, double strexp, FullMatrix<double> SCHMID_TENSOR1, unsigned int n_slip_systems, unsigned int n_Tslip_systems, Vector<double> s_alpha_tau, FullMatrix<double> Dmat, FullMatrix<double> CE_tau_trial, Vector<double> W_kh_t1, Vector<double> W_kh_t2, double hb1, double hb2, double mb1, double mb2, double rb1, double rb2, double bb1, double bb2);
  /**
  * Global deformation gradient F
  */
  FullMatrix<double> F;
  /**
  * Deformation gradient in crystal plasticity formulation. By default F=F_tau
  */
  FullMatrix<double> F_tau;
  /**
  * Plastic deformation gradient in crystal plasticity formulation. F_tau=Fe_tau*Fp_tau
  */
  FullMatrix<double> FP_tau;
  /**
  * Elastic deformation gradient in crystal plasticity formulation. F_tau=Fe_tau*Fp_tau
  */
  FullMatrix<double> FE_tau;
  /**
  * Cauchy Stress T
  */
  FullMatrix<double> T;
  /**
  * First Piola-Kirchhoff stress
  */
  FullMatrix<double> P;
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
  /**
  * Tangent modulus dPK1/dF
  */
  Tensor<4,dim,double> dP_dF;
  double No_Elem, N_qpts,local_F_e,local_F_r,F_e,F_r,local_microvol,microvol,F_s,local_F_s,F_T;
  std::vector<std::vector<double> > outputQuadrature;
  /**
  * Stores original crystal orientations as rodrigues vectors by element number and quadratureID
  */
  std::vector<std::vector<  Vector<double> > >  rot_conv,rot_iter,rot;
  /**
  * Stores deformed crystal orientations as rodrigues vectors by element number and quadratureID
  */
  std::vector<std::vector<  Vector<double> > >  rotnew_conv,rotnew_iter;
  /**
  * Stores the additional voxel data by element number and quadratureID
  */
  std::vector<std::vector<  Vector<double> > >  VoxelData;
  //Store history variables
  /**
  * Stores Plastic deformation gradient by element number and quadratureID at each iteration
  */
  std::vector< std::vector< FullMatrix<double> > > Fp_iter;
  /**
  * Stores Plastic deformation gradient by element number and quadratureID at each increment
  */
  std::vector< std::vector< FullMatrix<double> > > Fp_conv;
  /**
  * Stores Elastic deformation gradient by element number and quadratureID at each iteration
  */
  std::vector< std::vector< FullMatrix<double> > >   Fe_iter;
  /**
  * Stores Elastic deformation gradient by element number and quadratureID at each increment
  */
  std::vector< std::vector< FullMatrix<double> > > Fe_conv;
  /**
  * Stores Cauchy Stress by element number and quadratureID at each increment
  */
  std::vector< std::vector< FullMatrix<double> > > CauchyStress;
  std::vector< std::vector< FullMatrix<double> > > F_lastIter_Global;
  std::vector< std::vector< FullMatrix<double> > > FirstPiolaStress;
  Vector<double> workDensityTotal1_Tr;
  std::vector< std::vector< FullMatrix<double> > > TinterStress;
  std::vector< std::vector< FullMatrix<double> > > TinterStress_diff;
  /**
  * Stores slip resistance by element number and quadratureID at each iteration
  */
  std::vector<std::vector<  Vector<double> > >  s_alpha_iter;
  std::vector<std::vector<  Vector<double> > >  s_alpha_conv;
  std::vector<std::vector<  Vector<double> > >  W_kh_conv, W_kh_iter;
  Vector<double> Wkh_tau;
  FullMatrix<double> T_inter ;
  /**
  * Stores state variables by element number and quadratureID
  */
  std::vector<std::vector<  Vector<double> > >  stateVar_conv,stateVar_iter;
  Vector<double> UserMatConstants;
  std::vector<std::vector<  std::vector<double> > >  twinfraction_iter, slipfraction_iter,twinfraction_conv, slipfraction_conv,TwinOutputfraction_iter,TwinOutputfraction_conv;
  std::vector<std::vector<std::vector<unsigned int> > >	TwinFlag_conv, ActiveTwinSystems_conv, TwinFlag_iter, ActiveTwinSystems_iter;
  std::vector<std::vector<unsigned int> > NumberOfTwinnedRegion_conv, TwinMaxFlag_iter, TwinMaxFlag_conv, NumberOfTwinnedRegion_iter;
  std::vector<std::vector<double> >  twin_ouput, TotaltwinvfK;
  std::vector<std::vector<unsigned int> >  twin_iter, twin_conv;
  std::vector<unsigned int>  n_slip_systems_MultiPhase, n_twin_systems_MultiPhase, n_Tslip_systems_MultiPhase,n_Tslip_systems_Real_MultiPhase, n_UserMatStateVar_MultiPhase;
  std::vector<std::vector<unsigned int>> phase;
  unsigned int n_slip_systems,n_Tslip_systems,n_twin_systems,n_slip_systems_SinglePhase,n_Tslip_systems_SinglePhase,n_twin_systems_SinglePhase, phaseMaterial, numberofPhases, n_UserMatStateVar, n_UserMatStateVar_SinglePhase; //No. of slip systems
  FullMatrix<double> m_alpha,n_alpha,m_alpha_SinglePhase,n_alpha_SinglePhase,m_alpha_MultiPhase,n_alpha_MultiPhase,q,q_phase1,q_phase2,q_phase3,q_phase4,sres,Dmat,Dmat_SinglePhase,Dmat_MultiPhase, eulerAngles2;
  double twinShear;
  Vector<double> initialHardeningModulus, saturationStress, powerLawExponent, initialHardeningModulusTwin, saturationStressTwin, powerLawExponentTwin;
  FullMatrix<double> elasticStiffnessMatrix;
  Vector<double> C_1, C_2; // Backstress
  bool enableTwinning;
  double twinThresholdFraction;
  double twinSaturationFactor;
  /**
  * slip resistance
  */
  Vector<double> sres_tau;
  bool initCalled;
  //orientatations data for each quadrature point
  /**
  * Stores grainID No. by element number and quadratureID
  */
  std::vector<unsigned int> cellOrientationMap;
  void loadOrientations();
};

#endif
