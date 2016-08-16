//implementation of the crystal plasticity material model
#ifndef CRYSTALPLASTICITY_H
#define CRYSTALPLASTICITY_H

//dealii headers
#include "../../../../include/ellipticBVP.h"
#include "../../../../src/utilityObjects/crystalOrientationsIO.cc"
#include <iostream>
#include <fstream>

typedef struct {
     FullMatrix<double> m_alpha,n_alpha;
} materialProperties;

//material model class for crystal plasticity
//derives from ellipticBVP base abstract class
template <int dim>
class crystalPlasticity : public ellipticBVP<dim>
{
public:
    crystalPlasticity();
#ifdef readExternalMeshes
#if readExternalMeshes==true
    void mesh();
#endif
#endif 
    void reorient();
    void tangent_modulus(FullMatrix<double> &F_trial, FullMatrix<double> &Fpn_inv, FullMatrix<double> &SCHMID_TENSOR1, FullMatrix<double> &A,FullMatrix<double> &A_PA,FullMatrix<double> &B,FullMatrix<double> &T_tau, FullMatrix<double> &PK1_Stiff, Vector<double> &active, Vector<double> &resolved_shear_tau_trial, Vector<double> &x_beta, Vector<double> &PA, int &n_PA, double &det_F_tau, double &det_FE_tau );
    void inactive_slip_removal(Vector<double> &active,Vector<double> &x_beta_old, Vector<double> &x_beta, int &n_PA, Vector<double> &PA, Vector<double> b,FullMatrix<double> A,FullMatrix<double> &A_PA);
    //material properties
    materialProperties properties;
    //orientation maps
    crystalOrientationsIO<dim> orientations;
private:
    void init(unsigned int num_quad_points);
    void setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value); 
    void calculatePlasticity(unsigned int cellID,
                             unsigned int quadPtID);
    void getElementalValues(FEValues<dim>& fe_values,
                            unsigned int dofs_per_cell,
                            unsigned int num_quad_points,
                            FullMatrix<double>& elementalJacobian,
                            Vector<double>&     elementalResidual);
    void updateAfterIncrement();
    void updateBeforeIteration();
    void updateBeforeIncrement();
    
    
    void odfpoint(FullMatrix <double> &OrientationMatrix,Vector<double> r);
    Vector<double> vecform(FullMatrix<double> A);
    void matform(FullMatrix<double> &A, Vector<double> Av);
    void right(FullMatrix<double> &Aright,FullMatrix<double> elm);
    void symmf(FullMatrix<double> &A,FullMatrix<double> elm);
    void left(FullMatrix<double> &Aleft,FullMatrix<double> elm);
    void ElasticProd(FullMatrix<double> &stress,FullMatrix<double> elm, FullMatrix<double> ElasticityTensor);
    void tracev(FullMatrix<double> &Atrace, FullMatrix<double> elm, FullMatrix<double> B);
    void Twin_image(double twin_pos,unsigned int cellID,unsigned int quadPtID);
    void rod2quat(Vector<double> &quat,Vector<double> &rod);
    void quatproduct(Vector<double> &quatp,Vector<double> &quat2,Vector<double> &quat1);
    void quat2rod(Vector<double> &quat,Vector<double> &rod);
    /**
     *calculates the matrix exponential of matrix A
     */
    
    FullMatrix<double> matrixExponential(FullMatrix<double> A);
    
    
    
    
    FullMatrix<double> F,F_tau,FP_tau,FE_tau,T,P,local_stress,local_strain,global_stress,global_strain;
    Tensor<4,dim,double> dP_dF;
    double No_Elem, N_qpts,local_F_e,local_F_r,F_e,F_r,local_microvol,microvol;
    
    //Store crystal orientations
    std::vector<std::vector<  Vector<double> > >  rot;
    std::vector<std::vector<  Vector<double> > >  rotnew;
    
    //Store history variables
    std::vector< std::vector< FullMatrix<double> > >   Fp_iter;
    std::vector< std::vector< FullMatrix<double> > > Fp_conv;
    std::vector< std::vector< FullMatrix<double> > >   Fe_iter;
    std::vector< std::vector< FullMatrix<double> > > Fe_conv;
    std::vector<std::vector<  Vector<double> > >  s_alpha_iter;
    std::vector<std::vector<  Vector<double> > >  s_alpha_conv;
    std::vector<std::vector<  vector<double> > >  twinfraction_iter, slipfraction_iter,twinfraction_conv, slipfraction_conv;
    std::vector<std::vector<double> >  twin;
    
    unsigned int n_slip_systems,n_twin_systems; //No. of slip systems
    FullMatrix<double> m_alpha,n_alpha,q,sres,Dmat;
    Vector<double> sres_tau;
    bool initCalled;
    
    //orientatations data for each quadrature point
    std::vector<std::vector<unsigned int> > quadratureOrientationsMap;  
    void loadOrientations();
};

//(these are source files, which will are temporarily treated as
//header files till library packaging scheme is finalized)
#include "model.cc"
#include "calculatePlasticity.cc"
#include "rotationOperations.cc"
#include "init.cc"
#include "matrixOperations.cc"
//#include "tangentModulus.cc"
#include "inactiveSlipRemoval.cc"
#include "reorient.cc"
#include "loadOrientations.cc"

#endif
