/*Class to statically condense and recover the following system:
  |K   G||d|=|F|
  |G^T M||s|=|H|
  where sizes are as follows: K(dof1,dof1), G(dof1,dof2), M(dof2,dof2), F(dof1), H(dof2)
 
  Call to staticCondense(), constructs K2, F2:
  K2=K-G*inv(M)*G^T
  F2=F-G*inv(M)*H

  And call to recover(), computes s:
  s = inv(M)*(H-G^T*d)
*/

/*
  NOTE: initialize the elemental level staticCondensation class objects by including the following line in the elasticity class definition (or whatever your main problem class is):
  std::vector<staticCondensation<dof1,dof2> > staticCondensationData;
  where dof1, dof2 are integer values defining the size of K,G,M,F,H. (For example - std::vector<staticCondensation<3,1> > staticCondensationData;)

  Then include the following line in setup_system(): 
  staticCondensationData.resize(triangulation.n_active_cells());

  Now you can access the elemental level staticCondensation object inside the assemble_system() function. For example, access the K matrix or print the F vector of 10th element as follows:
  staticCondensationData[9].K;
  staticCondensationData[9].F.print();
*/
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

using namespace dealii;

template <int _dof1, int _dof2>
class  staticCondensation{
public:
  /**
   *Class constructor.
   */
  staticCondensation();
  /**
   *There are two sets of degrees of freedom. dof1 is the number of dofs in the first set,
   *dof2 is the number in the second set.
   */
  int dof1, dof2;
  /**
   *K, G, and M are the submatrices of the full symmetric matrix. K is top left, M bottom right,
   *and G is the top right (its transpose is the bottom left block).
   *K2 is the resulting condensed matrix.
   */
  FullMatrix<double> K, G, M, K2; 
  /**
   *F and H are the subvectors (F top, H bottom). d is the vector of the first set of dofs,
   *s is the vector of the second set of dofs.
   *F2 is the resulting condensed vector.
   */
  Vector<double> F, H, d, s, F2;
  /**
   *Perform static condensation to compute K2 and F2.
   */
  void staticCondense();
  /**
   *Recover the second set of dofs, i.e. compute s after having solved for d.
   */
  void recover();
  /**
   *Reset all matrices and vectors to zero.
   */
  void reset();
};

template <int _dof1, int _dof2>
staticCondensation<_dof1,_dof2>::staticCondensation(): 
  dof1(_dof1), dof2(_dof2), K(_dof1, _dof1), G(_dof1, _dof2), M(_dof2, _dof2), K2(_dof1, _dof1),
  F(_dof1), H(_dof2), d(_dof1), s(_dof2), F2(_dof1) {
  //initialize all matrices and vectors to zero
  K=0; G=0; M=0; K2=0;
  F=0; H=0; d=0; s=0; F2=0;  
}

//initialize matrices (K,G,M) and vectors (F, H) before calling this function
//after a call to this function: K2, F2 are computed
template <int _dof1, int _dof2>
void staticCondensation<_dof1,_dof2>::staticCondense(){
  FullMatrix<double> invM(dof2,dof2), GinvM(dof1,dof2), GinvMG(dof1,dof1);
  Vector <double> GinvMH(dof1);
  invM.invert(M); //invM=inverse(M)
  G.mmult(GinvM,invM); //GinvM=G*inv(M)
  GinvM.mTmult(GinvMG, G); //GinvMG=G*inv(M)*G^T
  GinvM.vmult(GinvMH, H); //GinvMH=G*inv(M)*H
  //compute K2=K-G*inv(M)*G^T
  K2=K;
  K2.add(-1.0, GinvMG); 
  //compute F2=F-G*inv(M)*H
  F2=F;
  F2.add(-1.0,GinvMH);
} 

//initialize vector (d) before calling this function
//after a call to this function: s is computed
template <int _dof1, int _dof2>
void staticCondensation<_dof1,_dof2>::recover(){
  FullMatrix<double> invM(dof2,dof2), invMG(dof2,dof1); 
  Vector <double> invMH(dof2), invMGd(dof2);
  invM.invert(M); //invM=inverse(M)
  invM.vmult(invMH,H); //invMH=invM*H
  invM.mTmult(invMG,G); //invMG=inv(M)*G^T
  invMG.vmult(invMGd,d); //invMGd=invMG*d
  //compute s = inv(M)*(H-G^T*d)
  s=0;
  s.add(-1.0, invMH, -1.0, invMGd);
} 

//reset all vectors and matrices to zero
template <int _dof1, int _dof2>
void staticCondensation<_dof1,_dof2>::reset(){
  K=0; G=0; M=0; K2=0;
  F=0; H=0; d=0; s=0; F2=0;  
}

