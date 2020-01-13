#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::reorient2(Vector<double> &rnew, Vector<double> rold, FullMatrix<double> FE_tau, FullMatrix<double> FE_t) {
    //Update the history variables

    LAPACKFullMatrix<double> C_old(dim,dim),C_new(dim,dim);

    Vector<double> eigenvalues(dim);
    FullMatrix<double> C_old_temp(dim,dim),C_new_temp(dim,dim), Fe_old(dim, dim), Fe_new(dim, dim),  eigenvectors(dim,dim),Lambda(dim,dim),U_old(dim,dim),U_new(dim,dim),R_old(dim,dim),R_new(dim,dim),Omega(dim,dim),temp(dim,dim);
    Lambda=IdentityMatrix(dim);
    Omega=0.0;
    Vector<double> rot1(dim),Omega_vec(dim);
    FullMatrix<double> rotmat(dim,dim);



            C_old_temp=0.0;
            C_old=0.0;

            Fe_old= FE_t;
            Fe_new= FE_tau;
			FE_t.Tmmult(C_old_temp, FE_t);
            C_old=C_old_temp;
            C_old.compute_eigenvalues_symmetric(0.0,200000.0,1e-15,eigenvalues, eigenvectors);

            for(unsigned int k=0;k<dim;k++){
                Lambda(k,k)=sqrt(eigenvalues(k));
            }


            eigenvectors.mmult(U_old,Lambda);
            temp=U_old; temp.mTmult(U_old,eigenvectors);
            R_old.invert(U_old);
            temp=R_old; FE_t.mmult(R_old,temp);
			FE_tau.Tmmult(C_new_temp, FE_tau);
            C_new=C_new_temp;
            C_new.compute_eigenvalues_symmetric(0.0,200000.0,1e-15,eigenvalues, eigenvectors);

            for(unsigned int k=0;k<dim;k++){
                Lambda(k,k)=sqrt(eigenvalues(k));
            }

            eigenvectors.mmult(U_new,Lambda);
            temp=U_new; temp.mTmult(U_new,eigenvectors);
            R_new.invert(U_new);
            temp=R_new; FE_tau.mmult(R_new,temp);

            Omega=0.0; Omega.add(1.0,R_new); Omega.add(-1.0,R_old);
            temp=Omega; temp.mTmult(Omega,R_new);

            Omega_vec(0)=-0.5*(Omega(1,2)-Omega(2,1));Omega_vec(1)=0.5*(Omega(0,2)-Omega(2,0));Omega_vec(2)=-0.5*(Omega(0,1)-Omega(1,0));

            double dot;
            dot=Omega_vec(0)*rold(0)+Omega_vec(1)*rold(1)+Omega_vec(2)*rold(2);
            Vector<double> cross(dim),dot_term(dim);

            dot_term(0)=dot*rold(0);dot_term(1)=dot*rold(1);dot_term(2)=dot*rold(2);

            cross(0)=Omega_vec(1)*rold(2)-Omega_vec(2)*rold(1);
            cross(1)=Omega_vec(2)*rold(0)-Omega_vec(0)*rold(2);
            cross(2)=Omega_vec(0)*rold(1)-Omega_vec(1)*rold(0);

			Vector<double> dr(dim);

            dr=0.0;	dr.add(1.0, Omega_vec); dr.add(1.0,dot_term); dr.add(1.0,cross); dr.equ(0.5,dr);

            rnew=0.0; rnew.add(1.0,rold); rnew.add(1.0,dr);
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
