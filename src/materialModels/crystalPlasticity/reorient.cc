#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::reorient() {
    //Update the history variables

    LAPACKFullMatrix<double> C_old(dim,dim),C_new(dim,dim);

    Vector<double> eigenvalues(dim);
    FullMatrix<double> C_old_temp(dim,dim),C_new_temp(dim,dim),Fe_old(dim,dim), Fe_new(dim,dim), eigenvectors(dim,dim),Lambda(dim,dim),U_old(dim,dim),U_new(dim,dim),R_old(dim,dim),R_new(dim,dim),Omega(dim,dim),temp(dim,dim);
    Lambda=IdentityMatrix(dim);
    Omega=0.0;
    Vector<double> rot1(dim),Omega_vec(dim),rold(dim),dr(dim),rnew(dim);
    FullMatrix<double> rotmat(dim,dim);
    //int itgno;
    unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();

    for (unsigned int i=0; i<num_local_cells; ++i) {
        for(unsigned int j=0;j<N_qpts;j++){

            C_old_temp=0.0;
            C_old=0.0;

            Fe_old=Fe_conv[i][j];
            Fe_new=Fe_iter[i][j];
            
            /////In the case of twin reorientation, do not update the orientation for one loading increment
            if (twin_conv[i][j]!=twin_iter[i][j]){
              Fe_old=Fe_new;
            }
            ///////////////////////////////////////////////////////////////////////////////////
            
            Fe_old.Tmmult(C_old_temp,Fe_old);
            C_old=C_old_temp;
            C_old.compute_eigenvalues_symmetric(0.0,200000.0,1e-15,eigenvalues, eigenvectors);

            for(unsigned int k=0;k<dim;k++){
                Lambda(k,k)=sqrt(eigenvalues(k));
            }


            eigenvectors.mmult(U_old,Lambda);
            temp=U_old; temp.mTmult(U_old,eigenvectors);
            R_old.invert(U_old);
            temp=R_old; Fe_old.mmult(R_old,temp);
            Fe_new.Tmmult(C_new_temp,Fe_new);
            C_new=C_new_temp;
            C_new.compute_eigenvalues_symmetric(0.0,200000.0,1e-15,eigenvalues, eigenvectors);

            for(unsigned int k=0;k<dim;k++){
                Lambda(k,k)=sqrt(eigenvalues(k));
            }

            eigenvectors.mmult(U_new,Lambda);
            temp=U_new; temp.mTmult(U_new,eigenvectors);
            R_new.invert(U_new);
            temp=R_new; Fe_new.mmult(R_new,temp);

            Omega=0.0; Omega.add(1.0,R_new); Omega.add(-1.0,R_old);
            temp=Omega; temp.mTmult(Omega,R_new);

            rold=rotnew_conv[i][j];


            Omega_vec(0)=-0.5*(Omega(1,2)-Omega(2,1));Omega_vec(1)=0.5*(Omega(0,2)-Omega(2,0));Omega_vec(2)=-0.5*(Omega(0,1)-Omega(1,0));

            double dot;
            dot=Omega_vec(0)*rold(0)+Omega_vec(1)*rold(1)+Omega_vec(2)*rold(2);
            Vector<double> cross(dim),dot_term(dim);

            dot_term(0)=dot*rold(0);dot_term(1)=dot*rold(1);dot_term(2)=dot*rold(2);

            cross(0)=Omega_vec(1)*rold(2)-Omega_vec(2)*rold(1);
            cross(1)=Omega_vec(2)*rold(0)-Omega_vec(0)*rold(2);
            cross(2)=Omega_vec(0)*rold(1)-Omega_vec(1)*rold(0);

            dr=0.0;	dr.add(1.0, Omega_vec); dr.add(1.0,dot_term); dr.add(1.0,cross); dr.equ(0.5,dr);

            rnew=0.0; rnew.add(1.0,rold); rnew.add(1.0,dr);
            
            ///////Very large Rodrigues vector norm leads to Nan or Inf for updated rnew. Accordingly, we keep the maximum norm to max_rnew_Norm=10000.
            double rnew_Norm, max_rnew_Norm;
            max_rnew_Norm=10000;
            rnew_Norm=sqrt(rnew(0)*rnew(0)+rnew(1)*rnew(1)+rnew(2)*rnew(2));

            if (rnew_Norm>max_rnew_Norm){
                rnew(0)=rnew(0)*max_rnew_Norm/rnew_Norm;
                rnew(1)=rnew(1)*max_rnew_Norm/rnew_Norm;
                rnew(2)=rnew(2)*max_rnew_Norm/rnew_Norm;
            }
            ///////////////////////

            rotnew_conv[i][j]=rnew;


        }
    }


}

#include "../../../include/crystalPlasticity_template_instantiations.h"
