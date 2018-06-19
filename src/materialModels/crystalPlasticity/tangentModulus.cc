#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::tangent_modulus(FullMatrix<double> &F_trial, FullMatrix<double> &Fpn_inv, FullMatrix<double> &SCHMID_TENSOR1, FullMatrix<double> &A,FullMatrix<double> &A_PA,FullMatrix<double> &B,FullMatrix<double> &T_tau, FullMatrix<double> &PK1_Stiff, Vector<double> &active, Vector<double> &resolved_shear_tau_trial, Vector<double> &x_beta, Vector<double> &PA, int &n_PA, double &det_F_tau, double &det_FE_tau ) {

    int iter=0,iter1=0,iter2=0,iter3=0;
    FullMatrix<double> temp(dim,dim),temp1(dim,dim),temp2(dim,dim),temp3(dim,dim);
    FullMatrix<double> Dmat2(2*dim,2*dim);
    FullMatrix<double> ElasticityTensor(6,6), TM(9,9);
    Vector<double> vec1(6), vec2(9);

    Dmat2=Dmat;
    Dmat2(3,3)=properties.C44; 	Dmat2(4,4)=properties.C44; 	Dmat2(5,5)=properties.C44;
    vec1(0)=0;vec1(1)=5;vec1(2)=4;vec1(3)=1;vec1(4)=3;vec1(5)=2;
    vec2(0)=0;vec2(1)=5;vec2(2)=4;vec2(3)=5;vec2(4)=1;vec2(5)=3;vec2(6)=4;vec2(7)=3;vec2(8)=2;

    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
            ElasticityTensor[i][j]=Dmat2(vec1(i),vec1(j));
        }
    }


    for(unsigned int i=0;i<9;i++){
        for(unsigned int j=0;j<9;j++){
            TM[i][j]=Dmat2(vec2(i),vec2(j));
        }
    }


    FullMatrix<double> F_temp(dim,dim);
    F_temp=Fpn_inv;
    FullMatrix<double> scratch_1,scratch_2,EtF;
    right(scratch_1, F_temp);
    for(unsigned int i=0;i<dim;i++){
        for(unsigned int j=0;j<dim;j++){
            F_temp[i][j]=F_trial[j][i];
        }
    }

    left(scratch_2,F_temp);
    F_temp.reinit(9,9);
    scratch_1.mmult(F_temp,scratch_2);
    symmf(EtF,F_temp);




    FullMatrix<double> s(dim,dim),p(n_PA,dim*dim),dgammadEmat(dim,dim), dgammadEtrial(n_PA,dim*dim),dgammadF(n_PA,dim*dim);
    Vector<double> P_as_vec3(dim*dim),s1(dim*dim),s2(dim*dim),temps2(dim*dim);

    if (n_PA > 0){

        s2=0.0;
        for(unsigned int alpha=0;alpha<n_PA;alpha++){
            iter=active(alpha);
            iter2=0;
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    s(j,k)=SCHMID_TENSOR1(dim*iter+j,k);
                    P_as_vec3(iter2)=s(j,k);
                    iter2++;
                }
            }

            TM.Tvmult(s1,P_as_vec3);
            if(resolved_shear_tau_trial(iter)<0){
                s1.equ(-1.0,s1);
            }
            s2=0.0;

            for(unsigned int beta=0;beta<n_PA;beta++){
                iter3=active(beta);

                temp3.reinit(dim,dim);
                temp.reinit(dim,dim);
                for(unsigned int k=0;k<dim;k++){
                    for(unsigned int l=0;l<dim;l++){
                        temp3(k,l)=SCHMID_TENSOR1(dim*iter3+k,l);
                        temp(l,k)=temp3(k,l);
                    }
                }

                temp1.reinit(dim*dim,dim*dim); temp2.reinit(dim*dim,dim*dim);
                right(temp1,temp3); left(temp2,temp);
                temp1.add(1.0,temp2);
                TM.mmult(temp2,temp1); temp2.Tvmult(temps2,P_as_vec3);

                if((resolved_shear_tau_trial(iter)<0.0)^(resolved_shear_tau_trial(iter3)<0.0))
                    temps2.equ(-1.0*x_beta(iter3),temps2);
                else
                    temps2.equ(x_beta(iter3),temps2);

                s2.add(temps2);

            }

            for(unsigned int index1=0;index1<dim*dim;index1++){
                p(alpha,index1)=s1(index1)-s2(index1);
            }

        }

        A_PA.reinit(n_PA,n_PA);
        for(unsigned int i=0;i<(n_PA);i++){

            for(unsigned int j=0;j<n_PA;j++){
                A_PA[i][j]=A[PA(i)][PA(j)];
            }
        }

        temp.reinit(n_PA,n_PA); temp.invert(A_PA);
        temp.mmult(dgammadEtrial,p); dgammadEtrial.mmult(dgammadF,EtF);

    }

    s.reinit(dim*dim,dim*dim); s=0.0;

    int alpha=0;
    if(n_PA>0){
        for(unsigned int i=0;i<n_PA;i++){
            alpha = active(i);
            iter=0;
            temp3.reinit(dim,dim);
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    dgammadEmat[k][j]=dgammadEtrial[i][dim*j+k];
                    temp3[j][k]=B[dim*alpha+j][k];

                }
            }

            temp1.reinit(dim*dim,dim*dim); right(temp1,dgammadEmat);
            temp.reinit(dim,dim);
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    temp[j][k]=B[dim*alpha+j][k];
                }
            }

            temp3.reinit(dim,dim); temp3=0.0;
            ElasticProd(temp3,temp,ElasticityTensor);
            temp2.reinit(dim*dim,dim*dim); tracev(temp2,temp1,temp3);

            if(resolved_shear_tau_trial(alpha)>0){
                temp2.add(-1.5,temp2);
            }
            else{
                temp2.add(-0.5,temp2);
            }
            s.add(1.0,temp2);
        }
    }

    FullMatrix<double> smat1(dim*dim,dim*dim),smat2(dim*dim,dim*dim), smat3(dim*dim,dim*dim);

    smat1=0.0; smat1.add(1.0,s); smat1.add(1.0,TM);
    smat3=0.0;

    if(n_PA>0){

        for(unsigned int i=0;i<n_PA;i++){
            alpha = active(i);
            F_temp.reinit(dim,dim);
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    F_temp[j][k]=SCHMID_TENSOR1[dim*alpha+j][k];
                    temp3[k][j]=F_temp[j][k];
                }
            }

            right(temp1,F_temp);
            left(temp2,temp3);
            temp1.add(1.0,temp2);
            TM.mmult(temp2,temp1);
            if(resolved_shear_tau_trial(alpha)>0){
                temp2.equ(-1.0*x_beta(alpha),temp2);
            }
            else{
                temp2.equ(1.0*x_beta(alpha),temp2);
            }
            smat3.add(1.0,temp2);
        }
    }


    FullMatrix<double> tangent_moduli(dim*dim,dim*dim),dgammadFmat(dim,dim);

    tangent_moduli=0.0; tangent_moduli.add(1.0,smat1); tangent_moduli.add(1.0,smat3);
    smat2=0.0;

    if(n_PA>0){

        for(unsigned int i=0;i<n_PA;i++){
            alpha = active(i);
            temp2.reinit(dim,dim);
            for(unsigned int j=0;j<dim;j++){
                for(unsigned int k=0;k<dim;k++){
                    dgammadFmat[k][j]=dgammadF[i][dim*j+k];
                    temp2[j][k]=SCHMID_TENSOR1[dim*alpha+j][k];
                }
            }

            right(temp1,dgammadFmat);
            F_trial.mmult(temp3,temp2);
            tracev(temp2,temp1,temp3);
            if(resolved_shear_tau_trial(alpha)>0){
                temp2.equ(-1.0,temp2);
            }
            smat2.add(1.0,temp2);
        }
    }

    temp2.reinit(dim,dim);
    smat1.reinit(dim*dim,dim*dim);
    temp2.invert(FP_tau);
    right(smat1,temp2);


    FullMatrix<double> FeF(dim*dim,dim*dim);
    FeF=0.0; FeF.add(1.0,smat1); FeF.add(1.0,smat2);
    temp1.reinit(dim,dim); temp1.invert(FE_tau);
    F_temp.reinit(dim,dim);


    //term 1; -trace(dFe Fe_inv) T * (1/det(Fe)) = [Term1]{dF}
    F_temp=temp1;
    FullMatrix<double> C4(dim*dim,dim*dim);

    right(C4,F_temp);

    FullMatrix<double> Term1(dim*dim,dim*dim),Term2(dim*dim,dim*dim),Term3(dim*dim,dim*dim),Term24(dim*dim,dim*dim);
    tracev(Term3,C4,T_tau);
    Term3.mmult(Term1,FeF); Term1.equ(-1.0, Term1);

    //term 3; Fe * dT_bar * FeT * (1/det(Fe)) = [Term3]{dF}

    F_temp.reinit(dim,dim);
    F_temp=FE_tau;

    FullMatrix<double> cauchy_Stiff(dim*dim,dim*dim);
    scratch_1.reinit(dim*dim,dim*dim); left(scratch_1,F_temp);
    temp1.reinit(dim,dim); temp1=IdentityMatrix(dim);

    temp2.reinit(dim,dim); temp2=F_temp;
    temp2.Tmmult(F_temp,temp1);

    right(C4,F_temp); C4.mmult(scratch_2,scratch_1);
    scratch_1.reinit(dim*dim,dim*dim); scratch_2.mmult(scratch_1,tangent_moduli);
    Term3.reinit(dim*dim,dim*dim); scratch_1.mmult(Term3,EtF);
    Term3.equ(1.0/det_FE_tau, Term3);




    //term 2&4: (dFe * Fe_inv) * T + T * (dFe * Fe_inv)^T = 2Symm(dFe * Fe_inv * T) = [Term24]{dF}

    temp2.invert(FE_tau); temp2.mmult(F_temp,T_tau);
    right(C4,F_temp);
    temp3.reinit(dim*dim,dim*dim); temp3=C4;
    symmf(C4,temp3);
    Term24=0.0; C4.mmult(Term24,FeF); Term24.equ(2.0,Term24);
    cauchy_Stiff=0.0; cauchy_Stiff.add(1.0, Term1); cauchy_Stiff.add(1.0, Term24); cauchy_Stiff.add(1.0, Term3);
    temp1.reinit(dim,dim); temp1.invert(F_tau);

    C4.reinit(dim*dim,dim*dim); right(C4,temp1);
    Term1.reinit(dim*dim,dim*dim); tracev(Term1,C4,T_tau);

    scratch_1.reinit(dim*dim,dim*dim); right(scratch_1,F_temp);
    Term2.reinit(dim*dim,dim*dim); symmf(Term2,scratch_1);
    Term2.equ(2.0,Term2); Term2.add(-1.0,scratch_1); Term2.equ(-1.0,Term2);




    Term3.reinit(dim*dim,dim*dim); Term3=cauchy_Stiff;

    PK1_Stiff=0.0; PK1_Stiff.add(1.0,Term1); PK1_Stiff.add(1.0,Term2); PK1_Stiff.add(1.0,Term3);
    temp2.reinit(dim,dim); temp2=IdentityMatrix(dim);
    temp3.reinit(dim,dim); temp3=temp1; temp3.Tmmult(temp1,temp2);
    temp2.reinit(dim*dim,dim*dim); right(temp2,temp1);
    temp3.reinit(dim*dim,dim*dim); temp3=PK1_Stiff;
    temp3.mmult(PK1_Stiff,temp2);
    PK1_Stiff.equ(det_F_tau,PK1_Stiff);


}

#include "../../../include/crystalPlasticity_template_instantiations.h"
