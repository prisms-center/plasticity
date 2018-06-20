#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::inactive_slip_removal(Vector<double> &active, Vector<double> &x_beta_old, Vector<double> &x_beta, unsigned int &n_PA, Vector<double> &PA, Vector<double> b,FullMatrix<double> A,FullMatrix<double> &A_PA){

    Vector<double> inactive;

    FullMatrix<double> temp,temp2;
    LAPACKFullMatrix<double> temp7;
    temp7.reinit(n_PA,n_PA);
    int iter=0,iter1=0,iter2=0,iter3=0;
    Vector<double> x_beta1(n_PA), x_beta2(n_PA),b_PA(n_PA);
    double flag1=0;
    A_PA.reinit(n_PA,n_PA);

    for(unsigned int i=0;i<n_PA;i++){
        b_PA(i)=b(PA(i));
        for(unsigned int j=0;j<n_PA;j++){
            A_PA[i][j]=A[PA(i)][PA(j)];
            temp7(i,j)=A[PA(i)][PA(j)];

        }

    }
    temp.reinit(n_PA,n_PA); temp2.reinit(n_PA,n_PA);
    temp7.compute_inverse_svd(0.0);
    Vector<double> tempv3;
    temp7.vmult(x_beta1,b_PA);

    //Check for model tolerance and activate adaptive time-stepping, if required
    if(x_beta1.l2_norm()> this->userInputs.modelMaxPlasticSlipL2Norm){
#ifdef enableAdaptiveTimeStepping
#if enableAdaptiveTimeStepping==true
        char buffer[200];
        sprintf (buffer, "processor %u: time-step is very large. Consider reducing the time-step. current model norm: %12.6e, tolerance: %12.6e\n", this->triangulation.locally_owned_subdomain(), x_beta1.l2_norm(), modelMaxPlasticSlipL2Norm);
        std::cout <<buffer;
        this->loadFactorSetByModel*=adaptiveLoadStepFactor;
        this->resetIncrement=true;
        throw 0;
#endif
#endif
    }

    x_beta2=x_beta1; x_beta1.reinit(n_slip_systems);
    x_beta1=0;
    for(unsigned int i=0;i<n_PA;i++){
        x_beta1(PA(i))=x_beta2(i);


    }



    Vector<double> row,PA_new,inactive2;
    double n_IA_new;

    // Continue the process till removal of all inactive slip systems
    while (flag1==0){

        iter=0;
        for(unsigned int i=0;i<n_slip_systems;i++){
            if((x_beta_old(i)+x_beta1(i))<0){
                iter++;
            }
        }

        if(iter>0){

            row.reinit(iter);
            iter=0;
            for(unsigned int i=0;i<n_slip_systems;i++){
                if((x_beta_old(i)+x_beta1(i))<0){
                    row(iter)=i;
                    iter++;
                }
            }

            n_IA_new=iter;
            PA_new.reinit(n_PA-n_IA_new);
            iter2=0; iter=0;
            for(unsigned int i=0;i<n_PA;i++){
                for(unsigned int j=0;j<n_IA_new;j++){
                    if(PA(i) != row(j)){
                        iter++;
                    }
                }

                if(iter==n_IA_new){
                    PA_new(iter2)=PA(i);
                    iter2++;
                }
                iter=0;

            }

            PA.reinit(n_PA-n_IA_new); PA=PA_new;

            A_PA.reinit(n_PA-n_IA_new,n_PA-n_IA_new); A_PA=0.0;
            b_PA.reinit(n_PA-n_IA_new); b_PA=0.0;

            temp7.reinit((n_PA-n_IA_new),(n_PA-n_IA_new));
            for(unsigned int i=0;i<(n_PA-n_IA_new);i++){
                b_PA(i)=b(PA(i));
                for(unsigned int j=0;j<(n_PA-n_IA_new);j++){
                    A_PA[i][j]=A[PA(i)][PA(j)];
                    temp7(i,j)=A[PA(i)][PA(j)];

                }

            }


            temp7.compute_inverse_svd(0.0);

            x_beta1.reinit(n_PA-n_IA_new); temp7.vmult(x_beta1,b_PA);
            x_beta2.reinit(n_PA-n_IA_new); x_beta2=x_beta1;
            x_beta1.reinit(n_slip_systems);x_beta1=0.0;
            n_PA=n_PA-n_IA_new;

            for(unsigned int i=0;i<n_PA;i++){
                x_beta1(PA(i))=x_beta2(i);
            }

        }

        else
            flag1=1;

    }

    active.reinit(n_PA); active=PA;
    x_beta.reinit(n_slip_systems); x_beta=x_beta1;
    x_beta_old.add(1.0,x_beta);




}

#include "../../../include/crystalPlasticity_template_instantiations.h"
