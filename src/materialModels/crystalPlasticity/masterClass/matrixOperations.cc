#include "../../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::tracev(FullMatrix<double> &Atrace, FullMatrix<double> elm, FullMatrix<double> B) {

    Atrace.reinit(dim*dim,dim*dim);
    Vector<double> C(dim*dim);

    for(unsigned int i=0;i<9;i++){
        C(i)=elm(0,i)+elm(4,i)+elm(8,i);
    }

    for(unsigned int i=0;i<9;i++){
        for(unsigned int j=0;j<dim;j++){
            for(unsigned int k=0;k<dim;k++){
                Atrace(dim*j+k, i) =  B(j, k) * C(i);
            }
        }
    }

}




template <int dim>
void crystalPlasticity<dim>::ElasticProd(FullMatrix<double> &stress,FullMatrix<double> elm, FullMatrix<double> ElasticityTensor) {

    int index=0;
    stress.reinit(dim,dim);
    Vector<double> temp1(2*dim),temp2(2*dim);

    for(unsigned int i=0;i<dim;i++){

        for(unsigned int j=i;j<dim;j++){

            temp1(index) = elm(i, j);
            index = index+1;
        }
    }

    ElasticityTensor.vmult(temp2,temp1);

    index = 0;

    for(unsigned int i=0;i<dim;i++){

        for(unsigned int j=i;j<dim;j++){

            stress(i, j) = temp2(index);
            index = index+1;
        }
    }


    stress(1, 0) = stress(0, 1);
    stress(2, 0) = stress(0, 2);
    stress(2, 1) = stress(1, 2);

}






template <int dim>
void crystalPlasticity<dim>::right(FullMatrix<double> &Aright,FullMatrix<double> elm) {

    Aright.reinit(dim*dim,dim*dim);

    Aright=0.0;

    for(unsigned int i=0;i<dim;i++){

        for(unsigned int j=0;j<dim;j++){

            Aright[i][j]=elm(j,i) ;
            Aright[i+3][j+3]=elm(j,i) ;
            Aright[i+6][j+6]=elm(j,i) ;
        }
    }

}


template <int dim>
void crystalPlasticity<dim>::left(FullMatrix<double> &Aleft,FullMatrix<double> elm) {

    Aleft.reinit(dim*dim,dim*dim);

    Aleft=0.0;

    for(unsigned int i=0;i<dim;i++){

        for(unsigned int j=0;j<dim;j++){

            Aleft[dim*i][dim*j]=elm(i,j) ;
            Aleft[dim*i+1][dim*j+1]=elm(i,j) ;
            Aleft[dim*i+2][dim*j+2]=elm(i,j) ;
        }
    }

}


template <int dim>
void crystalPlasticity<dim>::symmf(FullMatrix<double> &A, FullMatrix<double> elm) {

    A.reinit(dim*dim,dim*dim);

    for(unsigned int i=0;i<9;i++){

        for(unsigned int j=0;j<9;j++){

            if (i == 1 || i == 3)
                A(i, j) = 0.5* (elm(1, j) + elm(3, j));
            else
                if (i == 2 || i == 6)
                    A(i, j) = 0.5* (elm(2, j) + elm(6, j));
                else
                    if( i == 5 || i == 7)
                        A(i, j) = 0.5* (elm(5, j) + elm(7, j));
                    else
                        A(i, j) = elm(i, j);
        }
    }

}



template <int dim>
Vector<double> crystalPlasticity<dim>::vecform(FullMatrix<double> A) {

    Vector<double> Av(6);

    Av(0) =A[0][0];
    Av(1) =A[1][1];
    Av(2) =A[2][2];
    Av(3) =A[1][2];
    Av(4) =A[0][2];
    Av(5) =A[0][1];

    return Av;
}

template <int dim>
void crystalPlasticity<dim>::matform(FullMatrix<double> &A, Vector<double> Av) {

    A.reinit(dim,dim);

    A[0][0]=Av(0) ;
    A[1][1]=Av(1) ;
    A[2][2]=Av(2) ;
    A[1][2]=Av(3) ;
    A[0][2]=Av(4) ;
    A[0][1]=Av(5) ;
    A[2][1]=Av(3) ;
    A[2][0]=Av(4) ;
    A[1][0]=Av(5) ;

}

template <int dim>
FullMatrix<double> crystalPlasticity<dim>::matrixExponential(FullMatrix<double> A) {

    FullMatrix<double> matExp(dim,dim),temp(dim,dim),temp2(dim,dim);
    matExp=IdentityMatrix(dim);
    temp=IdentityMatrix(dim);

    double count=1;

    while(temp.frobenius_norm()>1e-15){

        temp.mmult(temp2,A);
        temp2.equ(1/count,temp2);
        temp.equ(1.0,temp2);
        matExp.add(1.0,temp);
        count=count+1.0;


    }

    return matExp;

}

#include "../../../../include/crystalPlasticity_template_instantiations.h"
