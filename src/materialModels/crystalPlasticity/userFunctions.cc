#include "../../../include/ellipticBVP.h"


void trpose(FullMatrix<double> &Atrpose, FullMatrix<double> elm);

void trpose(FullMatrix<double> &Atrpose, FullMatrix<double> elm) {

    Atrpose.reinit(9,9);

        for(unsigned int k=0;k<9;k++){
                for(unsigned int i=0;i<3;i++){
                     for(unsigned int j=0;j<3;j++){

                                Atrpose(3*i+j,k) = elm(3*j+i,k);

                        }

	           }

	    }


}

void traceval(FullMatrix<double> &Atrace, FullMatrix<double> elm);

void traceval(FullMatrix<double> &Atrace, FullMatrix<double> elm){
	
	Atrace.reinit(9,9);
	Atrace=0.0;
	
	for(unsigned int i=0;i<3;i++){
                for(unsigned int j=0;j<3;j++){
					
					Atrace(3*i+j,0)=elm(i,j);
					Atrace(3*i+j,4)=elm(i,j);
					Atrace(3*i+j,8)=elm(i,j);
					
				}
	}
	
	
}

void vecform9(Vector<double> &Avec9, FullMatrix<double> elm);

void vecform9(Vector<double> &Avec9, FullMatrix<double> elm)
{
Avec9.reinit(9);

 for(unsigned int i=0;i<3;i++){
                for(unsigned int j=0;j<3;j++){

                                        Avec9(3*i+j)=elm(i,j);


                                }
        }


}

void matform9(FullMatrix<double> &Amat9, Vector<double> elm);

void matform9(FullMatrix<double> &Amat9, Vector<double> elm)
{
Amat9.reinit(3,3);

 for(unsigned int i=0;i<3;i++){
                for(unsigned int j=0;j<3;j++){

                                        Amat9(i,j)=elm(3*i+j);


                                }
        }


}


