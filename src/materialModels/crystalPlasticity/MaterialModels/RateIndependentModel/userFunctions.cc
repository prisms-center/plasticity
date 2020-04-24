#include "../../../include/ellipticBVP.h"

////////The required functions should be defined here for user defined material modesl///
//// As an example, the function rearrng is defined here which will be used in the rate-dependent model////////
void rearrng(FullMatrix<double> &Aarrange, FullMatrix<double> elm);

void rearrng(FullMatrix<double> &Aarrange, FullMatrix<double> elm) {

    Aarrange.reinit(9,9);

        for(unsigned int k=0;k<9;k++){
                for(unsigned int i=0;i<3;i++){
                     for(unsigned int j=0;j<3;j++){

                                Aarrange(3*i+j,k) = elm(3*j+i,k);

                        }

	           }

	    }


}
