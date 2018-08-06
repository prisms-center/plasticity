#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::odfpoint(FullMatrix <double> &OrientationMatrix,Vector<double> r) {


    //function OrientationMatrix = odfpoint(r)
    //%USAGE: [C] = odfpoint(R)
    //%TO OBTAIN ORIENTATION MATRICES FROM RODRIGUES FORM

    double rdotr = 0.0;

    for(unsigned int i=0;i<dim;i++){
        rdotr = rdotr + r(i)*r(i);
    }

    double term1 = 1.0 - (rdotr);
    double term2 = 1.0 + (rdotr);

    OrientationMatrix.reinit(dim,dim);OrientationMatrix = IdentityMatrix(dim);

    for(unsigned int i=0;i<dim;i++){

        OrientationMatrix[i][i]=OrientationMatrix[i][i]*term1;
    }

    for(unsigned int i=0;i<dim;i++){

        for(unsigned int j=0;j<dim;j++){
            OrientationMatrix[i][j] = OrientationMatrix[i][j] + 2.0*(r(i)*r(j));
        }
    }

    OrientationMatrix[0][1] = OrientationMatrix[0][1]-2.0*r(2);
    OrientationMatrix[0][2] =  OrientationMatrix[0][2]+2.0*r(1);
    OrientationMatrix[1][2] = OrientationMatrix[1][2]-2.0*r(0);
    OrientationMatrix[1][0] =  OrientationMatrix[1][0]+2.0*r(2);
    OrientationMatrix[2][0] = OrientationMatrix[2][0]-2.0*r(1);
    OrientationMatrix[2][1] =  OrientationMatrix[2][1]+2.0*r(0);


    for(unsigned int i=0;i<dim;i++){

        for(unsigned int j=0;j<dim;j++){
            OrientationMatrix[i][j] = OrientationMatrix[i][j]*1.0/term2;
        }
    }



}

#include "../../../include/crystalPlasticity_template_instantiations.h"
