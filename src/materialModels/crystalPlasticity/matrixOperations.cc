#include "../../../include/crystalPlasticity.h"


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

    while(temp.frobenius_norm()>1.0e-10){

        temp.mmult(temp2,A);
        temp2.equ(1/count,temp2);
        temp.equ(1.0,temp2);
        matExp.add(1.0,temp);
        count=count+1.0;


    }

    return matExp;

}

template <int dim>
void crystalPlasticity<dim>::elasticmoduli(FullMatrix<double> &Ar, FullMatrix<double> R, FullMatrix<double> Av) {

	Ar.reinit(2 * dim, 2 * dim);
	FullMatrix<double> newr(2 * dim, 2 * dim), temp(2 * dim, 2 * dim);

	newr[0][0] = R[0][0]*R[0][0];
	newr[0][1] = R[0][1]*R[0][1];
	newr[0][2] = R[0][2]*R[0][2];
	newr[0][3] = 2.0*R[0][1]*R[0][2];
	newr[0][4] = 2.0*R[0][2]*R[0][0];
	newr[0][5] = 2.0*R[0][0]*R[0][1];
	newr[1][0] = R[1][0]*R[1][0];
	newr[1][1] = R[1][1]*R[1][1];
	newr[1][2] = R[1][2]*R[1][2];
	newr[1][3] = 2.0*R[1][1]*R[1][2];
	newr[1][4] = 2.0*R[1][2]*R[1][0];
	newr[1][5] = 2.0*R[1][0]*R[1][1];
	newr[2][0] = R[2][0]*R[2][0];
	newr[2][1] = R[2][1]*R[2][1];
	newr[2][2] = R[2][2]*R[2][2];
	newr[2][3] = 2.0*R[2][1]*R[2][2];
	newr[2][4] = 2.0*R[2][2]*R[2][0];
	newr[2][5] = 2.0*R[2][0]*R[2][1];
	newr[3][0] = R[1][0]*R[2][0];
	newr[3][1] = R[1][1]*R[2][1];
	newr[3][2] = R[1][2]*R[2][2];
	newr[3][3] = R[1][1]*R[2][2] + R[1][2]*R[2][1];
	newr[3][4] = R[1][2]*R[2][0] + R[1][0]*R[2][2];
	newr[3][5] = R[1][0]*R[2][1] + R[1][1]*R[2][0];
	newr[4][0] = R[2][0]*R[0][0];
	newr[4][1] = R[2][1]*R[0][1];
	newr[4][2] = R[2][2]*R[0][2];
	newr[4][3] = R[0][1]*R[2][2] + R[0][2]*R[2][1];
	newr[4][4] = R[0][2]*R[2][0] + R[0][0]*R[2][2];
	newr[4][5] = R[0][0]*R[2][1] + R[0][1]*R[2][0];
	newr[5][0] = R[0][0]*R[1][0];
	newr[5][1] = R[0][1]*R[1][1];
	newr[5][2] = R[0][2]*R[1][2];
	newr[5][3] = R[0][1]*R[1][2] + R[0][2]*R[1][1];
	newr[5][4] = R[0][2]*R[1][0] + R[0][0]*R[1][2];
	newr[5][5] = R[0][0]*R[1][1] + R[0][1]*R[1][0];

	newr.mmult(temp, Av);
	temp.mTmult(Ar, newr);

}

template <int dim>
FullMatrix<double> crystalPlasticity<dim>::matrixExponential6(FullMatrix<double> A) {

	FullMatrix<double> matExp(2*dim,2*dim),temp(2*dim,2*dim),temp2(2*dim,2*dim);
	matExp=IdentityMatrix(2*dim);
	temp=IdentityMatrix(2*dim);
	double count=1;

	while(temp.frobenius_norm()>1e-10){
		temp.mmult(temp2,A);
		temp2.equ(1/count,temp2);
		temp.equ(1.0,temp2);
		matExp.add(1.0,temp);
		count=count+1.0;
	}

	return matExp;
}



template <int dim>
FullMatrix<double> crystalPlasticity<dim>::matrixExponentialGateauxDerivative(FullMatrix<double> A, FullMatrix<double> B) {

	FullMatrix<double> matExp(2*dim,2*dim),matExpGatDer(dim,dim);
	FullMatrix<double> C(2*dim,2*dim) ;
	C = 0 ;
	for(unsigned int i=0 ; i<3 ; i++){
		for(unsigned int j=0 ; j<3 ; j++){
			C[i][j] = A[i][j] ;
			C[i+3][j+3] = A[i][j] ;
			C[i][j+3] = B[i][j] ;
		}
	}

	matExp.equ(1.0,matrixExponential6(C)) ;

	for(unsigned int i=0 ; i<3 ; i++){
		for(unsigned int j=3 ; j<6 ; j++ ){
			matExpGatDer[i][j-3] = matExp[i][j] ;
		}
	}



	return matExpGatDer;
}

template <int dim>
FullMatrix<double> crystalPlasticity<dim>::matrixExponentialGateauxDerivative2(FullMatrix<double> A, FullMatrix<double> B){

	FullMatrix<double> temp1(dim,dim),temp2(dim,dim),temp3(dim,dim),matExpGatDer(dim,dim), temp4(dim,dim),temp5(dim,dim);
	double resdl = 1.0, count = 1.0, sctemp = 1.0 ;
	temp1 = 0;
	temp2.equ(1.0,B) ;
	matExpGatDer.equ(1.0,temp1) ;


	while(resdl > 1.0e-10){
		sctemp = sctemp/count ;
		matExpGatDer.add(sctemp,temp2) ;
		temp3.equ(1.0,temp2) ;
		temp1.mmult(temp4,A) ;
		temp4.equ(-1.0,temp4) ;
		temp4.add(1.0,temp2) ;
		A.mmult(temp5,temp4) ;
		temp2.mmult(temp4,A) ;
		temp4.add(1.0,temp5) ;
		temp2.equ(1.0,temp4) ;
		temp1.equ(1.0,temp3) ;

		temp5.equ(sctemp/(count + 1.0),temp2) ;
		resdl = temp5.frobenius_norm() ;

		count = count + 1.0 ;
    }

	return matExpGatDer ;


}



template <int dim>
void crystalPlasticity<dim>::lnsrch(Vector<double> &statenew, unsigned int n, Vector<double> stateold, double Fold, Vector<double> gradFold, Vector<double> srchdir,double delgam_ref, double strexp, FullMatrix<double> SCHMID_TENSOR1, unsigned int n_slip_systems, unsigned int n_Tslip_systems, Vector<double> s_alpha_tau, FullMatrix<double> Dmat, FullMatrix<double> CE_tau_trial, Vector<double> W_kh_t1, Vector<double> W_kh_t2, double hb1, double hb2, double mb1, double mb2, double rb1, double rb2, double bb1, double bb2)
{

double stpmax,tolx,alf,TestNan ;
unsigned int chck,errnum,itr1;

double magsrch,slope,tst,temp,lam_min,lam;

TestNan=0;
tolx = 1.0e-15 ;
alf = 1.0e-4 ;
chck = 0 ;

magsrch = srchdir.l2_norm();

stpmax = magsrch ;

if(magsrch > stpmax)
srchdir.equ(stpmax/magsrch,srchdir) ;


slope = 0.0 ;
for(unsigned int i=0 ; i<n ; i++){
	slope = slope + srchdir(i)*gradFold(i) ;
}

errnum = 0 ;

if(slope > 0.0)
	errnum = 1 ;

tst = 0.0 ;
for(unsigned int i=0 ; i<n ; i++){
	temp = fabs(srchdir(i))/std::max(fabs(stateold(i)),1.0) ;
	if(temp > tst)
		tst = temp ;
}


lam_min = tolx/tst ;
lam = 1.0 ;
itr1 = 0 ;

FullMatrix<double> Tstar(dim,dim),LP_acc(dim,dim),LP_acc2(dim,dim) ;
FullMatrix<double> temp1(dim,dim),temp2(dim,dim),temp3(dim,dim) ;
Vector<double> Tstar_vec(2*dim),bck_strs(2*n_Tslip_systems),delgam(n_Tslip_systems),W_kh_tau(n_Tslip_systems),resval(2*dim+2*n_Tslip_systems) ;
Vector<double> W_kh_tau1(n_Tslip_systems), W_kh_tau2(n_Tslip_systems) ;
Vector<double> vtmp1(2*dim), vtmp2(2*dim), vtmp3(2*dim);
double loc_res_mag,tmplam,rhs1,rhs2,f2,lam2,acf,bcf,disc ;
double sctmp1,sgnm ;

while(itr1 < 300)
{
	statenew = 0 ;
	statenew.add(1.0,stateold,lam,srchdir) ;
	for(unsigned int i=0 ; i<2*dim ; i++){
		Tstar_vec(i) = statenew(i) ;
	}
	for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
		W_kh_tau1(i) = statenew(i+2*dim) ;
		W_kh_tau2(i) = statenew(i+2*dim+n_Tslip_systems) ;
        W_kh_tau(i) = W_kh_tau1(i) + W_kh_tau2(i)  ;
	}

	matform(Tstar,Tstar_vec) ;
	LP_acc = 0 ;
    delgam = 0 ;

	for(unsigned int i=0 ; i<n_Tslip_systems ; i++){
		for (unsigned int j = 0;j < dim;j++) {
     		     for (unsigned int k = 0;k < dim;k++) {
   		         temp1[j][k]=SCHMID_TENSOR1[dim*i + j][k];

        	  	}
       		 }

	Tstar.mTmult(temp2,temp1);
    sctmp1=temp2.trace();

	if(i<n_slip_systems){ // For slip systems due to symmetry of slip
          if((sctmp1-W_kh_tau(i))<0)
          sgnm=-1;
          else
          sgnm=1 ;
        }
        else               // For twin systems due to asymmetry of slip
        {
          if((sctmp1-W_kh_tau(i))<=0)
          sgnm=0;
          else
          sgnm=1 ;
        }

        delgam(i) = delgam_ref*pow(fabs((sctmp1 - W_kh_tau(i))/s_alpha_tau(i)),1.0/strexp)*sgnm ;

		LP_acc.add(delgam(i),temp1) ;

		resval(i+2*dim) = W_kh_tau1(i) - W_kh_t1(i) - hb1*delgam(i) + rb1*pow(fabs(W_kh_tau1(i)/bb1),mb1)*W_kh_tau1(i)*fabs(delgam(i)) ;
        resval(i+2*dim+n_Tslip_systems) = W_kh_tau2(i) - W_kh_t2(i) - hb2*delgam(i) + rb2*pow(fabs(W_kh_tau2(i)/bb2),mb2)*W_kh_tau2(i)*fabs(delgam(i)) ;

	}
	LP_acc.equ(-1.0,LP_acc) ;
	LP_acc2.equ(1.0,matrixExponential(LP_acc)) ;
	LP_acc.equ(1.0,LP_acc2) ;



	vtmp1 = vecform(Tstar) ;
	CE_tau_trial.mmult(temp1,LP_acc) ;
	LP_acc.Tmmult(temp2,temp1) ;
	temp1 = IdentityMatrix(dim) ;
	   temp3 = 0 ;
	   temp3.add(0.5,temp2,-0.5,temp1) ;
	   Dmat.vmult(vtmp2,vecform(temp3)) ;
	   vtmp3 = 0 ;
	   vtmp3.add(1.0,vtmp1,-1.0,vtmp2) ;

	   for(unsigned int i=0; i<2*dim ; i++){
		   resval(i) = vtmp3(i) ;

	   }


     //loc_res_mag = resval.l2_norm() ;

    	//loc_res_mag = 0.5*loc_res_mag*loc_res_mag ;
     loc_res_mag=0;
     for(unsigned int i=0; i<2*dim+2*n_Tslip_systems ; i++){
       loc_res_mag = loc_res_mag+0.5*resval(i)*resval(i) ;
     }
     if (TestNan==1){
       break ;
     }
     if ((std::isnan(loc_res_mag))||(std::isinf(loc_res_mag))){
       lam=0.001*lam;
       TestNan=1;
     }
     else{
	// Lots of conditional statements involved here
	if(lam < lam_min){
		statenew.equ(1.0,stateold) ;
		chck = 1 ;
		break ;

	}
	else if((loc_res_mag - Fold - alf*lam*slope)<=0)
	{
		break ;
	}
	else{

	if(fabs(lam - 1.0) < 1.0e-10)
	{
		tmplam = -1.0*slope/(2.0*(loc_res_mag - Fold - slope )) ;


	}
	else{
		rhs1 = loc_res_mag - Fold - lam*slope ;
		rhs2 = f2 - Fold - lam2*slope ;
		acf = (rhs1/(lam*lam) - rhs2/(lam2*lam2))/(lam - lam2)  ;
		bcf = (-lam2*rhs1/(lam*lam) + lam*rhs2/(lam2*lam2))/(lam - lam2) ;

		if(fabs(acf)<1.0e-10){
			tmplam = -1.0*slope/(2.0*bcf) ;

		}

		else{
			disc = bcf*bcf - 3.0*acf*slope ;

			if(disc<0)
				tmplam = 0.5*lam ;
			else if(bcf<=0)
			tmplam = (-bcf + sqrt(disc))/(3.0*acf) ;
			else
				tmplam = -slope/(bcf  + sqrt(disc)) ;


		}

		if(tmplam>0.5*lam)
			tmplam = 0.5*lam ;
	}



	}


	lam2 = lam ;
	f2 = loc_res_mag ;
	lam = std::max(tmplam,0.1*lam) ;
}
	itr1 = itr1 + 1 ;



} //end while


} // end  function




#include "../../../include/crystalPlasticity_template_instantiations.h"
