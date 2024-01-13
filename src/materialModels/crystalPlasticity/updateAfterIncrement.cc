#include "../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

template <int dim>
void crystalPlasticity<dim>::updateAfterIncrement()
{
  this->updateAfterIncrementBase();
	local_F_r=0.0;
	local_F_s=0.0;
	local_F_e = 0.0;
	unsigned int CheckBufferRegion,dimBuffer;
	double lowerBuffer,upperBuffer,workDensity_Element1,workDensity_Element2,workDensity_Element1_Tr,workDensity_Element2_Tr,det_Fe,det_Fp;
	Point<dim> pnt2;
	Vector<double> userDefinedAverageOutput,local_userDefinedAverageOutput;
  FullMatrix<double> F_pl(dim,dim),F_el(dim,dim),C_pl(dim,dim),E_pl(dim,dim),temp1(dim,dim),CauchyStress_cell(dim,dim),P_LastIter(dim,dim),S_LastIter(dim,dim),E_pl_cell(dim,dim),E_cell(dim,dim),dE_cell(dim,dim);
  FullMatrix<double> F_lastIter(dim,dim),deltaF(dim,dim),deltaE(dim,dim);

	if (this->userInputs.flagUserDefinedAverageOutput){
		userDefinedAverageOutput.reinit(this->userInputs.numberUserDefinedAverageOutput);
		userDefinedAverageOutput=0;
		local_userDefinedAverageOutput.reinit(this->userInputs.numberUserDefinedAverageOutput);
		local_userDefinedAverageOutput=0;
	}

	QGauss<dim>  quadrature(this->userInputs.quadOrder);
	FEValues<dim> fe_values(this->FE, quadrature, update_quadrature_points | update_gradients | update_JxW_values);
	const unsigned int num_quad_points = quadrature.size();
	const unsigned int   dofs_per_cell = this->FE.dofs_per_cell;
	std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
  Vector<double> Ulocal(dofs_per_cell);
  //Vector<double> Ulocal_lastIter(dofs_per_cell), Rlocal (dofs_per_cell)
	if (this->userInputs.flagTaylorModel){
		if(initCalled == false){
			if(this->userInputs.enableAdvancedTwinModel){
				init2(num_quad_points);
			}
			else{
				init(num_quad_points);
			}
		}
	}
	//loop over elements
	unsigned int cellID = 0;
	typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			fe_values.reinit(cell);
			//loop over quadrature points
			cell->set_user_index(fe_values.get_cell()->user_index());
			cell->get_dof_indices(local_dof_indices);

			//////////Buffer layer feature/////////////
			if (this->userInputs.flagBufferLayer){
				pnt2=cell->center();
				dimBuffer=this->userInputs.dimBufferLayer;
				lowerBuffer=this->userInputs.lowerBufferLayer;
				upperBuffer=this->userInputs.upperBufferLayer;
				if ((pnt2[dimBuffer]>=lowerBuffer)&&(pnt2[dimBuffer]<=upperBuffer)){
					CheckBufferRegion=1;
				}
				else{
					CheckBufferRegion=0;
				}
			}
			else{
				CheckBufferRegion=1;
			}
			/////////////////////////////////////////////////////

      Ulocal = 0.0;
      //Rlocal = 0.0; Ulocal_lastIter = 0.0;
			if (!this->userInputs.flagTaylorModel){
				for (unsigned int i = 0; i < dofs_per_cell; i++) {
					Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
				}
			}
      CauchyStress_cell=0;E_pl_cell=0;E_cell=0;dE_cell=0;P_LastIter=0;S_LastIter=0;workDensity_Element1=0;workDensity_Element2=0;workDensity_Element1_Tr=0;workDensity_Element2_Tr=0;
      det_Fp=0;det_Fe=0;

			for (unsigned int q = 0; q < num_quad_points; ++q) {
				//Get deformation gradient
				F = 0.0;
				if (this->userInputs.flagTaylorModel){
					F =this->Fprev;
				}
				else{
					for (unsigned int d = 0; d < dofs_per_cell; ++d) {
						unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
						for (unsigned int j = 0; j < dim; ++j) {
							F[i][j] += Ulocal(d)*fe_values.shape_grad(d, q)[j]; // u_{i,j}= U(d)*N(d)_{,j}, where d is the DOF correonding to the i'th dimension
						}
					}
					for (unsigned int i = 0; i < dim; ++i) {
						F[i][i] += 1;
					}
				}

        F_lastIter = 0.0;deltaF=0.0;deltaE=0.0;
        F_lastIter=F_lastIter_Global[cellID][q];

        for (unsigned int i = 0; i < dim; ++i) {
          for (unsigned int j = 0; j < dim; ++j) {
            deltaF[i][j]=F[i][j]-F_lastIter[i][j];
          }
        }

				//Update strain, stress, and tangent for current time step/quadrature point
				calculatePlasticity(cellID, q, 0);

        /////////I assigned Fe to F just for postprocessing output for Aaron
        stateVar_iter[cellID][q][62]=deltaF[0][0];
        stateVar_iter[cellID][q][63]=deltaF[0][1];
        stateVar_iter[cellID][q][64]=deltaF[0][2];
        stateVar_iter[cellID][q][65]=deltaF[1][0];
        stateVar_iter[cellID][q][66]=deltaF[1][1];
        stateVar_iter[cellID][q][67]=deltaF[1][2];
        stateVar_iter[cellID][q][68]=deltaF[2][0];
        stateVar_iter[cellID][q][69]=deltaF[2][1];
        stateVar_iter[cellID][q][70]=deltaF[2][2];

        stateVar_iter[cellID][q][71]=F[0][0];
        stateVar_iter[cellID][q][72]=F[0][1];
        stateVar_iter[cellID][q][73]=F[0][2];
        stateVar_iter[cellID][q][74]=F[1][0];
        stateVar_iter[cellID][q][75]=F[1][1];
        stateVar_iter[cellID][q][76]=F[1][2];
        stateVar_iter[cellID][q][77]=F[2][0];
        stateVar_iter[cellID][q][78]=F[2][1];
        stateVar_iter[cellID][q][79]=F[2][2];
        /////////I assigned Fe to F just for postprocessing output for Aaron

				FullMatrix<double> temp,temp_lastIter,temp3,temp4, C_tau(dim, dim), E_tau(dim, dim), b_tau(dim, dim), C_tau_lastIter(dim, dim), E_tau_lastIter(dim, dim);
				Vector<double> temp2;
				temp.reinit(dim, dim); temp = 0.0;
        temp_lastIter.reinit(dim, dim); temp_lastIter = 0.0;
				temp2.reinit(dim); temp2 = 0.0;
				temp3.reinit(dim, dim); temp3 = 0.0;
				temp4.reinit(dim, dim); temp4 = 0.0;
				C_tau_lastIter = 0.0;
        E_tau_lastIter= 0.0;
        C_tau = 0.0;
        E_tau= 0.0;
				temp = F;
        temp_lastIter = F_lastIter;
        F.Tmmult(C_tau, temp);
				F_lastIter.Tmmult(C_tau_lastIter, temp_lastIter);
				F.mTmult(b_tau, temp);
				//E_tau = CE_tau;
				temp = IdentityMatrix(dim);
				for (unsigned int i = 0;i<dim;i++) {
					temp2[i] = 0.5*log(b_tau[i][i])*fe_values.JxW(q);
					for (unsigned int j = 0;j<dim;j++) {
						E_tau[i][j] = 0.5*(C_tau[i][j] - temp[i][j]);
            E_tau_lastIter[i][j] = 0.5*(C_tau_lastIter[i][j] - temp[i][j]);
						temp3[i][j] = T[i][j] * fe_values.JxW(q);
						temp4[i][j] = E_tau[i][j]*fe_values.JxW(q);
					}
				}

        for (unsigned int i = 0; i < dim; ++i) {
          for (unsigned int j = 0; j < dim; ++j) {
            deltaE[i][j]=E_tau[i][j]-E_tau_lastIter[i][j];
          }
        }
/////////This CauchyStress is First Piola now
				CauchyStress[cellID][q]=T;

///////////Calculation of plastic components of Green strain tensor
        F_pl=0;C_pl=0;E_pl=0;temp1=0;F_el=0;
        F_pl=Fp_iter[cellID][q];
        det_Fp += F_pl.determinant()/8;  //The cell value of det_Fp=1/8*Sum(det_Fp for all quadratures) where 8 is the number of quadrature we used for first order elements.
        F_pl.Tmmult(C_pl,F_pl);
        temp1= IdentityMatrix(dim) ;
        E_pl.add(0.5,C_pl,-0.5,temp1) ;

        F_el=Fe_iter[cellID][q];
        det_Fe += F_el.determinant()/8; //The cell value of det_Fe=1/8*Sum(det_Fe for all quadratures) where 8 is the number of quadrature we used for first order elements.

/////////////////////////////////////////////////

////////////////Calculation of element nodal force -> To be use in work density calculation////////
        P_LastIter=FirstPiolaStress[cellID][q];
        S_LastIter=SecondPiolaStress[cellID][q];
      //  for (unsigned int d=0; d<dofs_per_cell; ++d) {
      //    unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
      //    for (unsigned int j = 0; j < dim; j++){
      //      Rlocal(d) +=  fe_values.shape_grad(d, q)[j]*((P[i][j]+P_LastIter[i][j])/2)*fe_values.JxW(q);
      //    }
      //  }

        for (unsigned int i = 0; i < dim; ++i) {
          for (unsigned int j = 0; j < dim; ++j) {
////////////////////////Rectangular integration rule////////////
            workDensity_Element1+=P[i][j]*deltaF[i][j]*fe_values.JxW(q);
            workDensity_Element2+=S[i][j]*deltaE[i][j]*fe_values.JxW(q);
/////////////////////Trapezoidal integration rule///////////////
            workDensity_Element1_Tr+=((P[i][j]+P_LastIter[i][j])/2)*deltaF[i][j]*fe_values.JxW(q);
            workDensity_Element2_Tr+=((S[i][j]+S_LastIter[i][j])/2)*deltaE[i][j]*fe_values.JxW(q);
          }
        }

        FirstPiolaStress[cellID][q]=P;
        SecondPiolaStress[cellID][q]=S;
///////////////////////////////////////////

        //CauchyStress_cell.add(1,CauchyStress_cell,0.125,T) ;
        //E_pl_cell.add(1,E_pl_cell,0.125,E_pl) ;
        //E_cell.add(1,E_cell,0.125,E_tau) ;
        CauchyStress_cell.add(0.125,T) ;
        E_pl_cell.add(0.125,E_pl) ;
        E_cell.add(0.125,E_tau) ;
        dE_cell.add(0.125,deltaE) ;

				if (this->userInputs.enableAdvRateDepModel){
					for(unsigned int i=0; i<dim ; i++){
						for(unsigned int j=0 ; j<dim ; j++){
							TinterStress_diff[cellID][q][i][j] = T_inter[i][j] - TinterStress[cellID][q][i][j] ;
							TinterStress[cellID][q][i][j] = T_inter[i][j];
						}
					}
				}

				//calculate von-Mises stress and equivalent strain
				double traceE, traceT, vonmises, eqvstrain;
				FullMatrix<double> deve(dim, dim), devt(dim, dim);


				traceE = E_tau.trace();
				traceT = T.trace();
				temp = IdentityMatrix(3);
				temp.equ(traceE / 3, temp);

				deve = E_tau;
				deve.add(-1.0, temp);

				temp = IdentityMatrix(3);
				temp.equ(traceT / 3, temp);

				devt = T;
				devt.add(-1.0, temp);

				vonmises = devt.frobenius_norm();
				vonmises = sqrt(3.0 / 2.0)*vonmises;
				eqvstrain = deve.frobenius_norm();
				eqvstrain = sqrt(2.0 / 3.0)*eqvstrain;

				//fill in post processing field values
				if (!this->userInputs.enableAdvancedTwinModel){
					twin_ouput[cellID][q]=twin_iter[cellID][q];
				}
				else{
					if (TotaltwinvfK[cellID][q]>=this->userInputs.criteriaTwinVisual){
						twin_ouput[cellID][q]=1;
					}
					else{
						twin_ouput[cellID][q]=0;
					}
				}
				if (this->userInputs.writeOutput){
					this->postprocessValues(cellID, q, 0, 0) = vonmises;
					this->postprocessValues(cellID, q, 1, 0) = eqvstrain;
					this->postprocessValues(cellID, q, 2, 0) = twin_ouput[cellID][q];

////////User Defined Variables for visualization outputs (output_Var1 to output_Var24)////////
					this->postprocessValues(cellID, q, 3, 0) = CauchyStress[cellID][q][0][0];
					this->postprocessValues(cellID, q, 4, 0) = CauchyStress[cellID][q][1][1];
					this->postprocessValues(cellID, q, 5, 0) = CauchyStress[cellID][q][2][2];
					this->postprocessValues(cellID, q, 6, 0) = CauchyStress[cellID][q][1][2];
					this->postprocessValues(cellID, q, 7, 0) = CauchyStress[cellID][q][0][2];
					this->postprocessValues(cellID, q, 8, 0) = CauchyStress[cellID][q][0][1];
					this->postprocessValues(cellID, q, 9, 0) = E_tau[0][0];
					this->postprocessValues(cellID, q, 10, 0) = E_tau[1][1];
					this->postprocessValues(cellID, q, 11, 0) = E_tau[2][2];
					this->postprocessValues(cellID, q, 12, 0) = E_tau[1][2];
					this->postprocessValues(cellID, q, 13, 0) = E_tau[0][2];
					this->postprocessValues(cellID, q, 14, 0) = E_tau[0][1];
					this->postprocessValues(cellID, q, 15, 0) = 0;
					this->postprocessValues(cellID, q, 16, 0) = 0;
					this->postprocessValues(cellID, q, 17, 0) = 0;
					this->postprocessValues(cellID, q, 18, 0) = 0;
					this->postprocessValues(cellID, q, 19, 0) = 0;
					this->postprocessValues(cellID, q, 20, 0) = 0;
					this->postprocessValues(cellID, q, 21, 0) = 0;
					this->postprocessValues(cellID, q, 22, 0) = 0;
					this->postprocessValues(cellID, q, 23, 0) = 0;
					this->postprocessValues(cellID, q, 24, 0) = 0;
					this->postprocessValues(cellID, q, 25, 0) = 0;
					this->postprocessValues(cellID, q, 26, 0) = 0;
				}

				if (CheckBufferRegion==1) {
					local_strain.add(1.0, temp4);
					local_stress.add(1.0, temp3);
					local_microvol = local_microvol + fe_values.JxW(q);
					for(unsigned int i=0;i<this->userInputs.numTwinSystems1;i++){
						local_F_r=local_F_r+twinfraction_iter[cellID][q][i]*fe_values.JxW(q);
					}

					if (!this->userInputs.enableAdvancedTwinModel){
						local_F_e = local_F_e + twin_ouput[cellID][q] * fe_values.JxW(q);
					}
					else{
						local_F_e=local_F_e+ TotaltwinvfK[cellID][q]*fe_values.JxW(q);
					}

					for(unsigned int i=0;i<this->userInputs.numSlipSystems1;i++){
						local_F_s=local_F_s+slipfraction_iter[cellID][q][i]*fe_values.JxW(q);
					}

					if (this->userInputs.flagUserDefinedAverageOutput){
						////One should define their UserDefinedAverageOutput here by defining based on the state Variables
						///In the following equation, "+0" should be substituted by "+targetVariable", where targetVariable is the variable you want to plot the average value.
						//Also, the first equation should be copy and paste depending on the number of variables you want to output,
						//and the integer in paranthesis should be updated accordingly.
						//The maximum number of lines (output variables) are defined in the input file as "set Number of Output Userdefined Average Variable".
						local_userDefinedAverageOutput(0)=local_userDefinedAverageOutput(0)+0;
					}
				}

        F_lastIter_Global[cellID][q]=F;
			}
			if (this->userInputs.writeOutput){

//////////////////Calculation of work density for the cell//////////
      //  workDensity=0;
        //for (unsigned int i = 0; i < dofs_per_cell; i++) {
				//	workDensity[cellID]+=Rlocal[i]*(Ulocal[i]-Ulocal_lastIter[i]);
				//}

        workDensity1[cellID]=workDensity_Element1;
        workDensity2[cellID]=workDensity_Element2;

        workDensity1_Tr[cellID]=workDensity_Element1_Tr;
        workDensity2_Tr[cellID]=workDensity_Element2_Tr;

        workDensityTotal1[cellID]+=workDensity_Element1;
        workDensityTotal2[cellID]+=workDensity_Element2;

        workDensityTotal1_Tr[cellID]+=workDensity_Element1_Tr;
        workDensityTotal2_Tr[cellID]+=workDensity_Element2_Tr;
///////////////////////////////////////////////////////////////////


				this->postprocessValuesAtCellCenters(cellID,0)=cellOrientationMap[cellID];
////////User Defined Variables for visualization outputs for cell_centers (outputoutputCellCenters_Var1 to outputoutputCellCenters_Var24)////////

				//this->postprocessValuesAtCellCenters(cellID,1)= 1;
        this->postprocessValuesAtCellCenters(cellID,1)= CauchyStress_cell[0][0];
				this->postprocessValuesAtCellCenters(cellID,2)= CauchyStress_cell[1][1];
//        this->postprocessValuesAtCellCenters(cellID,2)= 2;
				this->postprocessValuesAtCellCenters(cellID,3)= CauchyStress_cell[2][2];
				this->postprocessValuesAtCellCenters(cellID,4)= CauchyStress_cell[1][2];
				this->postprocessValuesAtCellCenters(cellID,5)= CauchyStress_cell[0][2];
				this->postprocessValuesAtCellCenters(cellID,6)= CauchyStress_cell[0][1];
				this->postprocessValuesAtCellCenters(cellID,7)= E_cell[0][0];
				this->postprocessValuesAtCellCenters(cellID,8)= E_cell[1][1];
				this->postprocessValuesAtCellCenters(cellID,9)= E_cell[2][2];
				this->postprocessValuesAtCellCenters(cellID,10)=E_cell[1][2];
				this->postprocessValuesAtCellCenters(cellID,11)=E_cell[0][2];
				this->postprocessValuesAtCellCenters(cellID,12)=E_cell[0][1];
				this->postprocessValuesAtCellCenters(cellID,13)=E_pl_cell[0][0];
				this->postprocessValuesAtCellCenters(cellID,14)=E_pl_cell[1][1];
				this->postprocessValuesAtCellCenters(cellID,15)=E_pl_cell[2][2];
				this->postprocessValuesAtCellCenters(cellID,16)=E_pl_cell[1][2];
				this->postprocessValuesAtCellCenters(cellID,17)=E_pl_cell[0][2];
				this->postprocessValuesAtCellCenters(cellID,18)=E_pl_cell[0][1];
				this->postprocessValuesAtCellCenters(cellID,19)=workDensity1[cellID];
				this->postprocessValuesAtCellCenters(cellID,20)=workDensity2[cellID];
				this->postprocessValuesAtCellCenters(cellID,21)=workDensityTotal1[cellID];
				this->postprocessValuesAtCellCenters(cellID,22)=workDensityTotal2[cellID];
        this->postprocessValuesAtCellCenters(cellID,23)=workDensity1_Tr[cellID];
				this->postprocessValuesAtCellCenters(cellID,24)=workDensity2_Tr[cellID];
				this->postprocessValuesAtCellCenters(cellID,25)=workDensityTotal1_Tr[cellID];
				this->postprocessValuesAtCellCenters(cellID,26)=workDensityTotal2_Tr[cellID];
				this->postprocessValuesAtCellCenters(cellID,27)=det_Fe;
				this->postprocessValuesAtCellCenters(cellID,28)=det_Fp;
				//this->postprocessValuesAtCellCenters(cellID,29)=Ulocal(0);
				//this->postprocessValuesAtCellCenters(cellID,30)=Ulocal_lastIter(0);
			}

			cellID++;
		}
	}

	//In Case we have twinning
	rotnew_conv=rotnew_iter;

	if (!this->userInputs.enableAdvancedTwinModel){
		//reorient() updates the rotnew_conv.
		reorient();
	}

	//Updating rotnew_iter using rotnew_conv updated by reorient();
	rotnew_iter=rotnew_conv;

	//Update the history variables when convergence is reached for the current increment
	Fe_conv=Fe_iter;
	Fp_conv=Fp_iter;
	s_alpha_conv=s_alpha_iter;
	W_kh_conv = W_kh_iter;
	twinfraction_conv=twinfraction_iter;
	slipfraction_conv=slipfraction_iter;
	rot_conv=rot_iter;
	twin_conv=twin_iter;

	if (this->userInputs.enableUserMaterialModel){
		stateVar_conv=stateVar_iter;
	}

	if (this->userInputs.enableAdvancedTwinModel){
		TwinMaxFlag_conv = TwinMaxFlag_iter;
		NumberOfTwinnedRegion_conv = NumberOfTwinnedRegion_iter;
		ActiveTwinSystems_conv = ActiveTwinSystems_iter;
		TwinFlag_conv = TwinFlag_iter;
		TwinOutputfraction_conv=TwinOutputfraction_iter;
	}





	char buffer[200];

	//////////////////////TabularOutput Start///////////////
	std::vector<unsigned int> tabularTimeInputIncInt;
	std::vector<double> tabularTimeInputInc;
	if (this->userInputs.tabularOutput){

		tabularTimeInputInc=this->userInputs.tabularTimeOutput;
		for(unsigned int i=0;i<this->userInputs.tabularTimeOutput.size();i++){
			tabularTimeInputInc[i]=tabularTimeInputInc[i]/this->delT;
		}

		tabularTimeInputIncInt.resize(this->userInputs.tabularTimeOutput.size(),0);
		///Converting to an integer always rounds down, even if the fraction part is 0.99999999.
		//Hence, I add 0.1 to make sure we always get the correct integer.
		for(unsigned int i=0;i<this->userInputs.tabularTimeOutput.size();i++){
			tabularTimeInputIncInt[i]=int(tabularTimeInputInc[i]+0.1);
		}
	}
	//////////////////////TabularOutput Finish///////////////
	if (this->userInputs.writeQuadratureOutput) {
		if (((!this->userInputs.tabularOutput)&&((this->currentIncrement+1)%this->userInputs.skipQuadratureOutputSteps == 0))||((this->userInputs.tabularOutput)&& (std::count(tabularTimeInputIncInt.begin(), tabularTimeInputIncInt.end(), (this->currentIncrement+1))==1))){
			//copy rotnew to output
			outputQuadrature.clear();
			//loop over elements
			cellID=0;
			cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
			for (; cell!=endc; ++cell) {
				if (cell->is_locally_owned()){
					fe_values.reinit(cell);
					//loop over quadrature points
					for (unsigned int q=0; q<num_quad_points; ++q){
						std::vector<double> temp;
						temp.push_back(cellOrientationMap[cellID]);

						if (!this->userInputs.enableAdvancedTwinModel){
							temp.push_back(phase[cellID][q]);
						}

						temp.push_back(fe_values.JxW(q));

			//			temp.push_back(twin_ouput[cellID][q]);

						temp.push_back(fe_values.get_quadrature_points()[q][0]);
						temp.push_back(fe_values.get_quadrature_points()[q][1]);
						temp.push_back(fe_values.get_quadrature_points()[q][2]);

		//				temp.push_back(rotnew_conv[cellID][q][0]);
		//				temp.push_back(rotnew_conv[cellID][q][1]);
		//				temp.push_back(rotnew_conv[cellID][q][2]);
            temp.push_back(CauchyStress[cellID][q][0][0]);
            temp.push_back(CauchyStress[cellID][q][0][1]);
            temp.push_back(CauchyStress[cellID][q][0][2]);
            temp.push_back(CauchyStress[cellID][q][1][0]);
            temp.push_back(CauchyStress[cellID][q][1][1]);
            temp.push_back(CauchyStress[cellID][q][1][2]);
            temp.push_back(CauchyStress[cellID][q][2][0]);
            temp.push_back(CauchyStress[cellID][q][2][1]);
            temp.push_back(CauchyStress[cellID][q][2][2]);

            temp.push_back(stateVar_conv[cellID][q][62]);
						temp.push_back(stateVar_conv[cellID][q][63]);
						temp.push_back(stateVar_conv[cellID][q][64]);
						temp.push_back(stateVar_conv[cellID][q][65]);
						temp.push_back(stateVar_conv[cellID][q][66]);
						temp.push_back(stateVar_conv[cellID][q][67]);
						temp.push_back(stateVar_conv[cellID][q][68]);
						temp.push_back(stateVar_conv[cellID][q][69]);
						temp.push_back(stateVar_conv[cellID][q][70]);

            temp.push_back(stateVar_conv[cellID][q][71]);
						temp.push_back(stateVar_conv[cellID][q][72]);
						temp.push_back(stateVar_conv[cellID][q][73]);
						temp.push_back(stateVar_conv[cellID][q][74]);
						temp.push_back(stateVar_conv[cellID][q][75]);
						temp.push_back(stateVar_conv[cellID][q][76]);
						temp.push_back(stateVar_conv[cellID][q][77]);
						temp.push_back(stateVar_conv[cellID][q][78]);
						temp.push_back(stateVar_conv[cellID][q][79]);

			//			temp.push_back(Fp_conv[cellID][q][0][0]);
			//			temp.push_back(Fp_conv[cellID][q][1][1]);
			//			temp.push_back(Fp_conv[cellID][q][2][2]);
			//			temp.push_back(Fp_conv[cellID][q][0][1]);
			//			temp.push_back(Fp_conv[cellID][q][0][2]);
			//			temp.push_back(Fp_conv[cellID][q][1][0]);
			//			temp.push_back(Fp_conv[cellID][q][1][2]);
			//			temp.push_back(Fp_conv[cellID][q][2][0]);
			//			temp.push_back(Fp_conv[cellID][q][2][1]);




					//	if (this->userInputs.enableAdvRateDepModel){
			//				temp.push_back(TinterStress[cellID][q][0][0]);
			//				temp.push_back(TinterStress[cellID][q][1][1]);
			//				temp.push_back(TinterStress[cellID][q][2][2]);
			//				temp.push_back(TinterStress[cellID][q][0][1]);
			//				temp.push_back(TinterStress[cellID][q][0][2]);
			//				temp.push_back(TinterStress[cellID][q][1][0]);
			//				temp.push_back(TinterStress[cellID][q][1][1]);
			//				temp.push_back(TinterStress[cellID][q][2][0]);
			//				temp.push_back(TinterStress[cellID][q][2][1]);

			//				temp.push_back(TinterStress_diff[cellID][q][0][0]);
			//				temp.push_back(TinterStress_diff[cellID][q][1][1]);
			//				temp.push_back(TinterStress_diff[cellID][q][2][2]);
			//				temp.push_back(TinterStress_diff[cellID][q][0][1]);
			//				temp.push_back(TinterStress_diff[cellID][q][0][2]);
			//				temp.push_back(TinterStress_diff[cellID][q][1][0]);
			//				temp.push_back(TinterStress_diff[cellID][q][1][1]);
			//				temp.push_back(TinterStress_diff[cellID][q][2][0]);
			//				temp.push_back(TinterStress_diff[cellID][q][2][1]);
			//			}

			//		//	temp.push_back(slipfraction_conv[cellID][q][0]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][1]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][2]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][3]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][4]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][5]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][6]);
			//		//	temp.push_back(slipfraction_conv[cellID][q][7]);
					//	temp.push_back(slipfraction_conv[cellID][q][8]);
					//	temp.push_back(slipfraction_conv[cellID][q][9]);
					//	temp.push_back(slipfraction_conv[cellID][q][10]);
					//	temp.push_back(slipfraction_conv[cellID][q][11]);
					//	temp.push_back(slipfraction_conv[cellID][q][12]);
					//	temp.push_back(slipfraction_conv[cellID][q][13]);
					//	temp.push_back(slipfraction_conv[cellID][q][14]);
					//	temp.push_back(slipfraction_conv[cellID][q][15]);
					//	temp.push_back(slipfraction_conv[cellID][q][16]);
					//	temp.push_back(slipfraction_conv[cellID][q][17]);
					//	temp.push_back(slipfraction_conv[cellID][q][18]);
					//	temp.push_back(slipfraction_conv[cellID][q][19]);
					//	temp.push_back(slipfraction_conv[cellID][q][20]);
					//	temp.push_back(slipfraction_conv[cellID][q][21]);
					//	temp.push_back(slipfraction_conv[cellID][q][22]);
					//	temp.push_back(slipfraction_conv[cellID][q][23]);
					//	temp.push_back(slipfraction_conv[cellID][q][24]);
					//	temp.push_back(slipfraction_conv[cellID][q][25]);
					//	temp.push_back(slipfraction_conv[cellID][q][26]);
					//	temp.push_back(slipfraction_conv[cellID][q][27]);
					//	temp.push_back(slipfraction_conv[cellID][q][28]);
					//	temp.push_back(slipfraction_conv[cellID][q][29]);
					//	temp.push_back(slipfraction_conv[cellID][q][30]);
					//	temp.push_back(slipfraction_conv[cellID][q][31]);
					//	temp.push_back(slipfraction_conv[cellID][q][32]);
					//	temp.push_back(slipfraction_conv[cellID][q][33]);
					//	temp.push_back(slipfraction_conv[cellID][q][34]);
					//	temp.push_back(slipfraction_conv[cellID][q][35]);
					//	temp.push_back(slipfraction_conv[cellID][q][36]);
					//	temp.push_back(slipfraction_conv[cellID][q][37]);
					//	temp.push_back(slipfraction_conv[cellID][q][38]);
					//	temp.push_back(slipfraction_conv[cellID][q][39]);
					//	temp.push_back(slipfraction_conv[cellID][q][40]);
					//	temp.push_back(slipfraction_conv[cellID][q][41]);
					//	temp.push_back(slipfraction_conv[cellID][q][42]);
					//	temp.push_back(slipfraction_conv[cellID][q][43]);
					//	temp.push_back(slipfraction_conv[cellID][q][44]);
					//	temp.push_back(slipfraction_conv[cellID][q][45]);
					//	temp.push_back(slipfraction_conv[cellID][q][46]);
					//	temp.push_back(slipfraction_conv[cellID][q][47]);
					//	temp.push_back(slipfraction_conv[cellID][q][48]);
					//	temp.push_back(slipfraction_conv[cellID][q][49]);
					//	temp.push_back(slipfraction_conv[cellID][q][50]);
					//	temp.push_back(slipfraction_conv[cellID][q][51]);
					//	temp.push_back(slipfraction_conv[cellID][q][52]);
					//	temp.push_back(slipfraction_conv[cellID][q][53]);
					//	temp.push_back(slipfraction_conv[cellID][q][54]);
					//	temp.push_back(slipfraction_conv[cellID][q][55]);
					//	temp.push_back(slipfraction_conv[cellID][q][56]);
					//	temp.push_back(slipfraction_conv[cellID][q][57]);
					//	temp.push_back(slipfraction_conv[cellID][q][58]);
					//	temp.push_back(slipfraction_conv[cellID][q][59]);
					//	temp.push_back(slipfraction_conv[cellID][q][60]);
					//	temp.push_back(slipfraction_conv[cellID][q][61]);
					//	temp.push_back(slipfraction_conv[cellID][q][62]);
					//	temp.push_back(slipfraction_conv[cellID][q][63]);
					//	temp.push_back(slipfraction_conv[cellID][q][64]);
					//	temp.push_back(slipfraction_conv[cellID][q][65]);
					//	temp.push_back(slipfraction_conv[cellID][q][66]);
					//	temp.push_back(slipfraction_conv[cellID][q][67]);
					//	temp.push_back(slipfraction_conv[cellID][q][68]);
					//	temp.push_back(slipfraction_conv[cellID][q][69]);
					//	temp.push_back(slipfraction_conv[cellID][q][70]);
					//	temp.push_back(slipfraction_conv[cellID][q][71]);
					//	temp.push_back(slipfraction_conv[cellID][q][72]);
					//	temp.push_back(slipfraction_conv[cellID][q][73]);
					//	temp.push_back(slipfraction_conv[cellID][q][74]);
					//	temp.push_back(slipfraction_conv[cellID][q][75]);
					//	temp.push_back(slipfraction_conv[cellID][q][76]);
					//	temp.push_back(slipfraction_conv[cellID][q][77]);
					//	temp.push_back(slipfraction_conv[cellID][q][78]);
					//	temp.push_back(slipfraction_conv[cellID][q][79]);
					//	temp.push_back(slipfraction_conv[cellID][q][80]);
					//	temp.push_back(slipfraction_conv[cellID][q][81]);
					//	temp.push_back(slipfraction_conv[cellID][q][82]);
					//	temp.push_back(slipfraction_conv[cellID][q][83]);


					//	temp.push_back(twinfraction_conv[cellID][q][0]);
					//	temp.push_back(twinfraction_conv[cellID][q][1]);
					//	temp.push_back(twinfraction_conv[cellID][q][2]);
					//	temp.push_back(twinfraction_conv[cellID][q][3]);
					//	temp.push_back(twinfraction_conv[cellID][q][4]);
					//	temp.push_back(twinfraction_conv[cellID][q][5]);

					//	if (this->userInputs.enableAdvancedTwinModel){
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][0]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][1]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][2]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][3]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][4]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][5]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][6]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][7]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][8]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][9]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][10]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][11]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][12]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][13]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][14]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][15]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][16]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][17]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][18]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][19]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][20]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][21]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][22]);
					//		temp.push_back(TwinOutputfraction_conv[cellID][q][23]);
					//	}

				////		if (this->userInputs.enableUserMaterialModel){
		//					temp.push_back(stateVar_conv[cellID][q][0]);
		//					temp.push_back(stateVar_conv[cellID][q][1]);
		//					temp.push_back(stateVar_conv[cellID][q][2]);
		//					temp.push_back(stateVar_conv[cellID][q][3]);
		//					temp.push_back(stateVar_conv[cellID][q][4]);
		//					temp.push_back(stateVar_conv[cellID][q][5]);
		//					temp.push_back(stateVar_conv[cellID][q][6]);
		//					temp.push_back(stateVar_conv[cellID][q][7]);
		//					temp.push_back(stateVar_conv[cellID][q][8]);
		//					temp.push_back(stateVar_conv[cellID][q][9]);
		//					temp.push_back(stateVar_conv[cellID][q][10]);
		//					temp.push_back(stateVar_conv[cellID][q][11]);
		//					temp.push_back(stateVar_conv[cellID][q][12]);
		//					temp.push_back(stateVar_conv[cellID][q][13]);
		//					temp.push_back(stateVar_conv[cellID][q][14]);
		//					temp.push_back(stateVar_conv[cellID][q][15]);
		//					temp.push_back(stateVar_conv[cellID][q][16]);
		//					temp.push_back(stateVar_conv[cellID][q][17]);
		//					temp.push_back(stateVar_conv[cellID][q][18]);
		//					temp.push_back(stateVar_conv[cellID][q][19]);
		//					temp.push_back(stateVar_conv[cellID][q][20]);
		//					temp.push_back(stateVar_conv[cellID][q][21]);
		//					temp.push_back(stateVar_conv[cellID][q][22]);
		//					temp.push_back(stateVar_conv[cellID][q][23]);
		//					temp.push_back(stateVar_conv[cellID][q][24]);
		//					temp.push_back(stateVar_conv[cellID][q][25]);
		//					temp.push_back(stateVar_conv[cellID][q][26]);
		//					temp.push_back(stateVar_conv[cellID][q][27]);
		//					temp.push_back(stateVar_conv[cellID][q][28]);
		//					temp.push_back(stateVar_conv[cellID][q][29]);
		//					temp.push_back(stateVar_conv[cellID][q][30]);
		//					temp.push_back(stateVar_conv[cellID][q][31]);
		//					temp.push_back(stateVar_conv[cellID][q][32]);
		//					temp.push_back(stateVar_conv[cellID][q][33]);
		//					temp.push_back(stateVar_conv[cellID][q][34]);
		//					temp.push_back(stateVar_conv[cellID][q][35]);
		//					temp.push_back(stateVar_conv[cellID][q][36]);
		//					temp.push_back(stateVar_conv[cellID][q][37]);
		//					temp.push_back(stateVar_conv[cellID][q][38]);
		//					temp.push_back(stateVar_conv[cellID][q][39]);
		//					temp.push_back(stateVar_conv[cellID][q][40]);
		//					temp.push_back(stateVar_conv[cellID][q][41]);
		//					temp.push_back(stateVar_conv[cellID][q][42]);
		//					temp.push_back(stateVar_conv[cellID][q][43]);
		//					temp.push_back(stateVar_conv[cellID][q][44]);
		//					temp.push_back(stateVar_conv[cellID][q][45]);
		//					temp.push_back(stateVar_conv[cellID][q][46]);
		//					temp.push_back(stateVar_conv[cellID][q][47]);
		//					temp.push_back(stateVar_conv[cellID][q][48]);
		//					temp.push_back(stateVar_conv[cellID][q][49]);
		//					temp.push_back(stateVar_conv[cellID][q][50]);
		//				}

						addToQuadratureOutput(temp);

					}
					cellID++;
				}
			}

			writeQuadratureOutput(this->userInputs.outputDirectory, this->currentIncrement);
		}
	}

	if (this->userInputs.writeGrainAveragedOutput) {
		if (((!this->userInputs.tabularOutput)&&((this->currentIncrement+1)%this->userInputs.skipGrainAveragedOutputSteps == 0))||((this->userInputs.tabularOutput)&& (std::count(tabularTimeInputIncInt.begin(), tabularTimeInputIncInt.end(), (this->currentIncrement+1))==1))){
			//loop over elements
			unsigned int numberOfGrainAverageDataOutput=this->userInputs.numberOfGrainAverageDataOutput;
			cellID=0;
			FullMatrix<double> grainAveragaData(orientations.numberOfGrains,numberOfGrainAverageDataOutput);
			FullMatrix<double> global_grainAveragaData(orientations.numberOfGrains,numberOfGrainAverageDataOutput);

			grainAveragaData=0;
			global_grainAveragaData=0;
		  unsigned int grainID;
			cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
			for (; cell!=endc; ++cell) {
				if (cell->is_locally_owned()){
					fe_values.reinit(cell);

					if (this->userInputs.flagBufferLayer){
						pnt2=cell->center();
						dimBuffer=this->userInputs.dimBufferLayer;
						lowerBuffer=this->userInputs.lowerBufferLayer;
						upperBuffer=this->userInputs.upperBufferLayer;
						if ((pnt2[dimBuffer]>=lowerBuffer)&&(pnt2[dimBuffer]<=upperBuffer)){
							CheckBufferRegion=1;
						}
						else{
							CheckBufferRegion=0;
						}
					}
					else{
						CheckBufferRegion=1;
					}

					if (CheckBufferRegion==1) {
						grainID=cellOrientationMap[cellID]-1;

						//loop over quadrature points
						for (unsigned int q=0; q<num_quad_points; ++q){
							std::vector<double> temp;

							temp.push_back(fe_values.JxW(q));

							temp.push_back(CauchyStress[cellID][q][0][0]*fe_values.JxW(q));
							temp.push_back(CauchyStress[cellID][q][1][1]*fe_values.JxW(q));
							temp.push_back(CauchyStress[cellID][q][2][2]*fe_values.JxW(q));
							temp.push_back(CauchyStress[cellID][q][0][1]*fe_values.JxW(q));
							temp.push_back(CauchyStress[cellID][q][0][2]*fe_values.JxW(q));
							temp.push_back(CauchyStress[cellID][q][1][2]*fe_values.JxW(q));

							if (temp.size()!=numberOfGrainAverageDataOutput){
								this->pcout << "The size of numberOfGrainAverageDataOutput defined in prm.prm is not correct\n";
				      	exit(1);
							}

							for (unsigned int i=0;i<numberOfGrainAverageDataOutput;i++){
								grainAveragaData(grainID,i)=grainAveragaData(grainID,i)+temp[i];
							}
						}
					}
					cellID++;
				}
			}


			for(unsigned int i=0;i<orientations.numberOfGrains;i++){
				for(unsigned int j=0;j<numberOfGrainAverageDataOutput;j++){
					global_grainAveragaData[i][j]=Utilities::MPI::sum(grainAveragaData[i][j],this->mpi_communicator);
				}
			}

			for(unsigned int i=0;i<orientations.numberOfGrains;i++){
				if (fabs(global_grainAveragaData[i][0])>1e-15){
					for(unsigned int j=1;j<numberOfGrainAverageDataOutput;j++){
						global_grainAveragaData[i][j]=global_grainAveragaData[i][j]/global_grainAveragaData[i][0];
					}
				}
			}

			std::string dir_GrainAverage(this->userInputs.outputDirectory);

			if(Utilities::MPI::this_mpi_process(this->mpi_communicator)==0){
				dir_GrainAverage+="/";
				std::ofstream outputFile_GrainAverage;
				dir_GrainAverage += std::string("GrainAverage");
				dir_GrainAverage += std::to_string(this->currentIncrement);
				dir_GrainAverage += std::string(".csv");
				outputFile_GrainAverage.open(dir_GrainAverage.c_str(),std::fstream::app);
				for(unsigned int i=0;i<orientations.numberOfGrains;i++){
					for(unsigned int j=0;j<numberOfGrainAverageDataOutput-1;j++){
						outputFile_GrainAverage<< std::setprecision(9) <<global_grainAveragaData[i][j]<<",";
					}
					outputFile_GrainAverage<< std::setprecision(9) <<global_grainAveragaData[i][numberOfGrainAverageDataOutput-1];
					outputFile_GrainAverage <<'\n';
				}
				outputFile_GrainAverage.close();
			}
		}
	}


	microvol=Utilities::MPI::sum(local_microvol,this->mpi_communicator);

	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=0;j<dim;j++){
			global_strain[i][j]=Utilities::MPI::sum(local_strain[i][j],this->mpi_communicator)/microvol;
			global_stress[i][j]=Utilities::MPI::sum(local_stress[i][j],this->mpi_communicator)/microvol;
		}
	}
	if (!this->userInputs.enableMultiphase){
		F_e = Utilities::MPI::sum(local_F_e , this->mpi_communicator)/ microvol;
		F_r=Utilities::MPI::sum(local_F_r,this->mpi_communicator)/microvol;
		F_s=Utilities::MPI::sum(local_F_s,this->mpi_communicator)/microvol;
	}
	else {
		F_e=0;F_r=0;F_s=0;
	}

	if (this->userInputs.flagUserDefinedAverageOutput){
		for(unsigned int i=0;i<this->userInputs.numberUserDefinedAverageOutput;i++){
			userDefinedAverageOutput(i)=Utilities::MPI::sum(local_userDefinedAverageOutput(i),this->mpi_communicator)/microvol;
		}
	}

	//check whether to write stress and strain data to file
	//write stress and strain data to file
	std::string dir(this->userInputs.outputDirectory);
	if(Utilities::MPI::this_mpi_process(this->mpi_communicator)==0){
		dir+="/";
		std::ofstream outputFile;
		dir += std::string("stressstrain.txt");

		if(this->currentIncrement==0){
			outputFile.open(dir.c_str());
			outputFile << "Exx" << '\t' << "Eyy" << '\t' << "Ezz" << '\t' << "Eyz" << '\t' << "Exz" << '\t' << "Exy"
			<< '\t' << "Txx" << '\t' << "Tyy" << '\t' << "Tzz" << '\t' << "Tyz" << '\t' << "Txz" << '\t'
			<< "Txy" << '\t' << "TwinRealVF" << '\t' << "TwinMade" << '\t' << "SlipTotal";
			if (this->userInputs.enableIndentationBCs)
			outputFile << "\tInd_Load\tInd_U";
			if (this->userInputs.flagUserDefinedAverageOutput){
				for(unsigned int i=0;i<this->userInputs.numberUserDefinedAverageOutput;i++){
					outputFile <<'\t'<<"userDefined"<<i;
				}
			}
			outputFile <<'\n';
			outputFile.close();
		}
		outputFile.open(dir.c_str(),std::fstream::app);
		outputFile << global_strain[0][0] << '\t' << global_strain[1][1] << '\t' << global_strain[2][2] << '\t'
		<< global_strain[1][2] << '\t' << global_strain[0][2] << '\t' << global_strain[0][1] << '\t'
		<< global_stress[0][0] << '\t' << global_stress[1][1] << '\t' << global_stress[2][2] << '\t'
		<< global_stress[1][2] << '\t' << global_stress[0][2] << '\t' << global_stress[0][1] << '\t' << F_r
		<< '\t' << F_e << '\t' << F_s;
		if (this->userInputs.enableIndentationBCs) {
			double Ind_displacement;
			double Ind_load;
			Ind_displacement = this->currentIndentDisp;
			Ind_load = this->indenterLoad;
			std::cout<<"IndenterLoad = "<<Ind_load<<'\n';
			outputFile << '\t' << Ind_load << '\t' << Ind_displacement ;
		}
		if (this->userInputs.flagUserDefinedAverageOutput){
			for(unsigned int i=0;i<this->userInputs.numberUserDefinedAverageOutput;i++){
				outputFile <<'\t'<<userDefinedAverageOutput(i);
			}
		}
		outputFile <<'\n';
		outputFile.close();
	}

	//call base class project() function to project post processed fields
	ellipticBVP<dim>::projection();
}


//------------------------------------------------------------------------
template <int dim>
void crystalPlasticity<dim>::rod2quat(Vector<double> &quat,Vector<double> &rod)
{
	double dotrod = rod(0)*rod(0) + rod(1)*rod(1) + rod(2)*rod(2);
	double cphiby2   = cos(atan(sqrt(dotrod)));
	quat(0) = cphiby2;
	quat(1) = cphiby2* rod(0);
	quat(2) = cphiby2* rod(1);
	quat(3) = cphiby2* rod(2);
}
//------------------------------------------------------------------------
template <int dim>
void crystalPlasticity<dim>::quatproduct(Vector<double> &quatp,Vector<double> &quat2,Vector<double> &quat1)
//R(qp) = R(q2)R(q1)
{
	double a = quat2(0);
	double b = quat1(0);
	double dot1 = quat1(1)*quat2(1) + quat1(2)*quat2(2) + quat1(3)*quat2(3);
	quatp(0) = (a*b) - dot1;
	quatp(1) = a*quat1(1) + b*quat2(1)+ quat2(2)*quat1(3) - quat1(2)*quat2(3);
	quatp(2) = a*quat1(2) + b*quat2(2)- quat2(1)*quat1(3) + quat1(1)*quat2(3);
	quatp(3) = a*quat1(3) + b*quat2(3)+ quat2(1)*quat1(2) - quat1(1)*quat2(2);
	if (quatp(0) < 0) {
		quatp(0) = -quatp(0);
		quatp(1) = -quatp(1);
		quatp(2) = -quatp(2);
		quatp(3) = -quatp(3);
	}
}
//------------------------------------------------------------------------
template <int dim>
void crystalPlasticity<dim>::quat2rod(Vector<double> &quat,Vector<double> &rod)
{
	if (fabs(quat(0))<1e-4){
	   quat(0)=1e-4;
	}
	double invquat1 = 1/quat(0);

	for (int i = 0;i <= 2;i++)
	rod(i) = quat(i+1)*invquat1;

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
