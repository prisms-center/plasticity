#include "../../../include/continuumPlasticity.h"
#include <iostream>
#include <fstream>

//implementation of the getElementalValues method
template <int dim>
void continuumPlasticity<dim>::updateAfterIncrement()
{
 // this->pcout << "continuum::updateAfterIncrement \n";
  this->updateAfterIncrementBase();
	unsigned int CheckBufferRegion,dimBuffer;
	double lowerBuffer,upperBuffer;
	Point<dim> pnt2;
	Vector<double> userDefinedAverageOutput,local_userDefinedAverageOutput;

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
	std::vector<unsigned int> local_dof_indices(dofs_per_cell);
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

			Vector<double> Ulocal(dofs_per_cell);


			for (unsigned int i = 0; i < dofs_per_cell; i++) {
				Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
			}

      enhStrain.staticCondensationData[cellID].reset();
      //Initialize the enhanced strain object with this information.
      enhStrain.reinit(Ulocal, cell);


			for (unsigned int q = 0; q < num_quad_points; ++q) {
        //Get enhanced deformation gradient
        enhStrain.get_F_enh(q, F);


				//Update strain, stress, and tangent for current time step/quadrature point
				calculatePlasticity(cellID, q);

        T=0;
        T.add(1/F.determinant(), tau);
				FullMatrix<double> temp,temp3,temp4, C_tau(dim, dim), E_tau(dim, dim), b_tau(dim, dim);
				Vector<double> temp2;
				temp.reinit(dim, dim); temp = 0.0;
				temp2.reinit(dim); temp2 = 0.0;
				temp3.reinit(dim, dim); temp3 = 0.0;
				temp4.reinit(dim, dim); temp4 = 0.0;
				C_tau = 0.0;
				temp = F;
				F.Tmmult(C_tau, temp);
				F.mTmult(b_tau, temp);
				//E_tau = CE_tau;
				temp = IdentityMatrix(dim);
				for (unsigned int i = 0;i<dim;i++) {
					temp2[i] = 0.5*log(b_tau[i][i])*fe_values.JxW(q);
					for (unsigned int j = 0;j<dim;j++) {
						E_tau[i][j] = 0.5*(C_tau[i][j] - temp[i][j]);
						temp3[i][j] = T[i][j] * fe_values.JxW(q);
						temp4[i][j] = E_tau[i][j]*fe_values.JxW(q);
					}
				}

				CauchyStress[cellID][q]=T;
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
				if (this->userInputs.writeOutput){
					this->postprocessValues(cellID, q, 0, 0) = vonmises;
					this->postprocessValues(cellID, q, 1, 0) = eqvstrain;
					this->postprocessValues(cellID, q, 2, 0) = histAlpha_iter[cellID][q];

					////////User Defined Variables for visualization outputs (output_Var1 to output_Var24)////////
					this->postprocessValues(cellID, q, 3, 0) = 0;
					this->postprocessValues(cellID, q, 4, 0) = 0;
					this->postprocessValues(cellID, q, 5, 0) = 0;
					this->postprocessValues(cellID, q, 6, 0) = 0;
					this->postprocessValues(cellID, q, 7, 0) = 0;
					this->postprocessValues(cellID, q, 8, 0) = 0;
					this->postprocessValues(cellID, q, 9, 0) = 0;
					this->postprocessValues(cellID, q, 10, 0) = 0;
					this->postprocessValues(cellID, q, 11, 0) = 0;
					this->postprocessValues(cellID, q, 12, 0) = 0;
					this->postprocessValues(cellID, q, 13, 0) = 0;
					this->postprocessValues(cellID, q, 14, 0) = 0;
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

					if (this->userInputs.flagUserDefinedAverageOutput){
						////One should define their UserDefinedAverageOutput here by defining based on the state Variables
						///In the following equation, "+0" should be substituted by "+targetVariable", where targetVariable is the variable you want to plot the average value.
						//Also, the first equation should be copy and paste depending on the number of variables you want to output,
						//and the integer in paranthesis should be updated accordingly.
						//The maximum number of lines (output variables) are defined in the input file as "set Number of Output Userdefined Average Variable".
						local_userDefinedAverageOutput(0)=local_userDefinedAverageOutput(0)+0;
					}
				}

			}

			cellID++;
		}
	}

  //Update the history variables when convergence is reached for the current increment
  histInvCP_conv = histInvCP_iter;
  histAlpha_conv = histAlpha_iter;
  histXi_conv = histXi_iter;

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
						temp.push_back(fe_values.get_quadrature_points()[q][0]);
						temp.push_back(fe_values.get_quadrature_points()[q][1]);
						temp.push_back(fe_values.get_quadrature_points()[q][2]);

						temp.push_back(CauchyStress[cellID][q][0][0]);
						temp.push_back(CauchyStress[cellID][q][1][1]);
						temp.push_back(CauchyStress[cellID][q][2][2]);
						temp.push_back(CauchyStress[cellID][q][0][1]);
						temp.push_back(CauchyStress[cellID][q][0][2]);
						temp.push_back(CauchyStress[cellID][q][1][0]);
						temp.push_back(CauchyStress[cellID][q][1][2]);
						temp.push_back(CauchyStress[cellID][q][2][0]);
						temp.push_back(CauchyStress[cellID][q][2][1]);

						addToQuadratureOutput(temp);

					}
					cellID++;
				}
			}

			writeQuadratureOutput(this->userInputs.outputDirectory, this->currentIncrement);
		}
	}

	microvol=Utilities::MPI::sum(local_microvol,this->mpi_communicator);
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=0;j<dim;j++){
			global_strain[i][j]=Utilities::MPI::sum(local_strain[i][j],this->mpi_communicator)/microvol;
			global_stress[i][j]=Utilities::MPI::sum(local_stress[i][j],this->mpi_communicator)/microvol;
		}
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
			outputFile << "Exx"<<'\t'<<"Eyy"<<'\t'<<"Ezz"<<'\t'<<"Eyz"<<'\t'<<"Exz"<<'\t'<<"Exy"<<'\t'<<"Txx"<<'\t'<<"Tyy"<<'\t'<<"Tzz"<<'\t'<<"Tyz"<<'\t'<<"Txz"<<'\t'<<"Txy";
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
		outputFile << global_strain[0][0]<<'\t'<<global_strain[1][1]<<'\t'<<global_strain[2][2]<<'\t'<<global_strain[1][2]<<'\t'<<global_strain[0][2]<<'\t'<<global_strain[0][1]<<'\t'<<global_stress[0][0]<<'\t'<<global_stress[1][1]<<'\t'<<global_stress[2][2]<<'\t'<<global_stress[1][2]<<'\t'<<global_stress[0][2]<<'\t'<<global_stress[0][1];
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

#include "../../../include/continuumPlasticity_template_instantiations.h"
