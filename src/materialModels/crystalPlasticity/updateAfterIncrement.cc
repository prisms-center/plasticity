#include "../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

template <int dim>
void crystalPlasticity<dim>::updateAfterIncrement()
{
	local_F_r=0.0;
	local_F_s=0.0;
	local_F_e = 0.0;
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

			Vector<double> Ulocal(dofs_per_cell);

			for (unsigned int i = 0; i < dofs_per_cell; i++) {
				Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
			}
			for (unsigned int q = 0; q < num_quad_points; ++q) {
				//Get deformation gradient
				F = 0.0;
				for (unsigned int d = 0; d < dofs_per_cell; ++d) {
					unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
					for (unsigned int j = 0; j < dim; ++j) {
						F[i][j] += Ulocal(d)*fe_values.shape_grad(d, q)[j]; // u_{i,j}= U(d)*N(d)_{,j}, where d is the DOF correonding to the i'th dimension
					}
				}
				for (unsigned int i = 0; i < dim; ++i) {
					F[i][i] += 1;
				}
				//Update strain, stress, and tangent for current time step/quadrature point
				calculatePlasticity(cellID, q, 0);

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


				local_strain.add(1.0, temp4);
				local_stress.add(1.0, temp3);
				local_microvol = local_microvol + fe_values.JxW(q);

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
				twin_ouput[cellID][q]=twin_iter[cellID][q];
				this->postprocessValues(cellID, q, 0, 0) = vonmises;
				this->postprocessValues(cellID, q, 1, 0) = eqvstrain;
				this->postprocessValues(cellID, q, 2, 0) = twin_ouput[cellID][q];

				local_F_e = local_F_e + twin_iter[cellID][q] * fe_values.JxW(q);

				for(unsigned int i=0;i<this->userInputs.numTwinSystems1;i++){
					local_F_r=local_F_r+twinfraction_iter[cellID][q][i]*fe_values.JxW(q);
				}

				for(unsigned int i=0;i<this->userInputs.numSlipSystems1;i++){
					local_F_s=local_F_s+slipfraction_iter[cellID][q][i]*fe_values.JxW(q);
				}

			}
			this->postprocessValuesAtCellCenters(cellID,0)=cellOrientationMap[cellID];

			cellID++;
		}
	}

	//In Case we have twinning
	rotnew_conv=rotnew_iter;
	
	//reorient() updates the rotnew_conv.
	reorient();
	
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

	



	char buffer[200];

	if (this->userInputs.writeQuadratureOutput) {
		if (this->currentIncrement%this->userInputs.skipQuadratureOutputSteps == 0) {
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
						temp.push_back(phase[cellID][q]);
						temp.push_back(fe_values.JxW(q));

						temp.push_back(twin_conv[cellID][q]);

						temp.push_back(fe_values.get_quadrature_points()[q][0]);
						temp.push_back(fe_values.get_quadrature_points()[q][1]);
						temp.push_back(fe_values.get_quadrature_points()[q][2]);

						temp.push_back(rotnew_conv[cellID][q][0]);
						temp.push_back(rotnew_conv[cellID][q][1]);
						temp.push_back(rotnew_conv[cellID][q][2]);

						temp.push_back(Fe_conv[cellID][q][0][0]);
						temp.push_back(Fe_conv[cellID][q][1][1]);
						temp.push_back(Fe_conv[cellID][q][2][2]);
						temp.push_back(Fe_conv[cellID][q][0][1]);
						temp.push_back(Fe_conv[cellID][q][0][2]);
						temp.push_back(Fe_conv[cellID][q][1][0]);
						temp.push_back(Fe_conv[cellID][q][1][2]);
						temp.push_back(Fe_conv[cellID][q][2][0]);
						temp.push_back(Fe_conv[cellID][q][2][1]);

						temp.push_back(Fp_conv[cellID][q][0][0]);
						temp.push_back(Fp_conv[cellID][q][1][1]);
						temp.push_back(Fp_conv[cellID][q][2][2]);
						temp.push_back(Fp_conv[cellID][q][0][1]);
						temp.push_back(Fp_conv[cellID][q][0][2]);
						temp.push_back(Fp_conv[cellID][q][1][0]);
						temp.push_back(Fp_conv[cellID][q][1][2]);
						temp.push_back(Fp_conv[cellID][q][2][0]);
						temp.push_back(Fp_conv[cellID][q][2][1]);

						temp.push_back(CauchyStress[cellID][q][0][0]);
						temp.push_back(CauchyStress[cellID][q][1][1]);
						temp.push_back(CauchyStress[cellID][q][2][2]);
						temp.push_back(CauchyStress[cellID][q][0][1]);
						temp.push_back(CauchyStress[cellID][q][0][2]);
						temp.push_back(CauchyStress[cellID][q][1][0]);
						temp.push_back(CauchyStress[cellID][q][1][2]);
						temp.push_back(CauchyStress[cellID][q][2][0]);
						temp.push_back(CauchyStress[cellID][q][2][1]);

						temp.push_back(slipfraction_conv[cellID][q][0]);
						temp.push_back(slipfraction_conv[cellID][q][1]);
						temp.push_back(slipfraction_conv[cellID][q][2]);
						temp.push_back(slipfraction_conv[cellID][q][3]);
						temp.push_back(slipfraction_conv[cellID][q][4]);
						temp.push_back(slipfraction_conv[cellID][q][5]);
						temp.push_back(slipfraction_conv[cellID][q][6]);
						temp.push_back(slipfraction_conv[cellID][q][7]);
						temp.push_back(slipfraction_conv[cellID][q][8]);
						temp.push_back(slipfraction_conv[cellID][q][9]);
						temp.push_back(slipfraction_conv[cellID][q][10]);
						temp.push_back(slipfraction_conv[cellID][q][11]);
						temp.push_back(slipfraction_conv[cellID][q][12]);
						temp.push_back(slipfraction_conv[cellID][q][13]);
						temp.push_back(slipfraction_conv[cellID][q][14]);
						temp.push_back(slipfraction_conv[cellID][q][15]);
						temp.push_back(slipfraction_conv[cellID][q][16]);
						temp.push_back(slipfraction_conv[cellID][q][17]);
						temp.push_back(slipfraction_conv[cellID][q][18]);
						temp.push_back(slipfraction_conv[cellID][q][19]);
						temp.push_back(slipfraction_conv[cellID][q][20]);
						temp.push_back(slipfraction_conv[cellID][q][21]);
						temp.push_back(slipfraction_conv[cellID][q][22]);
						temp.push_back(slipfraction_conv[cellID][q][23]);
						temp.push_back(slipfraction_conv[cellID][q][24]);
						temp.push_back(slipfraction_conv[cellID][q][25]);
						temp.push_back(slipfraction_conv[cellID][q][26]);
						temp.push_back(slipfraction_conv[cellID][q][27]);
						temp.push_back(slipfraction_conv[cellID][q][28]);
						temp.push_back(slipfraction_conv[cellID][q][29]);
						temp.push_back(slipfraction_conv[cellID][q][30]);
						temp.push_back(slipfraction_conv[cellID][q][31]);
						temp.push_back(slipfraction_conv[cellID][q][32]);
						temp.push_back(slipfraction_conv[cellID][q][33]);
						temp.push_back(slipfraction_conv[cellID][q][34]);
						temp.push_back(slipfraction_conv[cellID][q][35]);
						temp.push_back(slipfraction_conv[cellID][q][36]);
						temp.push_back(slipfraction_conv[cellID][q][37]);
						temp.push_back(slipfraction_conv[cellID][q][38]);
						temp.push_back(slipfraction_conv[cellID][q][39]);
						temp.push_back(slipfraction_conv[cellID][q][40]);
						temp.push_back(slipfraction_conv[cellID][q][41]);
						temp.push_back(slipfraction_conv[cellID][q][42]);
						temp.push_back(slipfraction_conv[cellID][q][43]);
						temp.push_back(slipfraction_conv[cellID][q][44]);
						temp.push_back(slipfraction_conv[cellID][q][45]);
						temp.push_back(slipfraction_conv[cellID][q][46]);
						temp.push_back(slipfraction_conv[cellID][q][47]);

						temp.push_back(twinfraction_conv[cellID][q][0]);
						temp.push_back(twinfraction_conv[cellID][q][1]);
						temp.push_back(twinfraction_conv[cellID][q][2]);
						temp.push_back(twinfraction_conv[cellID][q][3]);
						temp.push_back(twinfraction_conv[cellID][q][4]);
						temp.push_back(twinfraction_conv[cellID][q][5]);
						temp.push_back(twinfraction_conv[cellID][q][6]);
						temp.push_back(twinfraction_conv[cellID][q][7]);
						temp.push_back(twinfraction_conv[cellID][q][8]);
						temp.push_back(twinfraction_conv[cellID][q][9]);
						temp.push_back(twinfraction_conv[cellID][q][10]);
						temp.push_back(twinfraction_conv[cellID][q][11]);

						if (this->userInputs.enableUserMaterialModel){
							temp.push_back(stateVar_conv[cellID][q][0]);
							temp.push_back(stateVar_conv[cellID][q][1]);
							temp.push_back(stateVar_conv[cellID][q][2]);
							temp.push_back(stateVar_conv[cellID][q][3]);
							temp.push_back(stateVar_conv[cellID][q][4]);
							temp.push_back(stateVar_conv[cellID][q][5]);
							temp.push_back(stateVar_conv[cellID][q][6]);
							temp.push_back(stateVar_conv[cellID][q][7]);
							temp.push_back(stateVar_conv[cellID][q][8]);
							temp.push_back(stateVar_conv[cellID][q][9]);
							temp.push_back(stateVar_conv[cellID][q][10]);
							temp.push_back(stateVar_conv[cellID][q][11]);
							temp.push_back(stateVar_conv[cellID][q][12]);
							temp.push_back(stateVar_conv[cellID][q][13]);
							temp.push_back(stateVar_conv[cellID][q][14]);
							temp.push_back(stateVar_conv[cellID][q][15]);
							temp.push_back(stateVar_conv[cellID][q][16]);
							temp.push_back(stateVar_conv[cellID][q][17]);
							temp.push_back(stateVar_conv[cellID][q][18]);
							temp.push_back(stateVar_conv[cellID][q][19]);
							temp.push_back(stateVar_conv[cellID][q][20]);
							temp.push_back(stateVar_conv[cellID][q][21]);
							temp.push_back(stateVar_conv[cellID][q][22]);
							temp.push_back(stateVar_conv[cellID][q][23]);
							temp.push_back(stateVar_conv[cellID][q][24]);
							temp.push_back(stateVar_conv[cellID][q][25]);
							temp.push_back(stateVar_conv[cellID][q][26]);
							temp.push_back(stateVar_conv[cellID][q][27]);
							temp.push_back(stateVar_conv[cellID][q][28]);
							temp.push_back(stateVar_conv[cellID][q][29]);
							temp.push_back(stateVar_conv[cellID][q][30]);
							temp.push_back(stateVar_conv[cellID][q][31]);
							temp.push_back(stateVar_conv[cellID][q][32]);
							temp.push_back(stateVar_conv[cellID][q][33]);
							temp.push_back(stateVar_conv[cellID][q][34]);
							temp.push_back(stateVar_conv[cellID][q][35]);
							temp.push_back(stateVar_conv[cellID][q][36]);
							temp.push_back(stateVar_conv[cellID][q][37]);
							temp.push_back(stateVar_conv[cellID][q][38]);
							temp.push_back(stateVar_conv[cellID][q][39]);
							temp.push_back(stateVar_conv[cellID][q][40]);
							temp.push_back(stateVar_conv[cellID][q][41]);
							temp.push_back(stateVar_conv[cellID][q][42]);
							temp.push_back(stateVar_conv[cellID][q][43]);
							temp.push_back(stateVar_conv[cellID][q][44]);
							temp.push_back(stateVar_conv[cellID][q][45]);
							temp.push_back(stateVar_conv[cellID][q][46]);
							temp.push_back(stateVar_conv[cellID][q][47]);
							temp.push_back(stateVar_conv[cellID][q][48]);
							temp.push_back(stateVar_conv[cellID][q][49]);
							temp.push_back(stateVar_conv[cellID][q][50]);
						}

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
			global_strain[i][j]=Utilities::MPI::sum(local_strain[i][j]/microvol,this->mpi_communicator);
			global_stress[i][j]=Utilities::MPI::sum(local_stress[i][j]/microvol,this->mpi_communicator);
		}
	}
	if (!this->userInputs.enableMultiphase){
		F_e = Utilities::MPI::sum(local_F_e / microvol, this->mpi_communicator);
		F_r=Utilities::MPI::sum(local_F_r/microvol,this->mpi_communicator);
		F_s=Utilities::MPI::sum(local_F_s/microvol,this->mpi_communicator);
	}
	else {
		F_e=0;F_r=0;F_s=0;
	}

	//check whether to write stress and strain data to file
	//write stress and strain data to file
	std::string dir(this->userInputs.outputDirectory);
	dir+="/";

	std::ofstream outputFile;
	if(this->currentIncrement==0){
		dir += std::string("stressstrain.txt");
		outputFile.open(dir.c_str());
		outputFile << "Exx"<<'\t'<<"Eyy"<<'\t'<<"Ezz"<<'\t'<<"Eyz"<<'\t'<<"Exz"<<'\t'<<"Exy"<<'\t'<<"Txx"<<'\t'<<"Tyy"<<'\t'<<"Tzz"<<'\t'<<"Tyz"<<'\t'<<"Txz"<<'\t'<<"Txy"<<'\t'<<"TwinRealVF"<<'\t'<<"TwinMade"<<'\t'<<"SlipTotal"<<'\n';
		outputFile.close();
	}
	else{
		dir += std::string("stressstrain.txt");
	}
	outputFile.open(dir.c_str(),std::fstream::app);
	if(Utilities::MPI::this_mpi_process(this->mpi_communicator)==0){
		outputFile << global_strain[0][0]<<'\t'<<global_strain[1][1]<<'\t'<<global_strain[2][2]<<'\t'<<global_strain[1][2]<<'\t'<<global_strain[0][2]<<'\t'<<global_strain[0][1]<<'\t'<<global_stress[0][0]<<'\t'<<global_stress[1][1]<<'\t'<<global_stress[2][2]<<'\t'<<global_stress[1][2]<<'\t'<<global_stress[0][2]<<'\t'<<global_stress[0][1]<<'\t'<<F_r<<'\t'<<F_e<<'\t'<<F_s<<'\n';
	}
	outputFile.close();


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
	double invquat1 = 1/quat(0);

	for (int i = 0;i <= 2;i++)
	rod(i) = quat(i+1)*invquat1;

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
