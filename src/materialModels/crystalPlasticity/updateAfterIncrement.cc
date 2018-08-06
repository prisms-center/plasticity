#include "../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

template <int dim>
void crystalPlasticity<dim>::updateAfterIncrement()
{
	QGauss<dim>  quadrature(this->userInputs.quadOrder);
	FEValues<dim> fe_values(this->FE, quadrature, update_quadrature_points | update_gradients | update_JxW_values);
	const unsigned int num_quad_points = quadrature.size();
	const unsigned int   dofs_per_cell = this->FE.dofs_per_cell;
	std::vector<unsigned int> local_dof_indices(dofs_per_cell);
	//loop over elements
	unsigned int cellID = 0;
	FullMatrix<double> rotmat(dim, dim);
	Vector<double> quat1(4), rod(3), quat2(4), quatprod(4);


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
				calculatePlasticity(cellID, q);
			}

			cellID++;
		}
	}

	reorient();

  twinfraction_conv=twinfraction_iter;
	twinfraction_conv_Twin = twinfraction_iter_Twin;
	slipfraction_conv=slipfraction_iter;
	local_F_r_Twin = 0.0; //adde by Reza as a modification
	local_F_r=0.0;
  local_F_s=0.0;

  char buffer[200];
  //copy rotnew to output
  orientations.outputOrientations.clear();
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
          temp.push_back(rotnew[cellID][q][0]);
          temp.push_back(rotnew[cellID][q][1]);
          temp.push_back(rotnew[cellID][q][2]);
          temp.push_back(fe_values.JxW(q));
					temp.push_back(cellOrientationMap[cellID]);

	        twin[cellID][q] = 0.0;

          orientations.addToOutputOrientations(temp);
          for(unsigned int i=0;i<this->userInputs.numTwinSystems;i++){
              local_F_r=local_F_r+twinfraction_conv[cellID][q][i]*fe_values.JxW(q);
							local_F_r_Twin = local_F_r_Twin + twinfraction_conv_Twin[cellID][q][i] * fe_values.JxW(q);
}

					for(unsigned int i=0;i<this->userInputs.numSlipSystems;i++){
              local_F_s=local_F_s+slipfraction_conv[cellID][q][i]*fe_values.JxW(q);
          }
      }
      cellID++;
    }
  }
  orientations.writeOutputOrientations(this->userInputs.writeOutput,this->userInputs.outputDirectory);

  //Update the history variables when convergence is reached for the current increment
  Fe_conv=Fe_iter;
  Fp_conv=Fp_iter;
  s_alpha_conv=s_alpha_iter;
	if (this->currentIncrement==0) {
		F_T = this->userInputs.twinThresholdFraction;
		local_F_e = 0.0;
	}

  microvol=Utilities::MPI::sum(local_microvol,this->mpi_communicator);

  for(unsigned int i=0;i<dim;i++){
    for(unsigned int j=0;j<dim;j++){
        global_strain[i][j]=Utilities::MPI::sum(local_strain[i][j]/microvol,this->mpi_communicator);
        global_stress[i][j]=Utilities::MPI::sum(local_stress[i][j]/microvol,this->mpi_communicator);
    }
  }

	F_r_Twin = Utilities::MPI::sum(local_F_r_Twin / microvol, this->mpi_communicator);
	F_r=Utilities::MPI::sum(local_F_r/microvol,this->mpi_communicator);
  F_s=Utilities::MPI::sum(local_F_s/microvol,this->mpi_communicator);

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

	if(this->userInputs.enableCyclicLoading){
  //Adding backstress term during loading reversal
  if (this->currentIncrement == 0) {
 	 signstress = 1;
 	 backstressflag = 0;
  }

  if (backstressflag>10) {
		if (signstress*global_stress[this->userInputs.cyclicLoadingFace/2-1][this->userInputs.cyclicLoadingDOF-1]<0) {
	 	  signstress = global_stress[this->userInputs.cyclicLoadingFace/2-1][this->userInputs.cyclicLoadingDOF-1];
	 	  if (signstress<0) {
	 		  unsigned int cellID = 0;
	 		  typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
	 		  for (; cell != endc; ++cell) {
	 			  if (cell->is_locally_owned()) {
	 				  fe_values.reinit(cell);
	 	 				  for (unsigned int q = 0; q<num_quad_points; ++q)
		 					  for (unsigned int i = 0;i<n_slip_systems;i++) {
			 	 					if(this->userInputs.backstressFactor>1.0e-5)
		 							s_alpha_conv[cellID][q][i] = s_alpha_conv[cellID][q][i] - this->userInputs.backstressFactor*s_alpha_conv[cellID][q][i];
				 				}

		 			backstressflag = 0;
	 				cellID++;
	 			 	}
	 	 		}
	 	 	}
	  }
  }
  backstressflag = backstressflag + 1;
 }
  cellID=0;
  cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
  for (; cell!=endc; ++cell) {
	  if (cell->is_locally_owned()){

      fe_values.reinit(cell);

			//loop over quadrature points
      for (unsigned int q = 0; q < num_quad_points; ++q){
        std::vector<double> local_twin;
        local_twin.resize(this->userInputs.numTwinSystems,0.0);
        local_twin=twinfraction_conv_Twin[cellID][q];
        std::vector<double>::iterator result;
        result = std::max_element(local_twin.begin(), local_twin.end());
        double twin_pos, twin_max;
        twin_pos= std::distance(local_twin.begin(), result);
        twin_max=local_twin[twin_pos];

        if(F_r_Twin>0){
          if(twin_max > F_T){
	          FullMatrix<double> FE_t(dim,dim), FP_t(dim,dim),Twin_T(dim,dim),temp(dim,dim);
	          FE_t=Fe_conv[cellID][q];
	          FP_t=Fp_conv[cellID][q];
						for(unsigned int i=0;i<this->userInputs.numTwinSystems;i++){
                   twinfraction_conv_Twin[cellID][q][i]=0;
//                   s_alpha_conv[cellID][q][numSlipSystems+twin_pos]=s_alpha_twin;
							}

						rod(0) = rot[cellID][q][0];rod(1) = rot[cellID][q][1];rod(2) = rot[cellID][q][2];

	    			odfpoint(rotmat, rod);


						temp = 0.0;
						//FE_t.mmult(temp, rotmat);
						//rotmat.Tmmult(FE_t, temp);

						//FP_t.mmult(temp, rotmat);
						//rotmat.Tmmult(FP_t, temp);


						rod2quat(quat2, rod);
						quat1(0) = 0;
			      quat1(1) = n_alpha[this->userInputs.numSlipSystems + twin_pos][0];
						quat1(2) = n_alpha[this->userInputs.numSlipSystems + twin_pos][1];
			      quat1(3) = n_alpha[this->userInputs.numSlipSystems + twin_pos][2];

						quatproduct(quatprod, quat2, quat1);


						quat2rod(quatprod, rod);

	      		odfpoint(rotmat, rod);

						rot[cellID][q][0] = rod(0);rot[cellID][q][1] = rod(1);rot[cellID][q][2] = rod(2);
						rotnew[cellID][q][0] = rod(0);rotnew[cellID][q][1] = rod(1);rotnew[cellID][q][2] = rod(2);


						//FE_t.mTmult(temp, rotmat);
						//rotmat.mmult(FE_t, temp);

						//FP_t.mTmult(temp, rotmat);
						//rotmat.mmult(FP_t, temp);

	          Fe_conv[cellID][q]=FE_t;
	          Fp_conv[cellID][q]=FP_t;


						//loop over elements
						twin[cellID][q] = 1.0;


	        }
        }
      }
      cellID++;
    }
  }

	cellID = 0;
	cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
	for (; cell != endc; ++cell) {
		if (cell->is_locally_owned()) {
			fe_values.reinit(cell);
			//loop over quadrature points
			for (unsigned int q = 0; q<num_quad_points; ++q) {
				local_F_e = local_F_e + twin[cellID][q] * fe_values.JxW(q);
			}
			cellID++;
		}
	}

	F_e = Utilities::MPI::sum(local_F_e / microvol, this->mpi_communicator);
	if (F_r > 0) {
		F_T = this->userInputs.twinThresholdFraction + (this->userInputs.twinSaturationFactor*F_e / F_r_Twin);
	}
	else {
		F_T = this->userInputs.twinThresholdFraction;
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
    //if (quatp(0) < 0)
    //quatp.mult(-1);
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
