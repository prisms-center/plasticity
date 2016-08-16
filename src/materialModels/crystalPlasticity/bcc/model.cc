 //constructor
template <int dim>
crystalPlasticity<dim>::crystalPlasticity() :
ellipticBVP<dim>(),
F(dim,dim),
F_tau(dim,dim),
FP_tau(dim,dim),
FE_tau(dim,dim),
T(dim,dim),
P(dim,dim)
{
    initCalled = false;
    
    //post processing
    ellipticBVP<dim>::numPostProcessedFields=3;
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Grain_ID");
}
        
//implementation of the getElementalValues method
template <int dim>
void crystalPlasticity<dim>::getElementalValues(FEValues<dim>& fe_values,
						unsigned int dofs_per_cell,
						unsigned int num_quad_points,
						FullMatrix<double>& elementalJacobian,
						Vector<double>&     elementalResidual)
 {

     //Initialized history variables and pfunction variables if unititialized
     if(initCalled == false){
	 init(num_quad_points);
     }

     unsigned int cellID = fe_values.get_cell()->user_index();
     std::vector<unsigned int> local_dof_indices(dofs_per_cell);
     Vector<double> Ulocal(dofs_per_cell);

     typename DoFHandler<dim>::active_cell_iterator cell(& this->triangulation,
							 fe_values.get_cell()->level(),
							 fe_values.get_cell()->index(),
							 & this->dofHandler);
     cell->set_user_index(fe_values.get_cell()->user_index());
     cell->get_dof_indices (local_dof_indices);
     for(unsigned int i=0; i<dofs_per_cell; i++){
	 Ulocal[i] = this->solutionWithGhosts[local_dof_indices[i]];
     }

     //local data structures
     FullMatrix<double> K_local(dofs_per_cell,dofs_per_cell),CE_tau(dim,dim),E_tau(dim,dim),temp,temp2,temp3;
     Vector<double> Rlocal (dofs_per_cell);
     K_local = 0.0; Rlocal = 0.0;


     //loop over quadrature points
     for (unsigned int q=0; q<num_quad_points; ++q){
	 //Get deformation gradient
	 F=0.0;
	 for (unsigned int d=0; d<dofs_per_cell; ++d){
	     unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
	     for (unsigned int j=0; j<dim; ++j){
		 F[i][j]+=Ulocal(d)*fe_values.shape_grad(d, q)[j]; // u_{i,j}= U(d)*N(d)_{,j}, where d is the DOF correonding to the i'th dimension
	     }
	 }
	 for (unsigned int i=0; i<dim; ++i){
	     F[i][i]+=1;
	 }

     //this->pcout<<F[0][0]<<"\t"<<F[1][1]<<"\t"<<F[2][2]<<"\n";

	 //Update strain, stress, and tangent for current time step/quadrature point
	 calculatePlasticity(cellID, q);

     //this->pcout<<P[0][0]<<"\t"<<P[1][1]<<"\t"<<P[2][2]<<"\n";
         
	 //Fill local residual
	 for (unsigned int d=0; d<dofs_per_cell; ++d) {
	     unsigned int i = fe_values.get_fe().system_to_component_index(d).first;
	     for (unsigned int j = 0; j < dim; j++){
		 Rlocal(d) -=  fe_values.shape_grad(d, q)[j]*P[i][j]*fe_values.JxW(q);
	     }

	 }

	 temp.reinit(dim,dim); temp=0.0;
	 temp2.reinit(dim,dim); temp2=0.0;
	 temp3.reinit(dim,dim); temp3=0.0;
	 CE_tau=0.0;
	 temp=F;
	 F.Tmmult(CE_tau,temp);
	 E_tau=CE_tau;
	 temp=IdentityMatrix(dim);
	 for(unsigned int i=0;i<dim;i++){
	     for(unsigned int j=0;j<dim;j++){
		 E_tau[i][j] = 0.5*(E_tau[i][j]-temp[i][j]);
		 temp2[i][j]=E_tau[i][j]*fe_values.JxW(q);
		 temp3[i][j]=T[i][j]*fe_values.JxW(q);
	     }
	 }

	 local_strain.add(1.0,temp2);
	 local_stress.add(1.0,temp3);
	 local_microvol=local_microvol+fe_values.JxW(q);

         //calculate von-Mises stress and equivalent strain
         double traceE, traceT,vonmises,eqvstrain;
         FullMatrix<double> deve(dim,dim),devt(dim,dim);
         
         CE_tau=0.0;
         temp=F;
         F.Tmmult(CE_tau,temp);
         E_tau=CE_tau;
         temp=IdentityMatrix(dim);
         for(unsigned int i=0;i<dim;i++){
             for(unsigned int j=0;j<dim;j++){
                 E_tau[i][j] = 0.5*(E_tau[i][j]-temp[i][j]);
             }
         }
         
         
         traceE=E_tau.trace();
         traceT=T.trace();
         temp=IdentityMatrix(3);
         temp.equ(traceE/3,temp);
         
         deve=E_tau;
         deve.add(-1.0,temp);
         
         temp=IdentityMatrix(3);
         temp.equ(traceT/3,temp);
         
         devt=T;
         devt.add(-1.0,temp);
         
         vonmises= devt.frobenius_norm();
         vonmises=sqrt(3.0/2.0)*vonmises;
         eqvstrain=deve.frobenius_norm();
         eqvstrain=sqrt(2.0/3.0)*eqvstrain;
         
         
         //fill in post processing field values
         
         this->postprocessValues(cellID, q, 0, 0)=eqvstrain;
         this->postprocessValues(cellID, q, 1, 0)=vonmises;
         this->postprocessValues(cellID, q, 2, 0)=quadratureOrientationsMap[cellID][q];




	 //evaluate elemental stiffness matrix, K_{ij} = N_{i,k}*C_{mknl}*F_{im}*F{jn}*N_{j,l} + N_{i,k}*F_{kl}*N_{j,l}*del{ij} dV
	 for (unsigned int d1=0; d1<dofs_per_cell; ++d1) {
	     unsigned int i = fe_values.get_fe().system_to_component_index(d1).first;
	     for (unsigned int d2=0; d2<dofs_per_cell; ++d2) {
		 unsigned int j = fe_values.get_fe().system_to_component_index(d2).first;
		 for (unsigned int k = 0; k < dim; k++){
		     for (unsigned int l= 0; l< dim; l++){
			 K_local(d1,d2) +=  fe_values.shape_grad(d1, q)[k]*dP_dF[i][k][j][l]*fe_values.shape_grad(d2, q)[l]*fe_values.JxW(q);
		     }
		 }
	     }
	 }
     }
     elementalJacobian = K_local;
     elementalResidual = Rlocal;
     
     //this->pcout<<K_local[0][0]<<"\t"<<K_local[0][1]<<"\t"<<K_local[0][2]<<"\t"<<K_local[1][0]<<"\t"<<K_local[1][1]<<"\t"<<K_local[1][2]<<"\t"<<K_local[2][0]<<"\t"<<K_local[2][1]<<"\t"<<K_local[2][2]<<"\n";
    // this->pcout<<Rlocal[0]<<"\t"<<Rlocal[1]<<"\t"<<Rlocal[2]<<"\n";
 }


 template <int dim>
 void crystalPlasticity<dim>::updateBeforeIteration()
 {
     local_strain=0.0;
     local_stress=0.0;
     local_microvol=0.0;

     //call base class project() function to project post processed fields
     //ellipticBVP<dim>::project();
 }

 template <int dim>
 void crystalPlasticity<dim>::updateBeforeIncrement()
 {
     microvol=0.0;
     //call base class project() function to project post processed fields
     //ellipticBVP<dim>::project();
 }


 //implementation of the getElementalValues method
 template <int dim>
 void crystalPlasticity<dim>::updateAfterIncrement()
 {
     reorient();

     //copy rotnew to output
     orientations.outputOrientations.clear();
     QGauss<dim>  quadrature(quadOrder);
     const unsigned int num_quad_points = quadrature.size();
     FEValues<dim> fe_values (this->FE, quadrature, update_quadrature_points | update_JxW_values);
     //loop over elements
     unsigned int cellID=0;
     typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
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
		 temp.push_back(quadratureOrientationsMap[cellID][q]);
		 orientations.addToOutputOrientations(temp);

	     }
	     cellID++;
	 }
     }
     orientations.writeOutputOrientations();

     //Update the history variables when convergence is reached for the current increment
     Fe_conv=Fe_iter;
     Fp_conv=Fp_iter;
     s_alpha_conv=s_alpha_iter;

     microvol=Utilities::MPI::sum(local_microvol,this->mpi_communicator);

     for(unsigned int i=0;i<dim;i++){
	 for(unsigned int j=0;j<dim;j++){
	     global_strain[i][j]=Utilities::MPI::sum(local_strain[i][j]/microvol,this->mpi_communicator);
	     global_stress[i][j]=Utilities::MPI::sum(local_stress[i][j]/microvol,this->mpi_communicator);
	 }

     }

     //check whether to write stress and strain data to file
#ifdef writeOutput
  if (!writeOutput) return;
#endif
     //write stress and strain data to file
#ifdef outputDirectory
     std::string dir(outputDirectory);
     dir+="/";
#else
     std::string dir("./");
#endif
     ofstream outputFile;
     if(this->currentIncrement==0){
       dir += std::string("stressstrain.txt");
       outputFile.open(dir.c_str());
       outputFile << "Exx"<<'\t'<<"Eyy"<<'\t'<<"Ezz"<<'\t'<<"Eyz"<<'\t'<<"Exz"<<'\t'<<"Exy"<<'\t'<<"Txx"<<'\t'<<"Tyy"<<'\t'<<"Tzz"<<'\t'<<"Tyz"<<'\t'<<"Txz"<<'\t'<<"Txy"<<'\n';
	 outputFile.close();
     }
     else{
     dir += std::string("stressstrain.txt");
     }
     outputFile.open(dir.c_str(),ios::app);
     if(Utilities::MPI::this_mpi_process(this->mpi_communicator)==0){
       outputFile << global_strain[0][0]<<'\t'<<global_strain[1][1]<<'\t'<<global_strain[2][2]<<'\t'<<global_strain[1][2]<<'\t'<<global_strain[0][2]<<'\t'<<global_strain[0][1]<<'\t'<<global_stress[0][0]<<'\t'<<global_stress[1][1]<<'\t'<<global_stress[2][2]<<'\t'<<global_stress[1][2]<<'\t'<<global_stress[0][2]<<'\t'<<global_stress[0][1]<<'\n';
     }
     outputFile.close();
     global_strain=0.0;
     global_stress=0.0;

     //call base class project() function to project post processed fields
     ellipticBVP<dim>::project();
 }
