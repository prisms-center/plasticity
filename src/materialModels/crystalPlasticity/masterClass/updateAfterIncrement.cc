#include "../../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

//implementation of the updateAfterIncrement method
template <int dim>
void crystalPlasticity<dim>::updateAfterIncrement()
{
    reorient();

    //copy rotnew to output
    orientations.outputOrientations.clear();
    QGauss<dim>  quadrature(this->userInputs.quadOrder);
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
    std::ofstream outputFile;
    if(this->currentIncrement==0){
      dir += std::string("stressstrain.txt");
      outputFile.open(dir.c_str());
      outputFile << "Exx"<<'\t'<<"Eyy"<<'\t'<<"Ezz"<<'\t'<<"Eyz"<<'\t'<<"Exz"<<'\t'<<"Exy"<<'\t'<<"Txx"<<'\t'<<"Tyy"<<'\t'<<"Tzz"<<'\t'<<"Tyz"<<'\t'<<"Txz"<<'\t'<<"Txy"<<'\n';
  outputFile.close();
    }
    else{
    dir += std::string("stressstrain.txt");
    }
    outputFile.open(dir.c_str(),std::fstream::app);
    if(Utilities::MPI::this_mpi_process(this->mpi_communicator)==0){
      outputFile << global_strain[0][0]<<'\t'<<global_strain[1][1]<<'\t'<<global_strain[2][2]<<'\t'<<global_strain[1][2]<<'\t'<<global_strain[0][2]<<'\t'<<global_strain[0][1]<<'\t'<<global_stress[0][0]<<'\t'<<global_stress[1][1]<<'\t'<<global_stress[2][2]<<'\t'<<global_stress[1][2]<<'\t'<<global_stress[0][2]<<'\t'<<global_stress[0][1]<<'\n';
    }
    outputFile.close();



    //Adding backstress term during loading reversal

    if(this->currentIncrement==0){
        signstress=global_stress.trace();
    }

   //this->pcout<<signstress<<'\t'<<backstressFactor<<'\n';

    if(signstress*global_stress.trace()<0){
        signstress=global_stress.trace();
        //if(signstress<0){
        //loop over elements
        unsigned int cellID=0;
        typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
        for (; cell!=endc; ++cell) {
            if (cell->is_locally_owned()){
                fe_values.reinit(cell);
                //loop over quadrature points
                for (unsigned int q=0; q<num_quad_points; ++q){
                    for(unsigned int i=0;i<(this->userInputs.numSlipSystems);i++){

                        s_alpha_conv[cellID][q][i]=s_alpha_conv[cellID][q][i]-this->userInputs.backstressFactor*s_alpha_conv[cellID][q][i];
                    }
                }
                cellID++;
            }

        }

        // }

    }

    global_strain=0.0;
    global_stress=0.0;


    //call base class projection() function to project post processed fields
    ellipticBVP<dim>::projection();
}

#include "../../../../include/crystalPlasticity_template_instantiations.h"
