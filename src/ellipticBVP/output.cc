//output method for ellipticBVP class
#include "../../include/ellipticBVP.h"
#include <fstream>

//output results
template <int dim>
void ellipticBVP<dim>::output(){
  DataOut<dim> data_out, data_out_Scalar;
  data_out.attach_dof_handler (dofHandler);
  data_out_Scalar.attach_dof_handler (dofHandler_Scalar);

  //add displacement field
  data_out.add_data_vector (solutionWithGhosts,
			    nodal_solution_names,
			    DataOut<dim>::type_dof_data,
			    nodal_data_component_interpretation);

  //add postprocessing fields
  unsigned int numPostProcessedFieldsWritten=0;
  for (unsigned int field=0; field<numPostProcessedFields; field++){
    if(!userInputs.output_Eqv_strain)
        if (postprocessed_solution_names[field].compare(std::string("Eqv_strain"))==0) continue;
    if(!userInputs.output_Eqv_stress)
        if (postprocessed_solution_names[field].compare(std::string("Eqv_stress"))==0) continue;
    if (userInputs.continuum_Isotropic){
      if(!userInputs.output_alpha)
        if (postprocessed_solution_names[field].compare(std::string("alpha"))==0) continue;
    }
    else {
      if(!userInputs.output_Twin)
        if (postprocessed_solution_names[field].compare(std::string("Twin"))==0) continue;
    }


//pcout<<"field="<<field<<"step1\n";

  if(!userInputs.output_Var1)
    if (postprocessed_solution_names[field].compare(std::string("output_Var1"))==0) continue;
  if(!userInputs.output_Var2)
    if (postprocessed_solution_names[field].compare(std::string("output_Var2"))==0) continue;
  if(!userInputs.output_Var3)
    if (postprocessed_solution_names[field].compare(std::string("output_Var3"))==0) continue;
  if(!userInputs.output_Var4)
    if (postprocessed_solution_names[field].compare(std::string("output_Var4"))==0) continue;
  if(!userInputs.output_Var5)
    if (postprocessed_solution_names[field].compare(std::string("output_Var5"))==0) continue;
  if(!userInputs.output_Var6)
    if (postprocessed_solution_names[field].compare(std::string("output_Var6"))==0) continue;
  if(!userInputs.output_Var7)
    if (postprocessed_solution_names[field].compare(std::string("output_Var7"))==0) continue;
  if(!userInputs.output_Var8)
    if (postprocessed_solution_names[field].compare(std::string("output_Var8"))==0) continue;
  if(!userInputs.output_Var9)
    if (postprocessed_solution_names[field].compare(std::string("output_Var9"))==0) continue;
  if(!userInputs.output_Var10)
    if (postprocessed_solution_names[field].compare(std::string("output_Var10"))==0) continue;
  if(!userInputs.output_Var11)
    if (postprocessed_solution_names[field].compare(std::string("output_Var11"))==0) continue;
  if(!userInputs.output_Var12)
    if (postprocessed_solution_names[field].compare(std::string("output_Var12"))==0) continue;
  if(!userInputs.output_Var13)
    if (postprocessed_solution_names[field].compare(std::string("output_Var13"))==0) continue;
  if(!userInputs.output_Var14)
    if (postprocessed_solution_names[field].compare(std::string("output_Var14"))==0) continue;
  if(!userInputs.output_Var15)
    if (postprocessed_solution_names[field].compare(std::string("output_Var15"))==0) continue;
  if(!userInputs.output_Var16)
    if (postprocessed_solution_names[field].compare(std::string("output_Var16"))==0) continue;
  if(!userInputs.output_Var17)
    if (postprocessed_solution_names[field].compare(std::string("output_Var17"))==0) continue;
  if(!userInputs.output_Var18)
    if (postprocessed_solution_names[field].compare(std::string("output_Var18"))==0) continue;
  if(!userInputs.output_Var19)
    if (postprocessed_solution_names[field].compare(std::string("output_Var19"))==0) continue;
  if(!userInputs.output_Var20)
    if (postprocessed_solution_names[field].compare(std::string("output_Var20"))==0) continue;
  if(!userInputs.output_Var21)
    if (postprocessed_solution_names[field].compare(std::string("output_Var21"))==0) continue;
  if(!userInputs.output_Var22)
    if (postprocessed_solution_names[field].compare(std::string("output_Var22"))==0) continue;
  if(!userInputs.output_Var23)
    if (postprocessed_solution_names[field].compare(std::string("output_Var23"))==0) continue;
  if(!userInputs.output_Var24)
    if (postprocessed_solution_names[field].compare(std::string("output_Var24"))==0) continue;




  //if(!userInputs.output_alpha)
    //if (postprocessed_solution_names[field].compare(std::string("alpha"))==0) continue;
  //if(!userInputs.output_tau_vm)
    //if (postprocessed_solution_names[field].compare(std::string("tau_vm"))==0) continue;
    //
    data_out_Scalar.add_data_vector (*postFieldsWithGhosts[field],
				     postprocessed_solution_names[field].c_str());
    numPostProcessedFieldsWritten++;
  }


  //add material id to output file
  Vector<float> material (triangulation.n_active_cells());
  unsigned int matID=0;
  unsigned int cellID=0;
  typename parallel::distributed::Triangulation<dim>::active_cell_iterator cell = triangulation.begin_active(), endc = triangulation.end();
  if(userInputs.readExternalMesh){
    for (; cell!=endc; ++cell){
      material(matID) = cell->material_id(); matID++;}
  }
  else{
    for (; cell!=endc; ++cell){
      if(cell->is_locally_owned()){
        material(cell->active_cell_index())=postprocessValuesAtCellCenters(cellID,0);
        cellID++;
      }
    }
  }

  data_out.add_data_vector (material, "meshGrain_ID");

  Vector<float> FieldsAtCellCenters (triangulation.n_active_cells());
  for (unsigned int field=1; field<numPostProcessedFieldsAtCellCenters; field++){

    if(!userInputs.outputCellCenters_Var1)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var1"))==0) continue;
    if(!userInputs.outputCellCenters_Var2)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var2"))==0) continue;
    if(!userInputs.outputCellCenters_Var3)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var3"))==0) continue;
    if(!userInputs.outputCellCenters_Var4)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var4"))==0) continue;
    if(!userInputs.outputCellCenters_Var5)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var5"))==0) continue;
    if(!userInputs.outputCellCenters_Var6)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var6"))==0) continue;
    if(!userInputs.outputCellCenters_Var7)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var7"))==0) continue;
    if(!userInputs.outputCellCenters_Var8)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var8"))==0) continue;
    if(!userInputs.outputCellCenters_Var9)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var9"))==0) continue;
    if(!userInputs.outputCellCenters_Var10)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var10"))==0) continue;
    if(!userInputs.outputCellCenters_Var11)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var11"))==0) continue;
    if(!userInputs.outputCellCenters_Var12)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var12"))==0) continue;
    if(!userInputs.outputCellCenters_Var13)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var13"))==0) continue;
    if(!userInputs.outputCellCenters_Var14)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var14"))==0) continue;
    if(!userInputs.outputCellCenters_Var15)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var15"))==0) continue;
    if(!userInputs.outputCellCenters_Var16)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var16"))==0) continue;
    if(!userInputs.outputCellCenters_Var17)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var17"))==0) continue;
    if(!userInputs.outputCellCenters_Var18)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var18"))==0) continue;
    if(!userInputs.outputCellCenters_Var19)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var19"))==0) continue;
    if(!userInputs.outputCellCenters_Var20)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var20"))==0) continue;
    if(!userInputs.outputCellCenters_Var21)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var21"))==0) continue;
    if(!userInputs.outputCellCenters_Var22)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var22"))==0) continue;
    if(!userInputs.outputCellCenters_Var23)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var23"))==0) continue;
    if(!userInputs.outputCellCenters_Var24)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var24"))==0) continue;
    if(!userInputs.outputCellCenters_Var25)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var25"))==0) continue;
    if(!userInputs.outputCellCenters_Var26)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var26"))==0) continue;
    if(!userInputs.outputCellCenters_Var27)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var27"))==0) continue;
    if(!userInputs.outputCellCenters_Var28)
      if (postprocessedFieldsAtCellCenters_solution_names[field].compare(std::string("outputCellCenters_Var28"))==0) continue;

    cellID=0;
    cell = triangulation.begin_active(), endc = triangulation.end();
    for (; cell!=endc; ++cell){
      if(cell->is_locally_owned()){
        FieldsAtCellCenters(cell->active_cell_index())=postprocessValuesAtCellCenters(cellID,field);
        cellID++;
      }
    }

    data_out.add_data_vector (FieldsAtCellCenters,
				     postprocessedFieldsAtCellCenters_solution_names[field].c_str());
  }

  //add subdomain id to output file
  Vector<float> subdomain (triangulation.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");
  if (numPostProcessedFieldsWritten>0){
    data_out_Scalar.add_data_vector (subdomain, "subdomain");
  }

   data_out.build_patches ();

   if (numPostProcessedFieldsWritten>0){
     data_out_Scalar.add_data_vector (material, "meshGrain_ID");
     data_out_Scalar.build_patches ();
   }

  //write to results file
  std::string dir(userInputs.outputDirectory);
  dir+="/";

  //
  unsigned int incrementDigits= (totalIncrements<10000 ? 4 : std::ceil(std::log10(totalIncrements))+1);
  unsigned int domainDigits   = (Utilities::MPI::n_mpi_processes(mpi_communicator)<10000 ? 4 : std::ceil(std::log10(Utilities::MPI::n_mpi_processes(mpi_communicator)))+1);

  const std::string filename = (dir+"solution-" +
				Utilities::int_to_string (currentIncrement,incrementDigits) +
				"." +
				Utilities::int_to_string (triangulation.locally_owned_subdomain(),
							  domainDigits));
  std::ofstream outputFile ((filename + ".vtu").c_str());
  data_out.write_vtu (outputFile);
  //write projected fields, if any
  if (numPostProcessedFieldsWritten>0){
    const std::string filenameForProjectedFields = (dir+"projectedFields-" +
						    Utilities::int_to_string (currentIncrement,incrementDigits) +
						    "." +
						    Utilities::int_to_string (triangulation.locally_owned_subdomain(),
									      domainDigits));
    std::ofstream outputFileForProjectedFields ((filenameForProjectedFields + ".vtu").c_str());
    data_out_Scalar.write_vtu (outputFileForProjectedFields);
  }


  //create pvtu record
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
    std::vector<std::string> filenames, filenamesForProjectedFields;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i){
      filenames.push_back ("solution-" +
			   Utilities::int_to_string (currentIncrement, incrementDigits) +
			   "." +
			   Utilities::int_to_string (i, domainDigits) +
			   + ".vtu");
      if (numPostProcessedFieldsWritten>0){
	filenamesForProjectedFields.push_back ("projectedFields-" +
					       Utilities::int_to_string (currentIncrement, incrementDigits) +
					       "." +
					       Utilities::int_to_string (i, domainDigits) +
					       + ".vtu");
      }
    }
    const std::string filenamepvtu = (dir+"solution-" +
				      Utilities::int_to_string (currentIncrement,incrementDigits) +
				      ".pvtu");
    std::ofstream master_output (filenamepvtu.c_str());
    data_out.write_pvtu_record (master_output, filenames);
    pcout << "output written to: " << filenamepvtu.c_str();
    //
    if (numPostProcessedFieldsWritten>0){
      const std::string filenamepvtuForProjectedFields = (dir+"projectedFields-" +
							  Utilities::int_to_string (currentIncrement,incrementDigits) +
							  ".pvtu");
      std::ofstream master_outputForProjectedFields (filenamepvtuForProjectedFields.c_str());
      data_out_Scalar.write_pvtu_record (master_outputForProjectedFields, filenamesForProjectedFields);
      pcout << " and " << filenamepvtuForProjectedFields.c_str() ;
    }
    pcout << " \n\n";
  }
}
#include "../../include/ellipticBVP_template_instantiations.h"
