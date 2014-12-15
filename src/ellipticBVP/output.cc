//output method for ellipticBVP class

#ifndef OUTPUT_ELLIPTICBVP_H
#define OUTPUT_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//output results
template <int dim>
void ellipticBVP<dim>::output(){
  DataOut<dim> data_out;
  data_out.attach_dof_handler (dofHandler);
  data_out.add_data_vector (solutionLocal, 
			    nodal_solution_names, 
			    DataOut<dim>::type_dof_data, 
			    nodal_data_component_interpretation);
  Vector<float> subdomain (triangulation.n_active_cells());
  for (unsigned int i=0; i<subdomain.size(); ++i)
    subdomain(i) = triangulation.locally_owned_subdomain();
  data_out.add_data_vector (subdomain, "subdomain");
  data_out.build_patches ();

  //write to results file
  const std::string filename = ("solution-" +
				Utilities::int_to_string (currentIncrement, std::ceil(std::log10(totalIncrements))+1) + 
				"." +
				Utilities::int_to_string (triangulation.locally_owned_subdomain(), 
							  std::ceil(std::log10(Utilities::MPI::n_mpi_processes(mpi_communicator)))+1));
  std::ofstream outputFile ((filename + ".vtu").c_str());
  data_out.write_vtu (outputFile);
  
  //create pvtu record
  if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
    std::vector<std::string> filenames;
    for (unsigned int i=0;i<Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD); ++i)
      filenames.push_back ("solution-" +
			   Utilities::int_to_string (currentIncrement, std::ceil(std::log10(totalIncrements))+1) + 
			   "." +
			   Utilities::int_to_string (triangulation.locally_owned_subdomain(), 
						     std::ceil(std::log10(Utilities::MPI::n_mpi_processes(mpi_communicator)))+1)
			   + ".vtu");
    std::ofstream master_output ((filename + ".pvtu").c_str());
    data_out.write_pvtu_record (master_output, filenames);
    pcout << "output written to: " << (filename + ".pvtu").c_str() << "\n\n";
  }
}

#endif

