//constructor and destructor for ellipticBVP class

#ifndef ELLIPTICBVP_SRC_H
#define ELLIPTICBVP_SRC_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//constructor
template <int dim>
ellipticBVP<dim>::ellipticBVP ()
  :
  mpi_communicator (MPI_COMM_WORLD),
  triangulation (mpi_communicator,
		 typename Triangulation<dim>::MeshSmoothing
		 (Triangulation<dim>::smoothing_on_refinement |
		  Triangulation<dim>::smoothing_on_coarsening)),
  FE (FE_Q<dim>(feOrder), dim),
  FE_Scalar (FE_Q<dim>(feOrder), 1),
  dofHandler (triangulation),
  dofHandler_Scalar (triangulation),
  currentIteration(0),
  currentIncrement(0),
  totalIncrements(totalNumIncrements),
  resetIncrement(false),
  loadFactorSetByModel(1.0),
  totalLoadFactor(0.0),
  pcout (std::cout, Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0),
  computing_timer (pcout, TimerOutput::summary, TimerOutput::wall_times),
  numPostProcessedFields(0)
{
  //Nodal Solution names - this is for writing the output file
  for (unsigned int i=0; i<dim; ++i){
    nodal_solution_names.push_back("u");
    nodal_data_component_interpretation.push_back(DataComponentInterpretation::component_is_part_of_vector);
  }
}

//destructor
template <int dim>
ellipticBVP<dim>::~ellipticBVP ()
{
}

#endif
