//tension BVP
//general headers
#include <fstream>
#include <sstream>

//parameters file
#include "parameters.h"

//dealIIheaders
#include "../../../src/materialModels/continuumPlasticity/continuumPlasticity.h"

//generate or import mesh
template <int dim>
void continuumPlasticity<dim>::mesh(){
  //creating mesh
  this->pcout << "generating problem mesh\n";

  //Define the limits of the domain (this example is in 3D)
  double x_max = 10., y_max = 1., z_max = 1.;
  Point<dim,double> min(0.,0.,0.), max(x_max,y_max,z_max);

  //Define the mesh refinement
  std::vector<unsigned int> numberOfElements(dim,std::pow(2,meshRefineFactor));
  numberOfElements[0] *= 10;

  GridGenerator::subdivided_hyper_rectangle (this->triangulation, numberOfElements, min, max);

  //Output image of the mesh in eps format                                                                                                           
  if ((this->triangulation.n_global_active_cells()<1000) and (Utilities::MPI::n_mpi_processes(this->mpi_communicator)==1)){
    std::ofstream out ("mesh.eps");
    GridOut grid_out;
    grid_out.write_eps (this->triangulation, out);
    this->pcout << "writing mesh image to mesh.eps" << std::endl;
  }
}
 
//Mark boundaries for applying Dirichlet BC's
template <int dim>
void continuumPlasticity<dim>::markBoundaries(){
  typename DoFHandler<dim>::active_cell_iterator 
    cell = this->dofHandler.begin_active(), 
    endc = this->dofHandler.end();

  //All boundaries are by marked with flag '0' by default. 
  //To pick specific boundaries, one needs to mark them 
  //with integer flags and use those flags in apply_dirichlet_conditons()
  for (;cell!=endc; ++cell){
    if (cell->is_locally_owned()){
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f){
	if (cell->face(f)->at_boundary()){
	  const Point<dim> face_center = cell->face(f)->center();
	  if (face_center[0] == 0.0){
	    cell->face(f)->set_boundary_indicator (1); //back boundary
	  }
	  else if(face_center[0] == 10.0){
	    cell->face(f)->set_boundary_indicator (2); //front boundary
	  }
	}
      }
    }
  }
}


//Class to set Dirichlet BC values 
template <int dim>
class BCFunction : public Function<dim>{
public:
  BCFunction(): Function<dim> (dim){}
  void vector_value (const Point<dim>   &p, Vector<double>   &values) const{
    Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));    
    values[2]=-0.1/totalNumIncrements; // total displacement along Z-Direction divide by total increments
  }
};

//Apply Dirchlet BCs for simple bending BVP
template <int dim>
void continuumPlasticity<dim>::applyDirichletBCs(){
  this->constraints.clear ();
  this->constraints.reinit (this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (this->dofHandler, this->constraints);
  std::vector<bool> allComponents (dim, true);
  //u3 applied on X1=1,u1=u2=0
  if (this->currentIteration==0) {
    VectorTools:: interpolate_boundary_values (this->dofHandler,
					       2, 
					       BCFunction<dim>(), 
					       this->constraints,
					       allComponents);
  }
  //Don't apply further displacment simply for a new solver iteration.
  else {
    VectorTools:: interpolate_boundary_values (this->dofHandler,
					       2, 
					       ZeroFunction<dim>(dim),
					       this->constraints,
					       allComponents);
  }
  //u1=u2=u3=0 on X1=0
  VectorTools:: interpolate_boundary_values (this->dofHandler, 
					     1, 
					     ZeroFunction<dim>(dim), 
					     this->constraints, 
					     allComponents);

  this->constraints.close ();
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);
      continuumPlasticity<3> problem;

      //Read material parameters
      problem.properties.lambda = lame_lambda;
      problem.properties.mu = lame_mu;
      problem.properties.tau_y = yield_stress;
      problem.properties.K = strain_hardening;

      //Read pfunction names for strain energy density and yield functions
      problem.properties.strainEnergyModel = strain_energy_function;
      problem.properties.yieldModel = yield_function;

      problem.run ();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  
  return 0;
}

