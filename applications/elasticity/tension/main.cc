//tension BVP
//general headers
#include <fstream>
#include <sstream>

#define feOrder   1
#define quadOrder 2 
#define meshRefineFactor 2
#define writeOutput true
#define linearSolverType PETScWrappers::SolverCG
#define totalNumIncrements 1
#define maxLinearSolverIterations 5000
#define relLinearSolverTolerance  1.0e-8
#define maxNonLinearIterations 30
#define absNonLinearTolerance 1.0e-15
#define relNonLinearTolerance 1.0e-10
#define stopOnConvergenceFailure true

//material properties
#define lambda_val 1.0
#define mu_val 1.0

//dealIIheaders
#include "../../../src/materialModels/elasticity/elasticity.cc"
 
//Mark boundaries for applying Dirichlet BC's
template <int dim>
void elasticity<dim>::markBoundaries(){
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
	  if (face_center[0]==0.0){
	    cell->face(f)->set_boundary_indicator (1); //boundary at X=0.0 marked with flag '1'
	  }
	  else if (face_center[0]==1.0){
	    cell->face(f)->set_boundary_indicator (2); //boundary at X=1.0 marked with flag '2'
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
    values[0]=0.01; // displacement along X-Direction
  }
};

//Apply Dirchlet BCs for tension BVP
template <int dim>
void elasticity<dim>::applyDirichletBCs(){
  this->constraints.clear ();
  this->constraints.reinit (this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (this->dofHandler, this->constraints);
  std::vector<bool> allComponenents (dim, true); 
  std::vector<bool> xComponenent    (dim, false); xComponenent[0]=true;
  //u=0 along X=0
  VectorTools::interpolate_boundary_values (this->dofHandler,
					    1, 
					    ZeroFunction<dim>(dim),
					    this->constraints,
					    allComponenents);
  //u=0.01 along X=1.00
  VectorTools::interpolate_boundary_values (this->dofHandler,
					    2, 
					    BCFunction<dim>(),
					    this->constraints,
					    xComponenent);
  this->constraints.close ();
}

//main
int main (int argc, char **argv)
{
  Utilities::System::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  try
    {
      deallog.depth_console(0);
      elasticity<3> problem;
			
			//load material properties
			problem.properties.lambda=lambda_val;
			problem.properties.mu=mu_val;	 
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

