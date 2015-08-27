//tension BVP
//general headers
#include <fstream>
#include <sstream>

#define feOrder   1
#define quadOrder 2 
#define meshRefineFactor 3
#define writeOutput true
#define linearSolverType PETScWrappers::SolverCG
#define totalNumIncrements 400
#define maxLinearSolverIterations 5000
#define relLinearSolverTolerance  1.0e-12
#define maxNonLinearIterations 30
#define absNonLinearTolerance 1.0e-15
#define relNonLinearTolerance 1.0e-10
#define stopOnConvergenceFailure true

//Read json input 
#include "../../../utils/json/json_spirit.h"
#include "../../../utils/json/json_spirit_reader_template.h"

//dealIIheaders
#include "../../../src/materialModels/continuumPlasticity/continuumPlasticity.cc"

//generate or import mesh
template <int dim>
void continuumPlasticity<dim>::mesh(){
  //creating mesh
  this->pcout << "generating problem mesh\n";

  //Define the limits of the domain (this example is in 3D)
  double x_max = 5., y_max = 1., z_max = 1.;
  Point<dim,double> min(0.,0.,0.), max(x_max,y_max,z_max);

  //Define the mesh refinement - more refined near constrained end
  unsigned int meshSize = std::pow(2.,meshRefineFactor);
  double stepSize = 1./std::pow(2.,meshRefineFactor);
  std::vector< std::vector<double> > stepSizes(dim,std::vector<double>(meshSize,stepSize));
  stepSizes[0].resize(7*meshSize,stepSize);
  for(unsigned int i=0; i<2*meshSize; i++){
    stepSizes[0][i] = stepSize/4.;
  }
  for(unsigned int i=2*meshSize; i<3*meshSize; i++){
    stepSizes[0][i] = stepSize/2.;
  }

  GridGenerator::subdivided_hyper_rectangle (this->triangulation, stepSizes, min, max,false);

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
	  if (face_center[0]==0.0){
	    cell->face(f)->set_boundary_indicator (1); //boundary at X=0.0 marked with flag '1'
	  }
	  else if (face_center[0]==5.0){
	    cell->face(f)->set_boundary_indicator (2); //boundary at X=5.0 marked with flag '2'
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
    values[0]=0.1/totalNumIncrements; //total displacement along X-Direction divide by total increments
  }
};

//Apply Dirchlet BCs for constrained tension BVP
template <int dim>
void continuumPlasticity<dim>::applyDirichletBCs(){
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
  //u=0.5 along X=5.00
  if (this->currentIteration==0){
    VectorTools::interpolate_boundary_values (this->dofHandler,
					      2, 
					      BCFunction<dim>(),
					      this->constraints,
					      xComponenent);
  }
  else{
    VectorTools::interpolate_boundary_values (this->dofHandler,
					      2, 
					      ZeroFunction<dim>(dim),
					      this->constraints,
					      xComponenent);
  }
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
			
      //load material properties
      std::ifstream is("materialProperties.json");
      json_spirit::Value value;
      json_spirit::read_stream(is,value);
      std::vector< json_spirit::Pair > material, strainEnergy, yield;

      //Read material parameters
      material = value.get_obj()[0].value_.get_obj();
      problem.properties.lambda = material[0].value_.get_real();
      problem.properties.mu = material[1].value_.get_real();
      problem.properties.tau_y = material[2].value_.get_real();
      problem.properties.K = material[3].value_.get_real();

      //Read pfunction names for strain energy density and yield functions
      strainEnergy = material[5].value_.get_obj();
      problem.properties.strainEnergyModel = strainEnergy[1].value_.get_str();
      yield = material[6].value_.get_obj();
      problem.properties.yieldModel = yield[1].value_.get_str();

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

