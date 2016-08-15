//tension BVP
//general headers
#include <fstream>
#include <sstream>
using namespace std;


//parameters
#include "parameters.h"


//dealIIheaders
#include "../../../../src/materialModels/crystalPlasticity/hcp/model.h"

//overload mesh() method to generate the required polycrystal geometry
template <int dim>
void crystalPlasticity<dim>::mesh(){
    //creating mesh
    this->pcout << "generating problem mesh\n";
    double spanX=1.0; //Span along x-axis
    double spanY=1.0; //Span along y-axis
    double spanZ=1.0; //Span along z-axis
    GridGenerator::hyper_rectangle (this->triangulation, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
    this->triangulation.refine_global (meshRefineFactor);
}

//Mark boundaries for applying Dirichlet BC's
template <int dim>
void crystalPlasticity<dim>::markBoundaries(){
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
                    else if(face_center[0] == 1.0){
                        cell->face(f)->set_boundary_indicator (2); //front boundary
                    }
                    else if(face_center[1] == 0.0){
                        cell->face(f)->set_boundary_indicator (3); //left boundary
                    }
                    else if(face_center[2] == 0.0){
                        cell->face(f)->set_boundary_indicator (4); //bottom boundary
                    }
                }
            }
        }
    }
}



//Class to set Dirichlet BC values 
template <int dim>
class BCFunction : public Function<dim>{
  double loadFactor;
public:
  BCFunction(double _loadFactor): Function<dim> (dim){loadFactor=_loadFactor;}
  void vector_value (const Point<dim>   &p, Vector<double>   &values) const{
    Assert (values.size() == dim, ExcDimensionMismatch (values.size(), dim));    
    values[0]=0.001*loadFactor; // displacement along X-Direction
  }
};

//Apply Dirchlet BCs for simple tension BVP
template <int dim>
void crystalPlasticity<dim>::applyDirichletBCs(){
  this->constraints.clear ();
  this->constraints.reinit (this->locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints (this->dofHandler, this->constraints);
  std::vector<bool> mechanicsBoundary_Z1 (dim, false); mechanicsBoundary_Z1[0]=true;
  std::vector<bool> mechanicsBoundary_Z2 (dim, false); mechanicsBoundary_Z2[0]=true;
  std::vector<bool> mechanicsBoundary_Z3 (dim, false); mechanicsBoundary_Z3[1]=true;
  std::vector<bool> mechanicsBoundary_Z4 (dim, false); mechanicsBoundary_Z4[2]=true;
  //u1 applied on X1=1
  if (this->currentIteration==0) {
    VectorTools:: interpolate_boundary_values (this->dofHandler,
					       2, 
					       BCFunction<dim>(this->loadFactorSetByModel), 
					       this->constraints,
					       mechanicsBoundary_Z2);
  }
  else {
    VectorTools:: interpolate_boundary_values (this->dofHandler,
					       2, 
					       ZeroFunction<dim>(dim),
					       this->constraints,
					       mechanicsBoundary_Z2);
  }
  //u1=0 on X1=0
  VectorTools:: interpolate_boundary_values (this->dofHandler, 
					     1, 
					     ZeroFunction<dim>(dim), 
					     this->constraints, 
					     mechanicsBoundary_Z1);
  //u2=0 on X2=0
  VectorTools:: interpolate_boundary_values (this->dofHandler, 
					     3, 
					     ZeroFunction<dim>(dim), 
					     this->constraints, 
					     mechanicsBoundary_Z3);
  //u3=0 on X3=0
  VectorTools:: interpolate_boundary_values (this->dofHandler, 
					     4, 
					     ZeroFunction<dim>(dim), 
					     this->constraints, 
					     mechanicsBoundary_Z4);
  
  this->constraints.close ();
}


//main
int main (int argc, char **argv)
{
    Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
    try
    {
        deallog.depth_console(0);
        crystalPlasticity<3> problem;
        
        FullMatrix<double> m_alpha,n_alpha; // Slip directions and Slip Normals
        const unsigned int n_slip_systems=numSlipSystems; //No. of slip systems
        const unsigned int n_twin_systems=numTwinSystems; //No. of twin systems

        n_alpha.reinit(n_slip_systems,3);
        m_alpha.reinit(n_slip_systems,3);
        string line;
        ifstream myfile;
        
        //open data file to read slip normals
        ifstream slipNormalsDataFile("slipNormals.txt");
        //read data
        unsigned int id=0;
        if (slipNormalsDataFile.is_open()){
            cout << "reading slip Normals file\n";
            //read data
            while (getline (slipNormalsDataFile,line) && id<n_slip_systems){
                stringstream ss(line);
                ss >> n_alpha[id][0];
                ss >> n_alpha[id][1];
                ss >> n_alpha[id][2];
                //cout<<id<<'\t'<<n_alpha[id][0]<<'\t'<<n_alpha[id][1]<<'\t'<<n_alpha[id][2]<<'\n';
                id=id+1;
            }
        }
        else{
            cout << "Unable to open slipNormals.txt \n";
            exit(1);
        }
        
        //open data file to read slip directions
        ifstream slipDirectionsDataFile("slipDirections.txt");
        //read data
        id=0;
        if (slipDirectionsDataFile.is_open()){
            cout << "reading slip Directions file\n";
            //read data
            while (getline (slipDirectionsDataFile,line)&& id<n_slip_systems){
                stringstream ss(line);
                ss >> m_alpha[id][0];
                ss >> m_alpha[id][1];
                ss >> m_alpha[id][2];
                //cout<<id<<'\t'<<m_alpha[id][0]<<'\t'<<m_alpha[id][1]<<'\t'<<m_alpha[id][2]<<'\n';
                id=id+1;
            }
        }
        else{
            cout << "Unable to open slipDirections.txt \n";
            exit(1);
        }
        
        problem.properties.m_alpha=m_alpha;
        problem.properties.n_alpha=n_alpha;
        
        //reading materials atlas files
    
        double stencil[3]={1.0/(numPts[0]-1), 1.0/(numPts[1]-1), 1.0/(numPts[2]-1)}; // Dimensions of voxel
        problem.orientations.loadOrientations("grainID.txt",
                                              5,
                                              "orientations.txt",
                                              numPts,
                                              stencil);
        problem.orientations.loadOrientationVector("orientations.txt");
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
