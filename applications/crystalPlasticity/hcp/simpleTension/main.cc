//tension BVP
//general headers
#include <fstream>
#include <sstream>
using namespace std;


//parameters
#include "parameters.h"


//dealIIheaders
#include "../../../../src/materialModels/crystalPlasticity/hcp/model.h"

//Specify Dirichlet boundary conditions 
template <int dim>
void crystalPlasticity<dim>::setBoundaryValues(const Point<dim>& node, const unsigned int dof, bool& flag, double& value){
  //back boundary:   u_x=0
  if (node[0] == 0.0){
    if (dof==0) {flag=true; value=0.0;}
  }
  //front boundary:  u_x=0.001
  if (node[0] == 1.0){
    if (dof==0) {flag=true; value=0.001;}
  }
  //left boundary:   u_y=0
  if (node[1] == 0.0){
    if (dof==1) {flag=true; value=0.0;}
  }
  //bottom boundary: u_z=0
  if (node[2] == 0.0){
    if (dof==2) {flag=true; value=0.0;}
  }
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
