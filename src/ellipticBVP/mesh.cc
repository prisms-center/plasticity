//mesh generation/importing method for ellipticBVP class
#include <fstream>
#include "../../include/ellipticBVP.h"

//generate or import mesh
template <int dim>
void ellipticBVP<dim>::mesh(){
  if(userInputs.readExternalMesh){
    //reading external mesh
    pcout << "reading problem mesh\n";
    GridIn<dim> gridin;
    gridin.attach_triangulation(this->triangulation);
    //Read mesh in UCD format generated from Cubit
    std::ifstream extMesh(userInputs.externalMeshFileName);
    gridin.read_msh(extMesh);

    //Output image for viewing
    std::string dir(userInputs.outputDirectory);
    dir+="/";

    const std::string filename = (dir+"mesh.vtk");
    std::ofstream out ((filename).c_str());
    GridOut grid_out;
    grid_out.write_vtk (this->triangulation, out);
    pcout << "writing mesh image to mesh.vtk\n";
  }

  else{
    //creating mesh
    pcout << "generating problem mesh\n";

    GridGenerator::subdivided_hyper_rectangle (triangulation, userInputs.subdivisions, Point<dim>(), Point<dim>(userInputs.span[0],userInputs.span[1],userInputs.span[2]));
    triangulation.refine_global (userInputs.meshRefineFactor);

    //Output image of the mesh in eps format
    if(userInputs.writeMeshToEPS)
      if ((triangulation.n_global_active_cells()<10000) and (Utilities::MPI::n_mpi_processes(mpi_communicator)==1)){
        std::string dir(userInputs.outputDirectory);
        dir+="/";

        const std::string filename = (dir+"mesh.eps");
        std::ofstream out ((filename).c_str());
        GridOut grid_out;
        grid_out.write_eps (triangulation, out);
        pcout << "writing mesh image to mesh.eps" << std::endl;
      }
  }
}
#include "../../include/ellipticBVP_template_instantiations.h"
