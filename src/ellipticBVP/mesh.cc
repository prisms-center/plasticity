//mesh generation/importing method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//generate or import mesh
template <int dim>
void ellipticBVP<dim>::mesh(){
  //creating mesh
  pcout << "generating problem mesh\n";
  //
  GridGenerator::subdivided_hyper_rectangle (triangulation, userInputs.subdivisions, Point<dim>(), Point<dim>(userInputs.span));
  triangulation.refine_global (userInputs.meshRefineFactor);

  //Output image of the mesh in eps format
#ifdef  writeMeshToEPS
#if writeMeshToEPS==true
  if ((triangulation.n_global_active_cells()<10000) and (Utilities::MPI::n_mpi_processes(mpi_communicator)==1)){
    std::ofstream out ("mesh.eps");
    GridOut grid_out;
    grid_out.write_eps (triangulation, out);
    pcout << "writing mesh image to mesh.eps" << std::endl;
  }
#endif
#endif
}
