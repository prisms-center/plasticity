//mesh generation/importing method for ellipticBVP class
#include "../../include/ellipticBVP.h"

#ifndef MESH_ELLIPTICBVP_H
#define MESH_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//generate or import mesh
template <int dim>
void ellipticBVP<dim>::mesh(){
  //creating mesh
  pcout << "generating problem mesh\n";
  //
  std::vector<unsigned int> subdivisions;
  subdivisions.push_back(subdivisionsX);
  subdivisions.push_back(subdivisionsY);
  subdivisions.push_back(subdivisionsZ);
  GridGenerator::subdivided_hyper_rectangle (triangulation, subdivisions, Point<dim>(), Point<dim>(spanX,spanY,spanZ));
  triangulation.refine_global (meshRefineFactor);

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

#endif
