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

    GridGenerator::subdivided_hyper_rectangle (triangulation2, userInputs.subdivisions, Point<dim>(), Point<dim>(userInputs.span[0],userInputs.span[1],userInputs.span[2]), true);

   //////In the case of element deletion, this section will be active
    if(userInputs.enableElementDeletion){
      unsigned int gID;
      orientations_Mesh.loadOrientations(userInputs.grainIDFile,
              userInputs.headerLinesGrainIDFile,
              userInputs.grainOrientationsFile,
              userInputs.numPts,
              userInputs.span);
      typename Triangulation< dim, dim >::active_cell_iterator cell = triangulation2.begin_active(), endc = triangulation2.end();
      for (; cell!=endc; ++cell) {
              double pnt3[3];
              const Point<dim> pnt2=cell->center();
              for (unsigned int i=0; i<dim; ++i){
                  pnt3[i]=pnt2[i];
              }
              //Do you want cell centers or quadrature
              gID=orientations_Mesh.getMaterialID(pnt3);

              cellOrientationMap_Mesh.push_back(gID);
      }

      unsigned int materialID;
      unsigned int cell2=0;
      cell = triangulation2.begin_active(); endc = triangulation2.end();
      std::set<typename Triangulation< dim, dim >::active_cell_iterator> cells_to_remove;
      unsigned int NumberOwned=0;
      for (; cell!=endc; ++cell) {

        materialID=cellOrientationMap_Mesh[cell2];

          NumberOwned=NumberOwned+1;
          if (materialID ==userInputs.deletionGrainID){
            cells_to_remove.insert (cell);
          }
          cell2=cell2+1;
        //}
      }
      char buffer[200];
      unsigned int CellRemove=cells_to_remove.size();
      sprintf(buffer,
        "cells_to_remove %3u \n",
        CellRemove);
        pcout << buffer;
        sprintf(buffer,
          "NumberOwned %3u \n",
          NumberOwned);
          pcout << buffer;
       if (cells_to_remove.size()>0){
        GridGenerator::create_triangulation_with_removed_cells(triangulation2,cells_to_remove,triangulation);
       }
       else{
        triangulation.copy_triangulation(triangulation2);
       }
    }
    else {
      triangulation.copy_triangulation(triangulation2);
    }

  //In the case of Periodic BCs, it connects the periodic faces to each other,
  //and add those dofs as ghost cells of the other face.
    if(userInputs.enablePeriodicBCs){
      // Set which (if any) faces of the triangulation are periodic
      setPeriodicity();
    }



    triangulation.refine_global (userInputs.meshRefineFactor);

    if(userInputs.enableIndentationBCs){
          // Set which (if any) faces of the triangulation are indentation
          meshRefineIndentation();
    }




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
