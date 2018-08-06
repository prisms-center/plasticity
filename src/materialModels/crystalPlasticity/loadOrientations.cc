#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::loadOrientations(){
    QGauss<dim>  quadrature(this->userInputs.quadOrder);
    const unsigned int num_quad_points = quadrature.size();
    FEValues<dim> fe_values (this->FE, quadrature, update_quadrature_points);
    unsigned int gID;

    //loop over elements
    typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){

            fe_values.reinit(cell);
            double pnt3[3];
            const Point<dim> pnt2=cell->center();
            for (unsigned int i=0; i<dim; ++i){
                pnt3[i]=pnt2[i];
            }


            if(this->userInputs.readExternalMesh)
              gID=cell->material_id();
            else
            //Do you want cell centers or quadrature
              gID=orientations.getMaterialID(pnt3);

            cellOrientationMap.push_back(gID);

        }
    }
}

#include "../../../include/crystalPlasticity_template_instantiations.h"
