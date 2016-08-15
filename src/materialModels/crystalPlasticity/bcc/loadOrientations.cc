
template <int dim>
void crystalPlasticity<dim>::loadOrientations(){
    QGauss<dim>  quadrature(quadOrder);
    const unsigned int num_quad_points = quadrature.size();
    FEValues<dim> fe_values (this->FE, quadrature, update_quadrature_points);
    //loop over elements
    typename DoFHandler<dim>::active_cell_iterator cell = this->dofHandler.begin_active(), endc = this->dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){
            quadratureOrientationsMap.push_back(std::vector<unsigned int>(num_quad_points,0));

            
            fe_values.reinit(cell);
            double pnt3[3];
            const Point<dim> pnt2=cell->center();
            for (unsigned int i=0; i<dim; ++i){
                pnt3[i]=pnt2[i];
            }
            
            
            //loop over quadrature points
            for (unsigned int q=0; q<num_quad_points; ++q){
                double pnt[3];
                pnt[0]=fe_values.get_quadrature_points()[q][0];
                pnt[1]=fe_values.get_quadrature_points()[q][1];
                pnt[2]=fe_values.get_quadrature_points()[q][2];
                //get orientation ID and store it in quadratureOrientationsMap
                unsigned int gID=orientations.getMaterialID(pnt3);
                //gID=(gID%10)*10+gID/10;
                //pcout << gid << " ";
                quadratureOrientationsMap.back()[q]=gID;
                //now one can ac orientations.getMaterialID(pnt) << std::endlcess the oreintation id for each quadrature point using quadratureOrientationsMap[cellID][q]
            }
            
        }
    }
    
}
