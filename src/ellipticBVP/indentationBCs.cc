//solve method for ellipticBVP class
#include "../../include/ellipticBVP.h"
#include <fstream>
#include <iostream>

//template <int dim>
//void ellipticBVP<dim>::setIndentation(){
//    //This functions connects the periodic faces to each other, and add those dofs as ghost cells of the other face.
//    std::vector<GridTools::PeriodicFacePair<typename parallel::distributed::Triangulation<dim>::cell_iterator> > periodicity_vector;
//
/////Here, we assumed that in the face constraints lines, we always start with FXP, FYP, or FZP.
//    GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 1, /*b_id2*/ 0, /*direction*/ 0, periodicity_vector);
//    GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 3, /*b_id2*/ 2, /*direction*/ 1, periodicity_vector);
//    GridTools::collect_periodic_faces(triangulation, /*b_id1*/ 5, /*b_id2*/ 4, /*direction*/ 2, periodicity_vector);
//
//
//    triangulation.add_periodicity(periodicity_vector);
//    pcout << "periodic facepairs: " << periodicity_vector.size() << std::endl;
//}
template <int dim>
void ellipticBVP<dim>::meshRefineIndentation() {
    pcout <<"Mesh refinement for indentation \n";
    // Plan: use final indenter position to determine Primary Zone (some factor times size and some factor times displacement vector
    // for included plastic volume)
    // Looked to step 1 for mesh refinement basics
    //
    // In order to demonstrate how to write a loop over all cells, we will
    // refine the grid in five steps towards the inner circle of the domain:
    Point<dim> centerPlasticZone;
    for (unsigned int i = 0; i <dim; ++i)
    {
        centerPlasticZone(i)= userInputs.refinementCenter[i];
    }
    for (unsigned int radius_j = 0; radius_j < userInputs.refinementFactor.size(); ++radius_j)
    {
        for (unsigned int step = 0; (int)step < userInputs.refinementFactor[radius_j]; ++step)
        {
            for (auto &cell : triangulation.active_cell_iterators())
            {
              #if (DEAL_II_VERSION_MAJOR > 9)
                for (const auto v : cell->vertex_indices())
              #elseif ((DEAL_II_VERSION_MAJOR == 9)&&(DEAL_II_VERSION_MINOR >= 3))
                for (const auto v : cell->vertex_indices())
              #else
                for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
              #endif
                {
                   // REFINING SPHERICALLY LOCKS CELLS!!!
//                    const double distance_from_center =
//                            centerPlasticZone.distance(cell->vertex(v));
                    double distance_from_center = 0.0;
                    double component_distance = 0.0;
                    for (unsigned int i = 0; i <dim; ++i)
                        {
                            component_distance = std::abs(centerPlasticZone[i] - cell->vertex(v)[i]);
                            if (component_distance > distance_from_center)
                            {
                                distance_from_center = component_distance;
                            }
                        }

                    if (distance_from_center <= userInputs.refinementZoneSize[radius_j])
                    {
                        cell->set_refine_flag();
                        break;
                    }
                }
            }
            triangulation.execute_coarsening_and_refinement();
        }
    }

}

template <int dim>
void ellipticBVP<dim>::updateIndentPos() {
    double pos;
    unsigned int nextFrame, prevFrame;
    for (unsigned int j = 0; j < userInputs.indentationKeyFrames; ++j){
        if (userInputs.indentationKeyFrameTimeIncs[j] > (int)currentIncrement){
            nextFrame = j;
            prevFrame = j-1;
            break;
        }

    }
    pos = (currentIncrement + 1 - userInputs.indentationKeyFrameTimeIncs[prevFrame]);
    pos = pos/(userInputs.indentationKeyFrameTimeIncs[nextFrame] -
            userInputs.indentationKeyFrameTimeIncs[prevFrame]) ;
    for (int i = 0; i < dim; ++i) {
        prevPosIndenter[i] = currentPosIndenter[i];
        currentPosIndenter[i] = KeyPosIndenter[prevFrame](i)
                + pos * (KeyPosIndenter[nextFrame](i) -
                KeyPosIndenter[prevFrame](i));
    }
    pcout <<"previous Position " << pos << " (x,y,z) "<<prevPosIndenter[0]<<" "
          << prevPosIndenter[1]<< " "<<prevPosIndenter[2]<<"\n";
    pcout <<"current Position " << pos << " (x,y,z) "<<currentPosIndenter[0]<<" "
    << currentPosIndenter[1]<< " "<<currentPosIndenter[2]<<"\n";
    currentIndentDisp = KeyPosIndenter[prevFrame](indentDof)
                        + pos * (KeyPosIndenter[nextFrame](indentDof) -
                                 KeyPosIndenter[prevFrame](indentDof)) -
            KeyPosIndenter[0](indentDof);
}

template <int dim>
void ellipticBVP<dim>::displaceNode(const Point<dim> & p, const Point<dim> & u)
{
    displace_local.resize(dim);
    for (int i = 0; i < dim; ++i) {
        displace_local[i] = p(i) + u(i);
    }
}

//INDENTATION
template <int dim>
bool ellipticBVP<dim>::flagActiveSet(const Point<dim> & p)
{
    double dist2;
    double dist1;
    double tmp;
    dist2 = 0.0;
    dist1 = 0.0;
    for (unsigned int i = 0; i < dim; i=i+1)
    {
        switch (indenterShape) {
            case 0: {
                tmp = ((double)currentPosIndenter[i] - (double)p(i)) * ((double)currentPosIndenter[i] - (double)p(i));
                dist2 = (double)dist2 + (double)tmp;
//                pcout << tmp << " " << std::setprecision(6) << dist2 << "\n";
                break;
            }
            case 1: {
                dist2 = dist2 + fabs((double)currentPosIndenter[i] - (double)p(i));
                break;
            }
            case 2: {
                if (i < dim - 1)
                    dist2 = dist2 + ((double)currentPosIndenter[i] - (double)p(i)) * ((double)currentPosIndenter[i] - (double)p(i));
                else
                    dist1 = dist1 + fabs((double)currentPosIndenter[i] - (double)p(i));
                break;
            }
        }
    }
//    pcout << p(0)<<" "<<p(1)<<" "<<p(2)<<"\n";
//    pcout << currentPosIndenter[0]<<" "<<currentPosIndenter[1]<<" "<<currentPosIndenter[2]<<"\n";
//    pcout<<"shape "<<indenterShape<<" size "<< indenterSize << " dist2 " << dist2 << "\n";
    switch (indenterShape) {
        case 0:{
            if (dist2 - indenterTolerance < (double)indenterSize * (double)indenterSize)
                return true;
            else
                return false;
        }
        case 1:{
            if (dist2 - indenterTolerance< (double)indenterSize)
                return true;
            else
                return false;
        }
        case 2: {
            if (dist1 - indenterTolerance < indenterSize && dist2 - indenterTolerance < (double)indenterSize * (double)indenterSize)
                return true;
            else
                return false;
        }
    }
    return false;
}

template <int dim>
bool ellipticBVP<dim>::flagActiveSetLambda(const Point<dim> & p, double & criterion)
{
    double dist2;
    double dist1;
    double tmp;

    dist2 = 0.0;
    dist1 = 0.0;
    for (unsigned int i = 0; i < dim; i=i+1)
    {
        switch (indenterShape) {
            case 0: {
                tmp = ((double)currentPosIndenter[i] - (double)p(i)) * ((double)currentPosIndenter[i] - (double)p(i));
                dist2 = (double)dist2 + (double)tmp;
//                pcout << tmp << " " << std::setprecision(6) << dist2 << "\n";
                break;
            }
            case 1: {
                dist2 = dist2 + fabs((double)currentPosIndenter[i] - (double)p(i));
                break;
            }
            case 2: {
                if (i < dim - 1)
                    dist2 = dist2 + ((double)currentPosIndenter[i] - (double)p(i)) * ((double)currentPosIndenter[i] - (double)p(i));
                else
                    dist1 = dist1 + fabs((double)currentPosIndenter[i] - (double)p(i));
                break;
            }
        }
    }
//    pcout << p(0)<<" "<<p(1)<<" "<<p(2)<<"\n";
//    pcout << currentPosIndenter[0]<<" "<<currentPosIndenter[1]<<" "<<currentPosIndenter[2]<<"\n";
//    pcout<<"shape "<<indenterShape<<" size "<< indenterSize << " dist2 " << dist2 << "\n";
    switch (indenterShape) {
        case 0:{
            if (dist2 - indenterTolerance < (double)indenterSize * (double)indenterSize){
                std::cout<< "criterion "<<criterion<<std::endl;
                if (criterion >= 0)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
        case 1:{
            if (dist2 - indenterTolerance< (double)indenterSize){
                if (criterion >= 0)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
        case 2: {
            if (dist1 - indenterTolerance < indenterSize &&
                dist2 - indenterTolerance < (double) indenterSize * (double) indenterSize) {

            if (criterion >= 0)
                return true;
            else
                return false;
        }
            else
                return false;
        }
    }
    return false;
}

template <int dim>
bool ellipticBVP<dim>::flagActiveSetLambda2(const Point<dim> & p, double & criterion)
{
    double dist2;
    double dist1;
    double tmp;

    dist2 = 0.0;
    dist1 = 0.0;
    for (unsigned int i = 0; i < dim; i=i+1)
    {
        switch (indenterShape) {
            case 0: {
                tmp = ((double)currentPosIndenter[i] - (double)p(i)) * ((double)currentPosIndenter[i] - (double)p(i));
                if (i!=indentDof)
                    dist2 = (double)dist2 + (double)tmp;
//                pcout << tmp << " " << std::setprecision(6) << dist2 << "\n";
                break;
            }
            case 1: {
                if (i!=indentDof)
                    dist2 = dist2 + fabs((double)currentPosIndenter[i] - (double)p(i));
                break;
            }
            case 2: {
                if (i != indentDof)
                    dist2 = dist2 + ((double)currentPosIndenter[i] - (double)p(i)) * ((double)currentPosIndenter[i] - (double)p(i));
                else
                    dist1 = dist1 + fabs((double)currentPosIndenter[i] - (double)p(i));
                break;
            }
        }
    }
//    pcout << p(0)<<" "<<p(1)<<" "<<p(2)<<"\n";
//    pcout << currentPosIndenter[0]<<" "<<currentPosIndenter[1]<<" "<<currentPosIndenter[2]<<"\n";
//    pcout<<"shape "<<indenterShape<<" size "<< indenterSize << " dist2 " << dist2 << "\n";
    switch (indenterShape) {
        case 0:{
            if (dist2 - indenterTolerance < (double)indenterSize * (double)indenterSize){
//                std::cout<< "criterion "<<criterion<<std::endl;
                if (criterion >= 0)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
        case 1:{
            if (dist2 - indenterTolerance< (double)indenterSize){
                if (criterion >= 0)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
        case 2: {
            if (dist1 - indenterTolerance < indenterSize &&
                dist2 - indenterTolerance < (double) indenterSize * (double) indenterSize) {

                if (criterion >= 0)
                    return true;
                else
                    return false;
            }
            else
                return false;
        }
    }
    return false;
}


template <int dim>
double ellipticBVP<dim>::Obstacle(const Point<dim> & p, const unsigned int & component, const std::vector<double> & ind)
{
    //Currently f1 and f2 are free/fixed in place if indentDof is interfering with indenter
    // use current solutions for the displacements and use to calculate the delta u imposed by indenter
    //Future options...
    // mindist projection out of indenter for all components
    // calculation of Dual solution
    // calculation of dual solution with friction tangent to surface
    // how do we track movement parallel to surface if surface isn't constant? separate loop?
    unsigned int f1; //Other DOF 1
    unsigned int f2; //Other DOF 2
    if (indentDof == 2){
        f1 = 0;
        f2 = 1;
    }
    else if (indentDof == 1){
        f1 = 0;
        f2 = 2;
    }
    else {
        f1 = 1;
        f2 = 2;
    }

    if (component == f1)
        return 0;//p(0);
    else if (component == f2)
        return 0;//p(1);
    else if (component == indentDof)
    {
        switch (indenterShape) {
            case 0: {
                if ((p(f1) - ind[f1]) * (p(f1) - ind[f1]) + (p(f2) - ind[f2]) * (p(f2) - ind[f2]) < indenterSize*indenterSize){
                    return (ind[indentDof] - std::sqrt(indenterSize*indenterSize - (p(f1) - ind[f1]) * (p(f1) - ind[f1]) -
                                      (p(f2) - ind[f2]) * (p(f2) - ind[f2])) -
                           p(indentDof));
                }
                else
                    return 0;
            }
            case 1: {
                if ((p(f1) - ind[f1]) + (p(f2) - ind[f2])  < indenterSize){
                    return (ind[indentDof] -(indenterSize - (p(f1) - ind[f1])  -
                                       (p(f2) - ind[f2]) ) - p(indentDof));
            }
                else
                    return 0;
            }
            case 2: {
                if ((p(f1) - ind[f1]) * (p(f1) - ind[f1]) +
                    (p(f2) - ind[f2]) * (p(f2) - ind[f2]) < indenterSize * indenterSize) {
                    return (ind[indentDof] - indenterSize - p(indentDof));
                }
                else
                    return 0;
            }
        }
//        if ((p(0) - 0.5) * (p(0) - 0.5) + (p(1) - 0.5) * (p(1) - 0.5) < 0.36)
//            return (-std::sqrt(0.36 - (p(0) - 0.5) * (p(0) - 0.5) -
//                               (p(1) - 0.5) * (p(1) - 0.5)) +
//                    z_surface + 0.59);
//        else
//            return 1000;
    }

    Assert(false, ExcNotImplemented());
    return 1e9; // an unreasonable value; ignored in debug mode because of the
    // preceding Assert
}

template <int dim>
void ellipticBVP<dim>::assemble_mass_matrix_diagonal()
{
    QGaussLobatto<dim - 1> face_quadrature_formula(FE.degree + 1);

    FEFaceValues<dim> fe_values_face(FE,
                                     face_quadrature_formula,
                                     update_values | update_JxW_values);

    const unsigned int dofs_per_cell   = FE.n_dofs_per_cell();
    const unsigned int n_face_q_points = face_quadrature_formula.size();

    FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

    const FEValuesExtractors::Vector displacement(0);
    pcout << "inside assemble mass matrix diagonal "<< std::endl;
    for (const auto &cell : dofHandler.active_cell_iterators())
        if (cell->is_locally_owned()) {
            //std::cout << "cell locally owned! "<< std::endl;
            for (const auto &face: cell->face_iterators())
                if (face->at_boundary() && (face->boundary_id() == indenterFace || userInputs.readExternalMesh)) {
                    fe_values_face.reinit(cell, face);
                    cell_matrix = 0;

                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                         ++q_point)
                        for (unsigned int i = 0; i < dofs_per_cell; ++i){
//                            std::cout << "cell matrix (" <<i<<") disp "<<fe_values_face[displacement].value(i, q_point)
//                            << " JxW " << fe_values_face.JxW(q_point) << std::endl;
                            cell_matrix(i, i) +=
                                    (fe_values_face[displacement].value(i, q_point) *
                                     fe_values_face[displacement].value(i, q_point) *
                                     fe_values_face.JxW(q_point));
                        }

                    cell->get_dof_indices(local_dof_indices);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
//                        std::cout << "cell matrix (" <<i<<", "
//                        << local_dof_indices[i]<< "): " << cell_matrix(i, i) << std::endl;
                        massMatrix.add(local_dof_indices[i],
                                       local_dof_indices[i],
                                       cell_matrix(i, i));
                    }
                }
        }
    massMatrix.compress(VectorOperation::add);

}

/*template <int dim>
void ellipticBVP<dim>::updateActiveSet(){
//    pcout << "in updateActiveSet\n";
    solution.compress(VectorOperation::add);
    newton_rhs_uncondensed.compress(VectorOperation::add);
    //diag_mass_matrix_vector.compress(VectorOperation::unknown);
    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    vectorType distributed_solution(locally_owned_dofs,
                                        mpi_communicator);
//    pcout << "in updateActiveSet 1\n";
    distributed_solution = solution;
    vectorType lambda(locally_relevant_dofs,
                          mpi_communicator);
//    pcout << "in updateActiveSet 2\n";
    lambda = newton_rhs_uncondensed; //NEED THIS COMPUTED PRIOR TO CALL Done
    vectorType diag_mass_matrix_vector_relevant(
            locally_relevant_dofs, mpi_communicator);
//    pcout << "in updateActiveSet 3\n";
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector; //NEED THIS COMPUTED PRIOR TO CALL Done
//    pcout << "in updateActiveSet 4\n";
    indentation_constraints.reinit(locally_relevant_dofs);
    active_set.clear(); // Add reporting of active set for convergence test and incrementation holding
    Quadrature<dim - 1> face_quadrature(FE.get_unit_face_support_points());
    FEFaceValues<dim>   fe_values_face(FE,
                                       face_quadrature,
                                       update_quadrature_points);
    const unsigned int dofs_per_face   = FE.n_dofs_per_face();
    const unsigned int n_face_q_points = face_quadrature.size();
    std::vector<types::global_dof_index> dof_indices(dofs_per_face);
    for (const auto &cell : dofHandler.active_cell_iterators())
        if (!cell->is_artificial())
            for (const auto &face : cell->face_iterators())
                if (face->at_boundary() && face->boundary_id() == indenterFace)
                {
                    fe_values_face.reinit(cell, face);
                    face->get_dof_indices(dof_indices);
                    for (unsigned int q_point = 0; q_point < n_face_q_points;
                         ++q_point)
                    {
                        const unsigned int component =
                                FE.face_system_to_component_index(q_point).first;
                        const unsigned int index_z = dof_indices[q_point];
                        if ((component == indentDof) && (dof_touched[index_z] == false))
                        {
                            dof_touched[index_z] = true;
                            const Point<dim> this_support_point =
                                    fe_values_face.quadrature_point(q_point);
//                            const double obstacle_value =
//                                    obstacle->value(this_support_point, 2);
                            const double solution_here = solution(index_z);
                            const double undeformed_gap = Obstacle(this_support_point, indentDof, currentPosIndenter);
                                    //obstacle_value - this_support_point(2);
                            const double c = 100.0 * userInputs.elasticStiffness1[0][0]; //e_modulus;

                            if ((lambda(index_z) /
                                 diag_mass_matrix_vector_relevant(index_z) +
                                 c * (solution_here - undeformed_gap) >
                                 0)   //&& !constraints_hanging_nodes.is_constrained(index_z)
                                 )

                            {
                                if (undeformed_gap != 0)
                                {

                                pcout << "dof "<<index_z<< " sol: "<<solution_here<<" gap: "<<undeformed_gap <<" lambda:"
                                      <<lambda(index_z)<<" mass "<<diag_mass_matrix_vector_relevant(index_z) <<" crit(>0?): "
                                      << (lambda(index_z) / diag_mass_matrix_vector_relevant(index_z) +
                                      c * (solution_here - undeformed_gap)) <<"\n";
                                indentation_constraints.add_line(index_z);
                                indentation_constraints.set_inhomogeneity(index_z,
                                                                  undeformed_gap);
                                distributed_solution(index_z) = undeformed_gap;
                                active_set.add_index(index_z);
                                }
                            }
                        }
                    }
                }
}*/

//INDENTATION
template <int dim>
void ellipticBVP<dim>::setIndentation(const Point<dim>& node,
                                      const unsigned int dof,
                                      bool& flag, double& value) {
    flag = flagActiveSet(node);
    if (flag) {
        value = Obstacle(node, dof, currentPosIndenter);// - Obstacle(node, dof, prevPosIndenter);

    }
    else value = 0.0;
    return;
}

template <int dim>
void ellipticBVP<dim>::setIndentation2(const Point<dim>& node,
                                      const unsigned int dof,
                                      bool& flag, double& value, double& criterion) {
    flag = flagActiveSetLambda2(node, criterion);
//    if (flag)
        value = Obstacle(node, dof, currentPosIndenter);// - Obstacle(node, dof, prevPosIndenter);
//    else value = 0.0;
    return;
}

template <int dim>
void ellipticBVP<dim>::setActiveSet(){
    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim> fe_values (FE, QGauss<dim>(userInputs.quadOrder), update_values);
    FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs.quadOrder), update_values);
    Quadrature<dim - 1> face_quadrature(FE.get_unit_face_support_points());
    IndexSet own_dofs = dofHandler.locally_owned_dofs();
    Vector<double> Ulocal(dofs_per_cell);
    const unsigned int n_face_q_points = face_quadrature.size();
    //updateIndentationSet();
    //parallel loop over all elements
    //pcout << "in setActiveSet\n";
    indentation_constraints.clear ();
    indentation_constraints.reinit (locally_relevant_dofs);
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){
            cell->get_dof_indices (local_dof_indices);
            fe_values.reinit (cell);
            for (unsigned int faceID=0; faceID<GeometryInfo<dim>::faces_per_cell; faceID++){ //(const auto &face: cell->face_iterators()){
                if (cell->face(faceID)->at_boundary() && (cell->face(faceID)->boundary_id() == indenterFace )){ //|| userInputs.readExternalMesh && cell->face(faceID)->boundary_id()==indenterFace){ //(face->at_boundary() && face->boundary_id()==indenterFace) {
                    fe_face_values.reinit(cell, faceID); //face);
                    //cell->get_dof_indices (local_dof_indices);
                    //std::cout<<"cell "<<cell<<" boundary_id "<<face->boundary_id()<<"\n";
//                for (unsigned int i = 0; i < dofs_per_cell; i++) {
//                    Ulocal[i] = 0;
//                    Ulocal[i] = solutionWithGhosts[local_dof_indices[i]];
//                }

                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                    if (fe_face_values.shape_value(i, 0) != 0) {
                        const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
                        unsigned int globalDOF = local_dof_indices[i];
                        unsigned int index_z;
                        bool flag = false;
                        double value = 0;
                        node = supportPoints[globalDOF];
                        for (unsigned int i2 = 0; i2 < dofs_per_cell; ++i2) {
                            if (fe_face_values.shape_value(i2, 0) != 0) {
                                const unsigned int dof2 = fe_values.get_fe().system_to_component_index(i2).first;
                                unsigned int globalDOF2 = local_dof_indices[i2];
                                if (supportPoints[globalDOF2] == supportPoints[globalDOF]){
                                    nodeU2(dof2) = solutionWithGhosts[globalDOF2];
                                    index_z = globalDOF2;
                                }
                            }
                        }
                        for (unsigned int i2 = 0; i2 < dim; ++i2) {
                            nodeU(i2) = node(i2) + nodeU2(i2);
                        }
                        //setIndentation(nodeU, dof, flag, value);

                        if ((dof == indentDof) || (roughIndenter == true)) {
                            if (dof_touched[globalDOF] == false) {
                                dof_touched[globalDOF] = true;
                                //std::cout<<"globalDOF "<<globalDOF<<" nodeU "<<nodeU<<"\n";
                                setIndentation(nodeU, dof, flag, value);
                            }
                            else
                                flag = false;
                            if (flag) {
                                if (userInputs.debugIndentation && own_dofs.is_element(globalDOF)) {
                                    std::cout << "dof# " << globalDOF << " value: " << value << " nodeU: " << nodeU
                                              << " soln:" << solutionWithGhosts[globalDOF] << " next? "
                                              << value + nodeU[dof] << "\n";
                                }
                                //active_set.add_line(globalDOF);
                                //active_set.set_inhomogeneity(globalDOF, value);
                                if (dof == indentDof) active_set.add_index(globalDOF);
                                indentation_constraints.add_line(globalDOF);
                                indentation_constraints.set_inhomogeneity(globalDOF, value);
                            }
                        }
                    }
                }
            }
        }
    }
    }
    indentation_constraints.close();

}

template <int dim>
void ellipticBVP<dim>::measureIndentationLoad(){
    vectorType lambda2(locally_owned_dofs,locally_relevant_ghost_dofs,mpi_communicator);
    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    //std::vector<types::global_dof_index> own_dofs(dofHandler.n_locally_owned_dofs());
    IndexSet own_dofs = dofHandler.locally_owned_dofs();
    FEValues<dim> fe_values (FE, QGauss<dim>(userInputs.quadOrder), update_values);
    FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs.quadOrder), update_values);
    Quadrature<dim - 1> face_quadrature(FE.get_unit_face_support_points());
    const unsigned int n_face_q_points = face_quadrature.size();
    solution.compress(VectorOperation::add);
    newton_rhs_uncondensed_inc.compress(VectorOperation::add);
    vectorType distributed_solution(locally_owned_dofs,
                                    mpi_communicator);
    distributed_solution = solution;
    lambda2 = newton_rhs_uncondensed_inc; //NEED THIS COMPUTED PRIOR TO CALL Done
    vectorType diag_mass_matrix_vector_relevant(locally_owned_dofs,
                                                locally_relevant_ghost_dofs, mpi_communicator); // this is needed for petsc to work (not needed in step-42 trilinos)
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;
    active_set.clear();
    indenterLoad = 0.;
    //own_dofs = dofHandler.locally_owned_dofs();
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){
            cell->get_dof_indices (local_dof_indices);
            fe_values.reinit (cell);
            for (unsigned int faceID=0; faceID<GeometryInfo<dim>::faces_per_cell; faceID++){ //(const auto &face: cell->face_iterators()){
                if (cell->face(faceID)->at_boundary() && (cell->face(faceID)->boundary_id()==loadFace)){ //|| userInputs.readExternalMesh //(face->at_boundary() && face->boundary_id()==indenterFace) {
                    fe_face_values.reinit(cell, faceID); //face);
                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        if (fe_face_values.shape_value(i, 0) != 0) {
                            const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
                            unsigned int globalDOF = local_dof_indices[i];
                            //setIndentation(nodeU, dof, flag, value);
                            if (dof == indentDof) {
                                if (!dof_touched[globalDOF]) {
                                    dof_touched[globalDOF] = true;

                                    if (own_dofs.is_element(globalDOF)){
                                        indenterLoad = indenterLoad + lambda2(globalDOF);
                                    }



                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

template <int dim>
void ellipticBVP<dim>::setActiveSet2(){
    vectorType lambda2(locally_owned_dofs,locally_relevant_ghost_dofs,mpi_communicator);
    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    bool active_set_empty = true;
    const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim> fe_values (FE, QGauss<dim>(userInputs.quadOrder), update_values);
    FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs.quadOrder), update_values);
    Quadrature<dim - 1> face_quadrature(FE.get_unit_face_support_points());
    IndexSet own_dofs = dofHandler.locally_owned_dofs();
    Vector<double> Ulocal(dofs_per_cell);
    const unsigned int n_face_q_points = face_quadrature.size();
    solution.compress(VectorOperation::add);
    newton_rhs_uncondensed_inc.compress(VectorOperation::add);
    //diag_mass_matrix_vector.compress(VectorOperation::unknown);
//    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    vectorType distributed_solution(locally_owned_dofs,
                                    mpi_communicator);
    //pcout << "in updateActiveSet2 1\n";
    distributed_solution = solution;

    //pcout << "in updateActiveSet2 2\n";
    lambda2 = newton_rhs_uncondensed_inc; //NEED THIS COMPUTED PRIOR TO CALL Done
    vectorType diag_mass_matrix_vector_relevant(locally_owned_dofs,
                                                locally_relevant_ghost_dofs, mpi_communicator); // this is needed for petsc to work (not needed in step-42 trilinos)
    //pcout << "in updateActiveSet2 3\n";
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;
    active_set.clear();
    indenterLoad = 0.;
    indentation_constraints.clear ();
    indentation_constraints.reinit (locally_relevant_dofs);
    double criterion;
    double solve_for_c;
    double stiffness;
    if (userInputs.continuum_Isotropic)
        stiffness = (userInputs.lame_lambda + 4/3 * userInputs.lame_mu);
    else
        stiffness = userInputs.elasticStiffness1[0][0];
    const double c = userInputs.activeSetCriterionCoefficient*stiffness;
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()){
            cell->get_dof_indices (local_dof_indices);
            fe_values.reinit (cell);
            for (unsigned int faceID=0; faceID<GeometryInfo<dim>::faces_per_cell; faceID++){ //(const auto &face: cell->face_iterators()){
                if (cell->face(faceID)->at_boundary() && (cell->face(faceID)->boundary_id()==indenterFace )){ //|| userInputs.readExternalMesh(face->at_boundary() && face->boundary_id()==indenterFace) {
                    fe_face_values.reinit(cell, faceID); //face);
                    //cell->get_dof_indices (local_dof_indices);
                    //fe_face_values.reinit(cell, indenterFace);
//                for (unsigned int i = 0; i < dofs_per_cell; i++) {
//                    Ulocal[i] = 0;
//                    Ulocal[i] = solutionWithGhosts[local_dof_indices[i]];
//                }
                for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                    if (fe_face_values.shape_value(i, 0) != 0) {
                        const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
                        unsigned int globalDOF = local_dof_indices[i];
                        unsigned int index_z;
                        bool flag = false;
                        double value = 0;
                        node = supportPoints[globalDOF];
                        for (unsigned int i2 = 0; i2 < dofs_per_cell; ++i2) {
                            if (fe_face_values.shape_value(i2, 0) != 0) {
                                const unsigned int dof2 = fe_values.get_fe().system_to_component_index(i2).first;
                                unsigned int globalDOF2 = local_dof_indices[i2];
                                if (supportPoints[globalDOF2] == supportPoints[globalDOF]){
                                    nodeU2(dof2) = solutionWithGhosts[globalDOF2];
                                    if (dof2 == indentDof)
                                        index_z = globalDOF2;
                                }
                            }
                        }
                        for (unsigned int i2 = 0; i2 < dim; ++i2) {
                            nodeU(i2) = node(i2) + nodeU2(i2);
                            if (i2 == indentDof)
                                nodedU(i2) = node(i2);
                            else
                                nodedU(i2) = nodeU(i2);
                        }
                        //setIndentation(nodeU, dof, flag, value);
                        if ((dof == indentDof) || roughIndenter) {

                            criterion = (lambda2(index_z) /
                                         diag_mass_matrix_vector_relevant(index_z) + userInputs.activeSetLambdaTolerance +
                                         c * (solutionWithGhosts[index_z] - Obstacle(nodedU, indentDof, currentPosIndenter)));// + indenterTolerance
//                            solve_for_c = (0.1* diag_mass_matrix_vector_relevant(index_z) - lambda2(index_z)) /
//                                    (diag_mass_matrix_vector_relevant(index_z) *
//                                    (solutionWithGhosts[index_z] - Obstacle(nodedU, indentDof, currentPosIndenter)));

                            if (!dof_touched[globalDOF]) {
                                dof_touched[globalDOF] = true;
                                setIndentation2(nodeU, dof, flag, value, criterion);
                                if (userInputs.debugIndentation && own_dofs.is_element(globalDOF))
                                {
                                    std::vector<unsigned int> arr1={globalDOF};
                                    if (std::includes(debug_set.begin(), debug_set.end(), arr1.begin(), arr1.end())) {
                                        std::cout << "globalDOF: " << globalDOF <<
                                                  " criterion: " << criterion
                                                  << " solution: " << solutionWithGhosts[index_z]
                                                  //<< " c->crit==0.1: " << solve_for_c
                                                  //<< " gap1: " << solutionWithGhosts[index_z] - Obstacle(nodedU, indentDof, currentPosIndenter)
                                                  << " gap: " << Obstacle(nodeU, indentDof, currentPosIndenter)
                                                  //<< " obs(node): " << Obstacle(node, indentDof, currentPosIndenter)
                                                  //<< " obs(nodeU): " << Obstacle(nodeU, indentDof, currentPosIndenter)
                                                  << " lambda/mass: "
                                                  << lambda2(index_z) / diag_mass_matrix_vector_relevant(index_z) <<
                                                  " Node: " << node << "\n";
                                    }
                                }
//                                if ((criterion > -10) && ((dof == indentDof) && (solutionWithGhosts[index_z] < 0))) //currently prints many nodes away from indenter, due to lack of categorical value for missed nodes in Obstacle()
//                                    std::cout<<"dof# "<<globalDOF<<" value: "<<value<<" nodeU: "<<nodeU<<" soln:"<<
//                                             solutionWithGhosts[globalDOF]<<" next? "<<value+nodeU[dof]<<" crit: "<<criterion<<"\n";
//                                if (solutionWithGhosts[globalDOF]<-1e-3)
//                                    std::cout<<"dof# "<<globalDOF<<" value: "<<value<<" nodeU: "<<nodeU<<" soln:"<<solutionWithGhosts[globalDOF]<<" gap "<<Obstacle(nodedU, indentDof, currentPosIndenter)<<"\n";
                            }
                            else
                                flag = false;

                            if (flag) {

                                //active_set.add_line(globalDOF);
                                //active_set.set_inhomogeneity(globalDOF, value);
                                if (dof == indentDof) {
                                    if (active_set_empty) active_set_empty = false;
                                    //indenterLoad = lambda2(index_z);


                                        //indenterLoad = indenterLoad + lambda2(index_z);
                                    active_set.add_index(globalDOF);
                                    //std::cout<<"load of indenter += "<<lambda2(index_z)<<"\n";
                                    //std::cout<<"indenterLoad = "<<indenterLoad<<"\n";
                                }
                                indentation_constraints.add_line(globalDOF);
                                indentation_constraints.set_inhomogeneity(globalDOF, value);
                            }
                        }
                    }
                }
            }
        }
    }
    }
    //pcout << "in updateActiveSet2 6\n";
    indentation_constraints.close();

}


template <int dim>
void ellipticBVP<dim>::setFrozenSet(){
    vectorType lambda2(locally_owned_dofs,locally_relevant_ghost_dofs,mpi_communicator);
    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    bool active_set_empty = true;
    const unsigned int   dofs_per_cell   = FE.dofs_per_cell;
    std::vector<types::global_dof_index> local_dof_indices (dofs_per_cell);
    FEValues<dim> fe_values (FE, QGauss<dim>(userInputs.quadOrder), update_values);
    FEFaceValues<dim> fe_face_values (FE, QGauss<dim-1>(userInputs.quadOrder), update_values);
    Quadrature<dim - 1> face_quadrature(FE.get_unit_face_support_points());
    Vector<double> Ulocal(dofs_per_cell);
    const unsigned int n_face_q_points = face_quadrature.size();
    solution.compress(VectorOperation::add);
    newton_rhs_uncondensed_inc.compress(VectorOperation::add);
    //diag_mass_matrix_vector.compress(VectorOperation::unknown);
//    std::vector<bool> dof_touched(dofHandler.n_dofs(), false);
    vectorType distributed_solution(locally_owned_dofs,
                                    mpi_communicator);
//    pcout << "in updateActiveSet 1\n";
    distributed_solution = solution;

//    pcout << "in updateActiveSet 2\n";
    lambda2 = newton_rhs_uncondensed_inc; //NEED THIS COMPUTED PRIOR TO CALL Done
    vectorType diag_mass_matrix_vector_relevant(locally_owned_dofs,
                                                locally_relevant_ghost_dofs, mpi_communicator); // this is needed for petsc to work (not needed in step-42 trilinos)
//    pcout << "in updateActiveSet 3\n";
    diag_mass_matrix_vector_relevant = diag_mass_matrix_vector;
    active_set.clear();
    indenterLoad = 0.;
    //updateIndentationSet();
    //parallel loop over all elements
    //pcout << "in setIndentationConstraints\n";
    indentation_constraints.clear ();
    indentation_constraints.reinit (locally_relevant_dofs);
    double criterion;
    double stiffness;
    if (userInputs.continuum_Isotropic)
        stiffness = (userInputs.lame_lambda + 4/3 * userInputs.lame_mu);
    else
        stiffness = userInputs.elasticStiffness1[0][0];
    const double c = userInputs.activeSetCriterionCoefficient*stiffness;
    typename DoFHandler<dim>::active_cell_iterator cell = dofHandler.begin_active(), endc = dofHandler.end();
    for (; cell!=endc; ++cell) {
        if (cell->is_locally_owned()) {
            cell->get_dof_indices(local_dof_indices);
            fe_values.reinit(cell);
            for (unsigned int faceID = 0; faceID < GeometryInfo<dim>::faces_per_cell; faceID++) { //(const auto &face: cell->face_iterators()){
                if (cell->face(faceID)->at_boundary() && (cell->face(faceID)->boundary_id()==indenterFace )){ //|| userInputs.readExternalMesh && cell->face(faceID)->boundary_id() ==indenterFace) { //(face->at_boundary() && face->boundary_id()==indenterFace) {
                    fe_face_values.reinit(cell, faceID); //face);

                    for (unsigned int i = 0; i < dofs_per_cell; ++i) {
                        if (fe_face_values.shape_value(i, 0) != 0) {
                            const unsigned int dof = fe_values.get_fe().system_to_component_index(i).first;
                            unsigned int globalDOF = local_dof_indices[i];
                            unsigned int index_z;
                            bool flag = false;
                            double value = 0;
                            node = supportPoints[globalDOF];
                            for (unsigned int i2 = 0; i2 < dofs_per_cell; ++i2) {
                                if (fe_face_values.shape_value(i2, 0) != 0) {
                                    const unsigned int dof2 = fe_values.get_fe().system_to_component_index(i2).first;
                                    unsigned int globalDOF2 = local_dof_indices[i2];
                                    if (supportPoints[globalDOF2] == supportPoints[globalDOF]) {
                                        nodeU2(dof2) = solutionWithGhosts[globalDOF2];
                                        if (dof2 == indentDof)
                                            index_z = globalDOF2;
                                    }
                                }
                            }
                            for (unsigned int i2 = 0; i2 < dim; ++i2) {
                                nodeU(i2) = node(i2) + nodeU2(i2);
                                if (i2 == indentDof)
                                    nodedU(i2) = node(i2);
                                else
                                    nodedU(i2) = nodeU(i2);
                            }
                            //setIndentation(nodeU, dof, flag, value);

                            if ((dof == indentDof) || (roughIndenter == true)) {
                                criterion = (lambda2(index_z) /
                                             diag_mass_matrix_vector_relevant(index_z) +
                                             userInputs.activeSetLambdaTolerance +
                                             c * (solutionWithGhosts[index_z] - Obstacle(nodedU, indentDof,
                                                                                         currentPosIndenter)));// + indenterTolerance
                                if (dof_touched[globalDOF] == false) {
                                    dof_touched[globalDOF] = true;
                                    std::vector<unsigned int> arr1 = {globalDOF};
                                    if (std::includes(frozen_set.begin(), frozen_set.end(), arr1.begin(), arr1.end())) {
                                        flag = true;
                                        double arbitrary = 100;
                                        setIndentation2(nodeU, dof, flag, value, arbitrary);
                                    }
                                } else
                                    flag = false;
                                if (flag) {

                                    //active_set.add_line(globalDOF);
                                    //active_set.set_inhomogeneity(globalDOF, value);
                                    if (dof == indentDof) {
                                        if (active_set_empty) active_set_empty = false;
                                        active_set.add_index(globalDOF);
                                    }
                                    indentation_constraints.add_line(globalDOF);
                                    indentation_constraints.set_inhomogeneity(globalDOF, value);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    indentation_constraints.close();

}

//methods to apply dirichlet BC's
template <int dim>
void ellipticBVP<dim>::setIndentationConstraints(){
    if (currentIteration == 0) {
        setActiveSet();
        debug_set = active_set;
        //frozen_set = active_set;
    } else {
        //pcout << "debug set indentation constraints 1, freeze = "<< freeze_out_iterations <<"\n";
        if (freeze_out_iterations > 0)
        {
            pcout << "Active Set Frozen to encourage convergence"<<freeze_out_iterations;
            freeze_out_iterations-=1;
            setFrozenSet();
        }
        else
        {
            setActiveSet2();
            frozen_set = active_set;
        }
        //pcout << "debug set indentation constraints 2\n";
        indenterLoad = Utilities::MPI::sum(indenterLoad, mpi_communicator); // Zero before the measurement elsewhere
        //updateActiveSet();
    }
    active_set_size = Utilities::MPI::sum((active_set & locally_owned_dofs).n_elements(),
                                          mpi_communicator);

    pcout << "         Size of active set: "
          << active_set_size
          << std::endl;

}

#include "../../include/ellipticBVP_template_instantiations.h"
