//methods to apply initial conditions 

#ifndef INITIALCONDITIONS_H
#define INITIALCONDITIONS_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//methods to apply initial conditions
template <int dim>
void ellipticBVP<dim>::applyInitialConditions(){
  //pcout << "applying the default zero initial condition\n";
  //default method to apply zero initial conditions on all fields
  VectorTools::interpolate (dofHandler,
			    ZeroFunction<dim>(dim),
			    solution);
}
#endif
