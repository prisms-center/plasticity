//methods to allow for pre/post increment level updates
#include "../../include/ellipticBVP.h"

//method called before each increment
template <int dim>
void ellipticBVP<dim>::updateBeforeIncrement(){
  //default method does nothing
}

//method called after each increment
template <int dim>
void ellipticBVP<dim>::updateAfterIncrement(){
  //default method does nothing
}
