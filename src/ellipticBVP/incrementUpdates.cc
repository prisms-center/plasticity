//methods to allow for pre/post increment level updates
#include "../../include/ellipticBVP.h"

#ifndef UPDATEINCREMENT_H
#define UPDATEINCREMENT_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

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

#endif
