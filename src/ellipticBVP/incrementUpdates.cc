//methods to allow for pre/post increment level updates
#include "../../include/ellipticBVP.h"

//method called before each increment
template <int dim>
void ellipticBVP<dim>::updateBeforeIncrement(){
  //default method does nothing
  //Overwritten in crystal plasticity
  //pcout << "eBVP::updateBeforeIncrement\n";
  if (userInputs.enableIndentationBCs){
  //    pcout << "enableIndentationBCs\n";
      updateIndentPos();
      //newton_rhs_uncondensed += n
  }
}

//method called after each increment
template <int dim>
void ellipticBVP<dim>::updateAfterIncrementBase(){
    if (userInputs.enableIndentationBCs){
        measureIndentationLoad();
        indenterLoad = Utilities::MPI::sum(indenterLoad, mpi_communicator);
        pcout << "         Indenter Load: "
              << indenterLoad
              << std::endl;
    }
  //default method does nothing
}

template <int dim>
void ellipticBVP<dim>::updateAfterIncrement(){

    //default method does nothing
}


#include "../../include/ellipticBVP_template_instantiations.h"
