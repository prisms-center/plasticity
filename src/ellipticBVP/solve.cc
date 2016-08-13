//solve method for ellipticBVP class

#ifndef SOLVE_ELLIPTICBVP_H
#define SOLVE_ELLIPTICBVP_H
//this source file is temporarily treated as a header file (hence
//#ifndef's) till library packaging scheme is finalized

//loop over increments and solve each increment 
template <int dim>
void ellipticBVP<dim>::solve(){
  pcout << "begin solve...\n\n";

  //increments
  for (;currentIncrement<totalIncrements; ++currentIncrement){
    pcout << "\nincrement: " 
	  << currentIncrement 
	  << std::endl;
    //call updateBeforeIncrement, if any
    updateBeforeIncrement();

    //solve time increment
    bool success=solveNonLinearSystem();
    
    //call updateAfterIncrement, if any
    if (success){
      updateAfterIncrement();

      //output results to file
      computing_timer.enter_section("postprocess");
      if (writeOutput) output();
      computing_timer.exit_section("postprocess");
    }
  }
}

#endif
