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
    pcout << "increment: " 
	  << currentIncrement 
	  << std::endl;
    
    //solve time increment
    computing_timer.enter_section("solve");
    solveIncrement();
    computing_timer.exit_section("solve");

    //output results to file
    computing_timer.enter_section("postprocess");
    if (writeOutput) output();
    computing_timer.exit_section("postprocess");
  }
}

#endif
