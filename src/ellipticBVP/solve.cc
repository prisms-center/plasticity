//solve method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//loop over increments and solve each increment
template <int dim>
void ellipticBVP<dim>::solve(){
  pcout << "begin solve...\n\n";

  //load increments
  unsigned int successiveIncs=0;
#ifdef enableAdaptiveTimeStepping
#if enableAdaptiveTimeStepping==true
  for (;totalLoadFactor<totalNumIncrements;){
    ++currentIncrement;
    loadFactorSetByModel=std::min(loadFactorSetByModel, totalNumIncrements-totalLoadFactor);
#else
  for (;currentIncrement<totalIncrements; ++currentIncrement){
#endif
#else
  for (;currentIncrement<totalIncrements; ++currentIncrement){
#endif
    pcout << "\nincrement: "  << currentIncrement << std::endl;
#ifdef enableAdaptiveTimeStepping
#if enableAdaptiveTimeStepping==true
    char buffer[100];
    sprintf(buffer, "current load factor: %12.6e\ntotal load factor:   %12.6e\n", loadFactorSetByModel, totalLoadFactor);
    pcout << buffer;
#endif
#endif

    //call updateBeforeIncrement, if any
    updateBeforeIncrement();

    //solve time increment
    bool success=solveNonLinearSystem();

    //call updateAfterIncrement, if any
    if (success){
      updateAfterIncrement();

      //update totalLoadFactor
      totalLoadFactor+=loadFactorSetByModel;

      //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
      successiveIncs++;
#ifdef enableAdaptiveTimeStepping
#if enableAdaptiveTimeStepping==true
#ifdef succesiveIncForIncreasingTimeStep
      if (successiveIncs>=succesiveIncForIncreasingTimeStep){
#ifdef adaptiveLoadIncreaseFactor
	loadFactorSetByModel*=adaptiveLoadIncreaseFactor;
#else
	loadFactorSetByModel*=2.0;
#endif
	char buffer1[100];
	sprintf(buffer1, "current increment increased. Restarting increment with loadFactorSetByModel: %12.6e\n", loadFactorSetByModel);
	pcout << buffer1;
      }
#endif
#endif
#endif
      //output results to file
      computing_timer.enter_section("postprocess");
#ifndef skipOutputSteps
      unsigned int skipSteps=1;
#elif skipOutputSteps<=0
      unsigned int skipSteps=1;
#else
      unsigned int skipSteps=skipOutputSteps;
#endif
      if (currentIncrement%skipSteps==0){
#ifdef writeOutput
      if (writeOutput) output();
#else
      output();
#endif
      }
      computing_timer.exit_section("postprocess");
    }
    else{
      successiveIncs=0;
    }
  }
#ifdef enableAdaptiveTimeStepping
#if enableAdaptiveTimeStepping==true
  char buffer[100];
  sprintf(buffer, "\nfinal load factor  : %12.6e\n", totalLoadFactor);
  pcout << buffer;
#endif
#endif
}
#include "../../include/ellipticBVP_template_instantiations.h"
