//solve method for ellipticBVP class
#include "../../include/ellipticBVP.h"

//loop over increments and solve each increment
template <int dim>
void ellipticBVP<dim>::solve(){
  pcout << "begin solve...\n\n";
  bool success;
  //load increments
  unsigned int successiveIncs=0;

  if(userInputs.enableAdaptiveTimeStepping){
    for (;totalLoadFactor<totalIncrements;){
      ++currentIncrement;
      loadFactorSetByModel=std::min(loadFactorSetByModel, totalIncrements-totalLoadFactor);
      pcout << "\nincrement: "  << currentIncrement << std::endl;
      char buffer[100];
      sprintf(buffer, "current load factor: %12.6e\ntotal load factor:   %12.6e\n", loadFactorSetByModel, totalLoadFactor);
      pcout << buffer;

      //call updateBeforeIncrement, if any
      updateBeforeIncrement();

      if (!userInputs.flagTaylorModel){
        //solve time increment
        success=solveNonLinearSystem();
      }

      //call updateAfterIncrement, if any
      if ((success)||(userInputs.flagTaylorModel)){
        updateAfterIncrement();

        //update totalLoadFactor
        totalLoadFactor+=loadFactorSetByModel;

        //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
        successiveIncs++;

        if (successiveIncs>=userInputs.succesiveIncForIncreasingTimeStep){
          loadFactorSetByModel*=userInputs.adaptiveLoadIncreaseFactor;
          char buffer1[100];
          sprintf(buffer1, "current increment increased. Restarting increment with loadFactorSetByModel: %12.6e\n", loadFactorSetByModel);
          pcout << buffer1;
        }
        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        computing_timer.enter_section("postprocess");
        #else
	    computing_timer.enter_subsection("postprocess");
        #endif

        if (currentIncrement%userInputs.skipOutputSteps==0)
        if (userInputs.writeOutput) output();

        #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
        computing_timer.exit_section("postprocess");
        #else
	    computing_timer.leave_subsection("postprocess");
        #endif

      }
      else{
        successiveIncs=0;
      }
    }
    char buffer[100];
    sprintf(buffer, "\nfinal load factor  : %12.6e\n", totalLoadFactor);
    pcout << buffer;
  }
  else
    for (;currentIncrement<totalIncrements; ++currentIncrement){
    pcout << "\nincrement: "  << currentIncrement << std::endl;
    if (userInputs.enableIndentationBCs){
        ellipticBVP<dim>::updateBeforeIncrement();
        if (!userInputs.continuum_Isotropic)
            updateBeforeIncrement();
    }
    else{
        updateBeforeIncrement();
    }
    //call updateBeforeIncrement, if any


    if (!userInputs.flagTaylorModel){
      //solve time increment
      success=solveNonLinearSystem();
    }

    //call updateAfterIncrement, if any
    if ((success)||(userInputs.flagTaylorModel)){
      updateAfterIncrement();

      //update totalLoadFactor
      totalLoadFactor+=loadFactorSetByModel;

      //increase loadFactorSetByModel, if succesiveIncForIncreasingTimeStep satisfied.
      successiveIncs++;
      //output results to file
      //
      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer.enter_section("postprocess");
      #else
      computing_timer.enter_subsection("postprocess");
      #endif

      //////////////////////TabularOutput Start///////////////
      std::vector<unsigned int> tabularTimeInputIncInt;
      std::vector<double> tabularTimeInputInc;
      if (userInputs.tabularOutput){

        tabularTimeInputInc=userInputs.tabularTimeOutput;
        for(unsigned int i=0;i<userInputs.tabularTimeOutput.size();i++){
          tabularTimeInputInc[i]=tabularTimeInputInc[i]/delT;
        }
        tabularTimeInputIncInt.resize(userInputs.tabularTimeOutput.size(),0);
        ///Converting to an integer always rounds down, even if the fraction part is 0.99999999.
        //Hence, I add 0.1 to make sure we always get the correct integer.
        for(unsigned int i=0;i<userInputs.tabularTimeOutput.size();i++){
          tabularTimeInputIncInt[i]=int(tabularTimeInputInc[i]+0.1);
        }
      }
      //////////////////////TabularOutput Finish///////////////
      if (((!userInputs.tabularOutput)&&((currentIncrement+1)%userInputs.skipOutputSteps==0))||((userInputs.tabularOutput)&& (std::count(tabularTimeInputIncInt.begin(), tabularTimeInputIncInt.end(), (currentIncrement+1))==1))){
        if (userInputs.writeOutput) output();
      }

      #if ((DEAL_II_VERSION_MAJOR < 9)||((DEAL_II_VERSION_MINOR < 3)&&(DEAL_II_VERSION_MAJOR==9)))
      computing_timer.exit_section("postprocess");
      #else
      computing_timer.leave_subsection("postprocess");
      #endif
    }
    else{
      successiveIncs=0;
    }
  }
}
#include "../../include/ellipticBVP_template_instantiations.h"
