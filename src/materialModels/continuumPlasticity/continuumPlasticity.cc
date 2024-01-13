#include "../../../include/continuumPlasticity.h"
//constructor
template <int dim>
continuumPlasticity<dim>::continuumPlasticity(userInputParameters & _userInputs)
:
ellipticBVP<dim>(_userInputs),
  enhStrain(this->FE,this->pcout, _userInputs)
{
  //initialize "initCalled"
  initCalled = false;
  //initCalled = false;

  //post processing
  ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("alpha");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var1");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var2");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var3");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var4");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var5");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var6");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var7");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var8");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var9");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var10");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var11");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var12");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var13");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var14");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var15");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var16");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var17");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var18");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var19");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var20");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var21");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var22");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var23");
  ellipticBVP<dim>::postprocessed_solution_names.push_back("output_Var24");
  ellipticBVP<dim>::numPostProcessedFields=27;
  ellipticBVP<dim>::numPostProcessedFieldsAtCellCenters=1; //grainID
}

#include "../../../include/continuumPlasticity_template_instantiations.h"
