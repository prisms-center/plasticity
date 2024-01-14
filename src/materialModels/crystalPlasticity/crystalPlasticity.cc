#include "../../../include/crystalPlasticity.h"
//constructor
template <int dim>
crystalPlasticity<dim>::crystalPlasticity(userInputParameters & _userInputs):
ellipticBVP<dim>(_userInputs),
F(dim,dim),
F_tau(dim,dim),
FP_tau(dim,dim),
FE_tau(dim,dim),
T(dim,dim),
P(dim,dim)
{
    initCalled = false;

    //post processing
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_stress");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Eqv_strain");
    ellipticBVP<dim>::postprocessed_solution_names.push_back("Twin");
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

    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("meshGrain_ID");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var1");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var2");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var3");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var4");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var5");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var6");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var7");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var8");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var9");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var10");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var11");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var12");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var13");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var14");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var15");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var16");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var17");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var18");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var19");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var20");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var21");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var22");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var23");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var24");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var25");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var26");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var27");
    ellipticBVP<dim>::postprocessedFieldsAtCellCenters_solution_names.push_back("outputCellCenters_Var28");
    ellipticBVP<dim>::numPostProcessedFieldsAtCellCenters=29; //grainID

}

#include "../../../include/crystalPlasticity_template_instantiations.h"
