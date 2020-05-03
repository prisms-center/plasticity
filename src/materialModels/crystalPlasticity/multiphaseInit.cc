#include "../../../include/crystalPlasticity.h"

template <int dim>
void crystalPlasticity<dim>::multiphaseInit(unsigned int cellID,
  unsigned int quadPtID) {
  elasticStiffnessMatrix.reinit(2 * dim, 2 * dim);
  if (!this->userInputs.enableMultiphase){
    enableTwinning=this->userInputs.enableTwinning1;
    twinThresholdFraction=this->userInputs.twinThresholdFraction1;
    twinSaturationFactor=this->userInputs.twinSaturationFactor1;
    n_Tslip_systems=n_Tslip_systems_SinglePhase;
    n_slip_systems=n_slip_systems_SinglePhase;
    n_twin_systems=n_twin_systems_SinglePhase;
    elasticStiffnessMatrix=Dmat_SinglePhase;
    C_1.reinit(n_Tslip_systems);
    C_2.reinit(n_Tslip_systems);
    if(this->userInputs.enableKinematicHardening1){
      for(unsigned int i=0;i<n_slip_systems;i++){
        C_1(i)=this->userInputs.C_1_slip1[i];
        C_2(i)=this->userInputs.C_2_slip1[i];
      }
      for(unsigned int i=0;i<n_twin_systems;i++){
        C_1(n_slip_systems+i)=this->userInputs.C_1_twin1[i];
        C_2(n_slip_systems+i)=this->userInputs.C_2_twin1[i];
      }
    }
    else{
      for(unsigned int i=0;i<n_Tslip_systems;i++){
        C_1(i)=0;
        C_2(i)=0;
      }
    }

    q.reinit(n_Tslip_systems,n_Tslip_systems);
    q=q_phase1;


    twinShear=this->userInputs.twinShear1;
    initialHardeningModulus.reinit(n_slip_systems);
    saturationStress.reinit(n_slip_systems);
    powerLawExponent.reinit(n_slip_systems);
    initialHardeningModulusTwin.reinit(n_twin_systems);
    saturationStressTwin.reinit(n_twin_systems);
    powerLawExponentTwin.reinit(n_twin_systems);

    for(unsigned int i=0;i<n_slip_systems;i++){
      initialHardeningModulus[i]=this->userInputs.initialHardeningModulus1[i];
      saturationStress[i]=this->userInputs.saturationStress1[i];
      powerLawExponent[i]=this->userInputs.powerLawExponent1[i];
    }
    for(unsigned int i=0;i<n_twin_systems;i++){
      initialHardeningModulusTwin[i]=this->userInputs.initialHardeningModulusTwin1[i];
      saturationStressTwin[i]=this->userInputs.saturationStressTwin1[i];
      powerLawExponentTwin[i]=this->userInputs.powerLawExponentTwin1[i];
    }

    m_alpha=m_alpha_SinglePhase;
    n_alpha=n_alpha_SinglePhase;

    if (this->userInputs.enableUserMaterialModel){
      UserMatConstants.reinit(this->userInputs.numberofUserMatConstants1);
      for(unsigned int i=0;i<this->userInputs.numberofUserMatConstants1;i++){
        UserMatConstants(i)=this->userInputs.UserMatConstants1[i];
      }
    }

  }
  if (this->userInputs.enableMultiphase){
    phaseMaterial=phase[cellID][quadPtID];
    if (phaseMaterial==1){
      enableTwinning=this->userInputs.enableTwinning1;
      twinThresholdFraction=this->userInputs.twinThresholdFraction1;
      n_Tslip_systems=n_Tslip_systems_MultiPhase[0];
      n_slip_systems=n_slip_systems_MultiPhase[0];
      n_twin_systems=n_twin_systems_MultiPhase[0];
      m_alpha.reinit(n_Tslip_systems,3);
      n_alpha.reinit(n_Tslip_systems,3);
      for(unsigned int i=0;i<n_Tslip_systems;i++){
        for(unsigned int j=0;j<dim;j++){
          m_alpha[i][j]=m_alpha_MultiPhase[i][j];
          n_alpha[i][j]=n_alpha_MultiPhase[i][j];
        }
      }



      for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
          elasticStiffnessMatrix[i][j] = Dmat_MultiPhase[i][j];
        }
      }

      C_1.reinit(n_Tslip_systems);
      C_2.reinit(n_Tslip_systems);
      if(this->userInputs.enableKinematicHardening1){
        for(unsigned int i=0;i<n_slip_systems;i++){
          C_1(i)=this->userInputs.C_1_slip1[i];
          C_2(i)=this->userInputs.C_2_slip1[i];
        }
        for(unsigned int i=0;i<n_twin_systems;i++){
          C_1(n_slip_systems+i)=this->userInputs.C_1_twin1[i];
          C_2(n_slip_systems+i)=this->userInputs.C_2_twin1[i];
        }
      }
      else{
        for(unsigned int i=0;i<n_Tslip_systems;i++){
          C_1(i)=0;
          C_2(i)=0;
        }
      }

      q.reinit(n_Tslip_systems,n_Tslip_systems);
      q=q_phase1;


      twinShear=this->userInputs.twinShear1;
      initialHardeningModulus.reinit(n_slip_systems);
      saturationStress.reinit(n_slip_systems);
      powerLawExponent.reinit(n_slip_systems);
      initialHardeningModulusTwin.reinit(n_twin_systems);
      saturationStressTwin.reinit(n_twin_systems);
      powerLawExponentTwin.reinit(n_twin_systems);
      for(unsigned int i=0;i<n_slip_systems;i++){
        initialHardeningModulus[i]=this->userInputs.initialHardeningModulus1[i];
        saturationStress[i]=this->userInputs.saturationStress1[i];
        powerLawExponent[i]=this->userInputs.powerLawExponent1[i];
      }
      for(unsigned int i=0;i<n_twin_systems;i++){
        initialHardeningModulusTwin[i]=this->userInputs.initialHardeningModulusTwin1[i];
        saturationStressTwin[i]=this->userInputs.saturationStressTwin1[i];
        powerLawExponentTwin[i]=this->userInputs.powerLawExponentTwin1[i];
      }

      if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel1)){
        UserMatConstants.reinit(this->userInputs.numberofUserMatConstants1);
        for(unsigned int i=0;i<this->userInputs.numberofUserMatConstants1;i++){
          UserMatConstants(i)=this->userInputs.UserMatConstants1[i];
        }
      }

    }
    else if (phaseMaterial==2){
      enableTwinning=this->userInputs.enableTwinning2;
      twinThresholdFraction=this->userInputs.twinThresholdFraction2;
      n_Tslip_systems=n_Tslip_systems_MultiPhase[1];
      n_slip_systems=n_slip_systems_MultiPhase[1];
      n_twin_systems=n_twin_systems_MultiPhase[1];

      m_alpha.reinit(n_Tslip_systems,3);
      n_alpha.reinit(n_Tslip_systems,3);

      for(unsigned int i=n_Tslip_systems_MultiPhase[0];i<n_Tslip_systems+n_Tslip_systems_MultiPhase[0];i++){
        for(unsigned int j=0;j<dim;j++){
          m_alpha[i-n_Tslip_systems_MultiPhase[0]][j]=m_alpha_MultiPhase[i][j];
          n_alpha[i-n_Tslip_systems_MultiPhase[0]][j]=n_alpha_MultiPhase[i][j];
        }
      }

      for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
          elasticStiffnessMatrix[i][j] = Dmat_MultiPhase[i+6][j];
        }
      }

      C_1.reinit(n_Tslip_systems);
      C_2.reinit(n_Tslip_systems);
      if(this->userInputs.enableKinematicHardening2){
        for(unsigned int i=0;i<n_slip_systems;i++){
          C_1(i)=this->userInputs.C_1_slip2[i];
          C_2(i)=this->userInputs.C_2_slip2[i];
        }
        for(unsigned int i=0;i<n_twin_systems;i++){
          C_1(n_slip_systems+i)=this->userInputs.C_1_twin2[i];
          C_2(n_slip_systems+i)=this->userInputs.C_2_twin2[i];
        }
      }
      else{
        for(unsigned int i=0;i<n_Tslip_systems;i++){
          C_1(i)=0;
          C_2(i)=0;
        }
      }

      q.reinit(n_Tslip_systems,n_Tslip_systems);
      q=q_phase2;


      twinShear=this->userInputs.twinShear2;
      initialHardeningModulus.reinit(n_slip_systems);
      saturationStress.reinit(n_slip_systems);
      powerLawExponent.reinit(n_slip_systems);
      initialHardeningModulusTwin.reinit(n_twin_systems);
      saturationStressTwin.reinit(n_twin_systems);
      powerLawExponentTwin.reinit(n_twin_systems);
      for(unsigned int i=0;i<n_slip_systems;i++){
        initialHardeningModulus[i]=this->userInputs.initialHardeningModulus2[i];
        saturationStress[i]=this->userInputs.saturationStress2[i];
        powerLawExponent[i]=this->userInputs.powerLawExponent2[i];
      }
      for(unsigned int i=0;i<n_twin_systems;i++){
        initialHardeningModulusTwin[i]=this->userInputs.initialHardeningModulusTwin2[i];
        saturationStressTwin[i]=this->userInputs.saturationStressTwin2[i];
        powerLawExponentTwin[i]=this->userInputs.powerLawExponentTwin2[i];
      }

      if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel2)){
        UserMatConstants.reinit(this->userInputs.numberofUserMatConstants2);
        for(unsigned int i=0;i<this->userInputs.numberofUserMatConstants2;i++){
          UserMatConstants(i)=this->userInputs.UserMatConstants2[i];
        }
      }


    }
    else if (phaseMaterial==3){
      enableTwinning=this->userInputs.enableTwinning3;
      twinThresholdFraction=this->userInputs.twinThresholdFraction3;
      n_Tslip_systems=n_Tslip_systems_MultiPhase[2];
      n_slip_systems=n_slip_systems_MultiPhase[2];
      n_twin_systems=n_twin_systems_MultiPhase[2];

      m_alpha.reinit(n_Tslip_systems,3);
      n_alpha.reinit(n_Tslip_systems,3);

      for(unsigned int i=n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1];i<n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems;i++){
        for(unsigned int j=0;j<dim;j++){
          m_alpha[i-n_Tslip_systems_MultiPhase[0]-n_Tslip_systems_MultiPhase[1]][j]=m_alpha_MultiPhase[i][j];
          n_alpha[i-n_Tslip_systems_MultiPhase[0]-n_Tslip_systems_MultiPhase[1]][j]=n_alpha_MultiPhase[i][j];
        }
      }

      for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
          elasticStiffnessMatrix[i][j] = Dmat_MultiPhase[i+6*2][j];
        }
      }

      C_1.reinit(n_Tslip_systems);
      C_2.reinit(n_Tslip_systems);
      if(this->userInputs.enableKinematicHardening3){
        for(unsigned int i=0;i<n_slip_systems;i++){
          C_1(i)=this->userInputs.C_1_slip3[i];
          C_2(i)=this->userInputs.C_2_slip3[i];
        }
        for(unsigned int i=0;i<n_twin_systems;i++){
          C_1(n_slip_systems+i)=this->userInputs.C_1_twin3[i];
          C_2(n_slip_systems+i)=this->userInputs.C_2_twin3[i];
        }
      }
      else{
        for(unsigned int i=0;i<n_Tslip_systems;i++){
          C_1(i)=0;
          C_2(i)=0;
        }
      }

      q.reinit(n_Tslip_systems,n_Tslip_systems);
      q=q_phase3;


      twinShear=this->userInputs.twinShear3;
      initialHardeningModulus.reinit(n_slip_systems);
      saturationStress.reinit(n_slip_systems);
      powerLawExponent.reinit(n_slip_systems);
      initialHardeningModulusTwin.reinit(n_twin_systems);
      saturationStressTwin.reinit(n_twin_systems);
      powerLawExponentTwin.reinit(n_twin_systems);
      for(unsigned int i=0;i<n_slip_systems;i++){
        initialHardeningModulus[i]=this->userInputs.initialHardeningModulus3[i];
        saturationStress[i]=this->userInputs.saturationStress3[i];
        powerLawExponent[i]=this->userInputs.powerLawExponent3[i];
      }
      for(unsigned int i=0;i<n_twin_systems;i++){
        initialHardeningModulusTwin[i]=this->userInputs.initialHardeningModulusTwin3[i];
        saturationStressTwin[i]=this->userInputs.saturationStressTwin3[i];
        powerLawExponentTwin[i]=this->userInputs.powerLawExponentTwin3[i];
      }

      if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel3)){
        UserMatConstants.reinit(this->userInputs.numberofUserMatConstants3);
        for(unsigned int i=0;i<this->userInputs.numberofUserMatConstants3;i++){
          UserMatConstants(i)=this->userInputs.UserMatConstants3[i];
        }
      }


    }
    else if (phaseMaterial==4){
      enableTwinning=this->userInputs.enableTwinning4;
      twinThresholdFraction=this->userInputs.twinThresholdFraction4;
      n_Tslip_systems=n_Tslip_systems_MultiPhase[3];
      n_slip_systems=n_slip_systems_MultiPhase[3];
      n_twin_systems=n_twin_systems_MultiPhase[3];

      m_alpha.reinit(n_Tslip_systems,3);
      n_alpha.reinit(n_Tslip_systems,3);

      for(unsigned int i=n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[2];i<n_Tslip_systems_MultiPhase[0]+n_Tslip_systems_MultiPhase[1]+n_Tslip_systems_MultiPhase[2]+n_Tslip_systems;i++){
        for(unsigned int j=0;j<dim;j++){
          m_alpha[i-n_Tslip_systems_MultiPhase[0]-n_Tslip_systems_MultiPhase[1]-n_Tslip_systems_MultiPhase[2]][j]=m_alpha_MultiPhase[i][j];
          n_alpha[i-n_Tslip_systems_MultiPhase[0]-n_Tslip_systems_MultiPhase[1]-n_Tslip_systems_MultiPhase[2]][j]=n_alpha_MultiPhase[i][j];
        }
      }

      for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
          elasticStiffnessMatrix[i][j] = Dmat_MultiPhase[i+6*3][j];
        }
      }

      C_1.reinit(n_Tslip_systems);
      C_2.reinit(n_Tslip_systems);
      if(this->userInputs.enableKinematicHardening4){
        for(unsigned int i=0;i<n_slip_systems;i++){
          C_1(i)=this->userInputs.C_1_slip4[i];
          C_2(i)=this->userInputs.C_2_slip4[i];
        }
        for(unsigned int i=0;i<n_twin_systems;i++){
          C_1(n_slip_systems+i)=this->userInputs.C_1_twin4[i];
          C_2(n_slip_systems+i)=this->userInputs.C_2_twin4[i];
        }
      }
      else{
        for(unsigned int i=0;i<n_Tslip_systems;i++){
          C_1(i)=0;
          C_2(i)=0;
        }
      }

      q.reinit(n_Tslip_systems,n_Tslip_systems);
      q=q_phase4;


      twinShear=this->userInputs.twinShear4;
      initialHardeningModulus.reinit(n_slip_systems);
      saturationStress.reinit(n_slip_systems);
      powerLawExponent.reinit(n_slip_systems);
      initialHardeningModulusTwin.reinit(n_twin_systems);
      saturationStressTwin.reinit(n_twin_systems);
      powerLawExponentTwin.reinit(n_twin_systems);
      for(unsigned int i=0;i<n_slip_systems;i++){
        initialHardeningModulus[i]=this->userInputs.initialHardeningModulus4[i];
        saturationStress[i]=this->userInputs.saturationStress4[i];
        powerLawExponent[i]=this->userInputs.powerLawExponent4[i];
      }
      for(unsigned int i=0;i<n_twin_systems;i++){
        initialHardeningModulusTwin[i]=this->userInputs.initialHardeningModulusTwin4[i];
        saturationStressTwin[i]=this->userInputs.saturationStressTwin4[i];
        powerLawExponentTwin[i]=this->userInputs.powerLawExponentTwin4[i];
      }

      if ((this->userInputs.enableUserMaterialModel)&&(this->userInputs.enableUserMaterialModel4)){
        UserMatConstants.reinit(this->userInputs.numberofUserMatConstants4);
        for(unsigned int i=0;i<this->userInputs.numberofUserMatConstants4;i++){
          UserMatConstants(i)=this->userInputs.UserMatConstants4[i];
        }
      }


    }
  }




}

#include "../../../include/crystalPlasticity_template_instantiations.h"
