#include "../../../include/crystalPlasticity.h"
#include <iostream>
#include <fstream>

template <int dim>
void crystalPlasticity<dim>::init(unsigned int num_quad_points)
{

    //call loadOrientations to load material orientations
    loadOrientations();

    local_Truestrain.reinit(dim);
    local_Truestrain = 0.0;
    global_Truestrain.reinit(dim);
    global_Truestrain = 0.0;
    local_strain.reinit(dim,dim);
    local_stress.reinit(dim,dim);
    global_strain.reinit(dim,dim);
    global_stress.reinit(dim,dim);
    local_strain=0.0;
    local_stress=0.0;
    global_strain=0.0;
    global_stress=0.0;

    unsigned int num_local_cells = this->triangulation.n_locally_owned_active_cells();
    F.reinit(dim, dim);

    n_slip_systems=this->userInputs.numSlipSystems;
    if(this->userInputs.enableTwinning){
      n_slip_systems+=this->userInputs.numTwinSystems;
      n_twin_systems=this->userInputs.numTwinSystems;
    }
    else{
      n_slip_systems+=1;
      n_twin_systems=1;
    }

    n_alpha.reinit(n_slip_systems,3);
      m_alpha.reinit(n_slip_systems,3);
      std::string line;

      //open data file to read slip normals
      std::ifstream slipNormalsDataFile(this->userInputs.slipNormalsFile);
      //read data
      unsigned int id=0;
      if (slipNormalsDataFile.is_open()){
	while (getline (slipNormalsDataFile,line) && id<this->userInputs.numSlipSystems){
	  std::stringstream ss(line);
	  ss >> n_alpha[id][0];
	  ss >> n_alpha[id][1];
	  ss >> n_alpha[id][2];
	  id=id+1;
	}
      }
      else{
	std::cout << "Unable to open slip normals file \n";
	exit(1);
      }

      //open data file to read slip directions
     std::ifstream slipDirectionsDataFile(this->userInputs.slipDirectionsFile);
      //read data
      id=0;
      if (slipDirectionsDataFile.is_open()){
	while (getline (slipDirectionsDataFile,line)&& id<this->userInputs.numSlipSystems){
	  std::stringstream ss(line);
	  ss >> m_alpha[id][0];
	  ss >> m_alpha[id][1];
	  ss >> m_alpha[id][2];
	  id=id+1;
	}
      }
      else{
	std::cout << "Unable to open slip directions file \n";
	exit(1);
      }

      if(this->userInputs.enableTwinning){
  	     //open data file to read twin normals
  	     std::ifstream twinNormalsDataFile(this->userInputs.twinNormalsFile);
         //read data
         id=this->userInputs.numSlipSystems;
         if (twinNormalsDataFile.is_open())
        	 while (getline (twinNormalsDataFile,line) && id<n_slip_systems){
          	  std::stringstream ss(line);
          	  ss >> n_alpha[id][0];
          	  ss >> n_alpha[id][1];
          	  ss >> n_alpha[id][2];
          	  id=id+1;
    	     }
         else{
  	        std::cout << "Unable to open twin normals file\n";
  	        exit(1);
         }

         //open data file to read twin directions
  	     std::ifstream twinDirectionsDataFile(this->userInputs.twinDirectionsFile);
         //read data
         id=this->userInputs.numSlipSystems;
         if (twinDirectionsDataFile.is_open())
  	        //read data
          	while (getline (twinDirectionsDataFile,line)&& id<n_slip_systems){
          	  std::stringstream ss(line);
          	  ss >> m_alpha[id][0];
          	  ss >> m_alpha[id][1];
          	  ss >> m_alpha[id][2];
          	  id=id+1;
         }
         else{
        	  std::cout << "Unable to open twin directions file \n";
        	  exit(1);
         }
      }
      else{
         n_alpha[id][0]=1;
         n_alpha[id][1]=0;
         n_alpha[id][2]=0;
         m_alpha[id][0]=0;
         m_alpha[id][1]=1;
         m_alpha[id][2]=0;
      }

    //q is a parameter in the hardening model
    q.reinit(n_slip_systems,n_slip_systems);
    for(unsigned int i=0;i<n_slip_systems;i++){
        for(unsigned int j=0;j<n_slip_systems;j++){
            q[i][j] = this->userInputs.latentHardeningRatio;
        }
    }

    for(unsigned int i=0;i<n_slip_systems;i++){
        q[i][i] = 1.0;
    }

    //Elastic Stiffness Matrix Dmat
    Dmat.reinit(6,6); Dmat=0.0;

    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=0;j<6;j++){
            Dmat[i][j] = this->userInputs.elasticStiffness[i][j];
        }
    }


    for(unsigned int i=0;i<6;i++){
        for(unsigned int j=3;j<6;j++){
            Dmat[i][j] = 2*Dmat[i][j];
        }
    }

    Vector<double> s0_init (n_slip_systems),rot_init(dim),rotnew_init(dim);
    std::vector<double> twin_init(this->userInputs.numTwinSystems),slip_init(this->userInputs.numSlipSystems);

    for (unsigned int i=0;i<this->userInputs.numSlipSystems;i++){
        s0_init(i)=this->userInputs.initialSlipResistance[i];
    }

    for (unsigned int i=0;i<this->userInputs.numTwinSystems;i++){
        s0_init(i+this->userInputs.numSlipSystems)=this->userInputs.initialSlipResistanceTwin[i];
    }


    for (unsigned int i=0;i<this->userInputs.numSlipSystems;i++){
        slip_init[i]=0.0;
    }

    for (unsigned int i=0;i<this->userInputs.numTwinSystems;i++){
        twin_init[i]=0.0;
    }

    for (unsigned int i=0;i<dim;i++){
        rot_init(i)=0.0;
        rotnew_init(i)=0.0;
    }


    //Resize the vectors of history variables
    Fp_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_conv.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_conv.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init));
    Fp_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    Fe_iter.resize(num_local_cells,std::vector<FullMatrix<double> >(num_quad_points,IdentityMatrix(dim)));
    s_alpha_iter.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,s0_init));
    twinfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
    twinfraction_iter_Twin.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
    slipfraction_iter.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,slip_init));
    twinfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,twin_init));
    twinfraction_conv_Twin.resize(num_local_cells, std::vector<std::vector<double> >(num_quad_points, twin_init));
    slipfraction_conv.resize(num_local_cells,std::vector<std::vector<double> >(num_quad_points,slip_init));
    rot.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rot_init));
    rotnew.resize(num_local_cells,std::vector<Vector<double> >(num_quad_points,rotnew_init));
    twin_ouput.resize(num_local_cells,std::vector<double>(num_quad_points,0.0));
    twin.resize(num_local_cells,std::vector<double>(num_quad_points,0.0));

    //load rot and rotnew
    for (unsigned int cell=0; cell<num_local_cells; cell++){
      unsigned int materialID=cellOrientationMap[cell];
        for (unsigned int q=0; q<num_quad_points; q++){
            for (unsigned int i=0; i<dim; i++){
                rot[cell][q][i]=orientations.eulerAngles[materialID][i];
                rotnew[cell][q][i]=orientations.eulerAngles[materialID][i];
            }
        }
    }
    N_qpts=num_quad_points;
    initCalled=true;


}

#include "../../../include/crystalPlasticity_template_instantiations.h"
