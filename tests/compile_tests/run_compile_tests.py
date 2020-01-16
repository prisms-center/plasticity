'''
Used for testing PRISMS-palsticity code after making changes to verify consistency with previous version.
run as: python run_compile_tests.py numprocs
numprocs is the number of processes for running test case
'''
import numpy as np
import sys
import os,subprocess,shutil,fnmatch
import re,itertools

np.seterr(divide='ignore', invalid='ignore')

if( len(sys.argv) - 1 > 0 ):
    numprocs = int(sys.argv[1])
else:
    numprocs = 1

Error = float(1e-2)

rootDirectory = os.getcwd()
applicationsDirectory = rootDirectory+'/../../applications/crystalPlasticity'
testDirectories = ['/fcc/simpleTension/'
                  ,'/bcc/simpleTension/'
                  ,'/fcc/compression/'
                  ,'/fcc/CyclicLoading/'
                  ,'/fcc/DeformationGraidientBCs/'
                  ,'/fcc/FCC_RandomOrientationBlock/'
                  ,'/fcc/simpleCompression_Multiphase/'
                  ,'/fcc/simpleTension_ExternalMesh/'
                  ,'/hcp/compression/'
                  ,'/hcp/simpleTension/'
                  ] 

for item in testDirectories:

    runtimeOutputFile = rootDirectory+item+'testing.out'
    os.chdir(applicationsDirectory+item)

    subprocess.call(['mpirun -np %d ../../main prm.in' %numprocs],stdout=open(runtimeOutputFile,'w'), stderr=subprocess.STDOUT,shell=True)

    old_stress_strain = np.loadtxt(rootDirectory+item+'correct_stressstrain.txt',skiprows=1)
    new_stress_strain = np.loadtxt(applicationsDirectory+item+'results/stressstrain.txt',skiprows=1)

    old_vonmisesStrain=np.sqrt(0.5*(np.power(old_stress_strain[:,0]-old_stress_strain[:,1],2)+np.power(old_stress_strain[:,1]-old_stress_strain[:,2],2)+np.power(old_stress_strain[:,2]-old_stress_strain[:,0],2)+6*(np.power(old_stress_strain[:,3],2)+np.power(old_stress_strain[:,4],2)+np.power(old_stress_strain[:,5],2))))
    old_vonmisesStress=np.sqrt(0.5*(np.power(old_stress_strain[:,6]-old_stress_strain[:,7],2)+np.power(old_stress_strain[:,7]-old_stress_strain[:,8],2)+np.power(old_stress_strain[:,8]-old_stress_strain[:,6],2)+6*(np.power(old_stress_strain[:,9],2)+np.power(old_stress_strain[:,10],2)+np.power(old_stress_strain[:,11],2))))
    new_vonmisesStrain=np.sqrt(0.5*(np.power(new_stress_strain[:,0]-new_stress_strain[:,1],2)+np.power(new_stress_strain[:,1]-new_stress_strain[:,2],2)+np.power(new_stress_strain[:,2]-new_stress_strain[:,0],2)+6*(np.power(new_stress_strain[:,3],2)+np.power(new_stress_strain[:,4],2)+np.power(new_stress_strain[:,5],2))))
    new_vonmisesStress=np.sqrt(0.5*(np.power(new_stress_strain[:,6]-new_stress_strain[:,7],2)+np.power(new_stress_strain[:,7]-new_stress_strain[:,8],2)+np.power(new_stress_strain[:,8]-new_stress_strain[:,6],2)+6*(np.power(new_stress_strain[:,9],2)+np.power(new_stress_strain[:,10],2)+np.power(new_stress_strain[:,11],2))))

    if (item=='/hcp/compression/' or item=='/hcp/simpleTension/'):
        CheckNumber=350
    else:
        CheckNumber=1000

    CheckStrainError=np.max(np.nan_to_num(np.abs((old_vonmisesStrain[0:CheckNumber]-new_vonmisesStrain[0:CheckNumber])/old_vonmisesStrain[0:CheckNumber])))
    CheckStressError=np.max(np.nan_to_num(np.abs((old_vonmisesStress[0:CheckNumber]-new_vonmisesStress[0:CheckNumber])/old_vonmisesStress[0:CheckNumber])))
    CheckTwinError=np.max(np.nan_to_num(np.abs((old_stress_strain[0:CheckNumber,13]-new_stress_strain[0:CheckNumber,13])/old_stress_strain[0:CheckNumber,13])))
    CheckSlipError=np.max(np.nan_to_num(np.abs((old_stress_strain[0:CheckNumber,14]-new_stress_strain[0:CheckNumber,14])/old_stress_strain[0:CheckNumber,14])))

    if(CheckStressError<Error and CheckStrainError<Error and CheckSlipError<Error):
        print('%s worked successfully\n' %item)
        print('%f Stress Relative Error\n' %CheckStressError)
        print('%f Strain Relative Error\n' %CheckStrainError)
        print('%f Twin Volume Relative Error\n' %CheckTwinError)
        print('%f Slip Relative Error\n' %CheckSlipError)
    else:
        print('%s failed\n' %item)
        print('%s check its stress-strain curve\n' %item)
        print('%f Stress Relative Error\n' %CheckStressError)
        print('%f Strain Relative Error\n' %CheckStrainError)
        print('%f Twin Volume Relative Error\n' %CheckTwinError)
        print('%f Slip Relative Error\n' %CheckSlipError)
