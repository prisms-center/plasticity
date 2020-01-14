'''
Used for testing PRISMS-palsticity code after making changes to verify consistency with previous version.
run as: python run_compile_tests.py numprocs
numprocs is the number of processes for running test case
'''
import numpy as np
import sys
import os,subprocess,shutil,fnmatch
import re,itertools

if( len(sys.argv) - 1 > 0 ):
    numprocs = int(sys.argv[1])
else:
    numprocs = 1

#ncols_totest = 12
strainError = 10^-7
stressError = 10^-5
slipError = 10^-6

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

    #runtimeOutputFile = applicationsDirectory+item+'testing.out'
    #os.chdir(applicationsDirectory+item)

    #subprocess.call(['mpirun -np %d ../../main prm.in' %numprocs],stdout=open(runtimeOutputFile,'w'), stderr=subprocess.STDOUT,shell=True)
    subprocess.call('mpirun -np %d ../../main prm.in' %numprocs, shell=True)
    old_stress_strain = np.loadtxt(rootDirectory+item+'correct_stressstrain.txt',skiprows=1)
    new_stress_strain = np.loadtxt(applicationsDirectory+item+'results/stressstrain.txt',skiprows=1)

    if((np.max(np.abs(old_stress_strain[0:999,:6] - new_stress_strain[0:999,:6]))<strainError) and (np.max(np.abs(old_stress_strain[0:999,6:12] - new_stress_strain[0:999,6:12]))<stressError) and (np.max(np.abs(old_stress_strain[0:999,12:15] - new_stress_strain[0:999,12:15]))<slipError)):
        print('%s worked successfully\n' %item)
    else:
        print('%s failed\n' %item)
        break

    #os.chdir('..')
