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
    numprocs = sys.argv[1]
else:
    numprocs = 1

ncols_totest = 12

rootDirectory = os.getcwd()
applicationsDirectory = rootDirectory+'/../../applications/crystalPlasticity'
testDirectories = ['/fcc/simpleTension/'
                  ,'/bcc/simpleTension/'
                  ]

for item in testDirectories:

    runtimeOutputFile = rootDirectory + item + 'testing.out'
    os.chdir(applicationsDirectory+item)

    subprocess.call(['mpirun -np %d ../../main prm.in' %numprocs],stdout=open(runtimeOutputFile,'w'), stderr=subprocess.STDOUT,shell=True)

    new_stress_strain = np.loadtxt('results/stressstrain.txt',skiprows=1)
    old_stress_strain = np.loadtxt(rootDirectory+item+'correct_stressstrain.txt',skiprows=1)

    if(np.max(np.abs(old_stress_strain[:,ncols_totest-1] - new_stress_strain[:,ncols_totest-1]))<1e-5):
        print('%s worked successfully\n' %item)
    else:
        print('%s failed\n' %item)

    #os.chdir('..')
