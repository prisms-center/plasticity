import numpy as np
import os,subprocess,shutil,fnmatch
import re,itertools



for root,_,_ in os.walk('./', topdown=False):
    if(not(root == './')):
        print(root)
        os.chdir(root)
        print(os.getcwd())

        newdir = 'testing'
        subprocess.call(['../../applications/crystalPlasticity/main prm.in'],stdout=open(newdir+'.out',append_write), stderr=subprocess.STDOUT,shell=True)

        os.chdir('..')
