# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import matplotlib.pyplot as plt
import os

# Locate the histogram toolkit
import sys
DL_MONTE_HOME = "/home/mc15800/Documents/GCMC_mW/htk"
sys.path.append(DL_MONTE_HOME)
import htk.util
import htk.histogram


max_mols = 0
total = 0
	
data = np.genfromtxt('../Results/multicanonical/864.5K/mu_5.37/run_0_0/data.dat', usecols=(2))
new_data = np.zeros(data.shape)
max_mols = max(data)
	
mols_hist = np.zeros((int(max_mols+1),3))
mols_hist_new = np.zeros((int(max_mols+1),3))
vol = pow(4*1.8*2.3925,3)

#Densities associated with molecules calculated
for entry in range(int(max_mols+1)):
    mols_hist[entry,2] = entry
    mols_hist[entry,0] = (entry*180.15)/(vol*6.022)
    mols_hist_new[entry,0] = (entry*180.15)/(vol*6.022)


new_data[:,0] = exp(-(beta_1(data[:,3]-mu_1*data[:,2])-beta_0*(data[:,3]-mu_0*data[:,2])))

#Generate histogram by loading data and adding 1 to each relevant bin (1 bin per molecule number.
#Could not use numpy histogram function as did not quite bin correctly (I assume a precision error).
#This is the biased distribution.

for entry in data:
	mols_hist[int(entry),1]+=1
	total+=1
for entry in new_data:
	mols_hist_new[int(entry),1]+=1


mols_hist[:,1]/=total
mols_hist_new[:,1]/=total



#results = htk.histogram.boltzmann_reweight_obs(mols_hist[:,1], 1, -5.37, -5.36, mols_hist[:,2])
print(results)