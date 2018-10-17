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
"""DL_MONTE_HOME = "~/Documents/GCMC_mW/htk"
sys.path.append(DL_MONTE_HOME)
import htk.util
import htk.histogram"""
colours = ['plum', 'lightblue', 'pink', 'lightgreen', 'khaki', 'tomato', 'silver', 'turquoise', 'orange']
marker_colours = ['purple','blue', 'violet', 'darkgreen', 'gold', 'crimson', 'grey','teal', 'coral']

def find_prob_dist(num_sims, file_path, volume, run):
    
    """ Loads data from multiple simulations (run under same parameters and conditions) and returns
        the probability distribution. """
    
    max_mols = 0
    
    #Find maximum number of molecules across sims
    for sim in range(num_sims):
        input_data = np.genfromtxt(file_path + 'run_' + str(run) + '_'  + str(sim) + '/data.dat', usecols=(2))
        if(max_mols<max(input_data)):
            max_mols = max(input_data)
        
    #Generate histogram arrays for later
    mols_hist = np.zeros((int(max_mols+1), 2))
    total = 0
    
    for sim in range(num_sims):
         #Load data from all sims
        input_data = np.genfromtxt(file_path + 'run_' + str(run) + '_'  + str(sim) + '/data.dat', usecols=(2))
         
        #Generate histogram by loading data and adding 1 to each relevant bin (1 bin per molecule number.
        #Could not use numpy histogram function as did not quite bin correctly (I assume a precision error).
        #This is the biased distribution.
         
        for entry in input_data:
            mols_hist[int(entry),1]+=1
            total+=1

    #Densities associated with molecules calculated
    for entry in range(int(max_mols+1)):
        mols_hist[entry,0] = (entry*180.15)/(volume*6.022)

    #Divide by the total number of entries to find the probability distribution
    mols_hist[:,1]/=total
    
    return max_mols, mols_hist
        
def reweight_mu(num_sims, file_path, volume, mus, mu_org, run):
    
    """ Reweight chemical potential only. Starts by reading in denisty data from a series of given input
        file and calculates the probability distribution. The new distribution is then calculated using the
        formula P_new = exp(-dmu*N)*P_old. Note that the chemical potential is defined as mu/kbT. The output
        is then plotted in a graph, with a maximum number of outputs of 8."""

    if(len(mus)>8):
        print("The maximum number of mus this program can reweight to is 8. Cut off will occur after the 8th my.")
    
    max_mols, mols_hist = find_prob_dist(num_sims,file_path,volume,run)
    shape_y, shape_x= mols_hist.shape

    new_cols = np.zeros((shape_y,len(mus)))
    mols_hist = np.append(mols_hist, new_cols, axis=1)
    
    #Generate new distributions
    for indx,mu in enumerate(mus):
        print("mu is {} indx is {}\n".format(mu,indx))
        dmu = mu_org - mu
        for entry in range(int(max_mols+1)):
            mols_hist[entry,indx+2] = np.exp(-1*dmu*entry)*mols_hist[entry,1]
    
    #Finally, plot the output on one graph and save to file.
    plt.figure()
    plt.plot(mols_hist[:,0],mols_hist[:,1], color=colours[0],marker = '^', markersize=2, markerfacecolor=marker_colours[0], markeredgecolor=marker_colours[0], label = r'$original\ \mu=' + str(mu_org) + '$')
    for indx,mu in enumerate(mus):
        plt.plot(mols_hist[:,0], mols_hist[:,indx+2], color = colours[indx+1],marker = 'o', markersize=2, markerfacecolor=marker_colours[indx+1], markeredgecolor=marker_colours[indx+1], label = r'$\mu=' + str(mu) + '$')

    plt.legend()
    plt.xlabel(r'$\rho (g/cm^3)$')
    plt.ylabel('Probability')
    plt.xlim(0, (max_mols*180.15)/(volume*6.022))
    plt.ylim(0, max(map(max,mols_hist[:,1:num_sims+2])))
    plt.savefig(file_path + '/hist_reweighting_mu.pdf')

def reweight_beta(num_sims, file_path, volume, mus, mu_org, run):
    
    """ Reweights the temperature of the system. Because in the simulation the parameters are defined as
        E/kbT and mu/kbT, the reweighting formula becomes P_new = exp(-(T_old/T_new - 1)*(E-muN))*P_old. """
        
    max_mols, mols_hist = find_prob_dist(num_sims,file_path,volume,run)
    shape_y, shape_x= mols_hist.shape

    new_cols = np.zeros((shape_y,len(mus)))
    mols_hist = np.append(mols_hist, new_cols, axis=1)
    
    #Generate new distributions
    for indx,mu in enumerate(mus):
        print("mu is {} indx is {}\n".format(mu,indx))
        dmu = mu_org - mu
        for entry in range(int(max_mols+1)):
            mols_hist[entry,indx+2] = np.exp(-1*(T_rat-1)*())*mols_hist[entry,1]
    
    #Finally, plot the output on one graph and save to file.
    plt.figure()
    plt.plot(mols_hist[:,0],mols_hist[:,1], color=colours[0],marker = '^', markersize=2, markerfacecolor=marker_colours[0], markeredgecolor=marker_colours[0], label = r'$original\ \mu=' + str(mu_org) + '$')
    for indx,mu in enumerate(mus):
        plt.plot(mols_hist[:,0], mols_hist[:,indx+2], color = colours[indx+1],marker = 'o', markersize=2, markerfacecolor=marker_colours[indx+1], markeredgecolor=marker_colours[indx+1], label = r'$\mu=' + str(mu) + '$')

    plt.legend()
    plt.xlabel(r'$\rho (g/cm^3)$')
    plt.ylabel('Probability')
    plt.xlim(0, (max_mols*180.15)/(volume*6.022))
    plt.ylim(0, max(map(max,mols_hist[:,1:num_sims+2])))
    plt.savefig(file_path + '/hist_reweighting.pdf')
            
reweight_mu(4, '../864.5K_multiprocessing/mu_5.37/', pow(4*1.8*2.3925,3), [-5.35,-5.355, -5.36, -5.365, -5.375, -5.38, -5.385], -5.37, 0)
