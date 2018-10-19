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

def find_max_mols(file_path, num_sims, run):
    
    """Returns the maximum number of molecules across all simulations in one run. """
    
    max_mols = 0
    for sim in range(num_sims):
        input_data = np.genfromtxt(file_path + 'run_' + str(run) + '_'  + str(sim) + '/data.dat', usecols=(2))
        if(max_mols<max(input_data)):
            max_mols = max(input_data)
    
    return max_mols

def find_prob_dist(num_sims, file_path, volume, run, length):
    
    """ Loads data from multiple simulations (run under same parameters and conditions) and returns
        the probability distribution. """
    
    max_mols = find_max_mols(file_path, num_sims, run)
    
    #Generate histogram arrays for later
    mols_hist = np.zeros((int(max_mols+1), length+2))

    for sim in range(num_sims):
         #Load data from all sims
        input_data = np.genfromtxt(file_path + 'run_' + str(run) + '_'  + str(sim) + '/data.dat', usecols=(2))
         
        #Generate histogram by loading data and adding 1 to each relevant bin (1 bin per molecule number.
        #Could not use numpy histogram function as did not quite bin correctly (I assume a precision error).
        #This is the biased distribution.
         
        for entry in input_data:
            mols_hist[int(entry),1]+=1
            
    total = np.sum(mols_hist[:,1])
    
    #Densities associated with molecules calculated
    for entry in range(int(max_mols+1)):
        mols_hist[entry,0] = (entry*180.15)/(volume*6.022)

    #Divide by the total number of entries to find the probability distribution
    mols_hist[:,1]/=total
    
    return mols_hist

def plot_histogram(mols_hist, leg_lab, values, org, out_name, max_mols, volume):
    
    """Plots a histogram of the original and reweighted data. """
    
    plt.figure()
    plt.plot(mols_hist[:,0],mols_hist[:,1], color=colours[0], marker = '^', markersize=2, markerfacecolor=marker_colours[0], markeredgecolor=marker_colours[0],  label = leg_lab + str(org))
    for indx,val in enumerate(values):
       plt.plot(mols_hist[:,0], mols_hist[:,indx+2], color = colours[indx+1],marker = 'o', markersize=2, markerfacecolor=marker_colours[indx+1], markeredgecolor=marker_colours[indx+1], label = leg_lab + str(val) )

    plt.legend()
    plt.xlabel(r'$\rho (g/cm^3)$')
    plt.ylabel('Probability Density Function')
    plt.xlim(0, (max_mols*180.15)/(volume*6.022))
    plt.ylim(0, 1.1*max(map(max,mols_hist[:,1:(len(values)+2)])))
    plt.savefig(out_name)
    

def calc_areas(mols_hist, mus, mu_org, volume):
    
    comparisons = np.zeros((len(mus)+1,5))
    
    comparisons[0,0] = mu_org
    comparisons[0,1] = np.sum(mols_hist[:,1]*mols_hist[:,0])
    break_point = np.argmin(mols_hist[:,1].cumsum() < comparisons[0,1])
    comparisons[0,1] = mols_hist[break_point,0]
    comparisons[0,2] = np.sum(((180.15)/(volume*6.022)*mols_hist[0:break_point,1]))
    comparisons[0,3] = np.sum(((180.15)/(volume*6.022)*mols_hist[break_point:len(mols_hist),1]))
    comparisons[0,4] = comparisons[0,3]-comparisons[0,2]
    
    for indx,mu in enumerate(mus):
        comparisons[indx+1,0] = mu
        comparisons[indx+1,1] = np.sum(mols_hist[:,indx+2]*mols_hist[:,0])
        break_point = np.argmin(mols_hist[:,indx+2].cumsum() < comparisons[indx+1,1])
        comparisons[indx+1,1] = mols_hist[break_point,0]
        comparisons[indx+1,2] = np.sum(((180.15)/(volume*6.022)*mols_hist[0:break_point,indx+2]))
        comparisons[indx+1,3] = np.sum(((180.15)/(volume*6.022)*mols_hist[break_point:len(mols_hist),indx+2]))
        comparisons[indx+1,4] = comparisons[indx+1,3]-comparisons[indx+1,2]
    
    print("{:11s} {:7s} {:7s} {:7s} {:7s}\n".format("mu" , "mean rho", "A_left", "A_right", "Delta A" ))
    for line in range(len(mus)+1):
        print("{:.8f} {:.5f} {:.5f} {:.5f} {:.5f}".format(comparisons[line,0], comparisons[line,1], comparisons[line,2], comparisons[line,3], comparisons[line,4]))
   

def reweight_mu(num_sims, file_path, volume, mus, mu_org, run):
    
    """ Reweight chemical potential only. Starts by reading in denisty data from a series of given input
        file and calculates the probability distribution. The new distribution is then calculated using the
        formula P_new = exp(-dmu*N)*P_old. Note that the chemical potential is defined as mu/kbT. The output
        is then plotted in a graph, with a maximum number of outputs of 8."""


    max_mols = find_max_mols(file_path, num_sims, run)
    mols_hist = find_prob_dist(num_sims, file_path, volume, run, len(mus))
    
    total = np.zeros(len(mus))
 
    for sim in range(num_sims):
        
        input_data = np.genfromtxt(file_path + 'run_' + str(run) + '_'  + str(sim) + '/data.dat', usecols=(2))
        i=0
        #Generate new distributions
        for indx,mu in enumerate(mus):
            print("mu is {} indx is {}\n".format(mu,indx))
            dmu = mu_org - mu
            print("dmu: {} -dmu is {}\n".format(dmu,-1*dmu))
            for entry in input_data:
                mols_hist[int(entry),indx+2] += np.exp(-1.0*dmu*entry)
                print("E: {} dmu: {} W: {}".format(entry, dmu, np.exp(-1.0*dmu*entry)))
                i+=1
                if i>10:
                    break
            total[indx] = np.sum(mols_hist[:,indx+2])
    
    
    for col in range(len(mus)):
        mols_hist[:,col+2]/=total[col]
        
    with open(file_path + "mols_profile", 'w') as dprof:
        for entry in range(0,int(max_mols)+1):
            dprof.write("{:.5f} {:.5f} ".format(mols_hist[entry,0],mols_hist[entry,1]))
            for mu in range(len(mus)):
                dprof.write("{:.5f} ".format(mols_hist[entry,mu+2]))
            dprof.write("\n")
		

    plot_histogram(mols_hist, r'$\mu = $', mus, mu_org, file_path+"reweight_mus_" + str(run) + ".pdf", max_mols, volume)
    calc_areas(mols_hist,mus,mu_org, volume)
    
     
def reweight_beta(num_sims, file_path, volume, T_org, Ts,mu_org, run):
    
    """ Reweights the temperature of the system. Because in the simulation the parameters are defined as
        E/kbT and mu/kbT, the reweighting formula becomes P_new = exp(-(T_old/T_new - 1)*(E-muN))*P_old. """
        
    max_mols = find_max_mols(file_path, num_sims, run)
    
    mols_hist = find_prob_dist(num_sims, file_path, volume, run, len(Ts))
    total = np.zeros(len(Ts))
    
    for sim in range(num_sims):
        input_data = np.genfromtxt(file_path + 'run_' + str(run) + '_'  + str(sim) + '/data.dat', usecols=(2,3))
        
        for entry in input_data[:,0]:
            mols_hist[int(entry),1] +=1
        
        total[0] = np.sum(mols_hist[:,1]*((180.15)/(volume*6.022)))
        
    #Generate new distributions
    for indx,T in enumerate(Ts):
        print("Temperature is {} indx is {}\n".format(T,indx))
        T_rat = T_org/Ts[indx]
        for entry in input_data[:,0]:
            mols_hist[int(entry),indx+2] += np.exp(-1.0*(T_rat-1.0)*(input_data[int(entry),1]-mu_org*input_data[int(entry),0]))
        total[indx] = np.sum(((180.15)/(volume*6.022))*mols_hist[:,indx+2])
    
    
    for col in range(len(Ts)):
        mols_hist[:,col+2]/=total[col]
        
    with open(file_path + "mols_profile_temp_reweight", 'w') as dprof:
        for entry in range(0,int(max_mols)+1):
            dprof.write("{:.5f} {:.5f} ".format(mols_hist[entry,0],mols_hist[entry,1]))
            for T in range(len(Ts)):
                dprof.write("{:.5f} ".format(mols_hist[entry,T+2]))
            dprof.write("\n")
		   
    #Finally, plot the output on one graph and save to file.
    plot_histogram(mols_hist, "T = ", Ts, T_org, file_path + "reweight_Ts_" + str(run) + ".pdf", max_mols, volume)
    
    calc_areas(mols_hist,Ts,T_org, volume)
            
reweight_mu(1, '../Results/multicanonical/864.5K/mu_5.37/', pow(4*1.8*2.3925,3), [-5.36], -5.37, 0)
#reweight_beta(4, '../Results/multicanonical/864.5K/mu_5.37/', pow(4*1.8*2.3925,3), 864.5, [850,855,860,865,870,875,880],-5.37, 0)