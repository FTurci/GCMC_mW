""" This program  is used to explore multicanonical weighting of the GCMC mW program.
    It can be used to run the program to generate data to get new weights (generate_weights)
    and to plot the biased/unbiased density profiles and time series. It also can calculate
    new weights and plot. 

    Written by Mary Coe, University of Bristol, July 2018.
    e-mail: m.k.coe@bristol.ac.uk
    Last updated October 2018.
"""


import matplotlib
matplotlib.use('Agg')
import numpy as np
import fileinput as fi
import subprocess
import os
import matplotlib.pyplot as plt
import multiprocessing as mp
import time 
import math

def plot_biased_unbiased_dens_profile(file_path, vol, num_sims, run):
	
	""" Finds the unbiased and biased density profiles for multiprocessing multicanonical
	    MC simulation (up to 4 sims) and plots. Unbiased probability found by multiplying
	    biased probability by weights (themselves unbiased probabilities from the previous
	    iteration). The unbiased weights are then rescaled. Notes on exact method behind
	    tis to follow.  """

	c = ['plum', 'lightblue', 'lightgreen','pink', 'khaki', 'salmon','crimson', 'wheat']

	max_mols = 0
	y_lim = 0
	total = np.zeros(4)
	
	#The array sizes are set by finding the maximum number of molecules across all simulations
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(2))
		if(max(dens_data)>max_mols):
			max_mols = max(dens_data)
	
	mols_hist = np.zeros((int(max_mols+1),2*num_sims+1))

	#Densities associated with molecules calculated
	for entry in xrange(int(max_mols+1)):
		mols_hist[entry,0] = (entry*180.15)/(vol*6.022)
	
	#Generate histogram by loading data and adding 1 to each relevant bin (1 bin per molecule number.
	#Could not use numpy histogram function as did not quite bin correctly (I assume a precision error).
	#This is the biased distribution.
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(2))
		
		for entry in dens_data:
			mols_hist[int(entry),2*sim+1]+=1
			total[sim]+=1
		
		if(max(mols_hist[:,sim+1])>y_lim):
			y_lim = max(mols_hist[:,2*sim+1])

	#Weights are loaded to be used in calculating unbiased probabilities.
	weights = np.genfromtxt(file_path + "run_" + str(run-1) + "_" + str(sim) + "/weights")
	max_weight = weights[len(weights)-1,1]
	
	#The unbiased probabilities are calculated by multiplying the biased with the weights (or
	#effectively the previous iterations probabilities) and rescaling.
	for sim in xrange(num_sims):
		mols_hist[:,2*sim+1]/=total[sim]
		total[sim] = 0
		for entry in xrange(int(max_mols)):
			if(entry<len(weights)):
				mols_hist[entry,2*(sim+1)] = weights[entry,1]*mols_hist[entry,2*sim+1]
			else:
				mols_hist[entry,2*(sim+1)] = max_weight*mols_hist[entry,2*sim+1]
			total[sim] += mols_hist[entry,2*(sim+1)]
		mols_hist[:,2*(sim+1)] /= total[sim]
	
	#1-4 figures generated in one column and saved to the Graphs folder. 
	fig = plt.figure(1)
	num_graphs = num_sims*100 +11
	for sim in xrange(num_sims):
		graph = num_graphs+sim
		ax1 = plt.subplot(graph)

		if(sim<num_sims-1):
			plt.setp(ax1.get_xticklabels(), visible=False)
		
		plt.semilogy(mols_hist[:-1,0],mols_hist[:-1,2*sim+1],color=c[2*sim],marker='o', label='Biased')
		plt.semilogy(mols_hist[:-1,0], mols_hist[:-1,2*(sim+1)],color=c[2*sim+1], marker = '^',label = 'Unbiased')
		#plt.locator_params(axis='y', nbins= 2)
		plt.xlim(0, (max_mols*180.15)/(vol*6.022))
		plt.legend(loc = 1)
	
	
	plt.xlabel(r'$ \rho\ (g/cm^3)$') 
	time.sleep(5)
	fig.text(0.06, 0.5, "Probability", ha='center', va='center', rotation='vertical')	

	plt.savefig(file_path+"Graphs/dens_profile_" + str(run) + ".pdf")

def plot_initial_dens_profile(file_path, vol, num_sims, run):
	
	""" The first simulation will produce only an unbiased probability distribution hence a
	    simpler function is used. Plots up to 4 graphs in one column. """

	c = ['plum', 'lightblue', 'lightgreen','pink']

	fig = plt.figure(1)
	max_mols = 0
	y_lim = 0
	total = np.zeros(4)
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(2))
		if(max(dens_data)>max_mols):
			max_mols = max(dens_data)
	
	mols_hist = np.zeros((int(max_mols+1),num_sims+1))
	for entry in xrange(int(max_mols+1)):
		mols_hist[entry,0] = (entry*180.15)/(vol*6.022)

	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(2))
		total[sim] = len(dens_data)
		
		for entry in dens_data:
			mols_hist[int(entry),sim+1]+=1
			
		
		if(max(mols_hist[:,sim+1])>y_lim):
			y_lim = max(mols_hist[:,sim+1])

	for sim in xrange(num_sims):
		mols_hist[:,sim+1]/=total[sim]


	num_graphs = num_sims*100 +11
	for sim in xrange(num_sims):
		graph = num_graphs+sim
		ax1 = plt.subplot(graph)

		if(sim<num_sims-1):
			plt.setp(ax1.get_xticklabels(), visible=False)
		
		plt.plot(mols_hist[:,0],mols_hist[:,sim+1],color=c[sim],marker='o')
		plt.locator_params(axis='y', nbins= 2)
		plt.xlim(0, (max_mols*180.15)/(vol*6.022))
	
	
	plt.xlabel(r'$ \rho\ (g/cm^3)$') 
	fig.text(0.06, 0.5, "Probability", ha='center', va='center', rotation='vertical')	

	plt.savefig(file_path+"Graphs/dens_profile_" + str(run) + ".pdf")

def extract_weights(file_path,vol, num_sims,run):

	""" Generates the weights for the next run. First reads in the biased probability distribution
	    and then the previous weights. These are multiplied to get the unbiased distribution. The
	    weights are then rescaled to prevent them getting too small and these become the new weights.
	    The output is written to file and plotted. """
	
	max_mols = 0
	
	#Finds the maximum number of molecules across all simulations to be used to size arrays.
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) + "/data.dat", usecols=(2), dtype=int)
		if(max_mols < max(dens_data)):
			max_mols = max(dens_data)

	weight_data = np.zeros((max_mols+1,2))
	total = 0

	#Biased histogram of density profiles generated.
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) + "/data.dat", usecols=(2), dtype=int)
		total += len(dens_data)
		for entry in dens_data:
			weight_data[entry,1] += 1
		
	weight_data[:,1]/=total
	total = 0

	#If weights were used in simulation, generate unbiased probabilities and rescale
	if run>0:	
		prev_weights = np.genfromtxt(file_path + "run_" + str(run-1) + "_0/weights", usecols=(1))
		max_weight = prev_weights[len(prev_weights)-1,]

		for entry in xrange(max_mols+1):
			weight_data[entry,0] = (entry*180.15)/(vol*6.022)
			if(entry<len(prev_weights)):
				weight_data[entry,1] *= (prev_weights[entry])
			else:
				weight_data[entry,1] *= (max_weight)
			total+=weight_data[entry,1]
		weight_data[:,1]/=total
	
	#else only calculate densities
	else:
		for entry in xrange(max_mols+1):
			weight_data[entry,0] = (entry*180.15)/(vol*6.022)
	
	#Save weights to file
	with open(file_path + "Weights/weights_" + str(run), 'w') as w:
		for entry in xrange(len(weight_data)):
			w.write("{} {}\n".format(weight_data[entry][0],weight_data[entry][1]))
	
	#Plot weights
	plt.figure(1)
	plt.plot(weight_data[:,0],weight_data[:,1],color='pink',marker='o')
	plt.xlabel(r'$ \rho\ (g/cm^3)$')
	plt.ylabel('Weight')
	plt.xlim(0,weight_data[-1,0])
	plt.savefig(file_path+"Weights/weight_trend_" + str(run) + ".pdf")


def plot_time_series(file_path, num_sims, run):

	""" The density of the system as time progresses in the simulation is plotted. This can be used
	    to show that the system is moving more freely between densities and hence is being sampled
	    at coexistence. """
	
	c = ['plum', 'lightblue', 'lightgreen', 'pink','crimson', 'salmon','yellow', 'beige', 'brown']
	fig = plt.figure(1)
	for sim in xrange(num_sims):
		graph = 411+sim
		ax1 = plt.subplot(graph)
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(0,1))
		plt.plot(dens_data[:,0],dens_data[:,1], color = c[sim])
		if(sim<num_sims-1):
			plt.setp(ax1.get_xticklabels(), visible=False)
		plt.locator_params(axis='y', nbins=2)

	plt.xlabel("Move Attempts")
	fig.text(0.06, 0.5, r'$ \rho\ (g/cm^3)$', ha='center', va='center', rotation='vertical')
	plt.savefig(file_path + "Graphs/time_series_" + str(run) + ".pdf")


def write_input_file(file_path, load, file_path_in, save_frq, num_gr, num_steps_eqb, num_steps, num_mols, num_cells_row, test):
	
	""" At the beginning of each simulation, the INPUT file must be written. This sets all parameters (bar
	    the chemical potential which must be done in program). """

	if not os.path.exists(file_path):
		os.mkdir(file_path)
		
	with open(file_path+"INPUT", 'w') as fi:
		if(load):
			fi.write("LOAD TRUE\n")
		else:
			fi.write("LOAD FALSE\n")
		fi.write("FILE_PATH {}\n".format(file_path))
		fi.write("IN_PATH {}\n".format(file_path_in))
		fi.write("SAVE_FREQUENCY {}\n".format(save_frq))
		fi.write("NUMBER_OF_GR {}\n".format(num_gr))
		fi.write("NUM_STEPS_EQUILIBRATION {}\n".format(num_steps_eqb))
		fi.write("NUM_STEPS {}\n".format(num_steps))
		fi.write("NUM_MOLS {}\n".format(num_mols))
		fi.write("NUM_CELLS_ROW {}\n".format(num_cells_row))
		if(test):
			fi.write("TEST TRUE\n")
		else:
			fi.write("TEST FALSE\n")

def run_sim(input_file_path):
	
	#Runs program
	subprocess.call(['./mW',input_file_path])
	return 0


def generate_weights(file_path, num_sims, num_processes, save_frq, num_gr, num_steps_eqb,num_steps,num_mols, num_cells_row,vol, run, load):
	
	""" Creates INPUT files and lists the locations of these files. Then calls the program to run.
	    This function can use multiprocessing, so multiple simulations may be run simultaneously. 
	    Finishes by plotting the density profiles and time series. """
	
	inputs = []
	d=mp.Queue()	#Set up queue for multiprocessing
	for sim in xrange(num_sims):
		
		#Generates INPUT files. Must leave time between to prevent program attempting to 
		#open multiple files and write to multiple files at once. This can lead to errors.
		time.sleep(2)		
		input_file_path = file_path+"run_" + str(run) + "_" + str(sim) + "/"
		prev_file_path = file_path+"run_" + str(run-1) + "_" + str(sim) + "/"
		if(load):			
			write_input_file(input_file_path,True,prev_file_path,save_frq, num_gr,num_steps_eqb, num_steps,num_mols,num_cells_row,False)
		else:
			write_input_file(input_file_path,False,file_path,save_frq, num_gr,num_steps_eqb, num_steps,num_mols,num_cells_row,False)
		inputs.append(input_file_path)
	
	#Simulations to be run loaded into a list.		
	processes = [mp.Process(target=run_sim, args=(inp,)) for inp in inputs]

	#Processes are run. The mW program uses the computer time to seed the random number generator
	#hence a time delay between initiated each process is needed.
	for p in processes:
		print "Starting sim{}\n".format(p)
		p.start()
		time.sleep(5)

	#Must wait for all simulations to end.
	for p in processes:
		p.join()
	
	#Plotting occurs here. Appears to be a problem with writing graphs to file. Require significant time 
	#delay between each plot, or run functions outside of this one. 
	plot_time_series(file_path, num_sims, run)
	time.sleep(10)
	if(run>0):
		plot_biased_unbiased_dens_profile(file_path, vol, num_sims, run)
	else:
		plot_initial_dens_profile(file_path, vol, num_sims, run)

#Call functions here.
	
#generate_weights('../October_2018/Multicanonical/450K_multiprocessing/mu_10.0/', 4, 4, 1000, 10, 1, 100000000, 1, 4, pow(4*1.8*2.3925,3), 5, True)
#time.sleep(10)
#extract_weights('../October_2018/Multicanonical/450K_multiprocessing/mu_10.0/',pow(4*1.8*2.3925,3),4,5)
#time.sleep(10)
plot_biased_unbiased_dens_profile('../../Old_Github/GCMC_mW_old/September_2018/Multicanonical_by_hand/864.5K_multiprocessing/mu_5.37/',pow(4*1.8*2.3925,3),4,1)
