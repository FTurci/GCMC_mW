""" This program uses the python multiprocessing package to generate a mu-rho
    plot of the mW GCMC simulation. All that must be changed to run this is 
    the filepaths below as well as the list of the range of mus the user
    wishes to sample. 

    Written by Mary Coe, University of Bristol, July 2018.
    e-mail: m.k.coe@bristol.ac.uk
    Last updated October 2018.
"""



import numpy as np
import fileinput as fi
import subprocess
import os
import matplotlib.pyplot as plt
import multiprocessing as mp
import time 
import datetime

def avg_dens(file_path):

	""" Function attempts to open the data.dat file and calculate the average
	    density of the system from the file. If the file is not found, the 
	    function attempts to rerun the simulation. This is repeated to a 
  	    maximum of 5 attempts. The function then returns the mean density. """
	
	dens_found = True
	attempt_counter=0
	
	while dens_found and attempt_counter<5:
		try: 
			dens_array = np.genfromtxt(file_path, usecols= 1)
			dens_found = False
		except IOError:
			attempt_counter+=1
			subprocess.call('./mW')
	if attempt_counter>4:
		print("Error, failed after 5 attempts")

	dens = np.mean(dens_array)
	print("average dens is {}\n".format(dens))
	return dens

def update_mu_infile(mu):
	
	""" This function updates the value of mu in the original C++ program and
	    recompiles it. """
	print("UPDATING MU TO {}\n".format(mu))
	
	for line in fi.input("mW.cpp", inplace=True):
		if line.startswith("#define mu "):
			print(f'#define mu {-1.0*mu}', end='\n')
		else:
			print(f'{line}', end='')

	subprocess.call('make')
	
	return 0

def update_input_file(mu, sim, file_path, mu_diff, method):

	
#    Before each simulation, the input file must be updated to reflect
#	the loading and saving file paths, as well as if the simulation 
#	will load from a previous simulation.
    
    
    for line in fi.input(file_path + "INPUT", inplace = True):
        if (line.startswith("LOAD")):
            print(f'LOAD TRUE', end='\n')
        elif line.startswith("FILE_PATH"):
            fp = file_path + 'mu_' + str(mu) + '_' + str(sim) + '/'
            print(f'FILE_PATH {fp}',end='\n')
        elif line.startswith("IN_PATH"):
            if(method=='decreasing'):
                fp = file_path + 'mu_' + str(mu-mu_diff) + '_' + str(sim) + '/'
                print(f'IN_PATH {fp}', end='\n')
            else:
                fp = file_path + 'mu_' + str(mu+mu_diff) + '_' + str(sim) + '/'
                print(f'IN_PATH {fp}', end='\n')
        else:
            print(f'{line}')

    os.mkdir(file_path + "mu_" + str(mu) + "_" + str(sim))
	
def write_input_file(file_path, load, file_path_in, save_frq, num_gr, num_steps_eqb, num_steps, num_mols, num_cells_row, test):
	""" Before the first simulation is run, the input file must be written to
	    reflect all paramaters of the system. If the file path does not exist 
	    yet, a folder is also made. """

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
	return 0

def sort_array(data, column):

	""" The density data from each mu is stored in an array. Some of these entries
	    may be in the wrong order or certain lines may need to be deleted, hence 
	    this function sorts the data into the correct format to be plotted well. """
	
	data_sorted = data[np.argsort(data[:,column])]
	END = True
	entry = 0
	while END:
	
		if(data[entry][0]<0.001):
			data_sorted = np.delete(data_sorted, entry, 0)
		else:
			entry+=1
	
		if(abs(entry-len(data_sorted))<0.0001):
			END=False
	
	return data_sorted

def plot_mu_dens(file_path, avg_densities_decreasing, avg_densities_increasing):
	
	""" This plotting function plots the increasing and decreasing lines for the
	    mu-rho hysteresis plot. """

	dens_dec = sort_array(avg_densities_decreasing, 0)
	dens_inc = sort_array(avg_densities_increasing, 0)
	
	plt.figure()
	plt.plot(-1*dens_dec[:,0],dens_dec[:,1], color='purple', marker='o', label = 'decreasing 'r'$\mu$')
	plt.plot(-1*dens_inc[:,0],dens_inc[:,1], color='blue', marker = '^',label = 'increasing 'r'$\mu$')
	plt.legend(loc=0)
	plt.ylabel(r'$ \rho\ (g/cm^3)$')
	plt.xlabel(r'$\mu$')
	
	plt.savefig(file_path + 'dens_mu.pdf')

def run_sim(input_file_path, output):
	
	""" This function runs the simulation and records the ouput (the average
	    density from the data.dat file). """
	
	subprocess.call(['./mW',input_file_path])
	dens_array = np.genfromtxt(input_file_path + "data.dat", usecols= 1)
	dens = np.mean(dens_array)
	print("dens is {}\n".format(dens))
	
	output.put(dens)

def run_curve(method, file_path, mus, num_processes, num_sims, save_frq, num_gr, num_steps_eqb,num_steps,num_mols, num_cells_row):
	
	""" This is the main function which manages the running of the curves. It
	    uses the python multiprocessing package to run a number of data points
	    simultaneously and then averages over them to get the final data point
	    for the graph. Because the program uses the computer time as the seed 
	    for the random number generator, for the first run, there is a delay 
	    between starting each simulation. This is not necessary for the 
	    subsequent simulations. """

	prev_mu=0
	inputs = []
	avg_densities = np.zeros([len(mus),2])
	
	for mu in mus:
		output = mp.Queue()
		update_mu_infile(mu)
		time.sleep(2)
		for sim in range(num_sims):
			input_file_path = file_path+"mu_" + str(mu) + "_" + str(sim) + "/"
			if(mus.index(mu)>0):
				if(method=='decreasing'):
					prev_file_path = file_path+"mu_" + str(prev_mu) + "_" + str(sim) + "/"
				else:
					prev_file_path = file_path+"mu_" + str(prev_mu) + "_" + str(sim) + "/"
				write_input_file(input_file_path,True,prev_file_path,save_frq, num_gr,num_steps_eqb, num_steps,num_mols,num_cells_row,False)
			else:
				write_input_file(input_file_path,False,file_path,save_frq, num_gr,num_steps_eqb, num_steps,num_mols,num_cells_row,False)
			inputs.append(input_file_path)
		
		processes = [mp.Process(target=run_sim, args=(inp, output)) for inp in inputs]
		for p in processes:
			p.start()
			if mu==mus[0]:
				time.sleep(5)
		for p in processes:
			p.join()
		density = [output.get() for p in processes]

		avg_densities[mus.index(mu)][0] = mu
		avg_densities[mus.index(mu)][1] = np.mean(density)
		prev_mu=mu
		inputs=[]
	
	return avg_densities

#Specify file path to results folder and parameters to run for. 

file_path = '../March_2019/Produce_mu_dens_graph/800K/decreasing/'
mus_dec = [3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9,5.0,5.1,5.2,5.3,5.4,5.45,5.5,5.55,\
           5.56,5.565,5.7,5.75,5.8,5.9,6.0,6.1,6.2,6.3,6.4,6.5,6.6,6.7,6.8,6.9,7.0,7.1,7.2,7.3,7.4,7.5]
num_processors = 2
num_sims = 2

if not os.path.exists(file_path):
	os.makedirs(file_path)

avg_dens_decreasing = run_curve('decreasing',file_path, mus_dec, num_processors, num_sims, 10, 2, 100,100,0,4) 

#When moving to the increasing mu curve, create new file path and reverse order of mus.
file_path = '../March_2019/Produce_mu_dens_graph/800K/increasing/'

mus_inc = sorted(mus_dec, reverse=True)
tuple(mus_inc)

if not os.path.exists(file_path):
	os.mkdir(file_path)

avg_dens_increasing = run_curve('increasing',file_path, mus_inc, num_processors, num_sims, 10, 2, 100,100,0,4) 

#Finally, for plots, output to general folder.
file_path = '../March_2019/Produce_mu_dens_graph/800K/'

plot_mu_dens(file_path,avg_dens_decreasing, avg_dens_increasing)

with open(file_path + "dens_profile", 'w') as dprof:
	for entry in range(0,len(avg_dens_decreasing)):
		dprof.write("{} {} {}\n".format(avg_dens_decreasing[entry][0],  avg_dens_decreasing[entry][1], avg_dens_increasing[-entry][1]))
		
