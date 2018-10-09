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

def plot_dens_trend_with_weights(file_path, vol, num_sims, run):
	
	c = ['plum', 'lightblue', 'lightgreen','khaki', 'crimson', 'salmon','yellow', 'beige', 'brown']
	fig = plt.figure(1)
	max_mols = 0
	y_lim = 0
	total = np.zeros(4)
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(2))
		if(max(dens_data)>max_mols):
			max_mols = max(dens_data)
	
	print "Max mols is {}\n".format(max_mols)
	mols_hist = np.zeros((int(max_mols+1),num_sims+1))
	for entry in xrange(int(max_mols+1)):
		mols_hist[entry,0] = (entry*180.15)/(vol*6.022)

	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) +"/data.dat", usecols=(2))
		
		for entry in dens_data:
			mols_hist[int(entry),sim+1]+=1
			total[sim]+=1
		
		if(max(mols_hist[:,sim+1])>y_lim):
			y_lim = max(mols_hist[:,sim+1])
			print y_lim
	for sim in xrange(num_sims):
		mols_hist[:,sim+1]/=total[sim]

	
	weights = np.genfromtxt(file_path + "run_" + str(run-1) + "_" + str(sim) + "/weights")
	print weights
	max_weight = weights[len(weights)-1,1]
	"""for sim in xrange(num_sims):
		for entry in xrange(int(max_mols)):
			if(entry<len(weights)):
				mols_hist[entry,sim+1] *= (weights[entry,1])
			else:
				mols_hist[entry,sim+1] *= (max_weight)"""

		
	"""div = 1
	powers = 0
	while y_lim>10:
		y_lim /= 10
		div *= 10
		powers += 1"""
	num_graphs = num_sims*100 +11
	for sim in xrange(num_sims):
		graph = num_graphs+sim
		ax1 = plt.subplot(graph)

		if(sim<num_sims-1):
			plt.setp(ax1.get_xticklabels(), visible=False)
		
		plt.plot(mols_hist[:,0],mols_hist[:,sim+1],color=c[sim],marker='o')
		plt.locator_params(axis='y', nbins= 2)
		plt.xlim(0, (max_mols*180.15)/(vol*6.022))
		#plt.ylim(0, math.ceil(y_lim))
		mols_hist[:,1]=0
	
	
	plt.xlabel(r'$ \rho\ (g/cm^3)$')
	#fig.text(0.06, 0.5, "Number of Occurances (x" + str(div) + ")", ha='center', va='center', rotation='vertical')
	fig.text(0.06, 0.5, "Probability", ha='center', va='center', rotation='vertical')	

	plt.savefig(file_path+"Graphs/current_trend_" + str(run) + ".pdf")

def extract_weights(file_path, vol, cut_off, num_sims, run):
	
	max_mols = 0
	
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) + "/data.dat", usecols=(2), dtype=int)
		if(max_mols < max(dens_data)):
			max_mols = max(dens_data)

	weight_data = np.zeros((max_mols+1,2))
	total = 0

	for  sim in xrange(num_sims):	
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) + "/data.dat", usecols=(1))
		max_mols = max(dens_data)
		max_mols*=vol*6.022
		max_mols/=180.15
		hist_data, bin_data = np.histogram(dens_data,bins=int(round(max_mols)))
		for entry in xrange(0,len(bin_data)):
			if(bin_data[entry]<cut_off):		
				weight_data[entry,1] += hist_data[entry]
				total = total + hist_data[entry]
	print weight_data
	weight_data[:,1]/=total
	with open(file_path + "Weights/weights_" + str(run), 'w') as w:
		for entry in xrange(len(weight_data)-1):
			w.write("{} {}\n".format(weight_data[entry][0],weight_data[entry][1]))


	plt.figure(1)
	plt.plot(weight_data[:-1,0],weight_data[:-1,1],color='pink',marker='o')
	plt.xlabel(r'$ \rho\ (g/cm^3)$')
	plt.ylabel('Weight')
	plt.savefig(file_path+"Weights/weight_trend_" + str(run) + ".pdf")
	#print weight_data

def extract_weights_full_dist_initial(file_path,vol, num_sims,run):
	
	max_mols = 0
	
	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) + "/data.dat", usecols=(2), dtype=int)
		if(max_mols < max(dens_data)):
			max_mols = max(dens_data)

	weight_data = np.zeros((max_mols+1,2))
	total = 0

	for sim in xrange(num_sims):
		dens_data = np.genfromtxt(file_path + "run_" + str(run) + "_" + str(sim) + "/data.dat", usecols=(2), dtype=int)
		total += len(dens_data)
		for entry in dens_data:
			weight_data[entry,1] += 1
			
	weight_data[:,1]/=total
	print weight_data

	prev_weights = np.genfromtxt(file_path + "run_" + str(run-1) + "_0/weights", usecols=(1))
	max_weight = prev_weights[len(prev_weights)-1,]
	total = 0

	for entry in xrange(max_mols+1):
		weight_data[entry,0] = (entry*180.15)/(vol*6.022)
		if(entry<len(prev_weights)):
			weight_data[entry,1] *= (prev_weights[entry])
		else:
			weight_data[entry,1] *= (max_weight)
		
	
	print weight_data
	with open(file_path + "Weights/weights_" + str(run), 'w') as w:
		for entry in xrange(len(weight_data)):
			w.write("{} {}\n".format(weight_data[entry][0],weight_data[entry][1]*100))
	
	
	plt.figure(1)
	plt.plot(weight_data[:,0],weight_data[:,1],color='pink',marker='o')
	plt.xlabel(r'$ \rho\ (g/cm^3)$')
	plt.ylabel('Weight')
	plt.xlim(0,weight_data[-1,0])
	plt.savefig(file_path+"Weights/weight_trend_" + str(run) + ".pdf")



def plot_time_series(file_path, num_sims, run):
	
	c = ['plum', 'lightblue', 'lightgreen', 'khaki','crimson', 'salmon','yellow', 'beige', 'brown']
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
	
	subprocess.call(['./mW',input_file_path])
	return 0


def generate_weights(file_path, num_sims, num_processes, save_frq, num_gr, num_steps_eqb,num_steps,num_mols, num_cells_row,vol, run, load):
	
	inputs = []
	d=mp.Queue()
	#for run in xrange(num_runs):
	for sim in xrange(num_sims):
		time.sleep(2)
		input_file_path = file_path+"run_" + str(run) + "_" + str(sim) + "/"
		prev_file_path = file_path+"run_" + str(run-1) + "_" + str(sim) + "/"
		if(load):			
			write_input_file(input_file_path,True,prev_file_path,save_frq, num_gr,num_steps_eqb, num_steps,num_mols,num_cells_row,False)
		else:
			write_input_file(input_file_path,False,file_path,save_frq, num_gr,num_steps_eqb, num_steps,num_mols,num_cells_row,False)
		inputs.append(input_file_path)
			
	processes = [mp.Process(target=run_sim, args=(inp,)) for inp in inputs]
	for p in processes:
		print "Starting sim{}\n".format(p)
		p.start()
		time.sleep(5)
	for p in processes:
		p.join()
	
	plot_time_series(file_path, num_sims, run)
	time.sleep(3)
	plot_dens_trend_with_weights(file_path, vol, num_sims, run)
	
generate_weights('../October_2018/multicanonical_by_hand/450K_multiprocessing/mu_10.0/', 4, 4, 1000, 10, 1, 100000000, 1, 4, pow(4*1.8*2.3925,3), 3, True)
#plot_time_series('../September_2018/Multicanonical_by_hand/864.5K_multiprocessing/', 4, 0)
#plot_dens_trend('../September_2018/Multicanonical_by_hand/864.5K_multiprocessing/', pow(4*1.8*2.3925,3), 4, 0)
#extract_weights('../September_2018/Multicanonical_by_hand/864.5K_multiprocessing/', pow(4*1.8*2.3925,3), 0.2,4,0)
#extract_weights_full_dist_initial('../October_2018/multicanonical_by_hand/450K_multiprocessing/mu_10.0/',pow(4*1.8*2.3925,3),4,2)
#plot_dens_trend_with_weights('../October_2018/multicanonical_by_hand/450K_multiprocessing/mu_10.0/', pow(4*1.8*2.3925,3), 4, 2)

