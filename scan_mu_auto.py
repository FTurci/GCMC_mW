import numpy as np
import fileinput as fi
import subprocess
import os

mu=5.0
while mu<10:
	if mu==5.0:
		subprocess.call('make')
		subprocess.call('./mW')
		"""for line in fi.input('mW.cpp', inplace = True):
			line = line.replace('#define mu '+str(mu), '#define mu '+str(-1*mu))
			print line,
		for line in fi.input('INPUT', inplace=True):
			line = line.replace('FILE_PATH ../March_2018/19_3_18/Scan_mu_298/mu_pos_'+str(mu)+'/','FILE_PATH ../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu)+'/')
			print line,
		os.makedirs('../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu))
		subprocess.call('make')
		subprocess.call('./mW')"""
		
	else:
		for line in fi.input('mW.cpp', inplace = True):
			line = line.replace('#define mu '+str(-1*mu+0.5), '#define mu '+str(-1*mu))
			print line,
		for line in fi.input('INPUT', inplace=True):
			line = line.replace('FILE_PATH ../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu-0.5)+'/','FILE_PATH ../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu)+'/')
			print line,
		os.makedirs('../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu))
		subprocess.call('make')
		subprocess.call('./mW')
		"""for line in fi.input('mW.cpp', inplace = True):
			line = line.replace('#define mu '+str(mu), '#define mu '+str(-1*mu))
			print line,
		for line in fi.input('INPUT', inplace=True):
			line = line.replace('FILE_PATH ../March_2018/19_3_18/Scan_mu_298/mu_pos_'+str(mu)+'/','FILE_PATH ../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu)+'/')
			print line,
		os.makedirs('../March_2018/19_3_18/Scan_mu_298/mu_neg_'+str(mu))
		subprocess.call('make')
		subprocess.call('./mW')"""
	
	mu+=0.5
	

