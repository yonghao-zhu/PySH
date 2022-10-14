#! /usr/bin/python3
# -*- conding=UTF-8 -*-
#  .--,       .--,
# ( (  \.---./  ) )
#  '.__/o   o\__.'
#     {=  ^  =}
#      >  -  <
#     /  Zhu  \
#    //  Yong \\
#   //|  Hao  |\\
#   "'\       /'"_.-~^`'-.
#      \  _  /--'         `
#    ___)( )(___
#   (((__) (__))) 

"""
	1. ref: J. Chem. Phys. 93, 1061(1990)
	        Phys. Rev. Lett. 95, 163001 (2005)
			J. Chem. Theory Comput. 2013, 9, 4595-4972
			https://github.com/QijingZheng/Hefei-NAMD
	2. current version: yonghao_zhu@163.com, 2022-02-20, python code
	3. Acknowledgement: Prof. Qijing Zheng and Jin Zhao
"""

import numpy as np
import pandas as pd
import os
import time

from Files       import ReadCSVFile, WriteCSVFile
from Hamiltonian import initC
from Propagation import Propagation
from Hopping     import runSH

"""
	1. Ensemble of independent trajectories have same coefficients Cj(t).
	2. Internal consistency condition Nj(t)∝Cj*(t)Cj(t)=ρjj(t).
	3. Hops from j to different k ‰ j are independent.
	4. Overall trajectory hops should be minimum.
"""
def NAMD(params):
	"""
		check params in main function
	"""
	#----------------------------step1: read dataframe----------------------------#
	# read nac and ene
	# note: 1. the unit of coupling is not eV.
	#       2. in the namd (no dish), NACoupling is float, no imaginary part
	#          However, this script contains that.

	if not os.path.isfile(params["csv"]):
		exit("No csv file: %s! Exiting..." %params["csv"])

	print("--------< step1: read dataframe >--------")
	olap = ReadCSVFile(csv_name = params["csv"], potim = params["POTIM"])
	#print("   --> olap.keys() =", olap.keys())

	# for testing, del imag of NAC
	olap["Dij"] = np.real( olap["Dij"] )

	#----------------------------step3: loop iconds----------------------------#
	# creat and write DataFrame for output 
	ndim     = params["BandMax"] - params["BandMin"] + 1
	col_name = ["Time", "energy"] + ["state-%s" %(i+1) for i in range(ndim)]
	Index_   = []
	for i in range(len(params["iconds"])):
		for j in range(params["NamdTime"]):
			Index_.append(1+i*params["NamdTime"]+j)
	DF_pop   = pd.DataFrame(index=Index_, columns=col_name)
	DF_C     = pd.DataFrame(index=Index_, columns=col_name)

	# main loop Nconds
	print("--------< step2: loop iconds >--------")
	start_time = time.time()
	for ic in range(params["Nconds"]):
		print("  --> main, ic = ", ic)

		# initiate KS dict, 
		print("  --> initiate KS matrix")
		iniband, namdtini = params["iconds"][ic, 1], params["iconds"][ic, 0]
		params["INIBAND"] = iniband; params["NAMDTINI"] = namdtini
		ks = initC(params=params, olap=olap)
	
		# Time propagation
		print("  --> Running time propagation")
		ks = Propagation(ks=ks, params=params)

		# Run surface hopping
		print("  --> Running surface hopping")
		ks = runSH(ks=ks, params=params)

		# ouput
		WriteCSVFile(ks, params, DF_pop, DF_C, ic)

	DF_pop.to_csv("out_pop.csv")
	if params["save_C"]:
		DF_C.to_csv("out_C.csv")
	print("--------< DONE! >--------")
	end_time = time.time()
	print("--> main loop time: %0.4s s!" %(end_time-start_time))

#------------------------------------------------------------------------
def main():
	
	params              = {}
	# control 
	params["BandMin"]     = 1     # bottom band index
	params["BandMax"]     = 2     # top band index
	params["Nbands"]      = 2     # number of bands
	params["NSW"]         = 2000  # number of ionic steps
	params["POTIM"]       = 1     # MD time step, fs
	params["Temp"]        = 300   # temperature, K
	params["Nconds"]      = 2     # number of initial conditions
	params["NamdTime"]    = 100   # time for NAMD run
	params["NELM"]        = 1000  # electron time step
	# Number of stochastic realizations for each initial condition
	params["num_sh_traj"] = 800   
	params["eleORhole"]   = "ele" # ele or hole
	# nac and energy --> EneNACs.csv           
	params["csv"]         = "../EneNACs.csv"
	# out
	params["save_C"]      = True # True or False, out_C.csv

	# initial condition: initial_time initial_band
	# length: Nconds
	iconds = np.zeros([params["Nconds"], 2], dtype=int)
	# initial band
	if params["eleORhole"] == "ele":
		iconds[:, 1] = params["BandMax"]
	elif pandas["eleORhole"] == "hole":
		iconds[:, 1] = params["BandMin"]
	# initial time
	iconds[:, 0] = np.array([i for i in range(1, params["Nconds"]+1)])

	# for test
	iconds[:, 0] = 2 # [2, 2]

	params["iconds"] = iconds

	# run na-md
	# output: out.csv 
	NAMD(params)

if __name__ == "__main__":
	main()