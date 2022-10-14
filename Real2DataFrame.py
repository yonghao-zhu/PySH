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
	1. convert real files to a file with DataFrame
	2. input: time-energy; real files
	3. author: yonghao_zhu@163.com
	4. date: 2022.02.18
"""

import pandas as pd
import numpy as np
import os

#====================================================================
# R2D
def Real2DataFrame(params):
	"""
	DataFrame:
		ene1 ene2 dij00-re dij00-im dij01-re dij01-im dij10-re dij10-im dij11-re dij11-im
	1   x    y    x1       x2       x3       x4       x5       x6       x7       x8      
	2
	3
	
	dij matrix (2D):
		x1 x2 x3 x4
		x5 x6 x7 x8

	dij.tolist() = [ [x1, x2, x3, x4], [x5, x6, x7, x8] ]

	dij = np.array( dij.tolist() )

	"""
	
	# read energy file
	energy = np.loadtxt(params["energy_file"])
	times, num_KS_orbitals = energy.shape

	# read real files
	if os.path.isfile("reals.npy"):
		reals = np.load("reals.npy")
	else:
		reals = np.zeros([times, num_KS_orbitals, num_KS_orbitals*2], dtype=float)

		for time in range(params["time_range"][0], params["time_range"][1]+1):
			real_name = params["real_path"] + "real%04d" % time
			# real file
			reals[time-params["time_range"][0], :] = np.loadtxt(real_name, dtype=float)

		np.save("reals.npy", reals)

	print("reals.shape =", reals.shape)

	# creat and write DataFrame
	energy_names = ["ene-%s" %i for i in range(num_KS_orbitals)]
	nacs_names = []
	for i in range(num_KS_orbitals):
		for j in range(num_KS_orbitals):
			for x in ["re", "im"]:
				nacs_names.append( "d%s%s-%s" %(i, j, x) )
	
	col_name = energy_names + nacs_names

	DF = pd.DataFrame( index = [time for time in range(params["time_range"][0], params["time_range"][1]+1)], 
		               columns = col_name )
	
	# write DF
	for time in range(times):
		ene = energy[time, :].tolist()
		nac = []
		for i in reals[time, :].tolist():
			nac += i
		tmp = ene + nac

		DF.loc[params["time_range"][0] + time] = tmp

	if not os.path.isfile(params["csv_file"]):
		print("Saving CSV file...")
		DF.to_csv(params["csv_file"], index_label="time")
	else:
		print("CSV file exist!")

	print(DF)

#====================================================================
# D2R
def DataFrame2Real(params):
	
	DF = pd.read_csv(params["csv_file"])

	times = DF.shape[0]
	num_KS_orbitals = int( (DF.shape[1] - 1) / 5)
	time_lit = []

	energy = np.zeros([times, num_KS_orbitals], dtype=float)
	reals = np.zeros([times, num_KS_orbitals, 2 * num_KS_orbitals])

	for time in range(times):
		tmp = [ i for i in DF.loc[time] ]
		time_lit.append(tmp[0])
		energy[time, :] = np.array(tmp[1:num_KS_orbitals+1])
		reals[time, :] = np.array(tmp[num_KS_orbitals+1:]).reshape(num_KS_orbitals, -1)

	print("energy.shape  =", energy.shape)
	print("reals.shape   =", reals.shape)
	print("len(time_lit) =", len(time_lit))

	if not os.path.isfile("reals.npy"):
		print("Saving reals.npy")
		np.save("reals.npy", reals)
	else:
		print("reals.npy exist! Exiting...")

	# save reals and energy
	# np.savetxt("energy", energy, fmt="%0.6f")
	# for time in range(times):
	# 	real_name = params["real_path"] + "real%04d" % time_list[time]
	# 	np.savetxt(real_name, reals[time, :], fmt="%0.6f")


#====================================================================

def main():

	params = {}
	params["type_"]       = "R2D" # R2D or D2R
	# including, and START 1
	params["time_range"]  = [1, 1900]     # R2D-->construct real**** files
	params["real_path"]   = "./vbm_cbm/"  # R2D-->input, D2R-->output
	# no time column
	params["energy_file"] = "./energy"    # R2D-->input, D2R-->output
	params["csv_file"]    = "EneNACs.csv" # R2D-->output, D2R-->input

	#========================================
	# check file for R2D type
	if params["type_"] == "R2D":
		"""
			times real**** files
		"""
		print("R2D type! File Checking...")

		real = np.loadtxt(params["real_path"]+"real%04d" % params["time_range"][1])
		print("real.shape   =", real.shape)
		energy = np.loadtxt(params["energy_file"])
		print("energy.shape =", energy.shape)
		
		if real.shape[0] != energy.shape[1]:
			print("Shape Error! Exiting...")
			exit(0)
		
		times = params["time_range"][1] - params["time_range"][0] + 1
		print("times        =", times)
		
		if energy.shape[0] != times:
			print("Times Error! Exiting...")
			exit(0)

	elif params["type_"] == "D2R":

		print("D2R type! File Checking...")

		if os.path.isfile(params["csv_file"]):

			DF = pd.read_csv(params["csv_file"])

			print("DF.shape    =", DF.shape)
			print("Times       =", DF.shape[0])
			print("KS Orbitals =", int((DF.shape[1]-1)/5) )
		else:
			print("No CSV file! Exiting...")
			exit(0)

		print("============================")

	else:
		print("ERROR type_! Exiting...")
		exit(0)

	#========================================
	# run
	if params["type_"] == "R2D":

		print("R2D1 type Running...")

		Real2DataFrame(params=params)

		print("R2D DONE!")

		print("============================")

	elif params["type_"] == "D2R":

		print("D2R type Running...")

		DataFrame2Real(params=params)

		print("D2R DONE!")

		print("============================")

	else:
		print("Type Error! Exiting...")
		exit(0)


if __name__ == "__main__":
	main()
