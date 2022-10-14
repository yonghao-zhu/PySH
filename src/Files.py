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

#--------------------------------------------------------------------------------------------
def ReadCSVFile(csv_name, potim):
	"""
	read CSV file
	
	return olap {a dict}
	olap["NBANDS"]: int, the number of ks orbitals
	olap["TSTEPS"]: int, the number of NACs
	olap["dt"]    : float, fs, POTIM for NACs calculating
	olap["Dij"]   : complex, NACs, dimensionless, Dij*hbar/(2*POTIM) = eV
				    hbar = 0.6582119281559802 eV*fs
	olap["Eig"]   : float, eig of Hamiltonian, eV

	DataFrame:
	time eig0 eig1 dij00-re dij00-im dij01-re dij01-im dij10-re dij10-im dij11-re dij11-im
	 1    x    y    x1       x2       x3       x4       x5       x6       x7       x8      
	 2
	 3
	
	dij matrix (2D):
		x1 x2 x3 x4
		x5 x6 x7 x8

	dij.tolist() = [ [x1, x2, x3, x4], [x5, x6, x7, x8] ]

	dij = np.array( dij.tolist() )

	"""
	
	DF              = pd.read_csv(csv_name)
	times           = DF.shape[0]
	num_KS_orbitals = int( (DF.shape[1] - 1) / 5)
	
	# initialization
	ene = np.zeros([times, num_KS_orbitals], dtype=float)
	nac = np.zeros([times, num_KS_orbitals, num_KS_orbitals], dtype=complex)

	# read data from dataframe
	for time in range(times):
		tmp          = [ i for i in DF.loc[time] ]
		ene[time, :] = np.array( tmp[1:num_KS_orbitals+1] )
		dij          = np.array( tmp[num_KS_orbitals+1:] ).reshape( num_KS_orbitals, -1 )
		nac[time, :] = dij[:, 0::2] + 1.0j*dij[:, 1::2]

	print("   --> ene.shape   =", ene.shape)
	print("   --> nac.shape   =", nac.shape)

	# write olap dict
	olap = {}
	olap["NBANDS"] = num_KS_orbitals
	olap["TSTEPS"] = times
	olap["dt"]     = potim
	olap["Dij"]    = nac
	olap["Eig"]    = ene

	return olap

#--------------------------------------------------------------------------------------------
def WriteCSVFile(ks, params, DF_out, DF_C, ic):
	"""
		output the results

		ic start 0
	"""

	for it in range(params["NamdTime"]):
		index_ = ic * params["NamdTime"] + it + 1
		time   = (it+1) * params["POTIM"]
		ene    = sum( ks["eigKs"][it, :] * ks["sh_pops"][it, :] )
		states = [ks["sh_pops"][it, i] for i in range(ks["ndim"])]	

		DF_out.loc[index_] = [time] + [ene] + states

	for it in range(params["NamdTime"]):
		index_ = ic * params["NamdTime"] + it + 1
		time   = (it+1) * params["POTIM"]
		ene    = sum( ks["eigKs"][it, :] * ks["sh_pops"][it, :] )
		C      = [ks["C_t"][it, i] for i in range(ks["ndim"])]

		DF_C.loc[index_] = [time] + [ene] + C
