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
import random

from Constant import BOLKEV

def calcprop(tion, ks, params, cstat):
	"""
		return ks
	"""

	# AKK = CONJG(C_k) * C_k
	Akk = ( ks["C_t"][tion, cstat].conjugate() * ks["C_t"][tion, cstat] ).real

	# Bkm = REAL(CONJG(Akm) * Ckm)
	ks["Bkm"] = 2. * ( ks["C_t"][tion, cstat].conjugate() * ks["C_t"][tion, :] 
		                * ks["NAcoup"][tion, cstat, :] ).real

	# P(k -> m) = 2 Re [ CONJG(C_k) * C_m * d_km ] / ( CONJG(C_k) * C_k )
	ks["sh_prop"][tion, :] = ks["Bkm"] / Akk * params["POTIM"]

	kbT = params["Temp"] * BOLKEV

	if params["eleORhole"] == "hole":
		for i in range(ks["ndim"]):
			dE = 0.5 * ( (ks["eigKs"][tion, cstat] + ks["eigKs"][tion+1, cstat]) - \
				         (ks["eigKs"][tion, i] + ks["eigKs"][tion+1, i]) )

			if dE > 0:
				ks["sh_prop"][tion, i] *= np.exp(-dE/kbT)

	elif params["eleORhole"] == "ele":
		for i in range(ks["ndim"]):
			dE = 0.5 * ( (ks["eigKs"][tion, i] + ks["eigKs"][tion+1, i]) - \
				         (ks["eigKs"][tion, cstat] + ks["eigKs"][tion+1, cstat]) )

			if dE > 0:
				ks["sh_prop"][tion, i] *= np.exp(-dE/kbT)

	for i in range(ks["ndim"]):
		if ks["sh_prop"][tion, i] < 0:
			ks["sh_prop"][tion, i] = 0

	return ks

def whichToHop(tion, ks):
	"""
		return which
	"""
	# creat a ramdom number production
	rand  = random.random()
	which = 0

	for i in range(ks["ndim"]):
		if i == 0:
			lower = 0
			upper = ks["sh_prop"][tion, i]
		else:
			lower = upper
			upper = upper + ks["sh_prop"][tion, i]

		if ( lower <= rand ) and ( rand < upper ):
			which = i + 1
			break

	return which

def runSH(ks, params):
	"""
		return ks
	"""

	istat = params["INIBAND"] - params["BandMin"] # note this index, start 0

	# Ensemble of independent trajectories have same coefficients Cj(t).
	for i in range(params["num_sh_traj"]):
		
		# in the first step, current step always equal initial step
		cstat = istat

		for tion in range(params["NamdTime"]-1):
			ks    = calcprop(tion=tion, ks=ks, params=params, cstat=cstat)
			which = whichToHop(tion=tion, ks=ks)

			if which > 0:
				cstat = which - 1

			ks["sh_pops"][tion, cstat] += 1

	ks["sh_pops"] = ks["sh_pops"] / params["num_sh_traj"]

	return ks
