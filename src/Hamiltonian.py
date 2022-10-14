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

from Constant import hbar

# initial TDKS
def initC(olap, params):
	"""
		TDKS initialization
		input: olap {dict} and params

		return ks {dict}
	"""
	ks = {}

	ndim        = params["BandMax"] - params["BandMin"] + 1
	ks["ndim"]  = ndim

	# [c,p,n] means current, previous, next
	# C: coefficient
	ks["C_c"] = np.zeros([ ndim ], dtype=complex) 
	ks["C_p"] = np.zeros([ ndim ], dtype=complex) 
	ks["C_n"] = np.zeros([ ndim ], dtype=complex) 
	ks["C_t"] = np.zeros([ params["NamdTime"], ndim ], dtype=complex)

	# the result of hamiltonian acting on a vector
	ks["HC"]  = np.zeros([ ndim ], dtype=complex)

	# population
	ks["pop_t"] = np.zeros([ params["NamdTime"], ndim ], dtype=float) 
	ks["norm"]  = np.zeros([ params["NamdTime"] ], dtype=float)
	ks["ham_c"] = np.zeros([ ndim, ndim ], dtype=complex)

	# KS eigenvalues
	ks["eigKs"] = np.zeros([ params["NamdTime"], ndim ], dtype=float)

	# Non-adiabatic couplings
	ks["NAcoup"]= np.zeros([ params["NamdTime"], ndim, ndim ], dtype=complex)

	# surface hopping related
	# Bkm = REAL(CONJG(Akm) * Ckm)
	ks["Bkm"]     = np.zeros([ ndim ], dtype=float)
	ks["sh_pops"] = np.zeros([ params["NamdTime"], ndim ], dtype=float)
	ks["sh_prop"] = np.zeros([ params["NamdTime"], ndim ], dtype=float)

	# initilization
	ks["C_c"][ params["INIBAND"]-params["BandMin"] ] = 1.0+0.0j

	for i in range(params["NamdTime"]):
		ks["eigKs"][i]  = olap["Eig"][params["NAMDTINI"]+i-1, :]
      	# Divide by 2 * POTIM here, because we didn't do this in the calculation
      	# of couplings
		ks["NAcoup"][i] = olap["Dij"][params["NAMDTINI"]+i-1, :] / (2*params["POTIM"])

	return ks

# constructing the hamiltonian
def make_hamil(ks, tele, params, tion):
	"""
		the hamiltonian contains two parts, 
		which are obtained by interpolation method between two ionic tims step:
		non-adiabatic coupling part and energy eigenvalue part

		return ks
	"""
	# The non-adiabatic coupling part
	# D_{ij}(t) + ( D_{ij}(t+1) - D_{ij}(t) ) * ( TELE / NELM )
	ks["ham_c"] = ks["NAcoup"][tion] + (ks["NAcoup"][tion+1] - ks["NAcoup"][tion]) * \
		                               ((tele+1)/params["NELM"])
	# multiply by -i * hbar
	ks["ham_c"] = -1.j * hbar * ks["ham_c"]

	# the energy eigenvalue part
	# \epsilon _{i}(t) + ( \epsilon _{i}(t+1) - \epsilon _{i}(t) ) * ( TELE / NELM )
	for i in range(ks["ndim"]):
		ks["ham_c"][i, i] = ks["eigKs"][tion, i] + (ks["eigKs"][tion+1, i] - ks["eigKs"][tion, i]) * \
			                                       ((tele+1)/params["NELM"])
	return ks

# Acting the hamiltonian on the state vector
def hamil_act(ks):
	"""
		return ks
	"""
	N = ks["ndim"]

	for i in range(N):
		tmp = complex(0.0, 0.0)
		for j in range(N):
			tmp = tmp + ks["ham_c"][j, i] * ks["C_c"][j]
		ks["HC"][i] = tmp

	return ks
	