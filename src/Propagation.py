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

from Hamiltonian import make_hamil, hamil_act
from Constant    import hbar

"""
	J. Chem. Theory Comput. 2013, 9, 4595-4972
	equation (11)
"""
def Propagation(ks, params):
	"""
		return ks {dict}
	"""

	edt = params["POTIM"] / params["NELM"]

	# the OUTER loop
	for tion in range(params["NamdTime"]-1):
		ks["pop_t"][tion, :] = ( ks["C_c"].conjugate() * ks["C_c"] ).real
		ks["norm"][tion]     = sum(ks["pop_t"][tion, :])
		ks["C_t"][tion, :]   = ks["C_c"].copy()

		# check the norm of the state
		# print( tion, ks%norm(ks["norm"][tion], ks["psi_c"][:] )

		# the INNER loop
		for tele in range(params["NELM"]-1):
			
			# construct hamiltonian matrix
			# get ks["ham_c"]
			ks = make_hamil(ks=ks, tele=tele, params=params, tion=tion)
			# apply hamiltonian to state vector
			# get ks["hpsi"]
			ks = hamil_act(ks=ks)

			if (tion == 0) and (tele == 0):
				# This is the very first step of the time propagation use first order difference
				# [c,n,p] meas current, next, previous respectively
				ks["C_n"] = ks["C_c"] - 1.j * edt * ks["HC"] / hbar
			else:
				# use second order difference
				# Verlet 
				ks["C_n"] = ks["C_p"] - 2.j * edt * ks["HC"] / hbar

			ks["C_p"] = ks["C_c"].copy()
			ks["C_c"] = ks["C_n"].copy()

	return ks
