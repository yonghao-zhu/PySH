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
	4. from: VaspBandUnfolding, author: Qijing Zheng
"""

'''
Physical constants used in VASP
'''

#  Some important Parameters, to convert to a.u.
#  - AUTOA     =  1. a.u. in Angstroem
#  - RYTOEV    =  1 Ry in Ev
#  - EVTOJ     =  1 eV in Joule
#  - AMTOKG    =  1 atomic mass unit ("proton mass") in kg
#  - BOLKEV    =  Boltzmanns constant in eV/K
#  - BOLK      =  Boltzmanns constant in Joule/K

AUTOA    = 0.529177249
RYTOEV   = 13.605826
CLIGHT   =  137.037          # speed of light in a.u.
EVTOJ    = 1.60217733E-19
AMTOKG   = 1.6605402E-27
BOLKEV   = 8.6173857E-5
BOLK     = BOLKEV * EVTOJ
EVTOKCAL = 23.06
hbar     = 0.6582119281559802

# FELECT    =  (the electronic charge)/(4*pi*the permittivity of free space)
#         in atomic units this is just e^2
# EDEPS    =  electron charge divided by the permittivity of free space
#         in atomic units this is just 4 pi e^2
# HSQDTM    =  (plancks CONSTANT/(2*PI))**2/(2*ELECTRON MASS)

PI     = 3.141592653589793238
TPI    = 2 * PI
CITPI  = 1j * TPI
FELECT = 2 * AUTOA * RYTOEV
EDEPS  = 4 * PI * 2 * RYTOEV * AUTOA
HSQDTM = RYTOEV * AUTOA * AUTOA

# vector field A times momentum times e/ (2 m_e c) is an energy
# magnetic moments are supplied in Bohr magnetons
# e / (2 m_e c) A(r) p(r)    =  energy
# e / (2 m_e c) m_s x ( r - r_s) / (r-r_s)^3 hbar nabla    = 
# e^2 hbar^2 / (2 m_e^2 c^2) 1/ lenght^3    =  energy
# conversion factor from magnetic moment to energy
# checked independently in SI by Gilles de Wijs

MAGMOMTOENERGY  = 1 / CLIGHT**2 * AUTOA**3 * RYTOEV

# dimensionless number connecting input and output magnetic moments
# AUTOA e^2 (2 m_e c^2)

MOMTOMOM   = AUTOA / CLIGHT / CLIGHT / 2
AUTOA2   = AUTOA * AUTOA
AUTOA3   = AUTOA2 * AUTOA
AUTOA4   = AUTOA2 * AUTOA2
AUTOA5   = AUTOA3 * AUTOA2

# dipole moment in atomic units to Debye
AUTDEBYE = 2.541746  

