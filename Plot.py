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

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os

out_pop = "out_pop.csv"
Nconds  = 2

if not os.path.isfile(out_pop):
	exit("No csv file, %s! Exiting..." %out_pop)

DF       = pd.read_csv(out_pop)
NamdTime = int(DF.shape[0] / 2)
NKS      = DF.shape[1] - 3

print("NamdTime = ", NamdTime)
print("NKS      = ", NKS)

# get pop and ene
pop = np.zeros([Nconds, NamdTime, NKS], dtype=np.float64)
ene = np.zeros([Nconds, NamdTime], dtype=np.float64)
for ic in range(Nconds):
	for it in range(NamdTime):
		pop[ic, it] = DF.loc[ic*NamdTime+it][3:]
		ene[ic, it] = DF.loc[ic*NamdTime+it][2]

# average for all iconds
pop = np.sum(pop, axis=0) / Nconds

# get time list
time = [i+1 for i in range(NamdTime)]

for ist in range(NKS-1):
	plt.plot(time, pop[:, ist], linewidth=2, linestyle='-')

plt.show()	
