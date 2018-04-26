#plots convergence for different mp grids with respect to gCut



import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import os
import os.path
import re

from util_2dPW_Interf import get_All_subDirs
from util_2dPW_Interf import read_gcut_probe_w90
from util_2dPW_Interf import get_pw_count


def plot_p0(ind, pwCount_kx4, pwCount_kx8, pwCount_kx16, wf_pol_kx4, wf_pol_kx8, wf_pol_kx16):
	title 	= 'x polarization'
	label 	= r'$p_x (p_Q)$'
	if ind is not 0:
		ind 	= 1
		title 	= 'y polarization'
		label 	= label = r'$p_y (p_Q)$'



	#kx4
	pol_kx4 = []
	for pol in wf_pol_kx4:
		pol_kx4.append(pol[ind])
	#kx8
	pol_kx8 = []
	for pol in wf_pol_kx4:
		pol_kx8.append(pol[ind])
	#kx16
	pol_kx16 = []
	for pol in wf_pol_kx16:
		pol_kx16.append(pol[ind])
	#kx32
	#pol_kx32 = []
	#for pol in wf_pol_kx32:
	#	pol_kx32.append(pol[ind])
	#kx64
	#pol_kx64 = []
	#for pol in wf_pol_kx64:
	#	pol_kx64.append(pol[ind])


	#color sheme defined by number of k points
	nMP_runs= 3

	colors = cm.viridis(np.linspace(.25, .75, nMP_runs))     # generate a bunch of colors

	#plot p_X
	fig_p0, ax_p0  = plt.subplots(1,1) 
	
	plt.plot(pwCount_kx4, pol_kx4, 			'+-',	color=colors[0],	label ='4x8')
	plt.plot(pwCount_kx8, pol_kx8, 			'+-',	color=colors[1],	label ='8x16')
	plt.plot(pwCount_kx16, pol_kx16, 		'+-',	color=colors[2],	label ='16x32')
	#plt.plot(gCut_kx32, pol_kx32, 		'+-',	color=colors[2],	label ='32x64')
	#plt.plot(gCut_kx64, pol_kx128, 		'+-',	color=colors[2],	label ='64x128')

	#set labels & descriptors
	plt.title(title)
	plt.xlabel('# basis functions')
	plt.ylabel(label)
	plt.legend(title='MP grid')
	
	#show plot
	plt.tight_layout()
	plt.show()















#constants
aUtoAngstrm 	= 0.52917721092
polQuantum = 1.144295347539959e-6
descriptor = 'gCut'





#set sys info
nAt 	= 2
nWfs 	= 6

#get Data
gCut_kx4, at_cent_kx4, mpGrid_kx4, wf_cent_kx4, wf_sprd_kx4, wf_shift_kx4, wf_pol_kx4 = read_gcut_probe_w90(nAt, nWfs, dirpath="./kx4/GcutTest/")
gCut_kx8, at_cent_kx8, mpGrid_kx8, wf_cent_kx8, wf_sprd_kx8, wf_shift_kx8, wf_pol_kx8 = read_gcut_probe_w90(nAt, nWfs, dirpath="./kx8/GcutTest/")
gCut_kx16, at_cent_kx16, mpGrid_kx16, wf_cent_kx16, wf_sprd_kx16, wf_shift_kx16, wf_pol_kx16 = read_gcut_probe_w90(nAt, nWfs, dirpath="./kx16/GcutTest/")


#get # plane waves used
pwCount_kx4 	= 	get_pw_count("./kx4")
pwCount_kx8 	=	get_pw_count("./kx8")
pwCount_kx16 	=	get_pw_count("./kx16")



#gCut_kx32, at_cent_kx32, mpGrid_kx32, wf_cent_kx32, wf_sprd_kx32, wf_shift_kx32, wf_pol_kx32 = read_gcut_probe_w90(nAt, nWfs, dirpath="./kx32ky64/GcutTest/")
#gCut_kx64, at_cent_kx64, mpGrid_kx64, wf_cent_kx64, wf_sprd_kx64, wf_shift_kx64, wf_pol_kx64 = read_gcut_probe_w90(nAt, nWfs, dirpath="./kx64ky128/GcutTest/")


#call plot routines
plot_p0(0, pwCount_kx4, pwCount_kx8, pwCount_kx16, wf_pol_kx4, wf_pol_kx8, wf_pol_kx16)
plot_p0(1, pwCount_kx4, pwCount_kx8, pwCount_kx16, wf_pol_kx4, wf_pol_kx8, wf_pol_kx16)









#spreads
#wf_sprd_sum_kx4 = []
#for sprd in wf_sprd_kx4:
#	wf_sprd_sum_kx4.append(sum(sprd))
#
#wf_sprd_sum_kx8 = []
#for sprd in wf_sprd_kx8:
#	wf_sprd_sum_kx8.append(sum(sprd))
#
#wf_sprd_sum_kx16 = []
#for sprd in wf_sprd_kx16:
#	wf_sprd_sum_kx16.append(sum(sprd))
#plot spread
#fig_sprd, ax_sprd  = plt.subplots(1,1) 
#plt.plot(gCut_kx4, wf_sprd_sum_kx4, label = 'nKx = 4')
#plt.plot(gCut_kx8, wf_sprd_sum_kx8, label = 'nKx = 8')
#plt.plot(gCut_kx16, wf_sprd_sum_kx16, label = 'nKx = 16')
#
#
#
#plt.xlabel('PW cutoff')
#plt.ylabel('total spread')
#
#plt.legend()
#plt.tight_layout()
#plt.show()

#read pol files
#gCut_kx4, p0_kx4, pf2_kx4, pf3_kx4, p0_interp, pf2_interp, pf3_interp	= get_All_subDirs(	descriptor,	"./kx4ky8/GcutTest/")
#
#gCut_kx8, p0_kx8, pf2_kx8, pf3_kx8, p0_interp, pf2_interp, pf3_interp		= get_All_subDirs(	descriptor,	"./kx8ky16/GcutTest/")
##gCut_kx12, p0_kx12, pf2_kx12, pf3_kx12, p0_interp, pf2_interp, pf3_interp	= get_All_subDirs(	descriptor,	"./kx12ky24/GcutTest/")
#gCut_kx16, p0_kx16, pf2_kx16, pf3_kx16, p0_interp, pf2_interp, pf3_interp	= get_All_subDirs(	descriptor,	"./kx16ky32/GcutTest/")
#
#
#aCell= np.array([10.0,5.0,0.0])
#
#aCell = aCell*aUtoAngstrm
#print('acell (ang)=',aCell)
#
#
#
#
#
##plot
#fig_p0, ax_p0  = plt.subplots(1,1) 
#
##plt.plot(gCut_kx4, p0_kx4[:,0],		'+-'	,label='nK =  32')
##plt.plot(gCut_kx8, p0_kx8[:,0],		'+-'	,label='nK = 128')
##plt.plot(gCut_kx12, p0_kx12[:,0],	'+-'	,label='nK = 288')
#plt.plot(gCut_kx16, p0_kx16[:,0]/(polQuantum),	'+-'	,label='nK = 512')
#
#plt.title('Wannier centers x-component')
#plt.ylabel('w_center (ang)')
#plt.xlabel('PW cutoff')
#
#plt.legend()
#
#plt.tight_layout()
#plt.show()
#
#
#
#
##plot
#fig_p0, ax_p0  = plt.subplots(1,1) 
#
#plt.plot(gCut_kx4, p0_kx4[:,1],		'+-'	,label='nK =  32')
#plt.plot(gCut_kx8, p0_kx8[:,1],		'+-'	,label='nK = 128')
##plt.plot(gCut_kx12, p0_kx12[:,1],	'+-'	,label='nK = 288')
#plt.plot(gCut_kx16, p0_kx16[:,1]/(polQuantum),	'+-'	,label='nK = 512')
#
#plt.title(' y-component')
#plt.ylabel('w_center (ang)')
#plt.xlabel('PW cutoff')
#
#plt.legend()
#
#plt.tight_layout()
#plt.show()
#p0_abs = np.sqrt(		p0[:,0]**2 + p0[:,1]**2 + p0[:,2]**2 )
#
#p0_max = []
#tmp_max = np.amax(p0_abs)
#for idx, p0 in enumerate(p0_abs):
#	p0_max.append(tmp_max)
#
#print('p0_max:',p0_max)
#Bz = 1.0
#
#pf2 	= Bz * pf2
#pf3 	= Bz * pf3
#
##p0_x 	= p0[:,0]
##p0_table =  np.concatenate(data,p0_x)
##print(p0_table)
#
#
##PLOT
#direction = ['x-pol', 'y-pol']
#for ind, string in enumerate(direction):
#	fig, ax  = plt.subplots(1,1) 
#	#ax.set_xlim(min(data),max(data))
#	ax.set_xlim(0.0,2.0)
#
#	#plt.plot(data,  p0[:,ind],	marker='+',		color='black',	label='p0'		)
#	#plt.plot(data,	p0_max[:],		color='black',	label='p0_abs')
#	#plt.plot(data,	-1.0*p0_max[:],		color='black',	label='p0_abs')
#	#plt.plot(data, pf2[:,ind],	marker='+',		color='red',	label='p_f2'	)
#	#plt.plot(data, pf3[:,ind],	marker='+',		color='green',	label='p_f3'	)
#	plt.plot(data,pf2[:,ind]+pf3[:,ind], marker='+', color='blue', label='niu')
#	plt.plot(data,-pf2[:,ind]-pf3[:,ind], marker='+', color='orange', label='essin')
#
#
#	#plt.plot(data,pf2[:,ind]-pf3[:,ind], marker='+', color='orange', label='f2-f3')
#
#	#plt.plot(data, pf2_interp[:,ind],	 marker='+',	color='orange',	label='p_f2(I)')
#	#plt.plot(data, pf3_interp[:,ind],	 marker='+',	color='lightgreen',	label='p_f2(I)')
#
#	plt.title('OMP: '+string+' Bz='+str(Bz)+' T')
#	plt.ylabel(r'electric pol. ($\mu C / cm$)')
#	plt.xlabel(descriptor)
#
#	
#
#	plt.legend()
#	plt.show()

