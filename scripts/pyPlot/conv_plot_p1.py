import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

from util_2dPW_Interf import get_All_subDirs
from util_2dPW_Interf import get_pw_count



#constants
aUtoAngstrm 	= 0.52917721092
polQuantum = 1.144295347539959e-6


aX = 10.0 *aUtoAngstrm
aY = 5.0 *aUtoAngstrm

#read pol files
descriptor = 'gCut'



def plot_essin(ind ,aX ,aY , pwCount_kx4, pwCount_kx8, pwCount_kx16, pwCount_kx32, pwCount_kx64, pf2_kx4, pf3_kx4, pf2_kx8, pf3_kx8, pf2_kx16, pf3_kx16,pf2_kx32,pf3_kx32,pf2_kx64,pf3_kx64):
	fig, ax  = plt.subplots(1,1) 
	label 	= r'$p^{(1)}_x ( p_Q)$'
	a_ind	= aX

	if ind is not 0:
		ind 	= 1
		label 	= label = r'$p^{(1)}_y ( p_Q)$'
		a_ind	= aY

	


	#define colors
	nMP_runs= 5
	colors = cm.viridis(np.linspace(.25, .75, nMP_runs))     # generate a bunch of colors
	
	print('*')
	print('*')
	print('plot essin info:')
	print('nMP_runs = ',nMP_runs)
	print('a_ind (ang)=',a_ind)
	print('pol factor =',1.0/(polQuantum*a_ind)*100.0)


	if len(pwCount_kx4)>0:
		plt.plot(pwCount_kx4, 		-(	pf2_kx4[:,ind]	+	pf3_kx4[:,ind]	)		/	(polQuantum*a_ind),		'+-'	,label=' 4x8'	,color= colors[0]	) 
	if len(pwCount_kx8)>0:
		plt.plot(pwCount_kx8, 		-(	pf2_kx8[:,ind]	+	pf3_kx8[:,ind]	)		/	(polQuantum*a_ind),		'+-'	,label=' 8x16'	,color= colors[1]	)
	if len(pwCount_kx16)>0:
		plt.plot(pwCount_kx16, 		-(	pf2_kx16[:,ind]	+	pf3_kx16[:,ind]	)		/	(polQuantum*a_ind),		'+-'	,label='16x32'	,color= colors[2]	)
	if len(pwCount_kx32)>0:
		plt.plot(pwCount_kx32, 		-(	pf2_kx32[:,ind]	+	pf3_kx32[:,ind]	)		/	(polQuantum*a_ind),		'+-'	,label='32x64'	,color= colors[3]	)
	if len(pwCount_kx64)>0:
		plt.plot(pwCount_kx64,		-(	pf2_kx64[:,ind]	+	pf3_kx64[:,ind]	)		/	(polQuantum*a_ind),		'+-'	,label='64x128'	,color= colors[4]	)



	plt.xlabel('# basis functions')
	plt.ylabel(label)
	plt.title('magnetoelectric response')
	plt.legend(title='MP grid')

	plt.tight_layout()
	plt.show()




def read_kmesh_folder(path="."):
	pwCount = []
	p0 		= []
	pf2		= []
	pf3 	= []
	#
	if os.path.exists(path):
		gCut, p0, pf2, pf3, p0_interp, pf2_interp, pf3_interp	= 	get_All_subDirs(	descriptor,	path+"/GcutTest/")
		pwCount	=	get_pw_count(path)
	#	
	return pwCount, p0, pf2, pf3





pwCount_kx4, p0_kx4, pf2_kx4, pf3_kx4 		= read_kmesh_folder("./kx4")
pwCount_kx8, p0_kx8, pf2_kx8, pf3_kx8 		= read_kmesh_folder("./kx8")
pwCount_kx16, p0_kx16, pf2_kx16, pf3_kx16 	= read_kmesh_folder("./kx16")
pwCount_kx32, p0_kx32, pf2_kx32, pf3_kx32 	= read_kmesh_folder("./kx32")
pwCount_kx64, p0_kx64, pf2_kx64, pf3_kx64 	= read_kmesh_folder("./kx64")














plot_essin(0,aX, aY, pwCount_kx4, pwCount_kx8, pwCount_kx16, pwCount_kx16, pwCount_kx64, pf2_kx4, pf3_kx4, pf2_kx8, pf3_kx8, pf2_kx16, pf3_kx16, pf2_kx16, pf3_kx16, pf2_kx64, pf3_kx64)
plot_essin(1,aX, aY, pwCount_kx4, pwCount_kx8, pwCount_kx16, pwCount_kx16, pwCount_kx64, pf2_kx4, pf3_kx4, pf2_kx8, pf3_kx8, pf2_kx16, pf3_kx16, pf2_kx16, pf3_kx16, pf2_kx64, pf3_kx64)














#PLOT
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

#read data
#gCut_kx4, p0_kx4, pf2_kx4, pf3_kx4, p0_interp, pf2_interp, pf3_interp	= 	get_All_subDirs(	descriptor,	"./kx4/GcutTest/")
#gCut_kx8, p0_kx8, pf2_kx8, pf3_kx8, p0_interp, pf2_interp, pf3_interp		= get_All_subDirs(	descriptor,	"./kx8/GcutTest/")
#gCut_kx16, p0_kx16, pf2_kx16, pf3_kx16, p0_interp, pf2_interp, pf3_interp	= get_All_subDirs(	descriptor,	"./kx16/GcutTest/")
#gCut_kx32, p0_kx32, pf2_kx32, pf3_kx32, p0_interp, pf2_interp, pf3_interp	= get_All_subDirs(	descriptor,	"./kx32/GcutTest/")
#gCut_kx64, p0_kx64, pf2_kx64, pf3_kx64, p0_interp, pf2_interp, pf3_interp	= get_All_subDirs(	descriptor,	"./kx64/GcutTest/")
#
##get # plane waves used
#pwCount_kx4 	= 	get_pw_count("./kx4")
#pwCount_kx8 	=	get_pw_count("./kx8")
#pwCount_kx16 	=	get_pw_count("./kx16")
#pwCount_kx32 	=	get_pw_count("./kx32")
#pwCount_kx64 	=	get_pw_count("./kx64")
