#PLOTS THE BAND STRUCTURE

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D






#READ IN RAW DATA
f0     		= open("rawData/qpts.dat",'rb')
qpts   		= np.fromfile(f0,dtype='float64',count=-1)
f0.close()

#f1	 		= open("rawData/rpts.dat",'rb') #rb = Read Binary
#rpts 		= np.fromfile(f1,dtype='float64',count=-1)
#f1.close()

f2			= open("rawData/bandStruct.dat",'rb') #rb = Read Binary
rawData		= np.fromfile(f2,dtype='float64',count=-1)
f2.close()


f5			= open("rawData/EnInterP.dat",'rb') #rb = Read Binary
rawInterP	= np.fromfile(f5,dtype='float64',count=-1)
f5.close()

f4			= open("rawData/sysPara.dat",'rb')
rawSysP 	= np.fromfile(f4,dtype='int32',count=-1)
f4.close()

nAt			= rawSysP[0]
nG			= rawSysP[1]
nK			= rawSysP[2]
nKx			= rawSysP[3]
nKy			= rawSysP[4]
nR			= rawSysP[5]
nRx			= rawSysP[6]
nRy			= rawSysP[7]
nWfs		= rawSysP[8]


f5			= open("rawData/cellInfo.dat",'rb')
cellI 	= np.fromfile(f5,dtype='float64',count=-1)
f5.close()

aX			= cellI[0]
aY			= cellI[1]




#integer function getKindex(kx,ky)
#	integer,	intent(in)		:: kx, ky
#	!
#	getKindex = (ky-1) * nKx + kx
#	return
#end
def Kind(kx,ky):
	ibar = int(ky * nKx + kx)
	#print("nx="+str(kx)+", ny="+str(ky)+"ibar="+str(ibar))
	return	ibar




if( nKx != nKy):
	print("warning the k point spacing per dimension is differnent, this affects the path through k space")

#RESHAPE RAW DATA
qpts	= np.reshape(	qpts		,	(nK,2)		)
En		= np.reshape(	rawData		,	(nK,nG)		)    
EI 		= np.reshape(	rawInterP	,	(nK,nWfs)	)

#print(EW)

#PATH THROUGH K SPACE
nPath 	= 4
nKplot	= nPath * nKx #+ nKy + nKy
kPlot	= np.linspace(0,nKplot,nKplot)
EnPlot	= np.empty(nKplot)
EWPlot	= np.empty(nKplot)
EIPlot	= np.empty(nKplot)

xticks = np.arange(0,nKplot+nKx,nKx)				 #steps in kspace
xtickLabel = np.array([r'$X$',r'$\Gamma$',r'$M$',r'$Y$',r'$\Gamma$']) #symmetry points visited



fig, ax = plt.subplots(1,1)

n 		= 0
#BLOCH STATES
#for each n go along k point path once
for n in range(0,nG):
	#X to G
	offs	= 0
	for i in range(0,nKx):	
		ibar			= Kind(nKx-1-i,0)
		EnPlot[offs+i]	= En[ibar,n]
	#G to M
	offs	= nKx		
	for i in range(0,nKx):
		ibar			= Kind(i,i)
		EnPlot[offs+i]	= En[ibar,n]
	#M to Y
	offs	= 2 * nKx
	for i in range(0,nKy):
		ibar			= Kind(nKx-1-i,nKy-1)
		EnPlot[offs+i]	= En[ibar,n]

	#Y to G
	offs	= 3 * nKx
	for i in range(0,nKy):
		ibar			= Kind(0,nKy-1-i)
		EnPlot[offs+i]	= En[ibar,n]

	ax.plot(kPlot,EnPlot,color='k',linewidth=0.4)


#PROJECTED STATES
for n in range(0,nWfs):
	#X to G
	offs	= 0
	for i in range(0,nKx):	
		ibar			= Kind(nKx-1-i,0)
		EIPlot[offs+i]	= EI[ibar,n]
	#G to M
	offs	= nKx		
	for i in range(0,nKx):
		ibar			= Kind(i,i)
		EIPlot[offs+i]	= EI[ibar,n]
	#M to Y
	offs	= 2 * nKx
	for i in range(0,nKy):
		ibar			= Kind(nKx-1-i,nKy-1)
		EIPlot[offs+i]	= EI[ibar,n]
	#Y to G
	offs	= 3 * nKx
	for i in range(0,nKy):
		ibar			= Kind(0,nKy-1-i)
		EIPlot[offs+i]	= EI[ibar,n]
	ax.plot(kPlot,EIPlot,marker='+',color='r',linewidth=0.4)



ax.set_xlim([0,nKplot])
ax.set_xticks(xticks)	
ax.set_xticklabels(xtickLabel,fontsize=14)
ax.grid(b=None, axis='x',color='black',alpha=0.2)
ax.set_xlabel('k space', fontsize = 12)
plt.title('band structure', fontsize= 16)
plt.show()


