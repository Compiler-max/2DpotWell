#PLOTS THE BAND STRUCTURE

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D






#READ IN RAW DATA
f0     		= open("rawData/qpts.dat",'rb')
qpts   		= np.fromfile(f0,dtype='float64',count=-1)
f0.close()

f1     		= open("rawData/kpts.dat",'rb')
kpts   		= np.fromfile(f1,dtype='float64',count=-1)
f1.close()



#f1	 		= open("rawData/rpts.dat",'rb') #rb = Read Binary
#rpts 		= np.fromfile(f1,dtype='float64',count=-1)
#f1.close()

f2			= open("rawData/bandStruct.dat",'rb') #rb = Read Binary
rawData		= np.fromfile(f2,dtype='float64',count=-1)
f2.close()


f5			= open("rawData/EnInterP.dat",'rb') #rb = Read Binary
rawInterP	= np.fromfile(f5,dtype='float64',count=-1)
f5.close()

f6			= open("rawData/EnPei.dat",'rb') #rb = Read Binary
rawENpei	= np.fromfile(f6,dtype='float64',count=-1)
f6.close()

f4			= open("rawData/sysPara.dat",'rb')
rawSysP 	= np.fromfile(f4,dtype='int32',count=-1)
f4.close()

nAt			= rawSysP[0]
nG			= rawSysP[1]
nQ			= rawSysP[2]
nQx			= rawSysP[3]
nQy			= rawSysP[4]
nR			= rawSysP[5]
nRx			= rawSysP[6]
nRy			= rawSysP[7]
nWfs		= rawSysP[8]
nSC			= rawSysP[9]		
nSCx		= rawSysP[10]
nSCy		= rawSysP[11]
nK			= rawSysP[12]
nKx			= rawSysP[13]
nKy			= rawSysP[14]
nSolve		= rawSysP[15]


f5			= open("rawData/cellInfo.dat",'rb')
cellI 	= np.fromfile(f5,dtype='float64',count=-1)
f5.close()

aX			= cellI[0]
aY			= cellI[1]




def Qind(qx,qy):
	return int(qy * nQx + qx)


def Kind(kx,ky):
	return	int(ky * nKx + kx)

if( nQx != nQy):
	print("warning the k point spacing per dimension is differnent, this affects the path through k space")

#RESHAPE RAW DATA
qpts	= np.reshape(	qpts		,	(nQ,2)		)
kpts	= np.reshape(	kpts		,	(nK,2)		)
En		= np.reshape(	rawData		,	(nQ,nG)		)    

#rawInterP= np.ones(nK*nWfs)
#rawENpei = np.ones(nK*nWfs)
EI 		= np.reshape(	rawInterP	,	(nK,nWfs)	)
EnP 	= np.reshape(	rawENpei	,	(nK,nWfs)	)


#print(EW)

#PATH THROUGH K SPACE
nPath 	= 4
nQplot	= nPath * nQx #+ nQy + nQy
nKplot	= nPath * nKx
qPlot	= np.linspace(0,nQplot,nQplot)
kPlot	= np.linspace(0,nKplot,nKplot)

EnPlot	= np.empty(nQplot)
EIPlot	= np.empty(nKplot)
EnpPlot	= np.empty(nKplot)

xticks = np.arange(0,nQplot+nQx,nQx)				 #steps in kspace
xtickLabel = np.array([r'$X$',r'$\Gamma$',r'$M$',r'$Y$',r'$\Gamma$']) #symmetry points visited



fig, ax = plt.subplots(1,1)

n 		= 0
#BLOCH STATES
#for each n go along k point path once
for n in range(0,nSolve):
	#X to G
	offs	= 0
	for i in range(0,nQx):	
		ibar			= Qind(nQx-1-i,0)
		EnPlot[offs+i]	= En[ibar,n]
	#G to M
	offs	= nQx		
	for i in range(0,nQx):
		ibar			= Qind(i,i)
		EnPlot[offs+i]	= En[ibar,n]
	#M to Y
	offs	= 2 * nQx
	for i in range(0,nQy):
		ibar			= Qind(nQx-1-i,nQy-1)
		EnPlot[offs+i]	= En[ibar,n]

	#Y to G
	offs	= 3 * nQx
	for i in range(0,nQy):
		ibar			= Qind(0,nQy-1-i)
		EnPlot[offs+i]	= En[ibar,n]

	ax.plot(qPlot,EnPlot,color='k',marker='*',markersize=0.5,linewidth=0.6)
	#ax.plot(qPlot,EnPlot,color='k',marker='*',markersize=0.5,linewidth=0.6)


#PROJECTED STATES
for n in range(0,nWfs):
	#X to G
	offs	= 0
	for i in range(0,nKx):	
		ibar			= Kind(nKx-1-i,0)
		EIPlot[offs+i]	= EI[ibar,n]
		EnpPlot[offs+i] = EnP[ibar,n]
	#G to M
	offs	= nKx		
	for i in range(0,nKx):
		ibar			= Kind(i,i)
		EIPlot[offs+i]	= EI[ibar,n]
		EnpPlot[offs+i] = EnP[ibar,n]
	#M to Y
	offs	= 2 * nKx
	for i in range(0,nKy):
		ibar			= Kind(nKx-1-i,nKy-1)
		EIPlot[offs+i]	= EI[ibar,n]
		EnpPlot[offs+i] = EnP[ibar,n]
	#Y to G
	offs	= 3 * nKx
	for i in range(0,nKy):
		ibar			= Kind(0,nKy-1-i)
		EIPlot[offs+i]	= EI[ibar,n]
		EnpPlot[offs+i] = EnP[ibar,n]
	ax.plot(kPlot,EIPlot,marker='+',color='r',linewidth=0.4)
	ax.plot(kPlot,EnpPlot,marker='+',color='g',linewidth=0.3)



ax.set_xlim([0,nQplot])
#ax.set_ylim([-2,2])
ax.set_xticks(xticks)	
ax.set_xticklabels(xtickLabel,fontsize=14)
ax.grid(b=None, axis='x',color='black',alpha=0.2)
ax.set_xlabel('k space', fontsize = 12)
plt.title('band structure', fontsize= 16)
plt.show()


