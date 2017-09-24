#PLOTS THE WANNIER FUNCTIONS OVER THE REAL SPACE GRID

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D






#READ IN RAW DATA
f0     		= open("rawData/kpts.dat",'rb')
kpts   		= np.fromfile(f0,dtype='float64',count=-1)
f0.close()

f1	 		= open("rawData/rpts.dat",'rb') #rb = Read Binary
rpts 		= np.fromfile(f1,dtype='float64',count=-1)
f1.close()

f2			= open("rawData/wnfR.dat",'rb') #rb = Read Binary
rawDataR	= np.fromfile(f2,dtype='float64',count=-1)
f2.close()

f3			= open("rawData/wnfI.dat",'rb') #rb = Read Binary
rawDataI	= np.fromfile(f3,dtype='float64',count=-1)
f3.close()

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




	

#RESHAPE RAW DATA
rpts	= np.reshape(rpts,(nR,2))

wnfR	= np.reshape(rawDataR,(nWfs,nK, nR))   #wnF( 	nR, nSC, nWfs		)	
wnfI	= np.reshape(rawDataI,(nWfs,nK, nR))
wnf		= wnfR**2 + wnfI**2


print(wnfI)


#2D HEATMAP
R = 0				#unit cell
n = 0				#state			

for n in range(nWfs):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	bWfcont	=	np.reshape(wnf[n,R,:],(nRy, nRx))
	xpts	= np.linspace(0.0,aX*nKx,nRx)
	ypts	= np.linspace(0.0,aY*nKy,nRy)
	#X, Y = numpy.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
	ax.contour(xpts, ypts,bWfcont,cmap='coolwarm', linewidth=0.1)
	
	#ax.set_xlim(0, nKx*aX)
	#ax.set_ylim(0, nKy*aY)
	
	xticks 		= np.arange(0,aX*(nKx+1),  aX	)		
	xtickLabel	= np.arange(int(0),int(nKx+1)  )
	yticks 		= np.arange(0,aY*(nKy+1),  aY	)		
	ytickLabel	= np.arange(int(0),int(nKy+1)  )
	
	
	ax.set_xticks(xticks)	
	ax.set_xticklabels(xtickLabel,fontsize=12)
	ax.set_yticks(yticks)	
	ax.set_yticklabels(ytickLabel,fontsize=12)
	ax.grid(b=None, axis='both',color='black',alpha=0.2)
	ax.set_xlabel('a')
	ax.set_ylabel('b')
	
	ax.set_title('Wannier function n='+str(n),fontsize =18)
	
	plt.show()








