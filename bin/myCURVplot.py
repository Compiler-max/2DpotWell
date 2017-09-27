#PLOTS THE CONNECTION OVER THE K SPACE GRID

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D



#READ IN RAW DATA
f0     		= open("rawData/qpts.dat",'rb')
kpts   		= np.fromfile(f0,dtype='float64',count=-1)
f0.close()

#f1	 		= open("rawData/rpts.dat",'rb') #rb = Read Binary
#rpts 		= np.fromfile(f1,dtype='float64',count=-1)
#f1.close()


f3			= open("rawData/FcurvR.dat",'rb') #rb = Read Binary
rawDataCurv	= np.fromfile(f3,dtype='float64',count=-1)
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
kpts	= np.reshape(kpts,(nK,2))

Fcurv	= np.reshape(rawDataCurv,(nK, nWfs, 3)) 		#	Fcurv(		3	,	nK	, nWfs	)
print(Fcurv)

#2D HEATMAP
n = 0				#which state to plot
k = 0				#k point index to plot				

for n in range(nWfs):
	fig, ax = plt.subplots(1)#, sharex='col', sharey='row')
	
	Fx		=	np.reshape(Fcurv[:,n,0],(nKy, nKx))
	Fy		=	np.reshape(Fcurv[:,n,1],(nKy, nKx))
	Fz		=	np.reshape(Fcurv[:,n,2],(nKy, nKx))
	
	
	xpts	= 	np.linspace(	-np.pi/aX	,	np.pi/aX	,	nKx		)
	ypts	= 	np.linspace(	-np.pi/aY	,	np.pi/aY	,	nKy		)
	#X, Y = numpy.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
	#strm1	= ax1.streamplot(xpts, ypts, Fx, Fy, color=Fx, linewidth=0.5, cmap=plt.cm.autumn)
	#strm4	= ax1.streamplot(xpts, ypts, Fx, Fz, color=Fz, linewidth=0.5, cmap=plt.cm.autumn)
	#fig.colorbar(strm1.lines)
	
	#ax.contour(xpts, ypts,bWfcont,cmap='coolwarm', linewidth=0.1)
	#ax1.contourf(xpts,ypts, Fx, cmap='coolwarm', linewidth=0.1)
	#ax1.contour( xpts,ypts, Fx, cmap='coolwarm', linewidth=0.1)
	#ax2.contourf(xpts,ypts, Fy, cmap='coolwarm', linewidth=0.1)
	#ax2.contour( xpts,ypts, Fy, cmap='coolwarm', linewidth=0.1)
	

	ax.contour( xpts,ypts, Fz, cmap='coolwarm', linewidth=0.1)
	CS = ax.contourf(xpts,ypts, Fz, cmap='coolwarm', linewidth=0.1)
	cbar = plt.colorbar(CS)
	#ax.set_xlim(0, nKx*aX)
	#ax.set_ylim(0, nKy*aY)
	#xticks 		= np.arange(0,aX*(nKx+1),  aX	)		
	#xtickLabel	= np.arange(int(0),int(nKx+1)  )
	#yticks 		= np.arange(0,aY*(nKy+1),  aY	)		
	#ytickLabel	= np.arange(int(0),int(nKy+1)  )
	
	
	#ax.set_xticks(xticks)	
	#ax.set_xticklabels(xtickLabel,fontsize=12)
	#ax.set_yticks(yticks)	
	#ax.set_yticklabels(ytickLabel,fontsize=12)
	#ax.grid(b=None, axis='both',color='black',alpha=0.2)
	#ax3.set_xlabel(r'$k_x$')
	ax.set_ylabel(r'$k_y$')
	#ax3.set_ylabel(r'$k_y$')
	ax.set_xlabel(r'$k_x$')

	ax.set_title('F_z curvature n='+str(n),fontsize =12)
	
	#ax3.set_title(r'$F_z$',fontsize =12)
	#ax4.set_title(r'$F_z$',fontsize =12)

	
	plt.show()
