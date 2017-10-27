#PLOTS THE BLOCH WAVE FUNCTIONS OVER THE REAL SPACE GRID

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation







#READ IN RAW DATA
f0     		= open("rawData/qpts.dat",'rb')
qpts   		= np.fromfile(f0,dtype='float64',count=-1)
f0.close()

f1	 		= open("rawData/rpts.dat",'rb') #rb = Read Binary
rpts 		= np.fromfile(f1,dtype='float64',count=-1)
f1.close()

f2			= open("rawData/unkR.dat",'rb') #rb = Read Binary
rawDataR	= np.fromfile(f2,dtype='float64',count=-1)
f2.close()

f3			= open("rawData/unkI.dat",'rb') #rb = Read Binary
rawDataI	= np.fromfile(f3,dtype='float64',count=-1)
f3.close()

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





f5			= open("rawData/cellInfo.dat",'rb')
cellI 	= np.fromfile(f5,dtype='float64',count=-1)
f5.close()

aX			= cellI[0]
aY			= cellI[1]

f6			= open("rawData/atPos.dat",'rb')
atPos 	= np.fromfile(f6,dtype='float64',count=-1)
f6.close()

f7			= open("rawData/atPos.dat",'rb')
atR 	= np.fromfile(f7,dtype='float64',count=-1)
f7.close()


	

#RESHAPE RAW DATA
rpts	= np.reshape(rpts,(nR,2))
atPos	= np.reshape(atPos,(nAt,2))

qpts	= np.reshape(qpts,(nQ,2))
unkR	= np.reshape(rawDataR,(nQ, nWfs , nR))  #	unk(		nR		, 	nK		, nWfs	)		
unkI	= np.reshape(rawDataI,(nQ, nWfs , nR))
unk		= unkR**2 + unkI**2



k=1
#for n in range(nWfs):
#	print("k="+str(k)+" n="+str(n)+" oLap="  +str( np.sum(unk[k,n,:]) / float(nR) ))

for q in range(nQ):
	for n in range(1):
		print('q='+str(q)+', n=,'+str(n)+', oLap='+str(np.sum(unk[q,n,:])/float(nR)))

		xC 	= np.sum(unk[q,n,:]*rpts[:,0]) / float(nR)
		yC 	= np.sum(unk[q,n,:]*rpts[:,1])	/ float(nR)
		xC	= np.mod(xC,aX)
		yC	= np.mod(yC,aY)
		print('xC='+str(xC)+' yC='+str(yC))





#2D HEATMAP
n = 0				#which state to plot
k = 1				#k point index to plot				

for n in range(6):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	
	unkcont	=	np.reshape(unk[k,n,:],(nRy, nRx))
	xpts	= np.linspace(-25.0,25.0,nRx)
	ypts	= np.linspace(-25.0,25.0,nRy)
	#X, Y = numpy.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
	CS=ax.contourf(xpts, ypts,unkcont,cmap='magma', linewidth=0.1)
	#cbar = plt.colorbar(CS)
	#ax.set_xlim(0, nKx*aX)
	#ax.set_ylim(0, nKy*aY)
	
	xticks 		= np.arange(0,aX*(nQx+1),  aX	)		
	xtickLabel	= np.arange(int(0),int(nQx+1)  )
	yticks 		= np.arange(0,aY*(nQy+1),  aY	)		
	ytickLabel	= np.arange(int(0),int(nQy+1)  )
	
	#plot atoms
	atColor = 'k'
	for at in range(nAt):
		ax.plot(atPos[at,0],atPos[at,1],marker='+',color='k')
		ax.plot([atPos[at,0]-0.5,atPos[at,0]-0.5]	, [atPos[at,1]-0.5,atPos[at,1]+0.5],color=atColor,linewidth=0.2)
		ax.plot([atPos[at,0]+0.5,atPos[at,0]+0.5]	, [atPos[at,1]-0.5,atPos[at,1]+0.5],color=atColor,linewidth=0.2)
		ax.plot([atPos[at,0]-0.5,atPos[at,0]+0.5]	, [atPos[at,1]-0.5,atPos[at,1]-0.5],color=atColor,linewidth=0.2)
		ax.plot([atPos[at,0]-0.5,atPos[at,0]+0.5]	, [atPos[at,1]+0.5,atPos[at,1]+0.5],color=atColor,linewidth=0.2)
		#
		ax.plot([atPos[at,0]-10,atPos[at,0]+10]	, [atPos[at,1],atPos[at,1]],color=atColor,linewidth=0.2)
		ax.plot([atPos[at,0],atPos[at,0]]	, [atPos[at,1]-10,atPos[at,1]+10],color=atColor,linewidth=0.2)
	
	#ax.grid(color='white',linestyle='-',linewidth=1)


	ax.set_xticks(xticks)	
	ax.set_xticklabels(xtickLabel,fontsize=12)
	ax.set_yticks(yticks)	
	ax.set_yticklabels(ytickLabel,fontsize=12)
	ax.grid(b=None, axis='both',color='black',alpha=0.2)
	ax.set_xlabel('a')
	ax.set_ylabel('b')
	
	plt.title('n='+str(n))
		
	ax.set_xlim([0,2*aX])
	ax.set_ylim([0,2*aY])

	plt.show()






## To save the animation, use the command: line_ani.save('lines.mp4')
#fig2 = plt.figure()
#
#x = np.arange(-9, 10)
#y = np.arange(-9, 10).reshape(-1, 1)
#base = np.hypot(x, y)
#ims = []
#for add in np.arange(15):
# 	ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))
#
#im_ani = animation.ArtistAnimation(fig2, ims, interval=100, repeat_delay=3000,
#                                   blit=True)
## To save this second animation with some metadata, use the following command:
## im_ani.save('im.mp4', metadata={'artist':'Guido'})
#
#plt.show()


#3D PLOT
#fig = plt.figure()
#ax = fig.add_subplot(111, projection='3d')
#
##X, Y = numpy.meshgrid(x, y)  # `plot_surface` expects `x` and `y` data to be 2D
#ax.plot_trisurf(rpts[:,0], rpts[:,1],bWf[0,0,:],cmap='coolwarm', linewidth=0.1)
#
#ax.set_xlim(0, 1*aX)
#ax.set_ylim(0, 1*aY)
#
#ax.set_xlabel('X axis')
#ax.set_ylabel('Y axis')
#
#plt.show()
#

