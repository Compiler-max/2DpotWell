import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata



subTitleSize = 16



def plotNiuColor(searchDir,filename,plotState):
	fpath	= searchDir+'/'+filename
	
	#read header
	with open(fpath) as f:
		dummy		= f.readline()
		dummy		= f.readline()
		paraString	= f.readline()
		fieldString	= f.readline()
	
	nWfs, nKpts = np.fromstring(paraString,dtype=int,count=2,sep=" ")
	bz			= np.fromstring(fieldString,dtype=float,count=1,sep=" ")[0]
	
	print('detected f2 nWfs =',nWfs)
	print('detected f2 nKpts=',nKpts)
	print('B_z =',bz," (T)")
	
	#read body
	outFile     = open(fpath,'r')
	data        = np.genfromtxt(outFile, dtype=[('int'),('int'),('float64'),('float64'),('float64'),('float64'),('float64'),('float64')],skip_header=4)
	
	kptsX	= np.array([])
	kptsY	= np.array([])
	kptsZ	= np.array([])
	aNiu_x 	= np.array([])
	aNiu_y 	= np.array([])
	aNiu_z 	= np.array([])  
	
	for counter, line in enumerate(data):
		stat = line[0]
		kpt = line[1]
	
		if stat == plotState:
			kptsX	= np.append(kptsX,line[2])
			kptsY	= np.append(kptsY,line[3])
			kptsZ	= np.append(kptsZ,line[4])
	
			aNiu_x	= np.append(aNiu_x,line[5])
			aNiu_y	= np.append(aNiu_y,line[6])
			aNiu_z	= np.append(aNiu_z,line[7])
	
		
	# create x-y points to be used in heatmap
	xi = np.linspace(kptsX.min(),kptsX.max(),len(kptsX))
	yi = np.linspace(kptsY.min(),kptsY.max(),len(kptsY))
	
	# Z is a matrix of x-y values
	aNiuX_plot = griddata((kptsX,kptsY), aNiu_x, (xi[None,:], yi[:,None]), method='cubic')
	aNiuY_plot = griddata((kptsX,kptsY), aNiu_y, (xi[None,:], yi[:,None]), method='cubic')

	
	# I control the range of my colorbar by removing data 
	# outside of my range of interest
	#zmin = 3
	#zmax = 12
	#aNiu2x_zi[(aNiu2x_zi<zmin) | (aNiu2x_zi>zmax)] = None
	


	# Create the contour plot
	zmin 	= min(aNiu_x)
	tmp 	= min(aNiu_y)
	zmin	= min(zmin,tmp) 
	zmax 	= max(aNiu_x)
	tmp		= max(aNiu_y)
	zmax	= max(zmax,tmp)


	# Plot each slice as an independent subplot
	fig, axes = plt.subplots(nrows=1, ncols=2)
	# Make an axis for the colorbar on the right side
	cax = fig.add_axes([.95, 0.1, 0.03, 0.8])	# [left, bottom, width, height] 

	plt.subplot(121)
	CSx 	= plt.contourf(xi, yi, aNiuX_plot, 15, cmap=plt.cm.rainbow,vmax=zmax, vmin=zmin)
	plt.title('x response',fontsize=subTitleSize)
	plt.ylabel('ky (a.u.)')
	plt.xlabel('kx (a.u.)')
	#cbarX	= plt.colorbar(CSx)

	plt.subplot(122)
	CSy 	= plt.contourf(xi,yi, aNiuY_plot, 15, cmap=plt.cm.rainbow,vmax=zmax, vmin=zmin)
	plt.title('y response',fontsize=subTitleSize)
	plt.xlabel('kx (a.u.)')
	#plt.colorbar(CSx) 
	

	
	cbar = fig.colorbar(CSx, cax=cax)
	cbar.set_label('a (ang)')


	descriptor = filename[:-12]
	plt.suptitle(descriptor+': n='+str(plotState)+' Bz='+str(bz)+'T')
	plt.savefig(descriptor+'N'+str(plotState)+'Bz'+str(bz)+'.pdf')
	plt.show()



#test
searchDir 	= '.'
fileName 	= 'f2response.txt'
plotState	= 3
plotNiuColor(searchDir, fileName, plotState)