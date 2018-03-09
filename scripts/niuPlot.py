import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata




def plotNiuColor(searchDir,filename,plotState):
	fpath	= searchDir+'/'+fileName
	
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
	pf2_x 	= np.array([])
	pf2_y 	= np.array([])
	pf2_z 	= np.array([])  
	
	for counter, line in enumerate(data):
		stat = line[0]
		kpt = line[1]
	
		if stat == plotState:
			kptsX	= np.append(kptsX,line[2])
			kptsY	= np.append(kptsY,line[3])
			kptsZ	= np.append(kptsZ,line[4])
	
			pf2_x	= np.append(pf2_x,line[5])
			pf2_y	= np.append(pf2_y,line[6])
			pf2_z	= np.append(pf2_z,line[7])
	
		
	# create x-y points to be used in heatmap
	xi = np.linspace(kptsX.min(),kptsX.max(),len(kptsX))
	yi = np.linspace(kptsY.min(),kptsY.max(),len(kptsY))
	
	# Z is a matrix of x-y values
	pf2x_zi = griddata((kptsX,kptsY), pf2_x, (xi[None,:], yi[:,None]), method='cubic')
	
	# I control the range of my colorbar by removing data 
	# outside of my range of interest
	#zmin = 3
	#zmax = 12
	#pf2x_zi[(pf2x_zi<zmin) | (pf2x_zi>zmax)] = None
	
	# Create the contour plot
	CS = plt.contourf(xi, yi, pf2x_zi, 15, cmap=plt.cm.rainbow)#,
	 #                 vmax=zmax, vmin=zmin)
	
	descriptor = fileName[:-12]
	
	plt.title(descriptor+': n='+str(plotState)+' Bz='+str(bz)+'T')
	plt.ylabel('ky (a.u.)')
	plt.xlabel('kx (a.u.)')
	
	plt.colorbar()  
	plt.savefig(descriptor+'Bz'+str(bz)+'.pdf')
	plt.show()
	
	


#test
searchDir 	= '../bin/output/'
fileName 	= 'f2response.txt'
plotState	= 3
plotNiuColor(searchDir, fileName, plotState)