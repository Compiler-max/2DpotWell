import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter




def get_example(use_angstrom=False,use_eV=False):
	#get example data to play around
	aX	= 10
	aY	= 10
	#
	nAt		= 2
	relXpos = np.array(		[0.25,0.75]			) 
	relYpos = np.array(		[0.5,0.5]			) 
	atRx	= np.array(		[1.0,2.0]			)
	atRy	= np.array(		[0.5,3.0]			)
	atPot	= np.array(		[-1.0,-1.5]			)
	dVpot	= np.array(		[.0,.0]				)
	#
	#
	#HANDLE UNITS
	aUtoAngstrm = 0.52917721092
	aUtoEV		= 27.211385
	#
	if use_angstrom:
		aX 		= aUtoAngstrm * aX
		aY 		= aUtoAngstrm * aY
		atRx 	= aUtoAngstrm * atRx
		atRy 	= aUtoAngstrm * atRy
	#
	if use_eV:
		atPot = aUtoEV * atPot
		dVpot = aUtoEV * dVpot
	#
	return aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot



def read_line(line,dtype=float,count=1):
	#extracts all numbers from string
	#
	#extract value(s)
	string 	= line.split('=')[-1]
	raw		= np.fromstring(string,dtype=dtype,count=count,sep=" ")
	#
	#make scalars scalar
	if raw.size == 1:
		raw	= raw[-1]
	#
	return raw


def get_input_data(fpath,use_angstrom=False,use_eV=False):
	#traverse the input file and read the following vars
	nAt 	= 0
	aX		= 0
	aY		= 0
	relXpos = []
	relYpos = []
	atRx	= []
	atRy	= []
	atPot	= []
	dVpot	= []	
	#read file
	with open(fpath,'r') as f:
		for counter,line in enumerate(f):
			#lattice constants
			if 'aX' in line:
				aX	= read_line(line,dtype=float,count=1)
			if 'aY' in line:
				aY	= read_line(line,dtype=float,count=1)
			#
			#atom count
			if 'nAt' in line:
				nAt	= read_line(line,dtype=int,count=1)
			#
			#atom containers
			if nAt > 0:
				if 'relXpos' in line:
					relXpos 	= 			read_line(line,dtype=float,count=nAt)
				if 'relYpos' in line:
					relYpos 	= 			read_line(line,dtype=float,count=nAt)
				if 'atRx' in line:
					atRx 		=  			read_line(line,dtype=float,count=nAt)
				if 'atRy' in line:
					atRy 		=  			read_line(line,dtype=float,count=nAt)
				if 'atPot' in line:
					atPot		=			read_line(line,dtype=float,count=nAt)
				if 'dVpot' in line:
					dVpot		=			read_line(line,dtype=float,count=nAt)
	#
	#
	#HANDLE UNITS
	aUtoAngstrm = 0.52917721092
	aUtoEV		= 27.211385
	#
	if use_angstrom:
		aX 		= aUtoAngstrm * aX
		aY 		= aUtoAngstrm * aY
		atRx 	= aUtoAngstrm * atRx
		atRy 	= aUtoAngstrm * atRy
	#
	if use_eV:
		atPot = aUtoEV * atPot
		dVpot = aUtoEV * dVpot



	return aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot




def plot_2D(aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot,length_unit,energy_unit):
	fig2D = plt.figure()
	ax1 = fig2D.add_subplot(111, aspect='equal')
	ax1.set_xlim([0,aX])
	ax1.set_ylim([0,aY])
	#
	for at in range(nAt):
		#set rectangle
		xmin 	= relXpos[at]*aX - atRx[at]
		ymin 	= relYpos[at]*aY - atRy[at]
		dx		= atRx[at]*2.0
		dy		= atRy[at]*2.0
		Vpot	= atPot[at]
		#
		#plot rectangle
		ax1.add_patch(		patches.Rectangle( (xmin, ymin),  dx,  dy, fill=True, alpha=.4, 
													facecolor='blue', edgecolor='black',
													label='well #'+str(at),
													linestyle='solid', linewidth=2.0
											)
					)
		#create energy string
		pot_formatted = str(	"{:4.2f}".format(atPot[at]) )
		print( pot_formatted + energy_unit )
	
		#plot energy string
		ax1.text(xmin+dx*.3, ymin+dy*.4, pot_formatted+'\n'+energy_unit, fontsize=12 )
	
		#plot center indicator
		plt.plot([aX/2.0], [aY/2.0], marker='+', markersize=5, color="black")
	#
	ax1.set_xticks([0.0,aX/2.0,aX])
	ax1.set_yticks([0.0,aY/2.0,aY])
	#
	plt.tick_params(axis='both', which='major', labelsize=11)
	plt.xlabel(r'$x$' + length_unit,fontsize=14)
	plt.ylabel(r'$y$' + length_unit,fontsize=14)
	#
	#plt.legend()
	plt.show()



def plot_3D(aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot,length_unit,energy_unit):
	fig3D		= plt.figure()
	axes			= fig3D.gca(projection='3d')

	mesh_den	= 0.25
	cmap		= cm.viridis
	alpha		= .2

	X			= np.arange(0, aX, mesh_den)
	Y			= np.arange(0, aY, mesh_den)
	X, Y		= np.meshgrid(X, Y)



	xmin		= 1.0
	xmax 		= 4.0
	ymin		= 1.0
	ymax 		= 4.0

	Vpot		= -2.5
	Z1 			= Vpot		*		np.heaviside(X-xmin,1) * (1.0-np.heaviside(X-xmax,1))	*		np.heaviside(Y-ymin,1) * (1.0-np.heaviside(Y-ymax,1))


	xmin		= 6.0
	xmax 		= 8.0
	ymin		= 1.0
	ymax 		= 4.0
	Vpot		= -1.2
	Z2 			= Vpot		*		np.heaviside(X-xmin,1) * (1.0-np.heaviside(X-xmax,1))	*		np.heaviside(Y-ymin,1) * (1.0-np.heaviside(Y-ymax,1))
	for at in range(nAt):
		#set rectangle
		xmin 	= relXpos[at]*aX - atRx[at]
		ymin 	= relYpos[at]*aY - atRy[at]
		dx		= atRx[at]*2.0
		dy		= atRy[at]*2.0
		Vpot	= atPot[at]


		

	# Plot the surface.
	surf1 = axes.plot_surface(X, Y, Z1, color='blue',alpha=alpha,
	                       linewidth=0, antialiased=False)
#
	surf2 = axes.plot_surface(X, Y, Z2, color='blue',alpha=alpha,
	                       linewidth=0, antialiased=False)
	


	Vmin	= -2.5 
	surfA1 	= axes.contourf(X,Y,Z1+Z2,offset=Vmin*2.,cmap=cmap, zmin=Vmin, zmax=0.0)
	#surfA2 	= axes.contourf(X,Y,Z2,offset=Vmin*2.,cmap=cmap, zmin=Vmin, zmax=1.0)


	# Customize the z axis.
	axes.set_xlim(0.0,aX)
	axes.set_ylim(0.0,aY)

	axes.set_zlim(Vpot*2., 1.01)
	axes.zaxis.set_major_locator(LinearLocator(10))
	axes.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))

	axes.set_zticklabels([])



	axes.set_xlabel('x'+length_unit)
	axes.set_ylabel('y'+length_unit)
	
	# Add a color bar which maps values to colors.
	fig3D.colorbar(surfA1, shrink=0.5, aspect=5, label='E'+energy_unit)

	axes.view_init(25, 90+30)
	
	plt.show()

#
#
#------------------------------------------------------------------------------------------------------------------------------------------------

#options
use_angstrom 	= False
use_eV			= False

fpath='./input.orig'




#read input file
aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot = get_input_data(fpath, use_angstrom, use_eV)

#uncomment for example data
#aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot = get_example(fpath, use_angstrom, use_eV)

#set unit strings 
if use_angstrom:
	length_unit = r' ($\AA$)'
else:
	length_unit = ' (a.u.)'

if use_eV:
	energy_unit = ' (eV)'
else:
	energy_unit = r' ($\mathrm{E}_\mathrm{H}$)'



#
#plot_2D(aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot,length_unit,energy_unit)



plot_3D(aX, aY, nAt, relXpos, relYpos, atRx, atRy, atPot, dVpot,length_unit,energy_unit)


